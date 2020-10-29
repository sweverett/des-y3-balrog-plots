import numpy as np
import os
import healpy as hp
import fitsio
import pymaster as nmt
import matplotlib.pyplot as plt

def get_synthetic_map(bin_means,binned_map_values,map_values):
    # Linear interpolation ought to be enough.
    interp_map_values = np.interp(map_values,bin_means,binned_map_values)
    return interp_map_values



def get_property_map(map_key, nside=2048, mapdir = '/data/des81.a/data/severett/paper-plots/desy3-balrog-plots/systematics-maps-compare/maps'):
    maps = {'fwhm_g': f'y3a2_g_o.{nside}_t.32768_FWHM.WMEAN_EQU.fits.gz',
            'fwhm_r': f'y3a2_r_o.{nside}_t.32768_FWHM.WMEAN_EQU.fits.gz',
            'fwhm_i': f'y3a2_i_o.{nside}_t.32768_FWHM.WMEAN_EQU.fits.gz',
            'fwhm_z': f'y3a2_z_o.{nside}_t.32768_FWHM.WMEAN_EQU.fits.gz',
            'skybrite_g' : f'y3a2_g_o.{nside}_t.32768_SKYBRITE.WMEAN_EQU.fits.gz',
            'skybrite_r' : f'y3a2_r_o.{nside}_t.32768_SKYBRITE.WMEAN_EQU.fits.gz',
            'skybrite_i' : f'y3a2_i_o.{nside}_t.32768_SKYBRITE.WMEAN_EQU.fits.gz',
            'skybrite_z' : f'y3a2_z_o.{nside}_t.32768_SKYBRITE.WMEAN_EQU.fits.gz',
            'sig_zp_g' : f'y3a2_g_o.{nside}_t.32768_SIGMA_MAG_ZERO.QSUM_EQU.fits.gz',
            'sig_zp_r' : f'y3a2_r_o.{nside}_t.32768_SIGMA_MAG_ZERO.QSUM_EQU.fits.gz',
            'sig_zp_i' : f'y3a2_i_o.{nside}_t.32768_SIGMA_MAG_ZERO.QSUM_EQU.fits.gz',
            'sig_zp_z' : f'y3a2_z_o.{nside}_t.32768_SIGMA_MAG_ZERO.QSUM_EQU.fits.gz',
            'skyvar_g': f'y3a2_g_o.{nside}_t.32768_SKYVAR.UNCERTAINTY_EQU.fits.gz',
            'skyvar_r': f'y3a2_r_o.{nside}_t.32768_SKYVAR.UNCERTAINTY_EQU.fits.gz',
            'skyvar_i': f'y3a2_i_o.{nside}_t.32768_SKYVAR.UNCERTAINTY_EQU.fits.gz',
            'skyvar_z': f'y3a2_z_o.{nside}_t.32768_SKYVAR.UNCERTAINTY_EQU.fits.gz',
            'airmass_g': f'y3a2_g_o.{nside}_t.32768_AIRMASS.WMEAN_EQU.fits.gz',
            'airmass_r': f'y3a2_r_o.{nside}_t.32768_AIRMASS.WMEAN_EQU.fits.gz',
            'airmass_i': f'y3a2_i_o.{nside}_t.32768_AIRMASS.WMEAN_EQU.fits.gz',
            'airmass_z': f'y3a2_z_o.{nside}_t.32768_AIRMASS.WMEAN_EQU.fits.gz',
            'exp_time_g' : f'y3a2_g_o.{nside}_t.32768_EXPTIME.SUM_EQU.fits.gz',
            'exp_time_r' : f'y3a2_r_o.{nside}_t.32768_EXPTIME.SUM_EQU.fits.gz',
            'exp_time_i' : f'y3a2_i_o.{nside}_t.32768_EXPTIME.SUM_EQU.fits.gz',
            'exp_time_z' : f'y3a2_z_o.{nside}_t.32768_EXPTIME.SUM_EQU.fits.gz'}
    # Construct the filename. Make sure it's correct.
    map_filename = os.path.join(mapdir,f'{nside}',maps[map_key])
    valid_file = os.path.exists(map_filename)
    print(f'map: {map_filename}  exists: {valid_file}')
    if not valid_file:
        raise OSError("map file missing.")
    # Get the requested map.
    map_array = fitsio.read(map_filename)
    hmap = np.zeros(hp.nside2npix(nside))+hp.UNSEEN
    hmap[map_array['PIXEL']] = map_array['SIGNAL']

    return hmap

def estimate_bandpowers(hmap,maskval = hp.UNSEEN,nest=True,nside=2048,apodize = True,aposcale = 0.1,nbins=32,work=None,cache_maps=False,key=None):
    bins = nmt.NmtBin.from_nside_linear(nside, nbins)
    # First, check whether map is NEST or RING. If NEST, re-order.
    if nest:
        print("Converting power spectrum from NEST to RING for pseudo-C_\ell estimation.")
        hmap =hp.pixelfunc.reorder(hmap,inp='NEST',out='RING')
 
    mask = np.ones_like(hmap)
    mask[hmap == maskval] = 0
    # Apodize mask?? 
    if apodize:
        mask = nmt.mask_apodization(mask, aposcale)
    
    # Convert to a fluctuation field.
    hmap[mask > 0] = hmap[mask > 0] / np.average(hmap[mask>0],weights=mask[mask>0]) - 1.
    if cache_maps:
        fitsio.write(f'map-{key}.fits',hmap_ngal,clobber=True)
        fitsio.write(f'map-{key}.fits',mask,clobber=False)


    field = nmt.NmtField(mask,hmap[np.newaxis,:])
    if work == None:
        w = nmt.NmtWorkspace()
        w.compute_coupling_matrix(field,field,bins)
    else:
        w=work

    cl = nmt.compute_full_master(field,field,bins,workspace=w)
    l_eff = bins.get_effective_ells()
    return l_eff,cl,w


if __name__ == '__main__':
    cache_maps = True
    table_file = 'trend-tables-2048.npz'
    tables_loader = np.load(table_file, allow_pickle=True)
    # Never save dictionaries in npz files, or you have to do the following trick. The loader helpfully mangles them into zero-d arrays.
    bin_means = tables_loader['bin_mean'].reshape(-1)[0]
    bal_ratio = tables_loader['bal_ratio'].reshape(-1)[0]
    gold_ratio = tables_loader['gld_ratio'].reshape(-1)[0]
    work=None
    for i,key in enumerate(bin_means.keys()):
        hmap = get_property_map(key)
        print(f"estimating power spectrum for {key}")
        hmap_ngal = get_synthetic_map(bin_means[key],gold_ratio[key],hmap)
        hmap_nbal = get_synthetic_map(bin_means[key],bal_ratio[key],hmap)
        l,cl_hmap,work = estimate_bandpowers(hmap,work=work)
        l,cl_ngal,work = estimate_bandpowers(hmap_ngal,work=work,cache_maps=True,key=f"gold-{key}")
        l,cl_nbal,work = estimate_bandpowers(hmap_nbal,work=work,cache_maps=True,key=f"bal-{key}")
        fig,(ax1,ax2) = plt.subplots(nrows=1,ncols=2,figsize=(14,6))
        ax1.plot(l,l*(l+1)*cl_hmap[0,:]/2/np.pi,label=f'{key}')
        ax1.set_xlabel('$\\ell$', fontsize=16)
        ax1.set_ylabel('$(\\ell(\\ell+1)C_\\ell/2\\pi$', fontsize=16)
        ax1.set_title(f'{key} power spectrum')
        ax1.set_xscale('log')
        ax1.set_yscale('log')

        ax2.plot(l,l*(l+1)*cl_ngal[0,:]/2/np.pi,label=f'Gold response')
        ax2.plot(l,l*(l+1)*cl_nbal[0,:]/2/np.pi,label=f'Balrog response')
        ax2.set_xlabel('$\\ell$', fontsize=16)
        ax2.set_ylabel('$(\\ell(\\ell+1)C_\\ell/2\\pi$', fontsize=16)
        ax2.set_title(f'{key} Balrog response power spectrum')
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        plt.legend(loc='best')
        #ax1.set_ylim(0,1e-4)
        '''
        ax2.plot(l,l*(l+1)*cl_ngal[0,:]/2/np.pi)
        ax2.set_xscale('log')
        ax2.set_xlabel('$\\ell$', fontsize=16)
        ax2.set_ylabel('$(\\ell(\\ell+1)C_\\ell/2\\pi$', fontsize=16)
        ax2.set_title(f'{key}-response galaxy power \n (Balrog)')
        ax2.set_ylim(0,1e-4)
        '''
        fig.suptitle(key)
        fig.savefig(f'{key}-power.png')
        print(f"saved results to {key}-power.png")
        plt.close(fig)
