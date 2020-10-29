import numpy as np
import healpy as hp
import fitsio
import matplotlib.pyplot as plt
import pymaster as nmt
import glob
import os
import astropy.table
import ipdb

nside=2048
map_path = f'/data/des81.a/data/severett/paper-plots/desy3-balrog-plots/systematics-maps-compare/maps/{nside}'
wsp_path = None
def get_synthetic_map(bin_means,binned_map_values,map_values):
    interp_map_values = np.interp(map_values,bin_means,binned_map_values)
    return interp_map_values

def read_map(fname, nside=nside):
   data = fitsio.read(fname)
   hmap = np.zeros(12*nside**2)
   hmap[data['PIXEL']] = data['SIGNAL']
   hmap = hp.pixelfunc.reorder(hmap, inp='NEST', out='RING')
   return hmap

def get_delta_field(mask, hmap, debug=True):
   _aux_map = np.zeros(len(hmap))
   _aux_map[mask==1.0] = hmap[mask==1.0]/np.mean(hmap[mask==1.0]) - 1
   if debug:
       hp.mollview(hmap)
       hp.mollview(_aux_map)
       plt.show()
   f = nmt.NmtField(mask, [_aux_map], n_iter=0)
   return f

table_file = '/home/s1/emhuff/Projects/balrog/desy3-balrog-plots/systematics-maps-compare/trend-tables-2048.npz'
tables_loader = np.load(table_file, allow_pickle=True)
bin_means = tables_loader['bin_mean'].reshape(-1)[0]
bal_ratio = tables_loader['bal_ratio'].reshape(-1)[0]
gold_ratio = tables_loader['gld_ratio'].reshape(-1)[0]

if wsp_path is not None:
   wsp.read_from(wsp_path)
else:
   wsp = nmt.NmtWorkspace()
map_filelist = glob.glob(os.path.join(map_path, '*.fits.gz'))
cls_out = dict()
for i, fname in enumerate(map_filelist):
    sys_name = fname.split('.32768_')[1].split('.')[0]
    band_name = fname.split('y3a2_')[1].split('_')[0]
    print(sys_name, band_name)
    hmap = read_map(fname)
    if i==0:
        mask = 1.0*((hmap !=0) & (hmap != hp.UNSEEN))
        fsky = np.sum(mask)/len(mask)
    key = f'{sys_name.lower()}_{band_name}'
    if 'exptime' in key:
        key = key.replace('exptime', 'exp_time')
        ipdb.set_trace()
    if 'sig' in key:
        key = f'sig_zp_{band_name}'
    hmap_ngal = get_synthetic_map(bin_means[key],gold_ratio[key],hmap)
    hmap_nbal = get_synthetic_map(bin_means[key],bal_ratio[key],hmap)
    f_map = get_delta_field(mask, hmap, debug=False)
    f_ngal = get_delta_field(mask, hmap_ngal, debug=False)
    f_nbal = get_delta_field(mask, hmap_nbal, debug=False)
    if i==0:
        b = nmt.NmtBin(nside, nlb=int(1/fsky)) # This creates a bunch of bins but alas
        wsp.compute_coupling_matrix(f_map, f_map, b)
    _cls_ngal = wsp.decouple_cell(nmt.compute_coupled_cell(f_ngal, f_ngal))[0]
    _cls_nbal = wsp.decouple_cell(nmt.compute_coupled_cell(f_nbal, f_nbal))[0]
    _cls_map = wsp.decouple_cell(nmt.compute_coupled_cell(f_map, f_map))[0]
    if i==0:
        ells = b.get_effective_ells()
        cls_out['ells'] = ells
    cls_out[f'{band_name}_{sys_name}_gal'] = _cls_ngal
    cls_out[f'{band_name}_{sys_name}_bal'] = _cls_nbal
    cls_out[f'{band_name}_{sys_name}_map'] = _cls_map
    f, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=False, figsize=(10, 3))
    ax[0].plot(ells, _cls_map*ells*(ells+1)/(2*np.pi), 'o')
    ax[1].plot(ells, _cls_ngal*ells*(ells+1)/(2*np.pi), 'o', label='Gold')
    ax[1].plot(ells, _cls_nbal*ells*(ells+1)/(2*np.pi), 'o', label='Balrog')
    ax[1].legend(loc='best')
    ax[1].set_title(f'{band_name} {sys_name}')
    ax[1].set_xlabel(r'$\ell$', fontsize=16)
    ax[0].set_title(f'Map power. {band_name} {sys_name}')
    ax[0].set_xlabel(r'$\ell$', fontsize=16)
    ax[0].set_ylabel(r'$\ell(\ell+1) C_{\ell}/2\pi$', fontsize=16)
    ax[0].set_xscale('log')
    ax[0].set_xlim(None, 2*nside)
    plt.tight_layout() 
    f.savefig('{band_name}_{sys_name}.pdf')
    #plt.show()
    plt.clear(f)
tab  = astropy.table.Table(cls_out)
tab.write(f'Cls_sysmap_balrog_gold_{nside}.fits.gz', overwrite=True)
