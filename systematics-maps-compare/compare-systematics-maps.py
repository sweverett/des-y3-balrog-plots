import numpy as np
import fitsio
import os
import ntpath
from astropy.table import Table
import healpy as hp
import map_plots
import bin_info

import pdb

make_plots = ['density', 'trends']
#make_plots = ['trends']

cache_galaxy_maps = True
show_plots = False

run_name = 'maglim_sample'

error_type = 'poisson'
#error_type = 'bootstrap'
#Nsamples = 50
use_cached_bootstrap = False # In the future, will check if present. For now hardcoded

if error_type == 'poisson':
    Nsamples = None

bal_file = '/data/des81.a/data/severett/paper-plots/cats/gold-compare/balrog_sof_galaxy_compare.fits'
gld_cache_file = '/data/des81.a/data/severett/paper-plots/cats/systematics-maps-compare/y3_gold_2_2_systematics_and_zp_corr_collated.fits'

# Smaller, has no mag cols
#gld_cache_file = '/data/des81.a/data/severett/paper-plots/cats/systematics-maps-compare/y3_gold_2_2_galaxy_compare_healpy.fits'

vb = True
vb_iter = True
overwrite = False

remove_stars = True # Cut on EXTENDED_CLASS_SOF > 1

# This cuts on EXTENDED_CLASS_SOF==3 & 17.5 < i_mag < 21.5
use_maglim_sample = True 

NSIDE_MAP = 4096
NSIDE_OUT = 2048
nest = True
partial = True

maps = {
        'fwhm_g': 'y3a2_g_o.4096_t.32768_FWHM.WMEAN_EQU.fits.gz',
        'fwhm_r': 'y3a2_r_o.4096_t.32768_FWHM.WMEAN_EQU.fits.gz',
        'fwhm_i': 'y3a2_i_o.4096_t.32768_FWHM.WMEAN_EQU.fits.gz',
        'fwhm_z': 'y3a2_z_o.4096_t.32768_FWHM.WMEAN_EQU.fits.gz',
        'skybrite_g' : 'y3a2_g_o.4096_t.32768_SKYBRITE.WMEAN_EQU.fits.gz',
        'skybrite_r' : 'y3a2_r_o.4096_t.32768_SKYBRITE.WMEAN_EQU.fits.gz',
        'skybrite_i' : 'y3a2_i_o.4096_t.32768_SKYBRITE.WMEAN_EQU.fits.gz',
        'skybrite_z' : 'y3a2_z_o.4096_t.32768_SKYBRITE.WMEAN_EQU.fits.gz',
        'sig_zp_g' : 'y3a2_g_o.4096_t.32768_SIGMA_MAG_ZERO.QSUM_EQU.fits.gz',
        'sig_zp_r' : 'y3a2_r_o.4096_t.32768_SIGMA_MAG_ZERO.QSUM_EQU.fits.gz',
        'sig_zp_i' : 'y3a2_i_o.4096_t.32768_SIGMA_MAG_ZERO.QSUM_EQU.fits.gz',
        'sig_zp_z' : 'y3a2_z_o.4096_t.32768_SIGMA_MAG_ZERO.QSUM_EQU.fits.gz',
        'skyvar_g': 'y3a2_g_o.4096_t.32768_SKYVAR.UNCERTAINTY_EQU.fits.gz',
        'skyvar_r': 'y3a2_r_o.4096_t.32768_SKYVAR.UNCERTAINTY_EQU.fits.gz',
        'skyvar_i': 'y3a2_i_o.4096_t.32768_SKYVAR.UNCERTAINTY_EQU.fits.gz',
        'skyvar_z': 'y3a2_z_o.4096_t.32768_SKYVAR.UNCERTAINTY_EQU.fits.gz',
        'airmass_g': 'y3a2_g_o.4096_t.32768_AIRMASS.WMEAN_EQU.fits.gz',
        'airmass_r': 'y3a2_r_o.4096_t.32768_AIRMASS.WMEAN_EQU.fits.gz',
        'airmass_i': 'y3a2_i_o.4096_t.32768_AIRMASS.WMEAN_EQU.fits.gz',
        'airmass_z': 'y3a2_z_o.4096_t.32768_AIRMASS.WMEAN_EQU.fits.gz',
        'exp_time_g' : 'y3a2_g_o.4096_t.32768_EXPTIME.SUM_EQU.fits.gz',
        'exp_time_r' : 'y3a2_r_o.4096_t.32768_EXPTIME.SUM_EQU.fits.gz',
        'exp_time_i' : 'y3a2_i_o.4096_t.32768_EXPTIME.SUM_EQU.fits.gz',
        'exp_time_z' : 'y3a2_z_o.4096_t.32768_EXPTIME.SUM_EQU.fits.gz'
#         'det_frac_i' : 'y3a2_griz_o.4096_t.32768_coverfoot_EQU.fits.gz'
#         'stellar_density' : 'psf_stellar_density_fracdet_binned_1024_nside_4096_cel.fits.gz'
}

# TODO: Generalize eventually
mapdir = '/data/des81.a/data/severett/paper-plots/desy3-balrog-plots/systematics-maps-compare/maps'

## Modified functions from E Huff from E Suchyta
def hpRaDecToHEALPixel(ra, dec, nside=4096, nest=True):
    phi = ra * np.pi / 180.0
    theta = (90.0 - dec) * np.pi / 180.0
    hpInd = hp.ang2pix(nside, theta, phi, nest=nest)
    return hpInd

def make_map(catalog, hp_index_tag = 'INDEX', nside=4096):
    # takes a catalog and healpix map properties.
    # returns a healpix map of the galaxy density.

    number_of_pixels = hp.nside2npix(nside)
    pixcounts = np.bincount(catalog.as_array()[hp_index_tag],minlength = number_of_pixels)
    hpmap = np.zeros(number_of_pixels,dtype=pixcounts.dtype)
    hpmap = hpmap + pixcounts
    return hpmap



def rescale_hp_map(m, NSIDE_MAP, NSIDE_OUT):

    # Make a full NSIDE_MAP array with NaNs in the other pixels
    full_map = np.empty(12*NSIDE_MAP**2, dtype=m.dtype)
    full_map[:] = np.NaN
    full_map[m['PIXEL']] = m['SIGNAL']

    # Scale to NSIDE_OUT
    scaled_map = np.empty(12*NSIDE_OUT**2, dtype=m.dtype)
    scaled_map[:] = np.NaN
    scaled_map['PIXEL'] = np.arange(len(scaled_map), dtype=m.dtype[0])
    scaled_map['SIGNAL'] = hp.ud_grade(full_map['SIGNAL'], nside_out=NSIDE_OUT, pess=False,
                                       order_in='NEST', order_out=nest, dtype=m.dtype[1])

    # Get rid of bad pixels
    BAD_INT = -1637499996306027037830206717952
    scaled_map = scaled_map[(np.isnan(scaled_map['SIGNAL']) == False) & (scaled_map['SIGNAL'] != BAD_INT)]

    return scaled_map

def add_maps_to_catalog(catalog, mapfile_dict, ratag='ra', dectag='dec', map_path=None, nside_map=4096,
                        nside_out=2048, vb=False):
    names = [*(mapfile_dict.keys())]
    nmaps = len(names)
    
    k = 0
    for iname in names:
        k += 1
        print (f'Adding map {iname} ({k} of {nmaps})')
        fname = mapfile_dict[iname]
        
        try:
            hmap = fitsio.read(fname)

        except OSError as e:
            if NSIDE_OUT == 4096:
                raise OSError('That systematics map is not in the "maps" directory!')

            if nside_map != nside_out:
        	# Rescale map
                if vb is True:
                    print(f'Rescaling map to {nside_out}')
                if nest is True:
                    order_in = 'NEST'
                else:
                    order_in = 'RING'

                print('fname:',fname)
                fname = fname.replace(f'/{NSIDE_OUT}/', '/')
                fname = fname.replace(str(NSIDE_OUT), str(NSIDE_MAP))
                hmap = fitsio.read(fname)
                hmap = rescale_hp_map(hmap, nside_map, nside_out)
                
                # Now cache map for future use
                print('Caching map...')
                outdir = os.path.join('maps', str(NSIDE_OUT))
                if not os.path.exists(outdir):
                    os.mkdir(outdir)
                base = ntpath.basename(fname).replace(str(NSIDE_MAP), str(NSIDE_OUT))
                outfile = os.path.join(outdir, base)
                fitsio.write(outfile, hmap)

            else:
                raise e 
                
        hmap_big = np.zeros(hp.nside2npix(nside_out)) + hp.UNSEEN
        hmap_big[hmap['PIXEL']] = hmap['SIGNAL']
        hInd = hpRaDecToHEALPixel(catalog[ratag], catalog[dectag], nside=nside_out)
        catalog[iname] = hmap_big[hInd]
        catalog[f'{iname}_index'] = hInd

    return

mapfiles = {}
for mname, mfile in maps.items():
    if NSIDE_OUT == NSIDE_MAP:
        mapfiles[mname] = os.path.join(mapdir, mfile)
    else:
        mfile = mfile.replace(str(NSIDE_MAP), str(NSIDE_OUT))
        mapfiles[mname] = os.path.join(mapdir, str(NSIDE_OUT), mfile)

# Some col defs
bal_m_col = 'meas_cm_mag_deredden'
gld_m_col = 'SOF_CM_MAG_CORRECTED'
ext_col = 'EXTENDED_CLASS_SOF'

# maglim def
bindx = dict(zip('griz', range(4)))
b = 'i'
bi = bindx[b]
mmin = 17.5
mmax = 21.5

# Read in same catalogs as galaxy-compare:
print('Reading Balrog...')
bal_cols = ['meas_ra', 'meas_dec', 'meas_tilename', f'meas_{ext_col}']

if use_maglim_sample is True:
    bal_cols.append(bal_m_col)

bal = Table(fitsio.read(bal_file, columns=bal_cols))

print('Reading GOLD...')
try:
    gld_cols = ['RA', 'DEC', 'TILENAME', ext_col]

    if use_maglim_sample is True:
        gld_cols.append(gld_m_col+f'_{b.upper()}')

    gld = Table(fitsio.read(gld_cache_file, columns=gld_cols))

except OSError:
    print('Could not find cached GOLD file, creating new one from scratch (will take some time)')
    print('Reading in base GOLD catalog...')
    # Create GOLD file w/ Balrog mask
    gld_file = '/data/des81.a/data/severett/db_queries/y3_gold_2_2_galaxy_compare/s2n/y3_gold_2_2_galaxy_compare_basic_sof_s2n_collated.fits'
    gld = Table(fitsio.read(gld_file, columns=gld_cols))
    
    # Need to cut GOLD to same tile footprint as Balrog
    tiles = np.unique(bal['meas_tilename'])
    Ntiles = len(tiles)
    print(f'Found {Ntiles} unique tiles')

    print('Checking for objects in Balrog footprint...')
    in_bal_footprint = np.isin(gld['TILENAME'], tiles)

    print('Cutting objects...')
    gld = gld[in_bal_footprint]

    print('Writing...')
    gld.write(gld_cache_file, overwrite=overwrite)

    print('GOLD catalog w/ Balrog mask written!')

# Optional: Cut to maglim-like sample
if use_maglim_sample is True:
    bal = bal[
             (bal[f'meas_{ext_col}'] == 3) &
             (bal[bal_m_col][:,bi] >= mmin) &
             (bal[bal_m_col][:,bi] <= mmax)
    ]
    gld = gld[
             (gld[ext_col] == 3) &
             (gld[gld_m_col+f'_{b.upper()}'] >= mmin) &
             (gld[gld_m_col+f'_{b.upper()}'] <= mmax)
    ]
else:
    # Optional: Remove stars
    if remove_stars is True:
        bal = bal[bal['meas_EXTENDED_CLASS_SOF'] > 1]
        gld = gld[gld['EXTENDED_CLASS_SOF'] > 1]
    else:
        # Still want to remove -9's
        bal = bal[bal['meas_EXTENDED_CLASS_SOF'] >= 0]
        gld = gld[gld['EXTENDED_CLASS_SOF'] >= 0]

# Add maps to catalogs

print('Adding systematics maps to Balrog...')
add_maps_to_catalog(bal, mapfiles, nside_map=NSIDE_MAP,
                    nside_out=NSIDE_OUT,
                    ratag='meas_ra', dectag='meas_dec', vb=vb)

print('Adding systematics maps to GOLD...')
add_maps_to_catalog(gld, mapfiles, nside_map=NSIDE_MAP,
                    nside_out=NSIDE_OUT,
                    ratag='RA', dectag='DEC', vb=vb)

# Due to differences in number density and injection gaps on 
# the edge of tiles, there are still some differences in pixels.
# Cleanest to only compare pixels with both bal & gld present

bal_pix_indices = bal[list(mapfiles.keys())[0]+'_index']
gld_pix_indices = gld[list(mapfiles.keys())[0]+'_index']

gld_in_bal_pixels = np.isin(gld_pix_indices, bal_pix_indices)

Ngld_before = len(gld)
gld = gld[gld_in_bal_pixels]

# Now recompute GOLD pixels to match to Balrog
gld_pix_indices = gld[list(mapfiles.keys())[0]+'_index']

bal_in_gld_pixels = np.isin(bal_pix_indices, gld_pix_indices)

Nbal_before = len(bal)
bal = bal[bal_in_gld_pixels]

Nbal = len(bal)
Ngld = len(gld)

print(f'Balrog objects before & after selecting common pixels: {Nbal_before} -> {Nbal}')
print(f'GOLD objects before & after selecting common pixels: {Ngld_before} -> {Ngld}')

# -------------------------------------------------------------
# Bootstrap Samples

if error_type == 'bootstrap':
    if use_cached_bootstrap is False:
        map_plots.calc_bootstrap_samples(bal, gld, mapfiles, Nsamples,
                                         remove_stars=remove_stars,
                                         vb=vb, vb_iter=vb_iter)

    # In future, can auto detect this
    #else:
    #    pass

# -------------------------------------------------------------
# Plots

# Can adjust this for atypical runs
plot_outdir_base = 'plots'
plot_outdir = os.path.join(plot_outdir_base, str(NSIDE_OUT))

if run_name is not None:
    plot_outdir = os.path.join(plot_outdir, run_name)

# Cache maps of the galaxy density.
if cache_galaxy_maps:
    map_outfile = os.path.join(plot_outdir,'galaxy_maps.fits')
    print(f"Making maps of galaxy density for Gold and balrog. \n Will write outputs to: {map_outfile}")
    bal_map = make_map(bal, hp_index_tag = 'skybrite_r_index', nside=4096)
    gld_map = make_map(gld, hp_index_tag = 'skybrite_r_index', nside=4096)
    map_arr = np.empty(bal_map.size,dtype=[('balrog_map',np.int),('gold_map',np.int)])
    map_arr['balrog_map'] = bal_map
    map_arr['gold_map'] = gld_map
    fitsio.write(map_outfile,map_arr)


# Density plots

if 'density' in make_plots:
    print('Starting density plots')
    map_plots.plot_map_densities(mapfiles, bal, gld,
                                 outdir=plot_outdir, vb=vb, show=show_plots,
                                 nside=NSIDE_OUT, remove_stars=remove_stars)

# Trend Plots

if 'trends' in make_plots:
    print('Starting trend plots')
    map_plots.plot_map_trends(mapfiles, bal, gld,
                              outdir=plot_outdir, vb=vb, show=show_plots,
                              nside=NSIDE_OUT, remove_stars=remove_stars,
                              error_type=error_type, Nsamples=Nsamples)
