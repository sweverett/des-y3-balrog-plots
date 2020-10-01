import numpy as np
import fitsio
import os
import ntpath
from astropy.table import Table
import healpy as hp
import map_plots

import pdb

make_plots = ['density', 'trends']
#make_plots = ['trends']

show_plots = False

bal_file = '/data/des81.a/data/severett/paper-plots/cats/gold-compare/balrog_sof_galaxy_compare.fits'
gld_cache_file = '/data/des81.a/data/severett/paper-plots/cats/systematics-maps-compare/y3_gold_2_2_galaxy_compare_healpy.fits'

vb = True
overwrite = True

remove_stars = True # Cut on EXTENDED_CLASS_SOF > 1

NSIDE_MAP = 4096
NSIDE_OUT = 2048
nest = True
partial = True

maps = {
        'fwhm_g': 'y3a2_g_o.4096_t.32768_FWHM.WMEAN_EQU.fits.gz',
        'fwhm_r': 'y3a2_r_o.4096_t.32768_FWHM.WMEAN_EQU.fits.gz',
        'fwhm_i': 'y3a2_i_o.4096_t.32768_FWHM.WMEAN_EQU.fits.gz',
        'fwhm_z': 'y3a2_z_o.4096_t.32768_FWHM.WMEAN_EQU.fits.gz',
        'airmass_g': 'y3a2_g_o.4096_t.32768_AIRMASS.WMEAN_EQU.fits.gz',
        'airmass_r': 'y3a2_r_o.4096_t.32768_AIRMASS.WMEAN_EQU.fits.gz',
        'airmass_i': 'y3a2_i_o.4096_t.32768_AIRMASS.WMEAN_EQU.fits.gz',
        'airmass_z': 'y3a2_z_o.4096_t.32768_AIRMASS.WMEAN_EQU.fits.gz',
        'skyvar_i': 'y3a2_i_o.4096_t.32768_SKYVAR.UNCERTAINTY_EQU.fits.gz',
        'sig_zp_i' : 'y3a2_i_o.4096_t.32768_SIGMA_MAG_ZERO.QSUM_EQU.fits.gz',
        'skybrite_i' : 'y3a2_i_o.4096_t.32768_SKYBRITE.WMEAN_EQU.fits.gz',
        'exp_time_i' : 'y3a2_i_o.4096_t.32768_EXPTIME.SUM_EQU.fits.gz'
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

# Read in same catalogs as galaxy-compare:
print('Reading Balrog...')
bal_cols = ['meas_ra', 'meas_dec', 'meas_tilename', 'meas_EXTENDED_CLASS_SOF']
bal = Table(fitsio.read(bal_file, columns=bal_cols))

print('Reading GOLD...')
try:
    gld_cols = ['RA', 'DEC', 'TILENAME', 'EXTENDED_CLASS_SOF']
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

# Optional: Remove stars
if remove_stars is True:
    bal = bal[bal['meas_EXTENDED_CLASS_SOF'] > 1]
    gld = gld[gld['EXTENDED_CLASS_SOF'] > 1]

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

Nbal_after = len(bal)
Ngld_after = len(gld)

print(f'Balrog objects before & after selecting common pixels: {Nbal_before} -> {Nbal_after}')
print(f'GOLD objects before & after selecting common pixels: {Ngld_before} -> {Ngld_after}')

# -------------------------------------------------------------
# Plots

plot_outdir = os.path.join('plots', str(NSIDE_OUT))

# Density plots

density_xlim = {
    'fwhm_g' : [0.8, 1.5],
    'fwhm_r' : [0.8, 1.2],
    'fwhm_i' : [0.75, 1.2],
    'fwhm_z' : [0.75, 1.2],
    'airmass_g' : [1., 1.4],
    'airmass_r' : [1., 1.4],
    'airmass_i' : [1., 1.4],
    'airmass_z' : [1., 1.4],
    'skyvar_i' : [5, 15],
    'skybrite_i' : [2000, 5000],
    'det_frac_i' : [0.75, 1.],
    'sig_zp_i' : [.005, .015],
    'exp_time_i' : [0, 800]
}

density_dx = {
    'fwhm_g' : 0.025,
    'fwhm_r' : 0.0125,
    'fwhm_i' : 0.0125,
    'fwhm_z' : 0.0125,
    'airmass_g' : 0.025,
    'airmass_r' : 0.025,
    'airmass_i' : 0.025,
    'airmass_z' : 0.025,
    'skyvar_i' : 0.25,
    'skybrite_i' : 100,
    'det_frac_i' : .05,
    'sig_zp_i' : .0025,
    'exp_time_i' : 25
}

if 'density' in make_plots:
    print('Starting density plots')
    map_plots.plot_map_densities(mapfiles, bal, gld, xlim=density_xlim, dx=density_dx,
                                 outdir=plot_outdir, vb=vb, show=show_plots,
                                 nside=NSIDE_OUT, remove_stars=remove_stars)

# Trend Plots

trend_xlim = {
    'fwhm_g' : [0.8, 1.5],
    'fwhm_r' : [0.8, 1.2],
    'fwhm_i' : [0.8, 1.2],
    'fwhm_z' : [0.8, 1.2],
    'airmass_g' : [1., 1.4],
    'airmass_r' : [1., 1.4],
    'airmass_i' : [1., 1.4],
    'airmass_z' : [1., 1.4],
    'skyvar_i' : [5, 15],
    'skybrite_i' : [2000, 5000],
    'det_frac_i' : [0.75, 1.],
    'sig_zp_i' : [.005, .015],
    'exp_time_i' : [0, 800]
}

trend_dx = {
    'fwhm_g' : 0.1,
    'fwhm_r' : 0.1,
    'fwhm_i' : 0.1,
    'fwhm_z' : 0.1,
    'airmass_g' : 0.025,
    'airmass_r' : 0.025,
    'airmass_i' : 0.025,
    'airmass_z' : 0.025,
    'skyvar_i' : 2.5,
    'skybrite_i' : 500,
    'det_frac_i' : .05,
    'sig_zp_i' : .0025,
    'exp_time_i' : 200
}

if 'trends' in make_plots:
    print('Starting trend plots')
    map_plots.plot_map_trends(mapfiles, bal, gld, xlim=trend_xlim, dx=trend_dx,
                              outdir=plot_outdir, vb=vb, show=show_plots,
                              nside=NSIDE_OUT, remove_stars=remove_stars)

