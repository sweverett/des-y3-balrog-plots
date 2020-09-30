import numpy as np
import fitsio
from astropy.table import Table, join
import os
import time

import bin_info

xlim = bin_info.xlim
dx = bin_info.dx
cols = bin_info.cols
mag_colname = bin_info.mag_colname

def make_bal_compare_cat(vb=True):
    match_file = '/data/des81.a/data/severett/paper-plots/cats/balrog_matched_catalog_sof_y3-merged_v1.2.fits'
    det_file = '/data/des81.a/data/severett/paper-plots/cats/balrog_detection_catalog_sof_y3-merged_v1.2.fits'

    print('Making Balrog compare catalog...')
    match_cols = ['bal_id', 'true_id', 'meas_cm_mag_deredden', 'meas_cm_mag', 'true_bdf_mag_deredden', 'true_gap_flux_fwhm4asec', 'meas_gapflux',
                  'meas_cm_fracdev', 'true_bdf_fracdev', 'meas_cm_T', 'true_bdf_T', 'meas_cm_s2n_r', 'meas_cm_TdByTe',
                  'meas_cm_g', 'meas_cm_flux_s2n', 'meas_ra', 'meas_dec', 'meas_tilename']
    det_cols = ['bal_id', 'meas_FLAGS_GOLD_SOF_ONLY', 'flags_foreground', 'flags_footprint', 'flags_badregions', 'match_flag_1.5_asec',
                'meas_EXTENDED_CLASS_SOF']
    
    print('Loading matched file...')
    match = Table(fitsio.read(match_file, columns=match_cols))
    print('Loading det file...')
    det = Table(fitsio.read(det_file, columns=det_cols))
    print('joining...')
    sof = join(match, det, keys='bal_id', join_type='inner')
    assert len(match) == len(sof)
    
    # Match to Ian's DF classifier
    #df_file = 'cats/ugriz-mof02-JHK-extcorr_27May20_kNN_class.fits'
    #df = Table.read(df_file)
    #df.rename_column('id', 'true_id')
    #
    #print('Joining with df classifier...')
    #sof = join(sof, df, keys='true_id', join_type='left')
    
    cuts = np.where(
        (sof['meas_FLAGS_GOLD_SOF_ONLY'] < 2) &
        (sof['flags_foreground'] == 0 ) &
        (sof['flags_footprint'] == 1) &
        (sof['flags_badregions'] < 2) &
    #     (sof['meas_EXTENDED_CLASS_SOF'] > 1) &
        (sof['match_flag_1.5_asec'] < 2)
    )
    
    print('Making cuts...')
    sof = sof[cuts]
    
    outdir = '.'
    sof_file = os.path.join(outdir, 'balrog_sof_galaxy_compare.fits')
    print('Writing catalog...')
    sof.write(sof_file, overwrite=True)

    return

bindx = dict(zip('griz', range(4)))

vb = True # High level verbose
vb_iter = False # Verbose on iteration details

save_samples = False
overwrite = True

normed = False
save_histograms = False

Nsamples = 1000

gld_file = '/data/des81.a/data/severett/db_queries/y3_gold_2_2_galaxy_compare/y3_gold_2_2_galaxy_compare_s2n_collated.fits'
#bal_file = '/data/des41.b/data/severett/Balrog/y3-merged/stacked_catalogs/1.2/sof/balrog_matched_catalog_sof_y3-merged_v1.2.fits'
bal_file = './balrog_sof_galaxy_compare.fits'

outdir = f'./samples/{Nsamples}'

if not os.path.exists(outdir):
    os.mkdir(outdir)

if vb is True:
    print('Reading Balrog Header...')
    if not os.path.exists(bal_file):
        make_bal_compare_cat()
    h = fitsio.read_header(bal_file, ext=1)
    Nbal = h['NAXIS2']
    if vb is True:
        print('Balrog has {} objects'.format(Nbal))

if vb is True:
    print('Reading GOLD file...')
gld = fitsio.read(gld_file)

gld_hist = {}

for n in range(Nsamples):
    t0 = time.time()
    if vb is True:
        print('Generating random sample {} of {}'.format(n+1, Nsamples))
    s = gld[(np.random.rand(Nbal)*len(gld)).astype(int)]
    assert len(s) == Nbal
    
    k = 0

    if vb is True:
        print('Calculating histograms...')

    if save_histograms is True:
        hist_out = {}
    
    # Mags
    for b in 'griz':
        k += 1
        
        if vb_iter is True:
            print(b, k)

        col = mag_colname

        if 'deredden' in col:
            gcol = f'SOF_CM_MAG_CORRECTED_{b}'.upper()
        else:
            gcol = f'SOF_CM_MAG_{b}'.upper()
            
        if n == 0:
            gld_hist[f'{col}_{b}'] = []
            
        bi = bindx[b]
                    
        x = s[gcol]
                
        xbins = np.arange(xlim[col][0], xlim[col][1]+dx[col], dx[col])
        x_h, *_ = np.histogram(x, bins=xbins, density=normed)
        gld_hist[f'{col}_{b}'].append(x_h)

        if save_histograms is True:
            hist_out[col] = x_h
       
    # Colors
    for b1, b2 in zip('gr', 'ri'):
        k += 1

        if vb_iter is True:
            print(f'{b1}-{b2}, {k}')

        col = f'{b1}-{b2}'
            
        if 'deredden' in mag_colname:
            gcol = 'SOF_CM_MAG_CORRECTED'
        else:
            gcol = 'SOF_CM_MAG'
            
        if n == 0:
            gld_hist[col] = []

        b1i = bindx[b1]
        b2i = bindx[b2]

        x = s[f'{gcol}_{b1}'.upper()] - s[f'{gcol}_{b2}'.upper()]

        xbins = np.arange(xlim[col][0], xlim[col][1]+dx[col], dx[col])
        s_h, *_ = np.histogram(x, bins=xbins, density=normed)
        gld_hist[col].append(s_h)

        if save_histograms is True:
            hist_out[col] = s_h
        
    # Rest of cols
    for col in cols:
        k += 1

        if vb_iter is True:
            print(col, k)
            
        if n == 0:
            gld_hist[col] = []

        x = s[col.replace('meas', 'sof').upper()]
        
        xbins = np.arange(xlim[col][0], xlim[col][1]+dx[col], dx[col])
        x_h, *_ = np.histogram(x, bins=xbins, density=normed)
        gld_hist[col].append(x_h)

        if save_histograms is True:
            hist_out[col] = gld_h

    if save_histograms is True:
        outfile = os.path.join(outdir, 'y3_gold_2_2_galaxy_compare_no_s2n_sample_{}_hists.fits'.format(n))
        np.savez(outfile, **hist_out, overwrite=overwrite)

    if save_samples is True:
        outfile = os.path.join(outdir, 'y3_gold_2_2_galaxy_compare_no_s2n_sample_{}.fits'.format(n))
        fitsio.write(outfile, s, overwrite=overwrite)

    t1 = time.time()

    if vb_iter is True:
        print(f'Iteration {n} took {t1-t0:.2f}s')

gld_mean = {}
gld_std = {}
for col, vals in gld_hist.items():
    if vb is True:
        print('Calculating bin means...')
    gld_mean[col] = np.mean(vals, axis=0)

    if vb is True:
        print('Calculating bin stds...')
    gld_std[col] = np.std(vals, axis=0)

if vb is True:
    print('Saving histogram means...')
outfile = os.path.join(outdir, f'y3_gold_2_2_galaxy_compare_no_s2n_{Nsamples}_iter_hist_means.fits')
np.savez(outfile, **gld_mean)
    
if vb is True:
    print('Saving histogram stds...')
outfile = os.path.join(outdir, f'y3_gold_2_2_galaxy_compare_no_s2n_{Nsamples}_iter_hist_stds.fits')
np.savez(outfile, **gld_std)

if vb is True:
    print('Done!')


