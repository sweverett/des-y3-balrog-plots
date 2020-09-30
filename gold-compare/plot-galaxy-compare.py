import numpy as np
import fitsio
from astropy.table import Table, join, vstack
import os, sys
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
import pandas as pd
import h5py
import corner
#import mpl_scatter_density

# Make the norm object to define the image stretch
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
norm = ImageNormalize(vmin=0., vmax=1000, stretch=LogStretch())

import sys
sys.path.insert(0, '/global/homes/s/severett/repos/balutils/')
from balutils import stacked_catalogs as sc

import seaborn as sb
plt.style.use('seaborn')
sb.set_context("notebook", font_scale=1.5)

run = 'y3-merged'
ver = '1.2'

make_new_cats = False

#match_file = '/project/projectdirs/des/severett/Balrog/{}/stacked_catalogs/{}/sof/balrog_matched_catalog_sof_{}_v{}.fits'.format(run, ver, run, ver)
#det_file = '/project/projectdirs/des/severett/Balrog/{}/stacked_catalogs/{}/sof/balrog_detection_catalog_sof_{}_v{}.fits'.format(run, ver, run, ver)
#
#match_cols = ['bal_id', 'true_id', 'meas_cm_mag_deredden', 'meas_cm_mag', 'true_bdf_mag_deredden', 'true_gap_flux_fwhm4asec', 'meas_gapflux',
#              'meas_cm_fracdev', 'true_bdf_fracdev', 'meas_cm_T', 'true_bdf_T', 'meas_cm_s2n_r', 'meas_cm_TdByTe',
#              'meas_cm_g', 'meas_cm_flux_s2n']
#det_cols = ['bal_id', 'meas_FLAGS_GOLD_SOF_ONLY', 'flags_foreground', 'flags_footprint', 'flags_badregions', 'match_flag_1.5_asec',
#            'meas_EXTENDED_CLASS_SOF']
## sof = sc.BalrogMatchedCatalog(match_file, det_file, match_cols=match_cols, match_type='sof_only', vb=True)
#
#print('Loading matched file...')
#match = Table(fitsio.read(match_file, columns=match_cols))
#print('Loading det file...')
#det = Table(fitsio.read(det_file, columns=det_cols))
#print('joining...')
#sof = join(match, det, keys='bal_id', join_type='inner')
#assert len(match) == len(sof)
#
## Match to Ian's DF classifier
#df_file = 'cats/ugriz-mof02-JHK-extcorr_27May20_kNN_class.fits'
#df = Table.read(df_file)
#df.rename_column('id', 'true_id')
#
#print('Joining with df classifier...')
#sof = join(sof, df, keys='true_id', join_type='left')
#
#cuts = np.where(
#    (sof['meas_FLAGS_GOLD_SOF_ONLY'] < 2) &
#    (sof['flags_foreground'] == 0 ) &
#    (sof['flags_footprint'] == 1) &
#    (sof['flags_badregions'] < 2) &
##     (sof['meas_EXTENDED_CLASS_SOF'] > 1) &
#    (sof['match_flag_1.5_asec'] < 2)
#)
#
#print('Making cuts...')
#sof = sof[cuts]
#
#outdir = '/project/projectdirs/des/severett/Balrog/paper-plots/cats'
#sof_file = os.path.join(outdir, 'balrog_sof_galaxy_compare.fits')
#print('Writing catalog...')
#sof.write(sof_file, overwrite=True)

print('Reading Balrog...')
outdir = '/project/projectdirs/des/severett/Balrog/paper-plots/cats/'
sof_file = os.path.join(outdir, 'balrog_sof_galaxy_compare.fits')
sof = Table.read(sof_file)

def grab_gold_bin_stats(Nsamples, gold_dir=f'/data/des81.a/data/severett/paper-plots-gold-compare/samples/{Nsamples}', vb=True):
    means_file = os.path.join(gld_dir, 'y3_gold_2_2_galaxy_compare_no_s2n_{Nsamples}_iter_hist_means.fits.npz')
    stds_file = os.path.join(gld_dir, 'y3_gold_2_2_galaxy_compare_no_s2n_{Nsamples}_iter_hist_stds.fits.npz')
    
    if use_all is True:
        gld_file = os.path.join(gold_dir, 'y3_gold_2_2_galaxy_compare_no_s2n_collated.fits')
        
        with fitsio.read_header(f, ext=1)
        h = fitsio.read_header(f, ext=1)
        l = h['NAXIS2']
        rows = np.random.choice(Nbal, Nsamples, replace=False)
        f = '/project/projectdirs/des/severett/Balrog/paper-plots/cats/y3_gold_2_2_galaxy_compare_no_s2n_collated.fits'
        t = fitsio.read(f, rows=rows)
        gld = np.random.choice(gld, size=N, replace=False)

    else:
        gld_file = np.random.choice(gld_files, Nfiles, replace=True)

        cats = []
        k = 0
        for f in gld_file:
            k += 1

            if vb is True:
                print(f'Loading {k} of {Nfiles}...')

            cats.append(Table(fitsio.read(f)))

        print('Stacking...')
        gld = vstack(cats)

        print(f'Sampling {N} objects...')
        gld = np.random.choice(gld, size=N, replace=False)
    
    return gld

def calc_stats(bal, gld):
    bal_med = np.median(bal)
    bal_avg = np.mean(bal)
    gld_med = np.median(gld)
    gld_avg = np.mean(gld)
    
    return bal_med, bal_avg, gld_med, gld_avgo

# gld = gld[gld['EXTENDED_CLASS_SOF'] > 1]
# sof = sof[sof['meas_EXTENDED_CLASS_SOF'] > 1]

sb.set_style('whitegrid')

fig = plt.figure(constrained_layout=True)

Nrows, Ncols = 3, 4
w, h = 0, 2
outer = fig.add_gridspec(Nrows, Ncols, wspace=w, hspace=h)

vb = True
vb_iter = False

N_iterations = 5

histtype = 'step'

normed = False

bindx = dict(zip('griz', range(4)))

xlim = {
    'meas_cm_mag_deredden' : [16, 26],
    'meas_cm_mag' : [17, 26],
    'g-r' : [-2, 4],
    'r-i' : [-2, 4],
    'i-z' : [-2, 4],
    'meas_cm_g_1' : [-1, 1],
    'meas_cm_g_2' : [-1, 1],
    'meas_cm_T' : [-5, 100],
    'meas_cm_fracdev': [0, 1],
    'meas_cm_flux_s2n_i' : [0, 100],
    'meas_cm_TdByTe' : [0, 100]
}

dx = {
    'meas_cm_mag_deredden' : 0.25,
    'meas_cm_mag' : 0.25,
    'g-r' : .25,
    'r-i' : .25,
    'i-z' : .25,
    'meas_cm_g_1' : 0.1,
    'meas_cm_g_2' : 0.1,
    'meas_cm_T' : 5,
    'meas_cm_fracdev' : 0.05,
    'meas_cm_flux_s2n_i' : 5,
    'meas_cm_TdByTe' : 5
}

bal_c = 'tab:blue'
gld_c = 'tab:red'

cols = ['meas_cm_g_1', 'meas_cm_g_2', 'meas_cm_T', 'meas_cm_fracdev', 'meas_cm_TdByTe', 'meas_cm_flux_s2n_i']

# mag_colname = 'meas_cm_mag_deredden'
mag_colname = 'meas_cm_mag'

N_bal = len(sof)
# N_bal = int(1e5)
Nfiles = 2

gld_hist = {}

print('Starting GOLD iterations...')
if N_iterations is None:
    gld_file = '/project/projectdirs/des/severett/Balrog/paper-plots/cats/db_queries/y3_gold_2_2_galaxy_compare_s2n_000001.fits'
    gld = Table.read(gld_file)
    
    #...
    
else: 
    for n in range(N_iterations):
        print(f'Iteration {n+1} of {N_iterations}')
        
        gld = grab_gold_sample(N_bal, Nfiles=Nfiles)
    
        k = 0
        
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
                        
            gld_x = gld[gcol]
                    
            xbins = np.arange(xlim[col][0], xlim[col][1]+dx[col], dx[col])
            gld_h, *_ = np.histogram(gld_x, bins=xbins, density=normed)
            gld_hist[f'{col}_{b}'].append(gld_h)
           
        # Colors
        for b1, b2 in zip('gr', 'ri'):
            k += 1

            if vb_iter is True:
                print(f'{b1}-{b2}, {k}')

            col = f'{b1}-{b2}'
                
            if 'deredden' in col:
                gcol = 'SOF_CM_MAG_CORRECTED'
            else:
                gcol = 'SOF_CM_MAG'
                
            if n == 0:
                gld_hist[col] = []

            b1i = bindx[b1]
            b2i = bindx[b2]

            gld_x = gld[f'{gcol}_{b1}'.upper()] - gld[f'{gcol}_{b2}'.upper()]

            xbins = np.arange(xlim[col][0], xlim[col][1]+dx[col], dx[col])
            gld_h, *_ = np.histogram(gld_x, bins=xbins, density=normed)
            gld_hist[col].append(gld_h)
            
        # Rest of cols
        for col in cols:
            k += 1

            if vb_iter is True:
                print(col, k)
                
            if n == 0:
                gld_hist[col] = []

            gld_x = gld[col.replace('meas', 'sof').upper()]
            
            xbins = np.arange(xlim[col][0], xlim[col][1]+dx[col], dx[col])
            gld_h, *_ = np.histogram(gld_x, bins=xbins, density=normed)
            gld_hist[col].append(gld_h)

gld_mean = {}
gld_std = {}
for col, vals in gld_hist.items():
    gld_mean[col] = np.mean(vals, axis=0)
    gld_std[col] = np.std(vals, axis=0)
    
# Now do Balrog
print('Plotting...')

height_ratios = [3, 1]

k = 0
for b in 'griz':
    k += 1
    
    bi = bindx[b]
    if vb is True:
        print(b, k)

    col = mag_colname
                
    bal_x = sof[col][:,bi]

    inner = outer[k-1].subgridspec(2, 1, wspace=0, hspace=0, height_ratios=height_ratios)
    
    # Histogram plot
    fig.add_subplot(inner[0, 0])
    
    if k == 1:
        l1 = 'Balrog'
        l2 = 'Y3 GOLD'
    else:
        l1, l2 = None, None
    
    ax1 = plt.gca()
    xbins = np.arange(xlim[col][0], xlim[col][1]+dx[col], dx[col])
    bin_means = np.mean([xbins[:-1], xbins[1:]], axis=0)
    bal_h = sb.distplot(bal_x, bins=xbins, hist_kws={'label':l1, 'histtype':histtype, 'lw':3},
                        norm_hist=normed, kde=False, kde_kws={'linewidth':3}, color=bal_c, axlabel=False)
    plt.errorbar(bin_means, gld_mean[f'{col}_{b}'], gld_std[f'{col}_{b}'], fmt='o', c='k', label=l2)
    
    if (normed is False) and (k == 1):
        plt.ylabel('Counts')
    plt.xlim(xlim[col])
    plt.yscale('log')
    if k == 1:
        plt.legend(loc='upper left')
        
    # Residuals plot
    fig.add_subplot(inner[1, 0])

    ax2 = plt.gca()
    residuals = (np.histogram(bal_x, bins=xbins)[0] - gld_mean[f'{col}_{b}'])# / gld_mean[f'{col}_{b}']
    err = gld_std[f'{col}_{b}']
    plt.errorbar(bin_means, residuals, err, fmt='ok')
    plt.axhline(0, c='k', lw=2, ls='--')
    plt.xlabel(f'Meas {b}-mag (cm)')
    if k == 1:
        plt.ylabel('Residuals')
    
    ax1.get_shared_x_axes().join(ax1, ax2)
    ax1.set_xticklabels([])
    ax2.autoscale()
            
for b1, b2 in zip('gr', 'ri'):
    k += 1
    
    if vb is True:
        print(f'{b1}-{b2}, {k}')
        
    col = mag_colname

    b1i = bindx[b1]
    b2i = bindx[b2]
    
    bal_x = sof[col][:,b1i] - sof[col][:,b2i]
    
    inner = outer[k-1].subgridspec(2, 1, wspace=0, hspace=0, height_ratios=height_ratios)
    
    # Histogram plot
    fig.add_subplot(inner[0, 0])
    
    ax1 = plt.gca()
    xbins = np.arange(xlim[f'{b1}-{b2}'][0], xlim[f'{b1}-{b2}'][1]+dx[f'{b1}-{b2}'], dx[f'{b1}-{b2}'])
    bin_means = np.mean([xbins[:-1], xbins[1:]], axis=0)
    bal_j = sb.distplot(bal_x, bins=xbins, norm_hist=normed, hist_kws={'histtype':histtype, 'lw':3}, kde=False, kde_kws={'linewidth':2}, color=bal_c)
    plt.errorbar(np.mean([xbins[:-1], xbins[1:]], axis=0), gld_mean[f'{b1}-{b2}'], gld_std[f'{b1}-{b2}'], fmt='o', c='k')

    if (normed is False) and (k == 5):
        plt.ylabel('Counts')
    plt.xlim(xlim[f'{b1}-{b2}'])
    plt.yscale('log')

    # Residuals plot
    fig.add_subplot(inner[1, 0])
        
    ax2 = plt.gca()
    residuals = np.histogram(bal_x, bins=xbins)[0] - gld_mean[f'{b1}-{b2}']
    err = gld_std[f'{b1}-{b2}']
    plt.errorbar(bin_means, residuals, err, fmt='ok')
    plt.axhline(0, c='k', lw=2, ls='--')
    plt.xlabel(f'{b1}-{b2}')
    if k == 5:
        plt.ylabel('Residuals')
    
    ax1.get_shared_x_axes().join(ax1, ax2)
    ax1.set_xticklabels([])
    ax2.autoscale()
    
# For ease
sof['meas_cm_g_1'] = sof['meas_cm_g'][:,0]
sof['meas_cm_g_2'] = sof['meas_cm_g'][:,1]
sof['meas_cm_flux_s2n_i'] = sof['meas_cm_flux_s2n'][:,2]
for col in cols:
    k += 1
    
    if vb is True:
        print(col, k)
        
    bal_x = sof[col]
    
    inner = outer[k-1].subgridspec(2, 1, wspace=0, hspace=0, height_ratios=height_ratios)
    
    # Histogram plot
    fig.add_subplot(inner[0, 0])
    
    if 'fracdev' in col:
        kws = {'bw' : 0.01}
    else:
        kws = None
    
    ax1 = plt.gca()
    xbins = np.arange(xlim[col][0], xlim[col][1]+dx[col], dx[col])
    bin_means = np.mean([xbins[:-1], xbins[1:]], axis=0)
    bal_h = sb.distplot(bal_x, bins=xbins, hist_kws={'histtype':histtype, 'lw':3}, kde_kws=kws, norm_hist=normed, kde=False, color=bal_c)#linewidth=2)
    plt.errorbar(np.mean([xbins[:-1], xbins[1:]], axis=0), gld_mean[col], gld_std[col], fmt='o', c='k')

    if (normed is False) and (k == 9):
        plt.ylabel('Counts')
    plt.xlim(xlim[col])
    plt.yscale('log')
        
    # Residuals plot
    fig.add_subplot(inner[1, 0])
    
    ax2 = plt.gca()
    residuals = np.histogram(bal_x, bins=xbins)[0] - gld_mean[col]
    err = gld_std[col]
    plt.errorbar(bin_means, residuals, err, fmt='ok')
    plt.axhline(0, c='k', lw=2, ls='--')
    plt.xlabel(f'{col}')
    if k == 9:
        plt.ylabel('Residuals')
    
    ax1.get_shared_x_axes().join(ax1, ax2)
    ax1.set_xticklabels([])
    ax2.autoscale()

# plt.tight_layout()
plt.rcParams.update({'font.size': 24})
fig.set_size_inches(26, 20)

