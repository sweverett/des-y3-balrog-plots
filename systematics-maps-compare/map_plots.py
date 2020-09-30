import numpy as np
import healpy as hp
import os
import matplotlib.pyplot as plt

import seaborn as sb
plt.style.use('seaborn')
sb.set_context("notebook", font_scale=1.5)

def plot_map_densities(mapfiles, bal, gld, xlim=None, dx=None, 
                       w=0.1, h=0.25, s=[16, 10], show=True,
                       outdir='plots', outfile='systematics-density-compare.png'):
    sb.set_style('whitegrid')
    
    if xlim is None:
        xlim = {
            'fwhm_g' : [0.8, 1.5],
            'fwhm_r' : [0.8, 1.2],
            'fwhm_i' : [0.75, 1.2],
            'fwhm_z' : [0.75, 1.2],
            'airmass_i' : [1., 1.4],
            'skyvar_i' : [5, 15],
            'skybrite_i' : [2000, 5000],
            'det_frac_i' : [0.75, 1.],
            'sig_zp_i' : None,
            'exp_time_i' : [0, 800]
        }
        
    if dx is None:
        dx = {
            'fwhm_g' : 0.025,
            'fwhm_r' : 0.0125,
            'fwhm_i' : 0.0125,
            'fwhm_z' : 0.0125,
            'airmass_i' : 0.025,
            'skyvar_i' : 0.25,
            'skybrite_i' : 100,
            'det_frac_i' : .05,
            'sig_zp_i' : None,
            'exp_time_i' : 25
        }
    
    Nrows, Ncols = 2, 4
    
    fig, axes = plt.subplots(nrows=Nrows, ncols=Ncols)
    
    histtype = 'step'
    alpha = 0.75
    lw = 2
    
    Nmaps = len(mapfiles)
    
    axes = []
    
    k = 0
    for mname in mapfiles.keys():
        k += 1
        print(f'Plotting {mname} ({k} of {Nmaps})')
        
        if k == 1:
            lbal = 'Balrog'
            lgld = 'Y3 GOLD'
        else:
            l = None
        
        if k == 1:
            axes.append(plt.subplot(Nrows, Ncols, k))
        else:
            axes.append(plt.subplot(Nrows, Ncols, k))
            
        plt.setp(axes[k-1].get_yticklabels(), visible=False)
            
        if (xlim[mname] is None):
            bins = 30
        else:
            bins = np.arange(xlim[mname][0], xlim[mname][1]+dx[mname], dx[mname])
        plt.hist(bal[mname], bins=bins, label=lbal, density=True, histtype=histtype, alpha=alpha, lw=lw)
        plt.hist(gld[mname], bins=bins, label=lgld, density=True, histtype=histtype, alpha=alpha, lw=lw)
        
        plt.xlabel(mname)
       
        if (k-1) % 4 == 0:
            plt.ylabel('Density')
        
        if k == 1:
            plt.legend()
            
        axes[k-1].patch.set_edgecolor('black')  
        axes[k-1].patch.set_linewidth('2')
            
            
    plt.subplots_adjust(wspace=w, hspace=h)
    plt.gcf().set_size_inches(s[0], s[1])

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    plt.savefig(os.path.join(outdir, outfile))

    if show is True:
        plt.show()

    return

def plot_map_trends(mapfiles, bal, gld, xlim=None, dx=None,
                    w=0.3, h=0.25, s=[18, 10], show=True,
                    outdir='plots', outfile='systematics-trend-compare.png'):
    sb.set_style('whitegrid')

    if xlim is None:
        xlim = {
            'fwhm_g' : [0.8, 1.5],
            'fwhm_r' : [0.8, 1.2],
            'fwhm_i' : [0.8, 1.2],
            'fwhm_z' : [0.8, 1.2],
            'airmass_i' : [1., 1.4],
            'skyvar_i' : [5, 15],
            'skybrite_i' : [2000, 5000],
            'det_frac_i' : [0.75, 1.],
            'sig_zp_i' : None,
            'exp_time_i' : [0, 800]
        }

    if dx is None:
        dx = {
            'fwhm_g' : 0.1,
            'fwhm_r' : 0.1,
            'fwhm_i' : 0.1,
            'fwhm_z' : 0.1,
            'airmass_i' : 0.1,
            'skyvar_i' : 2.5,
            'skybrite_i' : 500,
            'det_frac_i' : .05,
            'sig_zp_i' : None,
            'exp_time_i' : 200
        }

    Nrows, Ncols = 2, 4
    
    fig, axes = plt.subplots(nrows=Nrows, ncols=Ncols)
    
    histtype = 'step'
    alpha = 0.75
    lw = 2
    
    Nmaps = len(mapfiles)
    
    axes = []
    
    Nmean_bal = {}
    Nmean_gld = {}
    
    bal_ratio = {}
    bal_err = {}
    gld_ratio = {}
    gld_err = {}
    bin_mean = {}
    
    k = 0
    for mname in mapfiles.keys():
        k += 1
        print(f'Computing {mname} ratio ({k} of {Nmaps})')
             
        if (xlim[mname] is None):
            bins = 30
        else:
            bins = np.arange(xlim[mname][0], xlim[mname][1]+dx[mname], dx[mname])
            
        # These should be the same for all mname, but let's double check
        Nmean_bal[mname] = len(bal) / len(np.unique(bal[f'{mname}_index']))
        Nmean_gld[mname] = len(gld) / len(np.unique(gld[f'{mname}_index']))
            
        N = len(bins)-1
        bal_ratio[mname] = np.zeros(N)
        bal_err[mname] = np.zeros(N)
        gld_ratio[mname] = np.zeros(N)
        gld_err[mname] = np.zeros(N)
        bin_mean[mname] = np.zeros(N)        
            
        j = 0
        for b1, b2 in zip(bins[:-1], bins[1:]):
            print(f'Bin {j+1} of {N}')
            bal_in_bin = np.where( (bal[mname] >= b1) & (bal[mname] < b2))
            gld_in_bin = np.where( (gld[mname] >= b1) & (gld[mname] < b2))
            
    #         print('Grabbing unique tiles...')
    #         Ntiles = len(np.unique(bal['meas_tilename'][bal_in_bin]))
    #         assert Ntiles == len(np.unique(gld['TILENAME'][gld_in_bin]))
    
            Npixels_bal = len(np.unique(bal[f'{mname}_index'][bal_in_bin]))
            Npixels_gld = len(np.unique(gld[f'{mname}_index'][gld_in_bin]))
            if vb is True:
                print(f'Npixels_bal = {Npixels_bal}')
                print(f'Npixels_gld = {Npixels_gld}')
            
            bal_ratio[mname][j] = len(bal[bal_in_bin]) / (Npixels_bal * Nmean_bal[mname])
            bal_err[mname][j] = np.sqrt(len(bal[bal_in_bin])) / (Npixels_bal * Nmean_bal[mname])
            gld_ratio[mname][j] = len(gld[gld_in_bin]) / (Npixels_gld * Nmean_gld[mname])
            gld_err[mname][j] = np.sqrt(len(gld[gld_in_bin])) / (Npixels_gld * Nmean_gld[mname])
            bin_mean[mname][j] = np.mean([b1, b2])
            
            j += 1
            
    k = 0
    for mname in mapfiles.keys():
        k += 1
        print(f'Plotting {mname} ({k} of {Nmaps})')
    
        if k == 1:
            lbal = 'Balrog'
            lgld = 'Y3 GOLD'
        else:
            l = None
        
        if k == 1:
            axes.append(plt.subplot(Nrows, Ncols, k))
        else:
            axes.append(plt.subplot(Nrows, Ncols, k))
        
        plt.errorbar(bin_mean[mname], bal_ratio[mname], bal_err[mname], label=lbal, lw=2, marker='o')
        plt.errorbar(bin_mean[mname], gld_ratio[mname], gld_err[mname], label=lgld, lw=2, marker='o')
        plt.axhline(1, lw=3, ls='--', c='k')
        
        plt.xlabel(mname)
       
        if (k-1) % 4 == 0:
            plt.ylabel('N / <N>')
        
        if k == 1:
            plt.legend()
            
        axes[k-1].patch.set_edgecolor('black')  
        axes[k-1].patch.set_linewidth('2')
        
    plt.subplots_adjust(wspace=w, hspace=h)
    plt.gcf().set_size_inches(s[0], s[1])

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    plt.savefig(os.path.join(outdir, outfile))

    if show is True:
        plt.show()

    return
