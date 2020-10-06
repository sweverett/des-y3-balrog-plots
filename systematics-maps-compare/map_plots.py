import numpy as np
import healpy as hp
import os
import time
import matplotlib.pyplot as plt

import bin_info

import seaborn as sb
plt.style.use('seaborn')
sb.set_context("notebook", font_scale=1.5)

def plot_map_densities(mapfiles, bal, gld, 
                       w=0.3, h=0.5, s=[16, 12], show=True, vb=False,
                       nside=None, remove_stars=False,
                       outdir='plots', outfile='systematics-density-compare.png'):
    sb.set_style('whitegrid')

    xlim = bin_info.density_xlim
    dx = bin_info.density_dx
    
    Nmaps = len(mapfiles)
    Nrows, Ncols = int(np.ceil(Nmaps / 4)), 4
    
    fig, axes = plt.subplots(nrows=Nrows, ncols=Ncols)
    
    histtype = 'step'
    alpha = 0.75
    lw = 2
    
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

        # Multiply some quantities by a scale factor
        if 'sig_zp' in mname:
            sf = 1000.
        else:
            sf = 1.
            
        if (xlim[mname] is None):
            bins = 30
        else:
            bins = np.arange(sf*xlim[mname][0], sf*xlim[mname][1]+sf*dx[mname], sf*dx[mname])

        plt.hist(sf*bal[mname], bins=bins, label=lbal, density=True, histtype=histtype, alpha=alpha, lw=lw)
        plt.hist(sf*gld[mname], bins=bins, label=lgld, density=True, histtype=histtype, alpha=alpha, lw=lw)
        
        if 'sig_zp' in mname:
            xl = f'{mname} (mmag)'
        else:
            xl = mname
        plt.xlabel(xl)
       
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

    title = ''
    if nside is not None:
        title = str(nside)
        outfile = outfile.replace('.png', f'_{nside}.png')
    if remove_stars is True:
        title += ', EXTENDED_CLASS_SOF > 1'
        outfile = outfile.replace('.png', '_no_stars.png')

    if title != '':
        plt.suptitle(title)

    plt.savefig(os.path.join(outdir, outfile))

    if show is True:
        plt.show()

    return

def plot_map_trends(mapfiles, bal, gld,
                    w=0.3, h=0.5, s=[18, 12], show=True, vb=False,
                    nside=None, remove_stars=False,
                    outdir='plots', outfile='systematics-trend-compare.png',
                    error_type='poisson', Nsamples=None):

    if remove_stars is False:
        stars_dir = ''
    else:
        stars_dir = 'no_stars'

    if error_type == 'bootstrap':
        if Nsamples is None:
            raise ValueError('If using bootstrap errors, must pass Nsamples!')

        bal_boot_file = os.path.join('./samples/', str(Nsamples), stars_dir,
                                     f'y3_gold_2_2_systematics_compare_{Nsamples}_bal_std.npz')
        gld_boot_file = os.path.join('./samples/', str(Nsamples), stars_dir,
                                     f'y3_gold_2_2_systematics_compare_{Nsamples}_gld_std.npz')

        bal_boot_std = np.load(bal_boot_file)
        gld_boot_std = np.load(gld_boot_file)

    sb.set_style('whitegrid')

    xlim = bin_info.trend_xlim
    dx = bin_info.trend_dx

    Nmaps = len(mapfiles)
    Nrows, Ncols = int(np.ceil(Nmaps / 4)), 4
    
    fig, axes = plt.subplots(nrows=Nrows, ncols=Ncols)
    
    histtype = 'step'
    alpha = 0.75
    lw = 2
    
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

        if k == 1:
            # Is the same for every map, given a sample
            Nmean_bal = len(bal) / len(np.unique(bal[f'{mname}_index']))
            Nmean_gld = len(gld) / len(np.unique(gld[f'{mname}_index']))
             
        if (xlim[mname] is None):
            bins = 30
        else:
            bins = np.arange(xlim[mname][0], xlim[mname][1]+dx[mname], dx[mname])
            
        N = len(bins)-1
        bal_ratio[mname] = np.zeros(N)
        gld_ratio[mname] = np.zeros(N)
        bin_mean[mname] = np.zeros(N)        

        if error_type == 'poisson':
            bal_err[mname] = np.zeros(N)
            gld_err[mname] = np.zeros(N)
            
        j = 0
        for b1, b2 in zip(bins[:-1], bins[1:]):
            print(f'Bin {j+1} of {N}')
            bal_in_bin = np.where( (bal[mname] >= b1) & (bal[mname] < b2))
            gld_in_bin = np.where( (gld[mname] >= b1) & (gld[mname] < b2))
            
            Npixels_bal = len(np.unique(bal[f'{mname}_index'][bal_in_bin]))
            Npixels_gld = len(np.unique(gld[f'{mname}_index'][gld_in_bin]))

            try:
                assert(Npixels_bal == Npixels_gld)
            except AssertionError as e:
                print(f'Npixels_bal = {Npixels_bal}')
                print(f'Npixels_gld = {Npixels_gld}')
                raise e

            bal_base_mean = (Npixels_bal * Nmean_bal)
            gld_base_mean = (Npixels_gld * Nmean_gld)
            
            bal_ratio[mname][j] = len(bal[bal_in_bin]) / bal_base_mean
            gld_ratio[mname][j] = len(gld[gld_in_bin]) / gld_base_mean
            bin_mean[mname][j] = np.mean([b1, b2])

            if error_type == 'poisson':
                bal_err[mname][j] = np.sqrt(len(bal[bal_in_bin])) / bal_base_mean
                gld_err[mname][j] = np.sqrt(len(gld[gld_in_bin])) / gld_base_mean
            
            j += 1

        if error_type == 'bootstrap':
            bal_err[mname] = bal_boot_std[mname]
            gld_err[mname] = gld_boot_std[mname]

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
        
        plt.errorbar(bin_mean[mname], bal_ratio[mname], bal_err[mname], label=lbal, lw=2, marker='o', zorder=10)
        plt.errorbar(bin_mean[mname], gld_ratio[mname], gld_err[mname], label=lgld, lw=2, marker='o', color='k', zorder=5)
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

    title = ''
    if nside is not None:
        title = str(nside)
        outfile = outfile.replace('.png', f'{nside}.png')
    if remove_stars is True:
        title += ', EXTENDED_CLASS_SOF > 1'
        outfile = outfile.replace('.png', '_no_stars.png')

    if title != '':
        plt.suptitle(title)

    plt.savefig(os.path.join(outdir, outfile))

    if show is True:
        plt.show()

    return

def calc_bootstrap_samples(bal, gld, mapfiles, Nsamples,
                           outdir='./samples/', remove_stars=None,
                           vb=True, vb_iter=False):

    if remove_stars is None:
        raise ValueError('Must pass `remove_stars`!')

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    outdir = os.path.join(outdir, str(Nsamples))

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if remove_stars is True:
        outdir = os.path.join(outdir, 'no_stars')

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    xlim = bin_info.trend_xlim
    dx = bin_info.trend_dx

    bal_ratios = {}
    gld_ratios = {}
    bin_means = {}

    Nmaps = len(mapfiles)

    # Indices are identical across maps
    bal_pix_indices = bal[list(mapfiles.keys())[0]+'_index']
    gld_pix_indices = gld[list(mapfiles.keys())[0]+'_index']

    # Bal & gold pixels should also also be identical at this stage
    try:
        unique_bal_pixels = np.unique(bal_pix_indices)
        unique_gld_pixels = np.unique(gld_pix_indices)
        assert (unique_bal_pixels == unique_gld_pixels).all()
        unique_pixels = unique_bal_pixels
        Npixels = len(unique_pixels)
    except AssertionError as e:
        print(f'Balrog has {len(bal_pix_indices)} unique pixels')
        print(f'GOLD has {len(gld_pix_indices)} unique pixels')
        raise e

    # Only used for sampling over galaxies
    #Nbal = len(bal)
    #Ngld = len(gld)

    for n in range(Nsamples):
        t0 = time.time()
        if vb is True:
            print('Generating random sample {} of {}'.format(n+1, Nsamples))

        sampled_pixels = unique_pixels[(np.random.rand(Npixels)*Npixels).astype(int)]
        
        if vb_iter is True:
            print('Resampling Balrog...')
        bal_in_sample = np.isin(bal_pix_indices, sampled_pixels)
        s_bal = bal[bal_in_sample]
        
        if vb_iter is True:
            print('Resampling GOLD...')
        gld_in_sample = np.isin(gld_pix_indices, sampled_pixels)
        s_gld = gld[gld_in_sample]

        # Use this to sample over galaxies, but this does not include cosmic variance
        #if vb_iter is True:
        #    print('Resampling Balrog...')
        #s_bal = bal[(np.random.rand(Nbal)*Nbal).astype(int)]
        #if vb_iter is True:
        #    print('Resampling GOLD...')
        #s_gld = gld[(np.random.rand(Ngld)*Ngld).astype(int)]

        k = 0
        for mname in mapfiles.keys():
            k += 1
            if vb is True:
                print(f'Computing {mname} ratio ({k} of {Nmaps})')

            if n == 0:
                bal_ratios[mname] = []
                gld_ratios[mname] = []
                bin_means[mname] = []

            if k == 1:
                # Is the same for every map, given a sample
                Nmean_bal = len(s_bal) / len(np.unique(s_bal[f'{mname}_index']))
                Nmean_gld = len(s_gld) / len(np.unique(s_gld[f'{mname}_index']))
                 
            if (xlim[mname] is None):
                bins = 10
            else:
                bins = np.arange(xlim[mname][0], xlim[mname][1]+dx[mname], dx[mname])
                
            Nbins = len(bins)-1
            bal_ratio = np.zeros(Nbins)
            gld_ratio = np.zeros(Nbins)
            bin_mean = np.zeros(Nbins)
                
            j = 0
            for b1, b2 in zip(bins[:-1], bins[1:]):
                if vb_iter is True:
                    print(f'Bin {j+1} of {Nbins}')
                bal_in_bin = np.where( (s_bal[mname] >= b1) & (s_bal[mname] < b2))
                gld_in_bin = np.where( (s_gld[mname] >= b1) & (s_gld[mname] < b2))
        
                Npixels_bal = len(np.unique(s_bal[f'{mname}_index'][bal_in_bin]))
                Npixels_gld = len(np.unique(s_gld[f'{mname}_index'][gld_in_bin]))

                assert(Npixels_bal == Npixels_gld)

                bal_base_mean = (Npixels_bal * Nmean_bal)
                gld_base_mean = (Npixels_gld * Nmean_gld)

                bal_ratio[j] = len(bal[bal_in_bin]) / bal_base_mean
                gld_ratio[j] = len(gld[gld_in_bin]) / gld_base_mean
                bin_mean[j] = np.mean([b1, b2])
                
                j += 1
                
            bal_ratios[mname].append(bal_ratio)
            gld_ratios[mname].append(gld_ratio)
            bin_means[mname].append(bin_mean)

        #if save_samples is True:
        #    outfile = os.path.join(outdir, 'y3_gold_2_2_systematics_compare_{}.fits'.format(n))
        #    fitsio.write(outfile, s, overwrite=overwrite)

        t1 = time.time()

        print(f'Iteration {n} took {t1-t0:.2f}s')

    bal_std = {}
    gld_std = {}

    for col, vals in bal_ratios.items():
        if vb is True:
            print('Calculating Balrog bin stds...')
        bal_std[col] = np.std(vals, axis=0)

    for col, vals in gld_ratios.items():
        if vb is True:
            print('Calculating GOLD bin stds...')
        gld_std[col] = np.std(vals, axis=0)
    
    if vb is True:
        print('Saving Balrog bin stds...')
    outfile = os.path.join(outdir, f'y3_gold_2_2_systematics_compare_{Nsamples}_bal_std.npz')
    np.savez(outfile, **bal_std)

    if vb is True:
        print('Saving GOLD bin stds...')
    outfile = os.path.join(outdir, f'y3_gold_2_2_systematics_compare_{Nsamples}_gld_std.npz')
    np.savez(outfile, **gld_std)

    if vb is True:
        print('Done!')

    return
