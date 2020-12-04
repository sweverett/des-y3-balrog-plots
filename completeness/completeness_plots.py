import numpy as np
import fitsio
import matplotlib.pyplot as plt
import os
import ipdb

stacked_catalog_path = "/data/des41.b/data/severett/Balrog/y3-merged/stacked_catalogs/1.2/sof/"
truecat_file = "balrog_detection_catalog_sof_y3-merged_v1.2.fits"
obscat_file = "balrog_matched_catalog_sof_y3-merged_v1.2.fits"

def completeness(catalog,comparison_type = 'sof_magnitude', band = 'i',nbins=25):
    bands = ['g','r','i','z']
    bin_edges = np.linspace(19,26,nbins+1)
    bin_centers = (bin_edges[1:] + bin_edges[:-1])/2.
    detfrac = np.zeros(nbins)
    good_detected = ( catalog['meas_FLAGS_GOLD_SOF_ONLY'] < 2 ) # & (catalog['meas_EXTENDED_CLASS_SOF'] > 1)
    if comparison_type is 'sof_magnitude':
        for i in range(nbins):
            all_these = (catalog['true_bdf_mag_deredden'][:,bands.index(band)] >  bin_edges[i]) &  \
                        (catalog['true_bdf_mag_deredden'][:,bands.index(band)] <= bin_edges[i+1])
            if np.sum(all_these) > 0:
                
                detected = np.sum(catalog[all_these & good_detected]['detected']) 
                detfrac[i] = np.sum(detected) * 1./np.sum(all_these)
        return bin_centers, detfrac

def select_catalog(pop = 'all', exclude_match = True):
    truth = fitsio.read(os.path.join(stacked_catalog_path,truecat_file))#,rows=np.arange(50000))
    tkeep = (truth['flags_footprint'] > 0) & (truth['flags_foreground'] == 0) & (truth['flags_badregions'] == 0) & (truth['run_name'] == 'run2')
    if exclude_match:
        tkeep = tkeep & ( truth['match_flag_2.0_asec'] < 2 )
    return truth[tkeep]


if __name__ == "__main__":

    cat = select_catalog(pop='galaxies')
    tilenames = np.unique(cat['meas_tilename'])
    np.random.shuffle(tilenames)
    njack = 50
    tile_lists = np.array_split(tilenames,njack)
    nbins = 25
    detfrac_array = np.zeros((njack,nbins))
    ipdb.set_trace()
    band_detfracs = np.zeros((4,nbins))
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(7,5))
    for i,iband in enumerate(['g','r','i','z']):
        bins, this_detfrac = completeness(cat, band = iband,nbins=25)
        band_detfracs[i,:] = this_detfrac
        ax.plot(bins,this_detfrac,label=iband)
        ax.set_ylim(0,1.1)
        ax.set_xlim(19,26)
        ax.set_ylabel('detected fraction')
        ax.axhline(1.0,color='grey',linestyle='--')
    plt.legend(loc='best')
    fig.savefig('balrog_completeness-multiband.png')
    plt.close(fig)


    for i,ilist in enumerate(tile_lists):
        these = np.in1d(cat['meas_tilename'], ilist)
        bins, this_detfrac = completeness(cat[these])
        detfrac_array[i,:] = this_detfrac

    detfrac_mean = np.average(detfrac_array,axis=0)
    detfrac_err = np.std(detfrac_array,axis=0)
    np.savez('balrog_multiband_completeness.npz',band_detfracs = band_detfracs,bin_centers=bins,detfrac_err = detfrac_err)
    xhsc, detfrac_hsc =np.loadtxt("HSC.dat",unpack=True,skiprows=1)
    xdeep,detfrac_deep = np.loadtxt("deep.dat",unpack=True,skiprows=1)

    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(7,5))
    ax.errorbar(bins,detfrac_mean,detfrac_err*np.sqrt(len(tilenames)*1./njack),label='Balrog')
    #ax.plot(xhsc,detfrac_hsc,label='HSC')
    ax.plot(xdeep,detfrac_deep,label='deep')
    plt.legend(loc='best')
    ax.set_ylim(0,1.1)
    ax.set_xlim(19,26)
    ax.set_xlabel('i-band truth magnitude')
    ax.set_ylabel('detected fraction')
    ax.axhline(1.0,color='grey',linestyle='--')
    fig.savefig('balrog_completeness-i.png')
    ipdb.set_trace()

    pass
