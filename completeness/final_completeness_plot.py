import numpy as np
import matplotlib.pyplot as plt




if __name__ == '__main__':

    data = np.load("./balrog_multiband_completeness.npz")

    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(7,5))
    ntiles = 1493
    njack = 50
    fA = np.sqrt(0.5/2.68) # sqrt ratio of tile area to deep field area

    # g:
    ax.errorbar(data['bin_centers'],data['band_detfracs'][0,:],data['detfrac_err']*np.sqrt(ntiles*1./njack)*fA,label='g',marker='o',color='blue',fmt='s')
    xdeep_g,detfrac_deep_g = np.loadtxt("deep_g.dat",unpack=True,skiprows=1)
    ax.plot(xdeep_g,detfrac_deep_g,color='blue')

    ax.errorbar(data['bin_centers'],data['band_detfracs'][1,:],data['detfrac_err']*np.sqrt(ntiles*1./njack)*fA,label='r',marker='o',color='green',fmt='s')
    xdeep_r,detfrac_deep_r = np.loadtxt("deep_r.dat",unpack=True,skiprows=1)
    ax.plot(xdeep_r,detfrac_deep_r,color='green')

    ax.errorbar(data['bin_centers'],data['band_detfracs'][2,:],data['detfrac_err']*np.sqrt(ntiles*1./njack)*fA,label='i',marker='o',color='orange',fmt='s')
    xdeep_i,detfrac_deep_i = np.loadtxt("deep_i.dat",unpack=True,skiprows=1)
    ax.plot(xdeep_i,detfrac_deep_i,color='orange')

    ax.errorbar(data['bin_centers'],data['band_detfracs'][3,:],data['detfrac_err']*np.sqrt(ntiles*1./njack)*fA,label='z',marker='o',color='red',fmt='s')
    xdeep_z,detfrac_deep_z = np.loadtxt("deep_z.dat",unpack=True,skiprows=1)
    ax.plot(xdeep_z,detfrac_deep_z,color='red')

    

    plt.legend(loc='best')
    ax.set_ylim(0,1.1)
    ax.set_xlim(19,26)
    ax.set_xlabel('truth magnitude')
    ax.set_ylabel('detected fraction')
    ax.axhline(1.0,color='grey',linestyle='--')
    fig.savefig('balrog_completeness.png')
    ipdb.set_trace()
