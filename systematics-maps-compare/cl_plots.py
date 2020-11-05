import os
import numpy as np
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')
rc('font', size=20)
import matplotlib.pyplot as plt
import fitsio
import camb
from camb import model
from camb.sources import GaussianSourceWindow, SplinedSourceWindow
from scipy.interpolate import interp1d
import ipdb

def get_galaxy_power(l,z=0.7,b=1.):
    lmax = np.max([np.max(l),2000])
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
    pars.InitPower.set_params(As=2e-9, ns=0.965)
    pars.set_for_lmax(lmax, lens_potential_accuracy=1)    
    pars.Want_CMB = False
    pars.NonLinear = model.NonLinear_both
    pars.SourceWindows = [GaussianSourceWindow(redshift=0.7, source_type='counts', sigma=0.1,dlog10Ndm=-0.2, bias = 1.0)]
    results = camb.get_results(pars)
    cls = results.get_source_cls_dict()
    l_model = np.arange(2,lmax+1)
    cl_model = cls['W1xW1'][2:lmax+1]
    # Interpolate this onto the provided grid.
    f = interp1d(l_model,cl_model,kind='cubic')
    return f(l)

datafile = './Cls_sysmap_balrog_gold_2048.fits.gz'
outdir = './power'
data = fitsio.read(datafile)
ell = data['ells']
D = data['ells']*(data['ells']+1)/(2*np.pi)

ell_chip = 3.5*180./(np.array([2048,4096])*(0.26/3600.))


filters = ['g','r','i','z']

keys = ['FWHM','SKYVAR','SKYBRITE','AIRMASS','SIGMA_MAG_ZERO','EXPTIME']

# Get a fiducial galaxy power spectrum.

llcl_gal = get_galaxy_power(ell.astype(int))

for i,ikey in enumerate(keys):
    figname = os.path.join(outdir,f"balrog_power_comparison-{ikey}.png")
    fig,ax = plt.subplots(nrows=len(filters),ncols=2,figsize=(21,6*len(filters)))
    for j,jfilt in enumerate(filters):
        fig1name = os.path.join(outdir,f"systematic_power-{jfilt}-{ikey}.png")
        fig1,ax1 = plt.subplots(nrows=1,ncols=2,figsize=(18,6))
        gldpower = data[f"{jfilt}_{ikey}_gal"]
        balpower = data[f"{jfilt}_{ikey}_bal"]
        mappower = data[f"{jfilt}_{ikey}_map"]
        
        ax[j,0].plot(ell,D*mappower,label=f"{ikey} map".replace('_',' '),color='blue')
        ax[j,0].plot(ell,llcl_gal,label='galaxy power spectrum \n (z=0.7, b=1)',color='orange')
        ax1[0].plot(ell,D*mappower,label=f"{ikey} map".replace('_',' '),color='blue')
        ax1[0].plot(ell,llcl_gal,label='galaxy power spectrum \n (z=0.7, b=1)',color='orange')
        #ax[j,0].set_xlabel(r'$\ell$')
        #ax[j,0].set_ylabel(r'$\frac{\ell(\ell+1)}{2\pi} C_{\ell}$')
        #ax[j,0].set_xscale('log')
        #ax[j,0].set_yscale('log')
        
        
        ax[j,0].plot(ell,D*gldpower,label='gold galaxies',color='purple')
        ax[j,0].plot(ell,D*balpower,label='balrog galaxies',color='green')
        ax[j,0].set_xlabel(r'$\ell$')
        ax[j,0].set_ylabel(r'$\frac{\ell(\ell+1)}{2\pi} C_{\ell}$')

        ax1[0].plot(ell,D*gldpower,label='gold galaxies',color='purple')
        ax1[0].plot(ell,D*balpower,label='balrog galaxies',color='green')
        ax1[0].set_xlabel(r'$\ell$')
        ax1[0].set_ylabel(r'$\frac{\ell(\ell+1)}{2\pi} C_{\ell}$')


        
        ax1[0].legend(loc='best')
        if j == 0:
            ax[j,0].legend(loc='best')
        ax[j,0].set_xscale('log')
        ax[j,0].set_yscale('log')
        ax1[0].set_xscale('log')
        ax1[0].set_yscale('log')

        
        ax[j,1].plot(ell,D*np.abs(balpower-gldpower)/llcl_gal,label='error in Balrog power',color='black')
        ax[j,1].axvline(ell_chip[0],color='grey',alpha=0.25,linestyle='--')
        ax[j,1].axvline(ell_chip[1],color='grey',alpha=0.25,linestyle='--')        
        ax[j,1].axhline(.01,color='red',linestyle='--',alpha=0.5)
        #ax[j,1].legend(loc='best')
        #ax[j,1].set_ylim(-0.1,0.1)
        ax[j,1].set_xscale('log')
        ax[j,1].set_yscale('log')        
        ax[j,1].set_xlabel(r'$\ell$')
        ax[j,1].set_ylabel(r'$\frac{|C_{\ell,{\rm Balrog}} - C_{\ell,{\rm Gold}}|}{C_{\ell,{\rm galaxy}}}$')

        ax1[1].plot(ell,D*np.abs(balpower-gldpower)/llcl_gal,label='error in Balrog power',color='black')
        ax1[1].axvline(ell_chip[0],color='grey',alpha=0.25,linestyle='--')
        ax1[1].axvline(ell_chip[1],color='grey',alpha=0.25,linestyle='--')        
        ax1[1].axhline(.01,color='red',linestyle='--',alpha=0.5)
        #ax[j,1].legend(loc='best')
        #ax[j,1].set_ylim(-0.1,0.1)
        ax1[1].set_xscale('log')
        ax1[1].set_yscale('log')        
        ax1[1].set_xlabel(r'$\ell$')
        ax1[1].set_ylabel(r'$\frac{|C_{\ell,{\rm Balrog}} - C_{\ell,{\rm Gold}}|}{C_{\ell,{\rm galaxy}}}$')        
        fig1.tight_layout()
        fig1.savefig(fig1name)
        plt.close(fig1)


        
    fig.tight_layout()
    fig.suptitle(f'{ikey} power spectra comparison'.replace('_',' '))
    fig.savefig(figname)
    plt.close(fig)
