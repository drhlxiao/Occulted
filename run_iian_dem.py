 #######################################
#display_aia_dem.py
# Erica Lastufka 15/03/2018  

#Description: Because OSX doesn't play well with XQuartz and IDL sucks
#######################################

#######################################
# Usage:

######################################

import numpy as np
import numpy.ma as ma

import os
from scipy.ndimage.filters import generic_filter as gf
import matplotlib.cm as cm
import matplotlib.colors as colors
import glob
from sunpy.net import Fido, attrs as a
import sunpy.map
from scipy.io import readsav
import astropy.units as u
from astropy.coordinates import SkyCoord
from datetime import datetime as dt
from scipy.interpolate import interp1d
from aia_dem_batch import bin_images
#how to import python version of demreg?

#global emcube
#global lgtaxis
#should store all this stuff in common block equivalents...how to deal with the map meta?

def gen_tresp_matrix(plot=False):
    '''from Iian's tutorial'''
    # Load in the SSWIDL generated response functions
    # Was produced by make_aiaresp_forpy.pro (can't escape sswidl that easily....)
    trin=readsav('/Users/wheatley/Documents/Solar/NuStar/demreg/python/aia_trespv9_en.dat')

    # Get rid of the b in the string name (byte vs utf stuff....)
    for i in np.arange(len(trin['channels'])):
        trin['channels'][i]=trin['channels'][i].decode("utf-8")
    print(trin['channels'])

    # Get the temperature response functions in the correct form for demreg
    tresp_logt=np.array(trin['logt'])
    nt=len(tresp_logt)
    nf=len(trin['tr'][:])
    trmatrix=np.zeros((nt,nf))
    for i in range(0,nf):
        trmatrix[:,i]=trin['tr'][i]
    if plot:
        # Setup some AIA colours
        clrs=['darkgreen','darkcyan','gold','sienna','indianred','darkslateblue']

        # Do the plot
        fig = plt.figure(figsize=(8, 7))
        for i in np.arange(6):
            plt.semilogy(tresp_logt,trmatrix[:,i],label=trin['channels'][i],color=clrs[i],lw=4)
        plt.xlabel('$\mathrm{\log_{10}T\;[K]}$')
        plt.ylabel('$\mathrm{AIA\;Response\;[DN\;s^{-1}\;px^{-1}\;cm^5]}$')
        plt.ylim([2e-29,5e-24])
        plt.xlim([5.2,7.6])
        plt.legend(ncol=2,prop={'size': 16})
        plt.rcParams.update({'font.size': 16})
        plt.grid(True,which='both',lw=0.5,color='gainsboro')
        plt.show()
    
    return nt,nf,trmatrix,tresp_logt

def make_simple_dfaia(submap,timerange=['20120912_09','20120912_'],folder='/Users/wheatley/Documents/Solar/NuStar/AIA',to_json=False):
    '''make the dataframe with just the unmasked data for given timerange'''
    wavs=[94,131,171,193,211,335]
    if type(submap[0]) !=SkyCoord:
        #convert to skycoord
        bl=SkyCoord()
        tr=SkyCoord() #fill in later, need coordframe, for now assume SkyCoord
    else:
        bl,tr=submap
    all_maps,pdicts=[],[]
    for w in wavs:
        afiles=glob.glob(folder+'/AIA'+timerange[0]+'*'+str(w)+'.fits')
        afiles.sort()
        amaps=[sunpy.map.Map(a).submap(bl,tr) for a in afiles]
        all_maps.append(amaps)

        pdict=track_region_box(amaps,circle=False,mask=False,plot=False) #currently might fail unless tr is 1 frame
        adf=pd.DataFrame(pdict)
        adf['filenames']=afiles
        adf['wavelength']=w
        
        pdicts.append(adf)

    dfaia=pd.concat(pdicts)
    dfaia['cutout_shape']=[np.product(np.shape(d)) for d in dfaia.data]
    #dfaia.apply() timestamps)
    dfaia.reset_index(inplace=True)
    if to_json != False:
        dfaia.to_json(to_json,default_handler=str)
    return dfaia

def get_datavec(dfaia,timestamp,waves=[94,131,171,193,211,335],plus=False,minus=False,total=True,data=False,filenames=False):
    datavec=[]
    for w in waves:
        wdf=dfaia.where(dfaia.wavelength==w).dropna(how='all')
        datavec.append(get_dataval(wdf,timestamp,plus=plus,minus=minus,total=total, data=data,filenames=filenames))
    return datavec

def get_dataval(wdf,timestamp,plus=False,minus=False,total=True,data=False,filenames=False):
    close_times= np.argmin(np.abs(wdf.timestamps - timestamp))
    #return close_times
    #print(close_times)#,close_times['flux_total_DNs-1'])
    if plus:
        val=wdf['flux_plus_DNs-1'][close_times]
    elif minus:
        val=wdf['flux_minus_DNs-1'][close_times]
    elif total:
        val=wdf['flux_total_DNs-1'][close_times]
    elif data:
        val=wdf.data[close_times]
    if filenames:
        val={'data':wdf.data[close_times],'filenames':wdf.filenames[close_times]}
    return val
    
def generate_errors(dn_in, npix=1):
    '''return error matrix/vector of same shape as input data'''
    serr_per=10.0
    #errors in dn/px/s
    #npix=4096.**2/(nx*ny)
    edata=np.zeros(len(dn_in))#([nx,ny,nf])
    gains=np.array([18.3,17.6,17.7,18.3,18.3,17.6])
    dn2ph=gains*[94,131,171,193,211,335]/3397.0
    rdnse=1.15*np.sqrt(npix)/npix
    drknse=0.17
    qntnse=0.288819*np.sqrt(npix)/npix
    for j in np.arange(nf):
        etemp=np.sqrt(rdnse**2.+drknse**2.+qntnse**2.+(dn2ph[j]*abs(dn_in[j]))/(npix*dn2ph[j]**2))
        esys=serr_per*dn_in[j]/100.
        edata[j]=np.sqrt(etemp**2. + esys**2.)
    return edata
    
def do_2D_dem(dfaia,timestamp, trmatrix,tresp_logt,temps,maxiter=20, binning=False,flat=True,plot=False):
    #another example -- 2D test
    aia= get_datavec(dfaia,timestamp,plus=False,minus=False,total=False,data=True)#get from dataframe
    nx,ny=np.shape(aia[0])
    #test that all dimensions are the same
    for i in aia[1:]:
        nxi,nyi=np.shape(i)
        if nxi !=nx or nyi !=ny:
            print('dimension mismatch in data!', nxi,nyi, ' do not match ',nx,ny)
            return None
    nf=len(aia)
    data=np.zeros([nx,ny,nf])
    #convert from our list to an array of data
    for j in np.arange(nf):
        data[:,:,j]=aia[j]
    data[data < 0]=0
    
    if binning:
        #print(np.shape(data))
        data=bin_images(data,n=binning)
        data=data.T
        nx,ny,_=np.shape(data)
        #print(np.shape(data.T))

    #calculate our dem_norm guess
    nt=len(temps)-1
    off=0.412
    gauss_stdev=12
    dem_norm0=np.zeros([nx,ny,nt])
    dem_norm_temp=np.convolve(np.exp(-(np.arange(nt)+1-(nt-2)*(off+0.1))**2/gauss_stdev),np.ones(3)/3)[1:-1]
    dem_norm0[:,:,:]=dem_norm_temp


    serr_per=10.0
    #errors in dn/px/s
    npix=4096.**2/(nx*ny)
    edata=np.zeros([nx,ny,nf])
    gains=np.array([18.3,17.6,17.7,18.3,18.3,17.6])
    dn2ph=gains*[94,131,171,193,211,335]/3397.0
    rdnse=1.15*np.sqrt(npix)/npix
    drknse=0.17
    qntnse=0.288819*np.sqrt(npix)/npix
    for j in np.arange(nf):
        etemp=np.sqrt(rdnse**2.+drknse**2.+qntnse**2.+(dn2ph[j]*abs(data[:,:,j]))/(npix*dn2ph[j]**2))
        esys=serr_per*data[:,:,j]/100.
        edata[:,:,j]=np.sqrt(etemp**2. + esys**2.)

    dem,edem,elogt,chisq,dn_reg=dn2dem_pos(data,edata,trmatrix,tresp_logt,temps,dem_norm0=dem_norm0,max_iter=maxiter)
    
    if plot:
        quick_plot(dem,np.log10(temps[:-1]))
        
    if flat:
        #flatten by tvec and turn into dataframe
        dfdict={}
        dfdict['timestamp']=pd.Series(timestamp)
        dfdict['chisq']=pd.Series([chisq])
        dfdict['chisq_mean']=pd.Series(np.mean(chisq[chisq != 0])) #drop zeros
        dfdict['chisq_std']=pd.Series(np.std(chisq[chisq != 0]))
        dfdict['dn_reg']=pd.Series([dn_reg])
        dfdict['maxiter']=pd.Series(maxiter)
        dfdict['nx']=pd.Series(nx)
        dfdict['ny']=pd.Series(ny)
        dfdict['fraction_nonzero']=pd.Series(fraction_nonzeros(dem[:,:,0]))

        if binning:
            dfdict['binning']=pd.Series(b)
        else:
            dfdict['binning']=pd.Series(None)
        
        for i,t in enumerate(np.log10(temps[:-1])):
            demkey='dem_'+str(t)
            edemkey='edem_'+str(t)
            elogtkey='elogt_'+str(t)
            dfdict[demkey]=pd.Series([dem[:,:,i]])
            dfdict[edemkey]=pd.Series([edem[:,:,i]])
            dfdict[elogtkey]=pd.Series([elogt[:,:,i]])
            dfdict[demkey+'_mean']=pd.Series(np.mean(dem[:,:,i][dem[:,:,i] !=0]))
            dfdict[edemkey+'_mean']=pd.Series(np.mean(edem[:,:,i][edem[:,:,i] !=0]))
            dfdict[elogtkey+'_mean']=pd.Series(np.mean(elogt[:,:,i][elogt[:,:,i] !=0]))
            dfdict[demkey+'_max']=pd.Series(np.max(dem[:,:,i][dem[:,:,i] !=0]))
            dfdict[edemkey+'_max']=pd.Series(np.max(edem[:,:,i][edem[:,:,i] !=0]))
            dfdict[elogtkey+'_max']=pd.Series(np.max(elogt[:,:,i][elogt[:,:,i] !=0]))
            #print(i,np.mean(np.nonzero(dem[:,:,i])))
            #calculate percent zeros cuz why not

    
        df=pd.DataFrame(dfdict)
    
        return df
        
    else:
        return dem,edem,elogt,chisq,dn_reg
    
def plot_dem(aia,dem,mlogt):
    fig=plt.figure(figsize=(6,10))
    xylabels=['Log T=' + str(np.round(m,1)) for m in mlogt]
    meta=aia[0].meta
    norm=colors.Normalize(vmin=np.min(dem),vmax=np.max(dem))
    for i in range(nt):
        ax=fig.add_subplot(7,3,i+1)
        demmap=sunpy.map.Map(dem[:,:,i],meta)
        cf=demmap.plot(axes=ax,cmap=cm.rainbow,norm=norm)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.title.set_visible(False)
        ax.annotate(xylabels[i],xy=(demmap.bottom_left_coord.Tx.value+5,demmap.bottom_left_coord.Ty.value+5),xytext=(demmap.bottom_left_coord.Tx.value+5,demmap.bottom_left_coord.Ty.value+5),color='w')#,textcoords='axes points')

        fig.tight_layout()
        fig.subplots_adjust(top=0.95)
        fig.suptitle('AIA DEM Analysis Results '+meta['date_obs'])
        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.875, 0.15, 0.05, 0.7])
        fig.colorbar(cf, cax=cbar_ax)
        fig.show()
        
def quick_plot(dem,mlogt):
    fig=plt.figure(figsize=(6,10))
    xylabels=['Log T=' + str(np.round(m,1)) for m in mlogt]
    norm=colors.Normalize(vmin=np.min(dem),vmax=np.max(dem))
    for i in range(len(mlogt)):
        ax=fig.add_subplot(5,4,i+1)
        cf=ax.imshow(dem[:,:,i],cmap=cm.rainbow,norm=norm)#demmap.plot(axes=ax,cmap=cm.rainbow,norm=norm)
        print(i,np.mean(np.nonzero(dem[:,:,i])),np.mean(dem[:,:,i]))
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.title.set_visible(False)
        #ax.annotate(xylabels[i],xy=(demmap.bottom_left_coord.Tx.value+5,demmap.bottom_left_coord.Ty.value+5),xytext=(demmap.bottom_left_coord.Tx.value+5,demmap.bottom_left_coord.Ty.value+5),color='w')#,textcoords='axes points')

    fig.tight_layout()
    fig.subplots_adjust(top=0.95)
    fig.suptitle('AIA DEM Analysis Results')
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.875, 0.15, 0.05, 0.7])
    fig.colorbar(cf, cax=cbar_ax)
    fig.show()
        
def percent_zeros(dem):
    '''quickly calculate % of zeros'''
    return ((np.product(np.shape(dem)) - np.count_nonzero(dem))/np.product(np.shape(dem)))*100.

def fraction_nonzeros(dem):
    '''quickly calculate fraction of nonzeros'''
    return np.count_nonzero(dem)/np.product(np.shape(dem))

def count_nans(dem):
    '''quickly calculate # of NaNs'''
    return np.count_nonzero(np.isnan(dem))
    

    
        

