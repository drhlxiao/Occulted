 #######################################
#make_timeseries_plot_flux.py
# Erica Lastufka 20/03/2018

#Description: Make flux timeseries plot for given region on a STEREO/AIA/whatever map
#######################################

#######################################
# Usage:

######################################
from __future__ import print_function, division
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import scipy.constants as sc
import os
from scipy import ndimage
from scipy.ndimage.filters import generic_filter as gf
from scipy.io import readsav
import matplotlib.cm as cm
import make_aligned_movie as mov
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from matplotlib.dates import date2num,num2date,DateFormatter
from datetime import datetime as dt
import glob
import sunpy.map
import astropy.units as u
from astropy.coordinates import SkyCoord
from datetime import datetime as dt
import pytz
import pickle


def get_submap(smap,bl,tr):
    bl = SkyCoord(bl[0] * u.arcsec, bl[1] * u.arcsec, frame=smap.coordinate_frame)
    tr = SkyCoord(tr[0] * u.arcsec, tr[1] * u.arcsec, frame=smap.coordinate_frame)
    ssmap = smap.submap(bl, tr)
    return ssmap

def check_coords(testmap,submap):
    '''check that coordinates in map match up with the ones of the desired submap'''
    tr=testmap.top_right_coord
    bl=testmap.bottom_left_coord
    print(int(bl.Tx.value),int(tr.Tx.value),int(bl.Ty.value),int(tr.Ty.value),submap)
    if int(bl.Tx.value) == submap[0][0] and int(tr.Tx.value) == submap[0][1] and int(bl.Ty.value) == submap[1][0] and int(tr.Ty.value) == submap[1][1]:
        return False
    else:
        return True

def mask_disk(submap,plot=False,fac=50.,filled=False,greater=False):
    x, y = np.meshgrid(*[np.arange(v.value) for v in submap.dimensions]) * u.pixel
    hpc_coords = submap.pixel_to_data(x, y)
    r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / (submap.rsun_obs+fac*u.arcsec)
    if not greater:
        mask = ma.masked_less_equal(r, 1)
    else:
        mask = ma.masked_greater_equal(r, 1)

    palette = submap.plot_settings['cmap']
    palette.set_bad('black')
    if filled:
        data_arr=np.ma.masked_array(submap.data,mask=mask.mask)
        filled_data=ma.filled(data_arr,0)
        scaled_map = sunpy.map.Map(filled_data, submap.meta)
    else:
        scaled_map = sunpy.map.Map(submap.data, submap.meta, mask=mask.mask)
    if plot:
        fig = plt.figure()
        fig.add_subplot(projection=scaled_map)
        scaled_map.plot(cmap=palette)
        scaled_map.draw_limb()
        fig.show()
    return scaled_map

def make_timeseries_plot_flux(ddir,mapcube=False,tags=False,n=False,submap=False,outfname=False,picklenewmaps=True,normalize=True,maskdisk=True,from_fits=False, background_subtract=True):
    ''' input for kwarg submap is [(bottom_left_tuple), (top_righ_tuple)]. should probably mpi this '''
    #check if any of the map pickles already exist, and if they do that they are the right submap
    cwd=os.getcwd()
    os.chdir(ddir)
    #reshape=True
    if tags:
        globparam='*'+tags
    else:
        globparam='*'
    all_fits=glob.glob(globparam+'.fts') #or .fits
    fitsnames=[af[:-4] for af in all_fits]
    #first search for pickles or pickle which holds cube...
    if mapcube != False:
        if type(mapcube) == list:
            all_maps=[]
            for mc in mapcube:
                all_maps.extend(pickle.load(open(mc,'rb')))
        else:
            all_maps=pickle.load(open(mapcube,'rb'))
    else:
        pickles=glob.glob(globparam+'*.p')
        if len(pickles) == 0 or from_fits == True:
            mapcube=sunpy.map.Map(all_fits,cube=True)
            if submap:
                all_maps=[mc.submap(SkyCoord(submap[0]*u.arcsec,submap[1]*u.arcsec,frame=mapcube[0].coordinate_frame)) for mc in mapcube]
            else:
                all_maps=mapcube
        else:
            all_maps=[pickle.load(open(p,'rb')) for p in pickles]
    if submap:
        reshape=check_coords(all_maps[0],submap)
    else:
        reshape=False
    if reshape: #reshape the maps.. from the fits files I guess. or full disk
        try:
            #print('test1')
            all_maps=[mc.submap(SkyCoord(submap[0]*u.arcsec,submap[1]*u.arcsec,frame=all_maps[0].coordinate_frame)) for mc in all_maps]
        except ValueError:
            #coordinates out of range
            fitsmaps=sunpy.map.Map(all_fits,cube=True)
            print(type(fitsmaps))
            all_maps=[mc.submap(SkyCoord(submap[0]*u.arcsec,submap[1]*u.arcsec,frame=fitsmaps[0].coordinate_frame)) for mc in fitsmaps]

    if picklenewmaps:
        pickle.dump(all_maps,open(picklenewmaps+'.p','wb'))

    #now sum the map data for each map, store in vector. make timestamp vector
    dtvec=[dt.strptime(m.meta['date-obs'][:-4],'%Y-%m-%dT%H:%M:%S') for m in all_maps if np.isfinite(np.sum(m.data))]
    if maskdisk:
        masked_maps=[mask_disk(m) for m in all_maps]
        sumdata0=[np.sum(m.data[np.where(m.mask == False)]) for m in masked_maps]
    else:
        sumdata0=[np.sum(m.data) for m in all_maps if np.isfinite(np.sum(m.data))]
    #subtract pre-flare image
    if background_subtract:
        try:
            bmap=sunpy.map.Map(background_subtract)
            if submap:
                bmap=bmap.submap(SkyCoord(submap[0]*u.arcsec,submap[1]*u.arcsec,frame=bmap.coordinate_frame))
                print('here',np.sum(bmap.data))
            if maskdisk:
                bmmap=mask_disk(bmap)
                sumdata=[s-np.sum(bmmap.data[np.where(bmmap.mask == False)]) for s in sumdata0]
            else:
                sumdata=[s-np.sum(bmap.data) for s in sumdata0]
        except ValueError:
            #print(sumdata[0])
            sumdata=[s-sumdata0[0] for s in sumdata0]
    #normalize
    if normalize:
        try:
            sumdata=sumdata/np.max(sumdata)
        except NameError:
            sumdata=sumdata0/np.max(sumdata0)

    dtvecnums=date2num(dtvec)
    fig,ax=plt.subplots()
    ax.plot_date(dtvecnums,sumdata,'b.-')
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    #etc
    fig.show()
    if outfname:
        plt.savefig(outfname)
    return dtvec,sumdata,sumdata0,all_maps

def make_big_lightcurve(ddir,mapcube=False,tags=False,n=False,submap=False,outfname=False,picklelc=True,normalize=True,maskdisk=True,from_fits=False, background_subtract=True):
    '''for when there are a lot of files, add intermediate step '''
    #check if any of the map pickles already exist, and if they do that they are the right submap
    cwd=os.getcwd()
    os.chdir(ddir)
    #reshape=True
    if tags:
        globparam='*'+tags
    else:
        globparam='*'
    all_fits=glob.glob(globparam+'.fts') #or .fits
    fitsnames=[af[:-4] for af in all_fits]
    #first search for pickles or pickle which holds cube...
    if mapcube != False:
        if type(mapcube) == list:
            all_maps=[]
            for mc in mapcube:
                all_maps.extend(pickle.load(open(mc,'rb')))
        else:
            all_maps=pickle.load(open(mapcube,'rb'))
    else:
        pickles=glob.glob(globparam+'*.p')
        if len(pickles) == 0 or from_fits == True:
            mapcube=sunpy.map.Map(all_fits,cube=True)
            if submap:
                all_maps=[mc.submap(SkyCoord(submap[0]*u.arcsec,submap[1]*u.arcsec,frame=mapcube[0].coordinate_frame)) for mc in mapcube]
            else:
                all_maps=mapcube
        else:
            all_maps=[pickle.load(open(p,'rb')) for p in pickles]

    #if picklelc:
    #    pickle.dump(all_maps,open(picklenewmaps+'.p','wb'))

    #now sum the map data for each map, store in vector. make timestamp vector
    #dtvec=[dt.strptime(m.meta['date-obs'][:-4],'%Y-%m-%dT%H:%M:%S') for m in all_maps if np.isfinite(np.sum(m.data))]
    sumdata,dtvec=[],[]
    if background_subtract:
        if type(background_subtract)==str:
            bmap=sunpy.map.Map(background_subtract)
        else:
            bmap=all_maps[0]
        bsumdata,foo=get_lc_from_map(bmap,submap=submap,maskdisk=maskdisk)
    else:
        bsumdata=False
    for m in all_maps:
        sd,dtv=get_lc_from_map(m,submap=submap,maskdisk=maskdisk,backgroundsubtract=bsumdata)
        sumdata.append(sd)
        dtvec.append(dtv)
    #normalize
    if normalize:
        try:
            maxv=np.max(sumdata)
            sumdata=sumdata/maxv
        except NameError:
            maxv=np.max(sumdata0)
            sumdata=sumdata0/np.max(sumdata0)
    else:
        maxv=np.max(sumdata)
    #pickle things
    if picklelc:
        pickle.dump([dtvec,sumdata,maxv],open('lightcurve.p','wb'))
    dtvecnums=date2num(dtvec)
    fig,ax=plt.subplots()
    ax.plot_date(dtvecnums,sumdata,'b.-')
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    #etc
    fig.show()
    if outfname:
        plt.savefig(outfname)
    return dtvec,sumdata,maxv

def get_background(bfilename,submap=False,maskdisk=True):
    bmap=sunpy.map.Map(bfilename)
    if submap:
        bmap=bmap.submap(SkyCoord(submap[0]*u.arcsec,submap[1]*u.arcsec,frame=bmap.coordinate_frame))
        #print('here',np.sum(bmap.data))
    if maskdisk:
        bmmap=mask_disk(bmap)
        bsumdata=np.sum(bmmap.data[np.where(bmmap.mask == False)])
    else:
        bsumdata=np.sum(bmap.data)
    return bsumdata

def nearest_time(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))

def calc_FE18_lc(lclist=['lightcurve94nobg.p','lightcurve171nobg.p','lightcurve211nobg.p'],num=200,plot=True,tnorm=False,yran=[-2,1],return_lc=False):
    from scipy.interpolate import interp1d
    dtv94,lc94,max94=pickle.load(open(lclist[0],'rb'))
    dtv171,lc171,max171=pickle.load(open(lclist[1],'rb'))
    dtv211,lc211,max211=pickle.load(open(lclist[2],'rb'))
    print(np.mean(lc94),max94)
    print(np.mean(lc211),max211)
    print(np.mean(lc171),max171)
    #dtvecnums=date2num(dtvec)
    if np.max(lc94) != max94:
        lc94=[l*max94 for l in lc94]
    if np.max(lc171) != max171:
        lc171=[l*max171 for l in lc171]
    if np.max(lc211) != max211:
        lc211=[l*max211 for l in lc211]

    #interpolate onto the same time grid
    dtnums=date2num(dtv171)
    x=np.linspace(dtnums[0],dtnums[-11],num=num)
    f94=interp1d(date2num(dtv94),lc94)
    f171=interp1d(date2num(dtv171),lc171)
    f211=interp1d(date2num(dtv211),lc211)
    fac211=f211(x)/120.
    fac171=f171(x)/450.
    #print(np.mean(lc94))
    #print(np.mean(lc211),np.mean(fac211))
    #print(np.mean(lc171),np.mean(fac171))

    fe18=f94(x)-f211(x)/120.-(f171(x)/450.)
    #print(np.mean(fe18),np.mean(fe18/np.max(fe18)))

    if tnorm: #normalize to this time
        sdate=dt.strftime(dtv94[0].date(),'%Y%m%d')
        dtnorm=dt.strptime(sdate+tnorm,'%Y%m%d%H:%M')
        adtnorm=dtnorm.replace(tzinfo=pytz.UTC)
        #print(adtnorm,adtnorm.tzinfo)
        xidx=[idx for idx,xx in enumerate(x) if num2date(xx) == nearest_time(num2date(x),adtnorm)][0]
        #print(nearest_time(num2date(x),adtnorm),xidx)
        n94idx=[idx for idx,xx in enumerate(dtv94) if xx == nearest_time(dtv94,dtnorm)][0]
        #print(nearest_time(dtv94,dtnorm),n94idx)
        n171idx=[idx for idx,xx in enumerate(dtv171) if xx == nearest_time(dtv171,dtnorm)][0]
        n211idx=[idx for idx,xx in enumerate(dtv211) if xx == nearest_time(dtv211,dtnorm)][0]
    else:
        xidx,n94idx,n171idx,n211idx=0,0,0,0
    #plot
    if plot:
        fig,ax=plt.subplots()
        ax.plot_date(x,(fe18-fe18[xidx])/np.max(fe18-fe18[xidx]),'r.-',label='Fe XVIII')
        ax.plot_date(dtv94,(lc94-lc94[n94idx])/np.max(lc94-lc94[n94idx]),'b-',label='94 A')
        ax.plot_date(dtv171,(lc171-lc171[n171idx])/np.max(lc171-lc171[n171idx]),'k-',label='171 A')
        ax.plot_date(dtv211,(lc211-lc211[n211idx])/np.max(lc211-lc211[n211idx]),'m-',label='211 A')
        myFmt = DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(myFmt)
        ax.set_ylim(yran)
        ax.legend(loc='upper right')

    #plt.gcf().autofmt_xdate()
    #etc
    fig.show()
    if return_lc:
        dtvFe18=x
        lcFe18=(fe18-fe18[xidx])/np.max(fe18-fe18[xidx])
        return dtvFe18,fe18,np.max(fe18)


def plot_lc_from_pickle(picklename,idx=False):
    dtvec,sumdata,maxv=pickle.load(open(picklename,'rb'))
    #dtvecnums=date2num(dtvec)
    fig,ax=plt.subplots()
    if not idx:
        ax.plot_date(dtvec,sumdata,'b.-')
        myFmt = DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(myFmt)
    else:
        ax.plot(range(0,len(sumdata)),sumdata,'b.-')
    #plt.gcf().autofmt_xdate()
    #etc
    fig.show()
    return x,fe18

def lc_from_idl(norpfile,mint=10,save=True):
    '''run the stuff in IDL and get the results here'''
    import pidly
    os.chdir('/Users/wheatley/Documents/Solar/occulted_flares/Nobeyama/NORP/Data')
    outfile=norpfile[:-4]+'_'+str(mint/10)+'s.sav'
    idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
    idl('ff',norpfile)
    idl('filename',outfile)
    idl('restore, ff')
    idl('mint',mint)
    idl('norp_mkint,mint,mvd,tim,fi ,fv,mvdav,timav, fiav,fvav')
    idl('foo=int2sec(timav)')
    idl('oof=anytim(foo,/vms)')
    idl('timav2=oof')
    if save:
        idl('save,fiav,timav2,mvdav,filename=filename')
        os.rename(outfile,'/Users/wheatley/Documents/Solar/occulted_flares/May1Flare/NORH/'+outfile)

    lcs=idl.fiav
    time=idl.timav2
    mvd=idl.mvdav
    idl.close()
    return lcs,time,mvd

def plot_norp_lc_may1(norpfile='norp130501.p',raw=False,from_idl=False,mint=10,plot=True,save=False):
    #from scipy.io import readsav
    #norpdict=readsav(norpfile,python_dict=True)
    #nonzero=norpdict['mvd']
    #time=norpdict['tim']
    #return lcs,time
    #lcs = norpdict['fi']
    if not raw:
        lcs,times_dt=pickle.load(open(norpfile,'rb'))
    elif from_idl:
        #set default file
        norpfile='norp130501.xdr'
        lc,time,mvd=lc_from_idl(norpfile, mint=mint)
    else:
        lc,time,mvd=pickle.load(open('norp130501raw.p','rb'))
        lcs,maxes,maxnorm=[],[],[]
        #remove and smooth where mvd==0
        for i in range(0,6):
            #badidx=np.where(mvd[:,i] == 0)
            #print(len(badidx))
            #normalize by continuum time interval of 00:00 - 00:25
            times_dt=[dt.strptime(t[:-4], '%d-%b-%Y %H:%M:%S') for t in time]
            cont_int=lc[:,i][times_dt.index(dt(2013,5,1,0,4,59)):times_dt.index(dt(2013,5,1,0,19,59))]
            cont=np.mean(cont_int)
            nlc=np.ma.array(lc[:,i],mask=np.logical_not(mvd[:,i]))
            print('continuum: ', cont, '    02:32:29: ', nlc[1995], '    02:34:19: ',nlc[2006],'     max near event: ', np.max(nlc[1800:2400]), '   ratio:', np.max(nlc[1800:2400])/cont )
            #mask bad values
            lcs.append(nlc/cont)
            maxes.append(np.max(nlc[1900:2100]))
            maxnorm.append(np.max(nlc[1900:2100])/cont)
            #lcs.append(lc[:,i]/cont)
        #normalize each curve

    fig,ax=plt.subplots(figsize=(12,9))
    ax.plot_date(times_dt,lcs[0],'g-',label='1')
    ax.plot_date(times_dt,lcs[1]+.1,'m-',label='2')
    ax.plot_date(times_dt,lcs[2]+.2,'k-',label='3.75')
    ax.plot_date(times_dt,lcs[3]+.3,'c-',label='9.4')
    ax.plot_date(times_dt,lcs[4]+.4,'r-',label='17')
    ax.plot_date(times_dt,lcs[5]+.5,'b-',label='35')
    #ax.plot_date(times_dt,lcs[6]+.6,'y-',label='80')
    ax.axvline(dt(2013,5,1,2,32,28))
    ax.axvline(dt(2013,5,1,2,34,15))
    #format dates
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    #ax.set_title('Off-limb emssion, STEREO B 195 '+dt.strftime(dtvec[0],'%Y-%m-%dT%H:%M:%S')[:-9])
    ax.set_xlabel('Start Time (01-May-2013 01:30:00 UTC)')
    ax.set_ylabel('Normalized Flux')
    #if not xran:
    #    xran=[dtvecnums[0],dtvecnums2[-1]]
    ax.set_xlim([dt(2013,5,1,1,30,0),dt(2013,5,1,7,0,0)])
    ax.set_ylim([.95,2])
    ax.legend(loc='upper left')#,fontsize='small')
    #etc
    fig.show()

    fig,ax=plt.subplots()
    xvals=np.array([1.,2.,3.75,9.4,17.,34.])
    ax.scatter(xvals,maxes)
    ax.scatter(xvals, maxnorm,marker='+',color='g')
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('Max Total Flux Near Event (sfu)')
    fig.show()
    fig,ax=plt.subplots()
    xvals=np.array([1.,2.,3.75,9.4,17.,34.])
    ax.scatter(xvals, maxnorm,marker='+',color='g')
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('Max Total Flux Near Event (sfu)')
    fig.show()
    #return lcs,times_dt,mvd

def gen_norh_lc2():
    '''using the IDL save file instead'''
    dates=pickle.load(open('lcdates.p','rb'))
    #norh17file='130501_17g10s/sfu17GHz_maps.sav'
    #norh34file='130501_34g10s/sfu34GHz_maps.sav'
    adata=readsav('norh_convol.sav',python_dict=True)
    f34=adata['f34_c'] #need to convert to sfu
    f17=adata['f17_c']

    lc17=[]#f17#readsav(norh17file)['f17_c']
    lc34=[]#f34#readsav(norh34file)['f34_c']

    for i,j in zip(f17,f34):
        sub17=i[80:185,0:125]
        sub34=j[80:185,0:125]
        #print(np.mean(sub17),np.mean(sub34))
        #get rid of the zeros
        #msub17=np.ma.masked_where(sub17 != 0.0,sub17)
        #msub34=np.ma.masked_where(sub34 != 0.0,sub34)
        #average
        #lc34.append(msub34.mean())
        #lc17.append(msub17.mean())
        lc34.append(np.mean(sub34))
        lc17.append(np.mean(sub17))

    print(len(lc17),len(lc34),np.mean(lc17),np.mean(lc34))
    fig,ax=plt.subplots()
    #dax=ax.twinx()
    #ax.plot(dates[1],lc34,'b',label='$\\alpha$')
    ax.plot(dates[1],lc17,'r',label='17 GHz')
    ax.plot(dates[1],lc34,'g',label='34 GHz')
    ax.set_ylabel('Flux (sfu/steradian)')
    ax.legend(loc='upper left')
    ax.set_xlim([dates[1][0],dates[1][-1]])
    #ax.set_ylabel('Mean Spectral Index $\\alpha$')
    ax.set_xlabel('Time on 01 May 2013')
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    #if with_lc:
    #    fig.savefig('alpha_and_lc.png')
    #else:
    #    fig.savefig('alpha.png')
    fig.show()
    return lc17,lc34


def gen_norh_lc(submap=False,plot=True,save=True,fac=10.,save_list=True):
    fits=glob.glob('ipz130501_02*')#.fits')
    fits.sort()
    maplist=[]
    lc=[]
    dtvec=[]
    for i,f in enumerate(fits): #fix headers and re-make maps
        try:
            sm=sunpy.map.Map(f)
            if f.startswith('ifa130430'):
                dtvec.append(dt.strptime('2013-04-30 '+f[-11:-5],'%Y-%m-%d %H%M%S'))
            else:
                dtvec.append(dt.strptime(f[3:],'%y%m%d_%H%M%S'))
        except IOError:
            continue
        data=sm.data
        meta=sm.meta
        meta['ctype1']='HPLN-TAN'
        meta['ctype2']='HPLT-TAN'
        nm=sunpy.map.Map(data,meta)
        maplist.append(nm)
    if save_list:
        pickle.dump(maplist,open('maplist.p','wb'))
    for m in maplist: #trim, mask, etc
        lcp=get_lc_from_map(m,submap=submap,maskdisk=True,fac=fac,backgroundsubtract=False,peek=False,dtv=False)
        lc.append(lcp)
        #dtvec.append(dt.strptime('2013-05-01 '+m.meta['time_obs'][-4],'%Y-%m-%d %H:%M:%S'))
    if plot:
        fig,ax=plt.subplots(figsize=(12,9))
        ax.plot_date(dtvec,lc,'g-')
        #ax.axvline(dt(2013,5,1,2,32,28))
        #ax.axvline(dt(2013,5,1,2,34,15))
        #format dates
        myFmt = DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(myFmt)
        plt.gcf().autofmt_xdate()
        #ax.set_title('Off-limb emssion, STEREO B 195 '+dt.strftime(dtvec[0],'%Y-%m-%dT%H:%M:%S')[:-9])
        ax.set_xlabel('Start Time (01-May-2013 01:30:00 UTC)')
        ax.set_ylabel('Normalized Flux')
        #if not xran:
        #    xran=[dtvecnums[0],dtvecnums2[-1]]
        #ax.set_xlim([dt(2013,5,1,1,30,0),dt(2013,5,1,7,0,0)])
        #ax.set_ylim([.95,2])
        ax.legend(loc='upper left')#,fontsize='small')
        #etc
        fig.show()
    return dtvec,lc

def gen_alpha_lc(submap=False,plot=True,save=True,fac=10.,mean=True):
    '''can't do it as a map so just mask the array'''
    data=readsav('spectral_index.sav',python_dict=True)
    alpha=data['alpha']
    bad=data['mvda']
    headers=data['indexal']
    lc=[]
    #load dtvec from previous lc
    dates17,dates34=pickle.load(open('lcdates.p','rb'))
    for m,a in zip(bad,alpha):#i in range(0,np.shape(alpha)[0]):
        masked_alpha=np.ma.masked_array(a,mask=m)
        good_alpha=np.ma.filled(masked_alpha,0.0)
        if np.sum(m) == 0:
            if not mean:
                lc.append(np.sum(good_alpha)/(256.*256.))
            else:
                lc.append(np.mean(good_alpha))
    if plot:
        fig,ax=plt.subplots(figsize=(12,9))
        ax.plot_date(dates17[:174],lc,'g-')
        #ax.axvline(dt(2013,5,1,2,32,28))
        #ax.axvline(dt(2013,5,1,2,34,15))
        #format dates
        myFmt = DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(myFmt)
        plt.gcf().autofmt_xdate()
        #ax.set_title('Off-limb emssion, STEREO B 195 '+dt.strftime(dtvec[0],'%Y-%m-%dT%H:%M:%S')[:-9])
        ax.set_xlabel('Start Time (01-May-2013 01:30:00 UTC)')
        ax.set_ylabel('Mean alpha')
        #if not xran:
        #    xran=[dtvecnums[0],dtvecnums2[-1]]
        #ax.set_xlim([dt(2013,5,1,1,30,0),dt(2013,5,1,7,0,0)])
        #ax.set_ylim([.95,2])
        ax.legend(loc='upper left')#,fontsize='small')
        #etc
        fig.show()
    #return dtvec,lc

def plot_norh_lc(norh17file='130501_17g10s/sfu17GHz.sav',norh34file='130501_34g10s/sfu34GHz.sav',get_times=False,bsub=True):
    if get_times:
        fits17=glob.glob('130501_17g10s/ifa130501_02*')
        fits34=glob.glob('130501_34g10s/ipz130501_02*')
        dates17,dates34=[],[]
        for f1 in fits17:
            m17=sunpy.map.Map(f1)
            dates17.append(dt.strptime(m17.meta['date-obs']+m17.meta['time-obs'][:-4],'%Y-%m-%d%H:%M:%S'))
        for f2 in fits34:
            m34=sunpy.map.Map(f2)
            dates34.append(dt.strptime(m34.meta['date-obs']+m34.meta['time-obs'][:-4],'%Y-%m-%d%H:%M:%S'))
    else:
        dates17,dates34=pickle.load(open('lcdates.p','rb'))
    lc17=readsav(norh17file)
    lc34=readsav(norh34file)
    print(len(dates17),np.shape(lc17['fi2']))

    if bsub:
        plot17=lc17['fi2']-lc17['fi2'][0]
        plot34=lc34['fi2']-lc34['fi2'][0]
    else:
        plot17=lc17['fi2']
        plot34=lc34['fi2']


    fig,ax=plt.subplots()
    ax.plot_date(dates17,plot17,'k',label='17 GHz')
    ax.plot_date(dates34,plot34,'m',label='34 GHz')
    ax.set_ylabel('SFU/steradian')
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    ax.legend(loc='upper left')
    fig.show()
    if get_times:
        pickle.dump([dates17,dates34],open('lcdates.p','wb'))

def plot_norp_lc(norpfile='norp130501.p',raw=False,from_idl=False,mint=10,norm_int=False,plot=True,save=False,write=False,peaktime=False,time_int=False):
    import matplotlib.pyplot as plt
    if from_idl:
        lc,time,mvd=lc_from_idl(norpfile, mint=mint)
    else:
        lc,time,mvd=pickle.load(open(norpfile,'rb'))
    lcs=[]
    #remove and smooth where mvd==0
    for i in range(0,6):
        times_dt=[dt.strptime(t[:-4], '%d-%b-%Y %H:%M:%S') for t in time]
        if not norm_int:
            cont=np.median(np.ma.array(lc[:,i],mask=np.logical_not(mvd[:,i])))#_int=lc[:,i][100:120]
        else:
            cont_int=lc[:,i][times_dt.index(dt(2013,5,1,0,4,59)):times_dt.index(dt(2013,5,1,0,19,59))]
            cont=np.mean(cont_int)
        #mask bad values
        lcs.append(np.ma.array(lc[:,i],mask=np.logical_not(mvd[:,i]))/cont)

    if not raw and not from_idl: #assume normalize already to desired inteval
        lcs,times_dt=pickle.load(open(norpfile,'rb'))

    if write:
        pickle.dump([lcs,times_dt,mvd],open(norpfile[:-4]+'_'+str(mint)+'s.p','wb'))

    fig,ax=plt.subplots(figsize=(12,9))
    ax.plot_date(times_dt,lcs[0],'g-',label='1')
    ax.plot_date(times_dt,lcs[1]+.1,'m-',label='2')
    ax.plot_date(times_dt,lcs[2]+.2,'k-',label='3.75')
    ax.plot_date(times_dt,lcs[3]+.3,'c-',label='9.4')
    ax.plot_date(times_dt,lcs[4]+.4,'r-',label='17')
    ax.plot_date(times_dt,lcs[5]+.5,'b-',label='35')
    if peaktime:
        ax.axvline(peaktime)
        pkstr=dt.strftime(peaktime,'%d %b %Y %H:%M:%S')
        #print(pkstr,pkstr[12:17])
    #ax.plot_date(times_dt,lcs[6]+.6,'y-',label='80')
    #format dates
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    #ax.set_title('Off-limb emssion, STEREO B 195 '+dt.strftime(dtvec[0],'%Y-%m-%dT%H:%M:%S')[:-9])
    ax.set_xlabel('Peak Time ('+ pkstr +' UTC)')
    ax.set_ylabel('Normalized Flux')
    #ax.set_xlim([dt(2013,5,1,1,30,0),dt(2013,5,1,7,0,0)])
    if time_int:
        ax.set_xlim(time_int)
    #ax.set_ylim([.95,2])
    ax.legend(loc='upper left')#,fontsize='small')
    #etc
    if plot:
        fig.show()
    if save:
        plt.savefig(norpfile[:-4]+'_'+pkstr[12:14]+'_'+pkstr[15:17]+'_lc.png')
    #return lcs,times_dt,mvd

def get_lc_from_map(imap,submap=False,maskdisk=False,fac=0.,backgroundsubtract=False,peek=False,dtv=True):
    '''extract single data point for lightcurve of given map'''
    reshape=False
    if submap:
        reshape=check_coords(imap,submap)
    if reshape: #reshape the maps.. from the fits files I guess. or full disk
        #try:
        #print('test1')
        imap=imap.submap(SkyCoord(submap[0]*u.arcsec,submap[1]*u.arcsec,frame=imap.coordinate_frame))
        #except ValueError:
        #coordinates out of range
        #fitsmaps=sunpy.map.Map(all_fits,cube=True)
        #print(type(fitsmaps))
        #all_maps=[mc.submap(SkyCoord(submap[0]*u.arcsec,submap[1]*u.arcsec,frame=fitsmaps[0].coordinate_frame)) for mc in fitsmaps]
    if maskdisk:
        masked_map=mask_disk(imap,fac=fac,plot=peek)
        sumdata=np.sum(masked_map.data[np.where(masked_map.mask == False)])
    else:
        sumdata=np.sum(imap.data)
    if type(backgroundsubtract)==np.float64:
        sumdata=sumdata-backgroundsubtract
    if np.isfinite(sumdata) and dtv:
        dtval=dt.strptime(imap.meta['date-obs'][:-4],'%Y-%m-%dT%H:%M:%S')
    if dtv:
        return sumdata,dtval
    else:
        return sumdata

def plot_from_pickle(picklename, picklename2=False,title=False,labels=['STEREO B 195','STEREO A 171'],yran=[.5,1.1],xran=False):
    dlist=pickle.load(open(picklename,'rb'))
    dtvec=dlist[0]
    sumdata=dlist[1]
    dtvecnums=date2num(dtvec)
    if picklename2:
        dlist2=pickle.load(open(picklename2,'rb'))
        dtvec2=dlist2[0]
        sumdata2=dlist2[1]
        dtvecnums2=date2num(dtvec2)

    fig,ax=plt.subplots()
    ax.plot_date(dtvecnums,sumdata,'c.-',label=labels[0])
    if picklename2:
        ax.plot_date(dtvecnums2,sumdata2,'m.-',label=labels[1])
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    if title:
        ax.set_title(title)
    else:
        ax.set_title('Off-limb emssion, STEREO B 195 '+dt.strftime(dtvec[0],'%Y-%m-%dT%H:%M:%S')[:-9])
    ax.set_xlabel('Time (UTC)')
    ax.set_ylabel('Normalized Flux')
    if not xran:
        xran=[dtvecnums[0],dtvecnums2[-1]]
        ax.set_xlim(xran)
    ax.set_ylim(yran)
    ax.legend()
    #etc
    fig.show()

def plot_xsi_rhessi_euv(xsi_rhessi='plot_curves.sav',stereoB='../STEREO/B/stereoB195curve.p',stereoA='../STEREO/A/171curvedata.p',background_subtract=True):
    '''make Sam's figure
        utplot,lc_tm.time,lc_tm.flux/max(lc_tm.flux(0:*)),timerange=['1-May-2013 01:30:00','1-May-2013 10:00:00'],xstyle=1,yrange=[-0.05,1.05],ystyle=1,psym=-1
        outplot,lc_be.time,lc_be.flux/max(lc_be.flux(0:*)),color=6,psym=-1,thick=2
        outplot,lc_tm.time,lc_tm.flux/max(lc_tm.flux(0:*)),color=1,psym=-1,thick=2
        ;same around CME escape time with RHESSI
        utplot,lc_tm.time,lc_tm.flux/max(lc_tm.flux(0:55)),timerange=['1-May-2013 01:30:00','1-May-2013 03:00:00'],xstyle=1,yrange=[-0.05,1.05],ystyle=1,psym=-1
        outplot,lc_be.time,lc_be.flux/max(lc_be.flux(0:18)),color=6,psym=-1,thick=2
        outplot,lc_tm.time,lc_tm.flux/max(lc_tm.flux(10:55)),color=1,psym=-1,thick=2
        outplot,hsi_time,hsi12/max(hsi12),psym=10,thick=2,color=2
        outplot,time_thermal,flux_thermal/max(flux_thermal),thick=2
    '''
    from scipy.io import readsav
    xrdict=readsav(xsi_rhessi,python_dict=True)
    hsi_time=xrdict['hsi_time']
    hsi12=xrdict['hsi12']/np.max(xrdict['hsi12'])
    time_thermal=xrdict['time_thermal']
    flux_thermal=xrdict['flux_thermal']/np.max(xrdict['flux_thermal'])
    lc_tm=xrdict['lc_tm'] #has .time and .flux. don't use .time
    lc_be=xrdict['lc_be']
    tm_flux=lc_tm['flux']/np.max(lc_tm['flux'][0][:55])
    be_flux=lc_be['flux']/np.max(lc_be['flux'][0][:18])
    print(np.shape(be_flux[0]),np.shape(tm_flux[0]))
    stereoBl=pickle.load(open(stereoB,'rb'))
    stereoAl=pickle.load(open(stereoA,'rb'))
    sB_time=stereoBl[0]
    sB_curve=stereoBl[1]
    sA_time=stereoAl[0]
    sA_curve=stereoAl[1]
    if background_subtract:
        sB_curve=sB_curve-sB_curve[0]
        sA_curve=sA_curve-sA_curve[0]

    #assume IDL anytim axes were converted via anytim(time,/vms) into format %d-%b%-Y %H:%M:%S
    hsi_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in hsi_time]
    thermal_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in time_thermal]
    tm_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in xrdict['lc_tm_time']]
    be_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in xrdict['lc_be_time']]
    print(np.shape(be_dt),np.shape(tm_dt))
    #plot
    fig,ax=plt.subplots(figsize=(12,9))
    ax.plot_date(sB_time,sB_curve,'g.-',label='STEREO-B 195')
    ax.plot_date(sA_time,sA_curve,'m.-',label='STEREO-A 195')
    #ax.plot_date(hsi_dt,hsi12,'b.-',label='RHESSI 10-30 keV')
    ax.step(hsi_dt,hsi12,'b-',label='RHESSI 10-30 keV')
    ax.plot_date(thermal_dt,flux_thermal,'k-',label='RHESSI 4-8 keV')
    ax.plot_date(tm_dt,tm_flux[0],'c+-',label='XSI TM')
    ax.plot_date(be_dt,be_flux[0],'r+-',label='XSI Be')
    #format dates
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    #ax.set_title('Off-limb emssion, STEREO B 195 '+dt.strftime(dtvec[0],'%Y-%m-%dT%H:%M:%S')[:-9])
    ax.set_xlabel('Start Time (01-May-2013 01:30:00 UTC)')
    ax.set_ylabel('Normalized Flux')
    #if not xran:
    #    xran=[dtvecnums[0],dtvecnums2[-1]]
    ax.set_xlim([tm_dt[22],tm_dt[43]])
    ax.set_ylim([-.1,1.1])
    ax.legend(loc='center left')#,fontsize='small')
    #etc
    fig.show()
    return tm_dt

def plot_rhessi_zoom(xsi_rhessi='plot_curves.sav',AIA='../AIA/high_cadence_cutout/lightcurve94nobg.p'):
    from scipy.io import readsav
    xrdict=readsav(xsi_rhessi,python_dict=True)
    hsi_time=xrdict['hsi_time']
    hsi12=xrdict['hsi12']/np.max(xrdict['hsi12'])
    time_thermal=xrdict['time_thermal']
    flux_thermal=xrdict['flux_thermal']/np.max(xrdict['flux_thermal'])
    dtv94,lc94,max94=pickle.load(open(AIA,'rb'))
    #mask bad value
    lc94[85]=np.mean([lc94[84],lc94[86]])
    os.chdir('../NORH')
    dtvNORH,lcNORH,aa=pickle.load(open('norh_lc.p','rb'))
    os.chdir('../')
    #assume IDL anytim axes were converted via anytim(time,/vms) into format %d-%b%-Y %H:%M:%S
    hsi_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in hsi_time]
    thermal_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in time_thermal]
    #plot
    fig,ax=plt.subplots(figsize=(6,6))
    #ax.plot_date(sB_time,sB_curve,'g.-',label='STEREO-B 195')
    #ax.plot_date(sA_time,sA_curve,'m.-',label='STEREO-A 195')
    #ax.plot_date(hsi_dt,hsi12,'b.-',label='RHESSI 10-30 keV')
    ax.step(hsi_dt,hsi12,'b-',label='RHESSI 10-30 keV',linewidth=2)
    ax.plot_date(thermal_dt,flux_thermal,'k-',label='RHESSI 4-8 keV',linewidth=2)
    ax.plot_date(dtv94,(lc94-lc94[70])/(max94-lc94[70]),'g.-',label='AIA 94 $\AA$')
    ax.plot_date(dtvNORH,(lcNORH-lcNORH[15])/np.max(lcNORH-lcNORH[15]),'m.-',label='NORH 17 GHz')
    #format dates
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    #ax.set_title('Off-limb emssion, STEREO B 195 '+dt.strftime(dtvec[0],'%Y-%m-%dT%H:%M:%S')[:-9])
    ax.set_xlabel('Start Time (01-May-2013 01:30:00 UTC)')
    ax.set_ylabel('Normalized Flux')
    #if not xran:
    #    xran=[dtvecnums[0],dtvecnums2[-1]]
    ax.set_xlim([thermal_dt[0],dtvNORH[-15]])
    ax.set_ylim([-.5,1.1])
    ax.legend(loc='lower right')#,fontsize='small')
    #etc
    fig.show()

def plot_xsi_euv(xsi_rhessi='plot_curves.sav',stereoB='../STEREO/B/stereoB195curve.p',stereoA='../STEREO/A/171curvedata.p'):
    '''make Sam's figure
        utplot,lc_tm.time,lc_tm.flux/max(lc_tm.flux(0:*)),timerange=['1-May-2013 01:30:00','1-May-2013 10:00:00'],xstyle=1,yrange=[-0.05,1.05],ystyle=1,psym=-1
        outplot,lc_be.time,lc_be.flux/max(lc_be.flux(0:*)),color=6,psym=-1,thick=2
        outplot,lc_tm.time,lc_tm.flux/max(lc_tm.flux(0:*)),color=1,psym=-1,thick=2
        ;same around CME escape time with RHESSI
        utplot,lc_tm.time,lc_tm.flux/max(lc_tm.flux(0:55)),timerange=['1-May-2013 01:30:00','1-May-2013 03:00:00'],xstyle=1,yrange=[-0.05,1.05],ystyle=1,psym=-1
        outplot,lc_be.time,lc_be.flux/max(lc_be.flux(0:18)),color=6,psym=-1,thick=2
        outplot,lc_tm.time,lc_tm.flux/max(lc_tm.flux(10:55)),color=1,psym=-1,thick=2
        outplot,hsi_time,hsi12/max(hsi12),psym=10,thick=2,color=2
        outplot,time_thermal,flux_thermal/max(flux_thermal),thick=2
    '''
    from scipy.io import readsav
    xrdict=readsav(xsi_rhessi,python_dict=True)
    #hsi_time=xrdict['hsi_time']
    #hsi12=xrdict['hsi12']/np.max(xrdict['hsi12'])
    #time_thermal=xrdict['time_thermal']
    #flux_thermal=xrdict['flux_thermal']/np.max(xrdict['flux_thermal'])
    lc_tm=xrdict['lc_tm'] #has .time and .flux. don't use .time
    lc_be=xrdict['lc_be']
    tm_flux=lc_tm['flux']/np.max(lc_tm['flux'][0])
    be_flux=lc_be['flux']/np.max(lc_be['flux'][0])
    #print(np.shape(be_flux[0]),np.shape(tm_flux[0]))
    stereoBl=pickle.load(open(stereoB,'rb'))
    stereoAl=pickle.load(open(stereoA,'rb'))
    sB_time=stereoBl[0]
    sB_curve=stereoBl[1]
    sA_time=stereoAl[0]
    sA_curve=stereoAl[1]

    #assume IDL anytim axes were converted via anytim(time,/vms) into format %d-%b%-Y %H:%M:%S
    #hsi_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in hsi_time]
    #thermal_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in time_thermal]
    tm_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in xrdict['lc_tm_time']]
    be_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in xrdict['lc_be_time']]
    print(np.shape(be_dt),np.shape(tm_dt))
    #plot
    fig,ax=plt.subplots(figsize=(12,9))
    ax.plot_date(sB_time,sB_curve,'g.-',label='STEREO-B 195')
    ax.plot_date(sA_time,sA_curve,'m.-',label='STEREO-A 195')
    #ax.plot_date(hsi_dt,hsi12,'b.-',label='RHESSI 10-30 keV')
    #ax.step(hsi_dt,hsi12,'b-',label='RHESSI 10-30 keV')
    #ax.plot_date(thermal_dt,flux_thermal,'k-',label='RHESSI 4-8 keV')
    ax.plot_date(tm_dt,tm_flux[0],'c+-',label='XSI TM')
    ax.plot_date(be_dt,be_flux[0],'r+-',label='XSI Be')
    #format dates
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    #ax.set_title('Off-limb emssion, STEREO B 195 '+dt.strftime(dtvec[0],'%Y-%m-%dT%H:%M:%S')[:-9])
    ax.set_xlabel('Start Time (01-May-2013 01:30:00 UTC)')
    ax.set_ylabel('Normalized Flux')
    #if not xran:
    #    xran=[dtvecnums[0],dtvecnums2[-1]]
    ax.set_xlim([tm_dt[22],tm_dt[-1]])
    ax.set_ylim([-.1,1.1])
    ax.legend(loc='lower right')#,fontsize='small')
    #etc
    fig.show()
    return tm_dt


def figure2(returnp=False):
    #read data for hot stuff
    from scipy.io import readsav
    xsi_rhessi='Xray/plot_curves.sav'
    xrdict=readsav(xsi_rhessi,python_dict=True)
    lc_tm=xrdict['lc_tm'] #has .time and .flux. don't use .time
    lc_be=xrdict['lc_be']
    tm_flux=lc_tm['flux']/np.max(lc_tm['flux'][0])
    #fix bad point
    tm_flux[0][3]=np.mean([tm_flux[0][2],tm_flux[0][4]])
    be_flux=lc_be['flux']/np.max(lc_be['flux'][0])
    #dtv94,lc94,aa=pickle.load(open('lightcurve_94_ext.p','rb'))#'lightcurve94.p','rb'))
    #dtvFe18,lcFe18=pickle.load(open('lightcurve_Fe18_ext.p','rb'))#'lightcurveFe18.p','rb'))
    #assume IDL anytim axes were converted via anytim(time,/vms) into format %d-%b%-Y %H:%M:%S
    tm_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in xrdict['lc_tm_time']]
    be_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in xrdict['lc_be_time']]

    #read data for cold stuff
    AIA=['lightcurve_131_ext.p','lightcurve_171_ext.p','lightcurve_193_ext.p']#['lightcurve131.p','lightcurve171.p','lightcurve193.p']
    #stereoB='../STEREO/B/stereoBreconstructed.p'
    #stereoA='../STEREO/A/195/stereoAreconstructed.p'
    os.chdir('AIA/low_cadence_cutout')
    dtv94,lc94,aa=pickle.load(open('lightcurve_94_ext.p','rb'))#'lightcurve94.p','rb'))
    #print(aa)
    dtvFe18,lcFe18,aa=pickle.load(open('lightcurve_Fe18_ext.p','rb'))#'lightcurveFe18.p','rb'))
    dtv131,lc131,aa=pickle.load(open(AIA[0],'rb'))
    dtv171,lc171,aa=pickle.load(open(AIA[1],'rb'))
    dtv193,lc193,aa=pickle.load(open(AIA[2],'rb'))
    stereoBl=pickle.load(open('lightcurve_sB195_ext.p','rb'))
    stereoAl=pickle.load(open('lightcurve_sA195_ext.p','rb'))
    os.chdir('../../NORH')
    dtvNORH,lcNORH,aa=pickle.load(open('norh_lc.p','rb'))
    #stereoBl=pickle.load(open(stereoB,'rb'))
    #stereoAl=pickle.load(open(stereoA,'rb'))
    os.chdir('../')
    sB_time=stereoBl[0]
    sB_curve=stereoBl[1]
    sA_time=stereoAl[0]
    sA_curve=stereoAl[1]

    #plot hot stuff on top, cold stuff on bottom
    fig,ax=plt.subplots(2,sharex=True,figsize=(12,9))
    ax[0].plot_date(dtv94,(lc94-lc94[0])/np.max((lc94-lc94[0])),'g.-',label='AIA 94')
    ax[0].plot_date(dtvFe18,(lcFe18-lcFe18[0])/np.max((lcFe18-lcFe18[0])),'m.-',label='AIA FeXVIII')
    ax[0].plot_date(tm_dt,tm_flux[0],'c+-',label='XSI TM')
    ax[0].plot_date(be_dt,be_flux[0],'r+-',label='XSI Be')
    #ax[0].plot_date(dtvNORH,(lcNORH-lcNORH[7])/np.max(lcNORH-lcNORH[7]),'k.-',label='NORH 17 GHz')
    ax[0].set_ylabel('Normalized Flux')
    ax[0].set_ylim([-.5,1.1])
    ax[0].legend(loc='upper left')#,fontsize='small')

    ax[1].plot_date(sB_time,sB_curve,'b.-',label='STEREO-B 195')
    ax[1].plot_date(sA_time,sA_curve,'m.-',label='STEREO-A 195')
    ax[1].plot_date(dtv131,(lc131-lc131[0])/np.max((lc131-lc131[0])),'c.-',label='AIA 131')
    ax[1].plot_date(dtv171,(lc171-lc171[0])/np.max((lc171-lc171[0])),'k.-',label='AIA 171')
    ax[1].plot_date(dtv193,(lc193-lc193[0])/np.max((lc193-lc193[0])),'g.-',label='AIA 193') #this really messes up the y-axis scale when normalized at pre-flare... maybe should normalize post-flare? Or pick an image way pre-flare (00:00) as the basis?
    ax[1].set_ylim([-3.5,1.1])
    ax[1].set_ylabel('Normalized Flux')
    ax[1].legend(loc='lower left')#,fontsize='small')

    #format dates
    myFmt = DateFormatter('%H:%M')
    ax[1].xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    ax[1].set_xlabel('Start Time (01-May-2013 01:30:00 UTC)')
    ax[1].set_xlim([dtv94[0],dtv94[-1]])

    #adjust subplots
    plt.subplots_adjust(top=.95,bottom=.05,right=.95,hspace=.05)
    #fig.show()
    if returnp:
        plotdata=[dtv94,lc94,dtvFe18,lcFe18,tm_dt,tm_flux,be_dt,be_flux,sB_time,sB_curve,sA_time,sA_curve,dtv131,lc131,dtv171,lc171,dtv193,lc193]
        return plotdata
    else:
        fig.show()

def animate_figure2(plotdata,timedata,moviename,framerate=12,show=False):
    '''Animate a bar sliding across the t-axis to indicate where it is snyced with the movie. timedata is a list of datetimes extracted from the file groups or synced maps'''
    import make_aligned_movie as mov
    #assign plotdata back to values
    dtv94,lc94,dtvFe18,lcFe18,tm_dt,tm_flux,be_dt,be_flux,sB_time,sB_curve,sA_time,sA_curve,dtv131,lc131,dtv171,lc171,dtv193,lc193=plotdata
    for i,t in enumerate(timedata):
        fig,ax=plt.subplots(2,sharex=True,figsize=(12,9))
        ax[0].plot_date(dtv94,(lc94-lc94[0])/np.max((lc94-lc94[0])),'g.-',label='AIA 94')
        ax[0].plot_date(dtvFe18,(lcFe18-lcFe18[0])/np.max((lcFe18-lcFe18[0])),'m.-',label='AIA FeXVIII')
        #ax[0].plot_date(dtvFe18,lcFe18,'m.-',label='AIA FeXVIII')
        ax[0].plot_date(tm_dt,tm_flux[0],'c+-',label='XSI TM')
        ax[0].plot_date(be_dt,be_flux[0],'r+-',label='XSI Be')
        ax[0].set_ylabel('Normalized Flux')
        ax[0].set_ylim([-.5,1.1])
        ax[0].legend(loc='upper left')#,fontsize='small')

        ax[1].plot_date(sB_time,sB_curve,'b.-',label='STEREO-B 195')
        ax[1].plot_date(sA_time,sA_curve,'m.-',label='STEREO-A 195')
        ax[1].plot_date(dtv131,(lc131-lc131[0])/np.max((lc131-lc131[0])),'c.-',label='AIA 131')
        ax[1].plot_date(dtv171,(lc171-lc171[0])/np.max((lc171-lc171[0])),'k.-',label='AIA 171')
        ax[1].plot_date(dtv193,(lc193-lc193[0])/np.max((lc193-lc193[0])),'g.-',label='AIA 193') #this really messes up the y-axis scale when normalized at pre-flare... maybe should normalize post-flare? Or pick an image way pre-flare (00:00) as the basis?
        ax[1].set_ylim([-3.5,1.1])
        ax[1].set_ylabel('Normalized Flux')
        ax[1].legend(loc='lower left')#,fontsize='small')

        #add the vlines
        ax[0].axvline(t,linestyle='dashed',color='k')
        ax[1].axvline(t,linestyle='dashed',color='k')

        #format dates
        myFmt = DateFormatter('%H:%M')
        ax[1].xaxis.set_major_formatter(myFmt)
        plt.gcf().autofmt_xdate()
        ax[1].set_xlabel('Start Time (01-May-2013 01:30:00 UTC)')
        ax[1].set_xlim([sB_time[0],sB_time[-47]])

        #adjust subplots
        plt.subplots_adjust(top=.95,bottom=.05,right=.95,hspace=.05)
        if show:
            fig.show()
        else:
            plt.savefig('fig2_'+'{0:03d}'.format(i)+'.png')
            plt.clf()

    #animate
    mov.run_ffmpeg('fig2_'+'%03d.png',moviename,framerate=framerate)

def extend_lc(preflare,original,wave):
    pf=pickle.load(open(preflare,'rb'))
    f=pickle.load(open(original,'rb'))
    fac=pf[1][-1]/f[1][0]
    pfadj=pf[1]/fac
    pf[0].extend(f[0]) #new timevector
    #nt=np.append(pf[0],f[0])
    lcfull=np.append(pfadj,f[1])
    pickle.dump([pf[0],lcfull,f[2]],open('lightcurve_'+wave+'_ext.p','wb'))

def plot_cold_stuff(AIA=['lightcurve131.p','lightcurve171.p','lightcurve193.p'],stereoB='../STEREO/B/stereoBreconstructed.p',stereoA='../STEREO/A/195/lc_bsub_stereoA.p'):
    #need to re-do stereo lightcurves without the funny normalization!
    os.chdir('AIA/low_cadence_cutout')
    dtv131,lc131,aa=pickle.load(open(AIA[0],'rb'))
    dtv171,lc171,aa=pickle.load(open(AIA[1],'rb'))
    dtv193,lc193,aa=pickle.load(open(AIA[2],'rb'))
    os.chdir('../')
    #print(np.shape(be_flux[0]),np.shape(tm_flux[0]))
    stereoBl=pickle.load(open(stereoB,'rb'))
    stereoAl=pickle.load(open(stereoA,'rb'))
    os.chdir('..')
    sB_time=stereoBl[0]
    sB_curve=stereoBl[1]
    sA_time=stereoAl[0]
    sA_curve=stereoAl[1]

    #plot
    fig,ax=plt.subplots(figsize=(12,9))
    ax.plot_date(sB_time,sB_curve,'b.-',label='STEREO-B 195')
    ax.plot_date(sA_time,sA_curve,'m.-',label='STEREO-A 195')
    ax.plot_date(dtv131,(lc131-lc131[0])/np.max((lc131-lc131[0])),'c.-',label='AIA 131')
    ax.plot_date(dtv171,(lc171-lc171[0])/np.max((lc171-lc171[0])),'k.-',label='AIA 171')
    #ax.plot_date(dtv193,(lc193-lc193[0])/np.max((lc193-lc193[0])),'g.-',label='AIA 193')
    #format dates
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    #ax.set_title('Off-limb emssion, STEREO B 195 '+dt.strftime(dtvec[0],'%Y-%m-%dT%H:%M:%S')[:-9])
    ax.set_xlabel('Start Time (01-May-2013 01:30:00 UTC)')
    ax.set_ylabel('Normalized Flux')
    #if not xran:
    #    xran=[dtvecnums[0],dtvecnums2[-1]]
    ax.set_xlim([sB_time[0],sB_time[-47]])
    ax.set_ylim([-1,1.1])
    ax.legend(loc='upper right')#,fontsize='small')
    #etc
    fig.show()
    #return tm_dt

def plot_hot_stuff(xsi_rhessi='Xray/plot_curves.sav',Fe18='lightcurveFe18.p',AIA94='lightcurve94.p'):
    '''make Sam's figure
        utplot,lc_tm.time,lc_tm.flux/max(lc_tm.flux(0:*)),timerange=['1-May-2013 01:30:00','1-May-2013 10:00:00'],xstyle=1,yrange=[-0.05,1.05],ystyle=1,psym=-1
        outplot,lc_be.time,lc_be.flux/max(lc_be.flux(0:*)),color=6,psym=-1,thick=2
        outplot,lc_tm.time,lc_tm.flux/max(lc_tm.flux(0:*)),color=1,psym=-1,thick=2
        ;same around CME escape time with RHESSI
        utplot,lc_tm.time,lc_tm.flux/max(lc_tm.flux(0:55)),timerange=['1-May-2013 01:30:00','1-May-2013 03:00:00'],xstyle=1,yrange=[-0.05,1.05],ystyle=1,psym=-1
        outplot,lc_be.time,lc_be.flux/max(lc_be.flux(0:18)),color=6,psym=-1,thick=2
        outplot,lc_tm.time,lc_tm.flux/max(lc_tm.flux(10:55)),color=1,psym=-1,thick=2
        outplot,hsi_time,hsi12/max(hsi12),psym=10,thick=2,color=2
        outplot,time_thermal,flux_thermal/max(flux_thermal),thick=2
    '''
    from scipy.io import readsav
    xrdict=readsav(xsi_rhessi,python_dict=True)
    #hsi_time=xrdict['hsi_time']
    #hsi12=xrdict['hsi12']/np.max(xrdict['hsi12'])
    #time_thermal=xrdict['time_thermal']
    #flux_thermal=xrdict['flux_thermal']/np.max(xrdict['flux_thermal'])
    lc_tm=xrdict['lc_tm'] #has .time and .flux. don't use .time
    lc_be=xrdict['lc_be']
    tm_flux=lc_tm['flux']/np.max(lc_tm['flux'][0])
    be_flux=lc_be['flux']/np.max(lc_be['flux'][0])
    #print(np.shape(be_flux[0]),np.shape(tm_flux[0]))
    dtv94,lc94,aa=pickle.load(open('lightcurve94.p','rb'))
    dtvFe18,lcFe18=pickle.load(open('lightcurveFe18.p','rb'))


    #assume IDL anytim axes were converted via anytim(time,/vms) into format %d-%b%-Y %H:%M:%S
    #hsi_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in hsi_time]
    #thermal_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in time_thermal]
    tm_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in xrdict['lc_tm_time']]
    be_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in xrdict['lc_be_time']]
    print(np.shape(be_dt),np.shape(tm_dt))
    #plot
    fig,ax=plt.subplots(figsize=(12,9))
    ax.plot_date(dtv94,(lc94-lc94[0])/np.max((lc94-lc94[0])),'g.-',label='AIA 94')
    ax.plot_date(dtvFe18,lcFe18,'m.-',label='AIA FeXVIII')
    ax.plot_date(tm_dt,tm_flux[0],'c+-',label='XSI TM')
    ax.plot_date(be_dt,be_flux[0],'r+-',label='XSI Be')
    #format dates
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    #ax.set_title('Off-limb emssion, STEREO B 195 '+dt.strftime(dtvec[0],'%Y-%m-%dT%H:%M:%S')[:-9])
    ax.set_xlabel('Start Time (01-May-2013 01:30:00 UTC)')
    ax.set_ylabel('Normalized Flux')
    #if not xran:
    #    xran=[dtvecnums[0],dtvecnums2[-1]]
    ax.set_xlim([tm_dt[22],tm_dt[-30]])
    ax.set_ylim([-.5,1.1])
    ax.legend(loc='lower right')#,fontsize='small')
    #etc
    fig.show()
    #return tm_dt


def reconstruct_lc(picklename,bg,newbgidx,newpicklename=False):
      dtv,lc,maxv=pickle.load(open(picklename,'rb'))
      #add the background,subtract new bg
      lcn=[(l*maxv)+bg for l in lc]
      nbg=lcn[newbgidx]
      lcn=[l-nbg for l in lcn]
      #re-normalize
      lcnn=lcn/np.max(lcn)
      if newpicklename:
            pickle.dump([dtv,lcnn,maxv],open(newpicklename,'wb'))
      return lcnn

def plot_AIA_euv(yran=[-1.1,1.1],tnorm=False,lpos='upper left'):
    '''plot all AIA lightcurvez'''
    dtv94,lc94,aa=pickle.load(open('lightcurve94.p','rb'))
    dtv131,lc131,aa=pickle.load(open('lightcurve131.p','rb'))
    dtv171,lc171,aa=pickle.load(open('lightcurve171.p','rb'))
    dtv193,lc193,aa=pickle.load(open('lightcurve193.p','rb'))
    dtv211,lc211,aa=pickle.load(open('lightcurve211.p','rb'))
    dtv335,lc335,aa=pickle.load(open('lightcurve335.p','rb'))

    if tnorm: #normalize to this time
        sdate=dt.strftime(dtv94[0].date(),'%Y%m%d')
        dtnorm=dt.strptime(sdate+tnorm,'%Y%m%d%H:%M')
        n94idx=[idx for idx,xx in enumerate(dtv94) if xx == nearest_time(dtv94,dtnorm)][0]
        n131idx=[idx for idx,xx in enumerate(dtv131) if xx == nearest_time(dtv131,dtnorm)][0]
        n171idx=[idx for idx,xx in enumerate(dtv171) if xx == nearest_time(dtv171,dtnorm)][0]
        n193idx=[idx for idx,xx in enumerate(dtv193) if xx == nearest_time(dtv193,dtnorm)][0]
        n211idx=[idx for idx,xx in enumerate(dtv211) if xx == nearest_time(dtv211,dtnorm)][0]
        n335idx=[idx for idx,xx in enumerate(dtv335) if xx == nearest_time(dtv335,dtnorm)][0]
    else:
        n94idx=np.where(lc94 == 0)[0]
        n131idx=np.where(lc131 == 0)[0]
        n171idx=np.where(lc171 == 0)[0]
        n193idx=np.where(lc193 == 0)[0]
        n211idx=np.where(lc211 == 0)[0]
        n335idx=np.where(lc335 == 0)[0]

    print(np.max(lc211-lc211[n211idx])/np.max(lc211))
    #plot
    fig,ax=plt.subplots(figsize=(12,9))
    ax.plot_date(dtv94,(lc94-lc94[n94idx])/np.max(lc94-lc94[n94idx]),'b-',label='94 A')
    ax.plot_date(dtv131,(lc131-lc131[n131idx])/np.max(lc131-lc131[n131idx]),'c-',label='131 A')
    ax.plot_date(dtv171,(lc171-lc171[n171idx])/np.max(lc171-lc171[n171idx]),'k-',label='171 A')
    ax.plot_date(dtv193,(lc193-lc193[n193idx])/np.max(lc193-lc193[n193idx]),'g-',label='193 A')
    ax.plot_date(dtv211,(lc211-lc211[n211idx])/np.max(lc211-lc211[n211idx]),'m-',label='211 A')
    ax.plot_date(dtv335,(lc335-lc335[n335idx])/np.max(lc335-lc335[n335idx]),'r-',label='335 A')
    #format dates
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    #ax.set_title('Off-limb emssion, STEREO B 195 '+dt.strftime(dtvec[0],'%Y-%m-%dT%H:%M:%S')[:-9])
    ax.set_xlabel('Start Time (01-May-2013 01:30:00 UTC)')
    ax.set_ylabel('Normalized Flux')
    #if not xran:
    #    xran=[dtvecnums[0],dtvecnums2[-1]]
    #ax.set_xlim([tm_dt[22],tm_dt[-1]])
    ax.set_ylim(yran)
    ax.legend(loc=lpos)#,fontsize='small')
    #etc
    fig.show()
    #return tm_dt
