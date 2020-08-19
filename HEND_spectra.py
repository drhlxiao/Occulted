 #######################################
#HEND_spectra.py
# Erica Lastufka 25/09/2018

#Description: Display HEND spectra
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
from scipy import ndimage as ndi
from scipy import interpolate
import astropy.units as u
from astropy.coordinates import SkyCoord
import datetime
from  datetime import datetime as dt
from datetime import timedelta as td
from matplotlib.dates import date2num,num2date,DateFormatter

#global emcube
#global lgtaxis
#should store all this stuff in common block equivalents...how to deal with the map meta?

def spec_from_sav(filename='HEND/hend_010513.sav'):
    spec_dict=readsav(filename,python_dict=True)
    return spec_dict

def fix_time(t):
    dd=[t[-1],t[-2],t[-3]]
    dd.extend(t[:3])
    ft = datetime.datetime(*map(int,dd))
    return ft

def plot_spec(spec_dict,bgsub=False,bgshow=False,yran=False,eran=False,rhessi=False,norm=True,add=False,err=False):
    '''what it sounds like. '''
    #xylabels=['Log T='+str(np.round(lg,2)) for lg in lgtaxis]

    #time axis
    taxis=spec_dict['hend_sc_out_time_ex']
    dtaxis=[fix_time(t) for t in taxis]

    #energy bin labels
    elabels=np.transpose(spec_dict['hend_sc_out_ch_kev'])
    lab=[str(int(e[0])) + ' - ' + str(int(e[1])) + ' keV' for e in elabels]

    #spectra
    if bgsub:
        spec=np.transpose(spec_dict['hend_sc_out_ch_sub_bgrd_countrate'])
    else:
        spec=np.transpose(spec_dict['hend_sc_out_ch_countrate']) #shape is 212,16
    bg=np.transpose(spec_dict['hend_sc_out_ch_bgrd_countrate'])

    fig,ax=plt.subplots(figsize=(8,8))

    if not eran: #given in index numbers
        ivec = range(0,16)
    else:
        ivec = range(eran[0],eran[1])

    if not add:
        for i in ivec: #16 channels
            if norm:
                fac=np.max(spec[i])
            else:
                fac=1.
            p=ax.step(dtaxis,spec[i]/fac,label=lab[i],where='mid')
            if err:
                ax.errorbar(dtaxis,spec[i]/fac,yerr=1./np.sqrt(spec[i]),fmt='none',ecolor=p[0].get_color())
    else:
        if norm:
            fac=np.max(np.sum(spec[eran[0]:eran[1]+1],axis=0))
        else:
            fac=1
        lab='HEND '+str(int(elabels[eran[0]][0])) + ' - ' + str(int(elabels[eran[1]][1])) + ' keV'
        ydata=np.sum(spec[eran[0]:eran[1]+1],axis=0)

        ax.step(dtaxis,ydata/fac,label=lab,where='mid',linewidth=2)
        if err:
            ax.errorbar(dtaxis,ydata/fac,yerr=1./np.sqrt(ydata*fac),fmt='none',ecolor='b') #is this the correct way to calculate/treat error?

    if rhessi: #overplot RHESSI
        xrdict=readsav('Xray/plot_curves.sav',python_dict=True)
        hsi_time=xrdict['hsi_time']
        hsi_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in hsi_time]
        time_thermal=xrdict['time_thermal']
        if norm:
            fac12=np.max(xrdict['hsi12'])
            facth=np.max(xrdict['flux_thermal'])
        else:
            fac12=1.
            facth=1.
        hsi12=xrdict['hsi12']
        flux_thermal=xrdict['flux_thermal']
        thermal_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in time_thermal]
        ax.step(hsi_dt,hsi12/fac12,c='darkorange',label='RHESSI 10-30 keV',linewidth=2,where='mid')
        #ax.step(thermal_dt,flux_thermal/facth,'k-',label='RHESSI 4-8 keV',linewidth=2,where='mid')
        plt.fill_between(thermal_dt,flux_thermal/facth,0,alpha=0.5,color='grey',step='mid',label='RHESSI 4-8 keV')
        #if err:
        #    ax.errorbar(hsi_dt,hsi12/fac12,yerr=1./np.sqrt(hsi12))
        #    ax.errorbar(thermal_dt,flux_thermal/fac12,yerr=1./np.sqrt(flux_thermal))
        ax.set_xlim([hsi_dt[0],hsi_dt[-1]])
    #ticks every minute
    ax.xaxis.set_ticks([dt.strptime('2013-05-01 02:25:00','%Y-%m-%d %H:%M:%S') + td(minutes=i) for i in range(0,20)])
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    #ax.set_ylim([0,1.1])
    ax.set_xlabel('Time on 1 May 2013')
    ax.set_xlim([dt.strptime('2013-05-01 02:25:00','%Y-%m-%d %H:%M:%S'),dt.strptime('2013-05-01 02:40:00','%Y-%m-%d %H:%M:%S')])
    #ax.set_xlim([hsi_dt[27],hsi_dt[-14]])
    if norm:
        ax.set_ylabel('Normalized Counts/s')
    else:
        ax.set_ylabel('Counts/s')

    if yran:
        loc='upper left'
        ax.set_ylim(yran)
    else:
        loc='center right'
    ax.legend(loc=loc,fontsize=14)
    fig.show()
    #return ydata
    #return mapcube

def new_fig3(spec_dict,bgsub=True,eran=[3,9],rhessi=True,norm=True,add=True):
    os.chdir('Xray')
    goes=pickle.load(open('goes_lc.p','rb'))
    os.chdir('../')
    xrsa=goes.data['xrsa']
    xrsb=goes.data['xrsb']

    xrsa_t=xrsa.keys()
    xrsb_t=xrsb.keys()

    l=1905
    u=2335

    taxis=spec_dict['hend_sc_out_time_ex']
    dtaxis=[fix_time(t) for t in taxis]

    #energy bin labels
    elabels=np.transpose(spec_dict['hend_sc_out_ch_kev'])
    lab=[str(int(e[0])) + ' - ' + str(int(e[1])) + ' keV' for e in elabels]

    #spectra
    spec=np.transpose(spec_dict['hend_sc_out_ch_sub_bgrd_countrate'])

    ivec = range(eran[0],eran[1])

    if norm:
        fac=np.max(np.sum(spec[eran[0]:eran[1]+1],axis=0))
        lab='HEND '+str(int(elabels[eran[0]][0])) + ' - ' + str(int(elabels[eran[1]][1])) + ' keV'
        ydata=np.sum(spec[eran[0]:eran[1]+1],axis=0)


    xrdict=readsav('Xray/plot_curves.sav',python_dict=True)
    hsi_time=xrdict['hsi_time']
    hsi_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in hsi_time]
    time_thermal=xrdict['time_thermal']
    if norm:
        fac12=np.max(xrdict['hsi12'])
        facth=np.max(xrdict['flux_thermal'])

    hsi12=xrdict['hsi12']
    flux_thermal=xrdict['flux_thermal']
    thermal_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in time_thermal]


    fig,ax=plt.subplots(3,sharex=True,figsize=(7,14),gridspec_kw={'height_ratios':[1,1,2]}) #aspect ratio
    ax[0].plot_date(xrsb_t[l:u],xrsb[l:u], 'r-',linewidth=2,label='GOES 1-8$\AA$')
    ax[0].plot_date(xrsa_t[l:u],xrsa[l:u], 'k-',linewidth=2,label='GOES .5-4$\AA$')
    ax[0].axvline(x=xrsb_t[2124],linestyle='dashed',color='k')
    ax[0].set_yscale('log')
    ax[0].set_ylabel('GOES Flux (W/m$^2$)')
    ax[0].legend(loc='center right',fontsize=12)
    #ax[1].plot_date(xrsa_t[l:u],xrsa[l:u],'k-') #do I want these to have a similar y-axis range for comparison?
    ax[1].plot_date(xrsb_t[l:u],xrsb[l:u]*1e7,'r-',linewidth=2)
    ax[1].axvline(x=xrsb_t[2124],linestyle='dashed',color='k')

#     ax[1].xaxis.set_ticks([dt.strptime('2013-05-01 02:25:00','%Y-%m-%d %H:%M:%S') + td(minutes=i) for i in range(0,16)])
#     myFmt = DateFormatter('%H:%M')
#     ax[1].xaxis.set_major_formatter(myFmt)
#     plt.gcf().autofmt_xdate()
#    ax[1].set_xlabel('Start Time (01-May-2013 01:30:00 UTC)')
    #ax[2].ticklabel_format(style='sci',axis='y')
    #ax[2].ticklabel_format(style='sci',axis='y')
    #ax[1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax[1].set_ylabel('GOES 1-8$\AA$ Flux x 10$^{-7}$ (W/m$^2$)')

    ax[0].set_ylim([1e-9,1.2e-6])
    #ax[1].set_ylim([8e-7,1.1e-6])

    ax[2].step(dtaxis,ydata/fac,label=lab,where='mid',linewidth=2)
    ax[2].step(hsi_dt,hsi12/fac12,c='darkorange',label='RHESSI 10-30 keV',linewidth=2,where='mid')
    plt.fill_between(thermal_dt,flux_thermal/facth,0,alpha=0.5,color='grey',step='mid',label='RHESSI 4-8 keV')
    ax[2].set_xlim([hsi_dt[0],hsi_dt[-1]])
    #ticks every minute
    ax[2].xaxis.set_ticks([dt.strptime('2013-05-01 02:25:00','%Y-%m-%d %H:%M:%S') + td(minutes=i) for i in range(0,20)])
    myFmt = DateFormatter('%H:%M')
    ax[2].xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    #ax.set_ylim([0,1.1])
    ax[2].set_xlabel('Time on 1 May 2013')
    ax[2].set_xlim([dt.strptime('2013-05-01 02:25:00','%Y-%m-%d %H:%M:%S'),dt.strptime('2013-05-01 02:40:00','%Y-%m-%d %H:%M:%S')])
    #ax.set_xlim([hsi_dt[27],hsi_dt[-14]])
    if norm:
        ax[2].set_ylabel('Normalized Counts/s')

    loc='upper left'
    ax[2].set_ylim([0,1.1])
    ax[2].legend(loc=loc,fontsize=12)
    plt.subplots_adjust(top=.95,bottom=.1,right=.95,hspace=.05)

    plt.savefig('plots/newfig3b.png')
    fig.show()


def plot_count_spec(spec_dict,bgsub=False,bgshow=False,yran=False,tran=False,rhessi=False,norm=True,add=False,err=False):
    '''what it sounds like. '''
    #xylabels=['Log T='+str(np.round(lg,2)) for lg in lgtaxis]

    #time axis
    taxis=spec_dict['hend_sc_out_time_ex']
    dtaxis=[fix_time(t) for t in taxis]

    #energy bin labels
    elabels=np.transpose(spec_dict['hend_sc_out_ch_kev'])
    lab=dt.strftime(dtaxis[tran[0]],'%H:%M:%S') + ' - ' + dt.strftime(dtaxis[tran[1]],'%H:%M:%S')
    eaxis=[e[0] for e in elabels]
    #eaxis.append(elabels[-1][1])

    #spectra
    if bgsub:
        spec=spec_dict['hend_sc_out_ch_sub_bgrd_countrate']
    else:
        spec=spec_dict['hend_sc_out_ch_countrate'] #shape is 212,16
    bg=spec_dict['hend_sc_out_ch_bgrd_countrate']

    fig,ax=plt.subplots(figsize=(8,10))

#     if not tran: #given in index numbers - peak is between 113 and 130. Max bin is 123.
#         ivec = range(0,16)
#     else:
#         ivec = range(eran[0],eran[1])
    ydata=np.sum(spec[tran[0]:tran[1]],axis=0)
    p=ax.step(eaxis,ydata,where='mid',label=lab)
    if err:
        ax.errorbar(eaxis,ydata, yerr=1./np.sqrt(ydata),fmt='none',ecolor=p[0].get_color())

    #now interpolate to smooth curve and plot with error ....
    ff=interpolate.interp1d(eaxis,ydata)
    xnew=np.arange(30,1014)
    ynew=ff(xnew)
    p=ax.step(xnew,ynew,where='mid')
    if err:
        ax.errorbar(xnew,ynew, yerr=1./np.sqrt(ynew),fmt='none',ecolor=p[0].get_color())




#     if rhessi: #overplot RHESSI
#         xrdict=readsav('Xray/plot_curves.sav',python_dict=True)
#         hsi_time=xrdict['hsi_time']
#         hsi_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in hsi_time]
#         time_thermal=xrdict['time_thermal']
#         if norm:
#             fac12=np.max(xrdict['hsi12'])
#             facth=np.max(xrdict['flux_thermal'])
#         else:
#             fac12=1.
#             facth=1.
#         hsi12=xrdict['hsi12']
#         flux_thermal=xrdict['flux_thermal']
#         thermal_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in time_thermal]
#         ax.step(hsi_dt,hsi12/fac12,c='darkorange',label='RHESSI 10-30 keV',linewidth=2,where='mid')
#         #ax.step(thermal_dt,flux_thermal/facth,'k-',label='RHESSI 4-8 keV',linewidth=2,where='mid')
#         #if err:
#         #    ax.errorbar(hsi_d,hsi12/fac12,yerr=1./np.sqrt(hsi12))
#         #    ax.errorbar(thermal_dt,flux_thermal/fac12,yerr=1./np.sqrt(flux_thermal))
#         ax.set_xlim([hsi_dt[0],hsi_dt[-1]])
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    ax.set_xlabel('Energy (keV)')
    ax.set_xlim([30,2000])
    ax.set_yscale('log')
    ax.set_xscale('log')
    if norm:
        ax.set_ylabel('Normalized Counts/s')
    else:
        ax.set_ylabel('Counts/s')

    if yran:
        loc='upper right'
        ax.set_ylim(yran)
    else:
        loc='center right'
        ax.set_ylim([0.1,1000])
    ax.legend(loc=loc,fontsize=10)
    fig.show()

def gamma_from_txt(f='../May1flare/HEND/ch_results.txt'):
    ''' get peak spectrum from text file.
    1 col: date
    2 col: time (at Earth)
    3 col: HXR spectral normalization factor (at 100 keV)
    4 col: error of HXR spectral normalization factor
    5 col: HXR spectral index
    6 col: HXR spectral index error
    '''
    from collections import OrderedDict
    hdict=OrderedDict()
    hdict['datetime']=[]
    hdict['norm']=[]
    hdict['norm_err']=[]
    hdict['gamma']=[]
    hdict['gamma_err']=[]
    with open(f,'r') as lines:
        txt=lines.readlines()

    for l in txt:
        line=l.replace('\t',' ')
        foo=line[:-1].split(' ')
        #combine the first 2 into datetime
        hdt=dt.strptime(foo[0]+'T'+foo[1],'%Y-%m-%dT%H:%M:%S')
        newfoo=[hdt]+foo[2:]
        for k,v in zip(hdict.keys(), newfoo):
            hdict[k].append(v)

    #convert numbers to float
    fnorm=[float(n) for n in hdict['norm']]
    fnorm_err=[float(n) for n in hdict['norm_err']]
    fgamma=[float(n) for n in hdict['gamma']]
    fgamma_err=[float(n) for n in hdict['gamma_err']]
    hdict['norm']=fnorm
    hdict['norm_err']=fnorm_err
    hdict['gamma']=fgamma
    hdict['gamma_err']=fgamma_err

    return hdict

def gamma_from_sav():
    from scipy.io import readsav
    #read all the relavent sav files and extract the gamma
    savs=glob.glob('fit_results*.sav')
    dlist=[]
    for s in savs:
        dlist.append(readsav(s, python_dict=True))
    return dlist

def spec_from_txt(f='spec_023206.txt'):
    from collections import OrderedDict
    hdict=OrderedDict()
    hdict['Energy']=[]
    hdict['Flux']=[]
    hdict['Error']=[]
    hdict['Threshold']=[]
    with open(f,'r') as lines:
        txt=lines.readlines()

    for l in txt[1:]:
        line=l.replace('\t',' ')
        foo=line[:-1].split(' ')
        for k,v in zip(hdict.keys(), foo):
            hdict[k].append(v)

    #convert numbers to float
    fenergy=[float(n) for n in hdict['Energy']]
    fflux=[float(n) for n in hdict['Flux']]
    ferror=[float(n) for n in hdict['Error']]
    fthresh=[float(n) for n in hdict['Threshold']]
    hdict['Energy']=fenergy
    hdict['Flux']=fflux
    hdict['Error']=ferror
    hdict['Threshold']=fthresh

    return hdict

def hend_and_rhessi_spec():
    shend=spec_from_txt()
    henergy=shend['Energy']
    hflux=shend['Flux']
    herr=shend['Error']
    thresh=shend['Threshold']
    from scipy.io import readsav
    srhessi=readsav('../Xray/fit_results_HEND_peak.sav',python_dict=True)
    pf=readsav('../Xray/peak_flux2.sav',python_dict=True)
    fitvals=readsav('../Xray/fit_vals_test.sav',python_dict=True)
    #fv=fitvals['foo']
    pflux=pf['ph_data']
    pflux[50]=np.mean([pflux[49],pflux[51]])
    rdict=srhessi['s']
    renergy=rdict.spex_summ_energy
    revec=[r[0] for r in renergy[0]]
    rct=rdict.spex_summ_ct_rate[0] # if it's count rate, I can bin it right? photons per cm^3 per keV

    rbg=rdict.spex_summ_bk_rate
    rcterr=rdict.spex_summ_ct_error
    rconv=rdict.spex_summ_conv
    #rarea=rdict.spex_summ_area
    rflux=rconv*(rct-rbg)
    rferr=rconv*rcterr #is this correct?
    rpflux=list(pflux)[:55]#*rconv[0]
    revec2=revec[:55]
    #rebin high energy counts - remember to divide by bin size
    for i in [10,15,20,25,30,35,40,45]:
        j=i-15
        rpflux.append(np.mean(pflux[55+j:55+i])/5.)
        revec2.append(revec[55+i])

    #fits
    ct_energy=fitvals['ct_energy']
    thfits=fitvals['vth']#fv.yvals[0][1,:]
    plawfits=fitvals['bpow']#fv.yvals[0][2,:]
    #background
    bgvals=readsav('../Xray/bk_photon_flux2.sav',python_dict=True)
    rbk=bgvals['bk_data']

    #sanity check plot
    #fig,ax=plt.subplots()
    #ax.step(revec,rct[0])
    #ax.step(revec,rbg[0])
    #ax.step(revec,rflux[0])
    #ax.plot(revec,rpflux)
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    #fig.show()

    #plot together
    fig,ax=plt.subplots()
    ax.plot(henergy,hflux,color='m',linewidth=2)
    ax.plot(henergy,thresh,color='m',alpha=0.5)
    ax.scatter(henergy,thresh,color='m',label='HEND background \n 02:32:04 - 02:32:24',alpha=0.5)
    ax.errorbar(henergy,hflux,herr,fmt='m.', label='HEND 02:32:04 - 02:32:24')
    #ax.step(revec,rflux[0],label='RHESSI',color='c') #max should be ~ 10^3
    #ax.errorbar(revec,rflux[0],rferr[0],fmt='c.')
    ax.step(revec2[10:-6],rpflux[10:-6],label='RHESSI 02:32:04-02:32:24',color='black',linewidth=4) #max should be ~ 10^3
    ax.step(revec[10:-30],rbk[10:-30],label='RHESSI background \n 02:20:00 - 02:20:30',color='gray',alpha=0.7,linewidth=2) #max should be ~ 10^3
    #ax.scatter(revec2[10:-13],rpflux[10:-13],color='c') #max should be ~ 10^3
    #ax.scatter(revec[10:-30],rbk[10:-30],color='k') #max should be ~ 10^3
    #ax.plot(revec,tfits,label='fit',color='r',linewidth=2) #max should be ~ 10^3
    ax.plot(ct_energy,thfits,label='vth fit',color='r',linestyle='--',linewidth=2) #max should be ~ 10^3
    ax.plot(ct_energy,plawfits,label='bpow fit',color='y',linestyle='--',linewidth=2) #max should be ~ 10^3
    #ax.text(10,30,'02:')
    #ax.errorbar(revec,rflux[0],rferr[0],fmt='c.')
    ax.set_xlabel('Energy (keV)')
    ax.set_xlim([5,1000])
    ax.set_ylim([.0001,100])
    ax.set_ylabel('Photon Flux per keV') #NOT flux (counts/s per cm^3). HOw to convert to flux? using the photon model?
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='upper right',fontsize='medium')
    fig.show()

#     def plot_rhessi_zoom(xsi_rhessi='plot_curves.sav',AIA='../AIA/high_cadence_cutout/lightcurve94nobg.p'):
#     from scipy.io import readsav
#     xrdict=readsav(xsi_rhessi,python_dict=True)
#     hsi_time=xrdict['hsi_time']
#     hsi12=xrdict['hsi12']/np.max(xrdict['hsi12'])
#     time_thermal=xrdict['time_thermal']
#     flux_thermal=xrdict['flux_thermal']/np.max(xrdict['flux_thermal'])
#     dtv94,lc94,max94=pickle.load(open(AIA,'rb'))
#     #mask bad value
#     lc94[85]=np.mean([lc94[84],lc94[86]])
#     os.chdir('../NORH')
#     dtvNORH,lcNORH,aa=pickle.load(open('norh_lc.p','rb'))
#     os.chdir('../')
#     #assume IDL anytim axes were converted via anytim(time,/vms) into format %d-%b%-Y %H:%M:%S
#     hsi_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in hsi_time]
#     thermal_dt=[dt.strptime(ht.lstrip()[:-4], '%d-%b-%Y %H:%M:%S') for ht in time_thermal]
#     #plot
#     fig,ax=plt.subplots(figsize=(6,6))
#     #ax.plot_date(sB_time,sB_curve,'g.-',label='STEREO-B 195')
#     #ax.plot_date(sA_time,sA_curve,'m.-',label='STEREO-A 195')
#     #ax.plot_date(hsi_dt,hsi12,'b.-',label='RHESSI 10-30 keV')
#     ax.step(hsi_dt,hsi12,'b-',label='RHESSI 10-30 keV',linewidth=2)
#     ax.plot_date(thermal_dt,flux_thermal,'k-',label='RHESSI 4-8 keV',linewidth=2)
#     ax.plot_date(dtv94,(lc94-lc94[70])/(max94-lc94[70]),'g.-',label='AIA 94 $\AA$')
#     ax.plot_date(dtvNORH,(lcNORH-lcNORH[15])/np.max(lcNORH-lcNORH[15]),'m.-',label='NORH 17 GHz')
#     #format dates
#     myFmt = DateFormatter('%H:%M')
#     ax.xaxis.set_major_formatter(myFmt)
#     plt.gcf().autofmt_xdate()
#     #ax.set_title('Off-limb emssion, STEREO B 195 '+dt.strftime(dtvec[0],'%Y-%m-%dT%H:%M:%S')[:-9])
#     ax.set_xlabel('Start Time (01-May-2013 01:30:00 UTC)')
#     ax.set_ylabel('Normalized Flux')
#     #if not xran:
#     #    xran=[dtvecnums[0],dtvecnums2[-1]]
#     ax.set_xlim([thermal_dt[0],dtvNORH[-15]])
#     ax.set_ylim([-.5,1.1])
#     ax.legend(loc='lower right')#,fontsize='small')
#     #etc
#     fig.show()

