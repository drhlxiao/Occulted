 #######################################
# messenger_angles.py
# Erica Lastufka 30/11/2016  

#Description: Check if Messenger can see an occulted flare
#######################################

#######################################
# Usage:

# for default output: python messenger_angles.py
######################################

import numpy as np
import scipy.constants as sc
import pandas as pd
from datetime import datetime
from datetime import timedelta as td
import os
import data_management as da
import pickle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import copy

def import_flare_list(instrument,filename='test.csv', ):
    '''reads in list of occulted flare positions and time intervals'''
    data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text

    dt = []
    for d,t in zip(data['Date'].tolist(),data['Start_Time'].tolist()):
        dt.append(d + ' '+ t)
    #print dt[0:5]
            
    flare_list = da.Data_Study(length=len(data['Date']))
    if instrument == 'Messenger':
        flare_list.Datetimes["Messenger_datetimes"] = dt
        flare_list.Flare_properties["Messenger_T"] = data['T1'].tolist()
        flare_list.Flare_properties["Messenger_EM1"] = data['EM1'].tolist()
    else:
        flare_list.Datetimes["RHESSI_datetimes"] = dt
        flare_list.Flare_properties['RHESSI_GOES'] = data['GOES'].tolist()

    flare_list.format_datetimes()    
    if 'Angle' in data.keys():
        flare_list.Angle=data['Angle'].tolist()
    if 'ID' in data.keys():
        flare_list.ID=data['ID'].tolist()
        
    return flare_list

def is_flare_observed(obj1,obj2, cutoff=7200):
    '''Checks if Messenger observed anything within the flare time intervals. Use Brian Dennis's flare list. They shouldn't be too different for an occulted flare? Cutoff time should be specified in seconds, default is 2 hours'''

    dt1=obj1.Datetimes['RHESSI_datetimes']
    dt2=obj2.Datetimes['Messenger_datetimes']
    id1=obj1.ID
    Rgoes=obj1.Flare_properties['RHESSI_GOES']
    
    IDs,newfoundlist, nMtime,nRtime,nMdt,nRdt=[],[],[],[],[],[]
    for dt,id in zip(dt1,id1):
        closest = map(lambda d: abs(d-dt), dt2) #this is a list
        for item,index in zip(closest, range(0,len(dt2))):
            if item.total_seconds() < cutoff:
                print dt, '    ',dt2[index]
                newfoundlist.append(str(dt.date()))
                #nMtime.append(str(dt2[index].time()))
                #nRtime.append(str(dt.time()))
                nMdt.append(dt2[index])
                nRdt.append(dt)
                IDs.append(id)

    new_list = da.Data_Study(length = len(IDs))
    new_list.ID = IDs
    new_list.Datetimes['RHESSI_datetimes'] = nRdt
    new_list.Datetimes['Messenger_datetimes'] = nMdt
    new_list.Flare_properties['RHESSI_GOES'] = Rgoes
    new_list.Flare_properties["Messenger_T"] = obj2.Flare_properties["Messenger_T"]
    new_list.Flare_properties["Messenger_EM1"] = obj2.Flare_properties["Messenger_EM1"]
    
    return new_list
       
def plot_angle_distribution(list_obj,ymax=10):
    '''make a histogram of the angle distributions'''
     
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    n, bins, patches = plt.hist(list_obj.Angle, 18, facecolor='green', alpha=0.75)
    
    plt.xlabel('Angle between Earth and Mercury (degrees)')
    plt.ylabel('Number of Events')
    ax1.set_ylim([0,ymax])
    ax1.set_xlim([0,180])

    ax1.plot()

    fig.show()

def convert_goes2flux(goes_class):
    '''Converts Goes class to flux value, for either a single value or list of values'''
    flux = -1
    if type(goes_class) == list:
        flux=[]
        for item in goes_class:
            try:
                val=item[0:1]
                if item.endswith('*'):
                    item = item[:-1]
                if val == 'A':
                    flux.append(float(item[1:])*10**-8)
                if val == 'B':
                    flux.append(float(item[1:])*10**-7)
                if val == 'C':
                    flux.append(float(item[1:])*10**-6)
                if val == 'M':
                    flux.append(float(item[1:])*10**-5)
                if val == 'X':
                    flux.append(float(item[1:])*10**-4)  
            except TypeError:
                pass
    else:
        try:
            val = goes_class[0:1]
            if goes_class.endswith('*'):
                goes_class = goes_class[:-1]
            if val == 'A':
                flux = float(goes_class[1:])*10**-8
            if val == 'B':
                flux = float(goes_class[1:])*10**-7
            if val == 'C':
                flux = float(goes_class[1:])*10**-6
            if val == 'M':
                flux = float(goes_class[1:])*10**-5
            if val == 'X':
                flux = float(goes_class[1:])*10**-4   
        except TypeError:
            pass
    return flux

def get_ratio(flare_list):
    '''Get ratio of Messenger goes v actual goes'''
    mvals,rvals=[],[]
    mc = flare_list.Flare_properties["Messenger_GOES"]
    gc=flare_list.Flare_properties["GOES_GOES"]
    for mval,rval in zip(mc,gc):
        if type(rval) != float and len(rval) > 3:#:type(rval) != float : #or np.isnan(rval) == False:          
            rf=convert_goes2flux(rval)
            mf=convert_goes2flux(mval)
        else:
            rf=-1
            mf=0
        mvals.append(mf)
        rvals.append(rf)
    ratio = np.array(mvals)/np.array(rvals)
    return ratio
   
    
def plot_goes_ratio(list_obj, title= "",ymin=0, ymax=3, labels=1,ylog=False, goes=False,mgoes=False, scatter = False,cc='GOES',save=False,show=True):
    '''make a plot of the GOEs ratio vs. angle'''
    ang,mvals,rvals,delta,colors,coordlabel,labelang,labelratio=[],[],[],[],[],[],[],[]
    gc = list_obj.Flare_properties["RHESSI_GOES"]
    mc=list_obj.Flare_properties["Messenger_GOES"]
    ylabel='Messenger_GOES/RHESSI_GOES'
    if goes:
            mc = list_obj.Flare_properties["GOES_GOES"]
            gc=list_obj.Flare_properties["RHESSI_GOES"]
            ylabel='Observed GOES/RHESSI_GOES'
    if mgoes:
            mc = list_obj.Flare_properties["Messenger_GOES"]
            gc=list_obj.Flare_properties["GOES_GOES"]
            ylabel='Messenger_GOES/Observed_GOES'
    for angle, mval,rval,chisq,ids,dts in zip(list_obj.Angle,mc,gc,list_obj.Notes,list_obj.ID,list_obj.Datetimes['Obs_start_time']):
        try:
            rval=float(rval)
            mval=float(mval)
            #ang.append(angle)
            if cc=='GOES':
                #print np.rint(-np.log10(rval))
                if np.rint(-np.log10(rval)) <= 4.0:colors.append('r')
                elif np.rint(-np.log10(rval)) == 5.0:colors.append('m')
                elif np.rint(-np.log10(rval)) == 6.0:colors.append('y')
                elif np.rint(-np.log10(rval)) == 7.0:colors.append('g')
                elif np.rint(-np.log10(rval)) == 8.0:colors.append('b')
                elif np.rint(-np.log10(rval)) == 9.0:colors.append('k')
            else: #color code by other one
                if np.rint(-np.log10(mval)) <= 4.0:colors.append('r')
                elif np.rint(-np.log10(mval)) == 5.0:colors.append('m')
                elif np.rint(-np.log10(mval)) == 6.0:colors.append('y')
                elif np.rint(-np.log10(mval)) == 7.0:colors.append('g')
                elif np.rint(-np.log10(mval)) == 8.0:colors.append('b')
                elif np.rint(-np.log10(mval)) == 9.0:colors.append('k')
            rf=rval                
            mf=mval
        except ValueError:
            if type(rval) != float and len(rval) > 3:#:type(rval) != float : #or np.isnan(rval) == False:
                #ang.append(angle)
                if cc=='GOES':
                    if rval.startswith('X'): colors.append('r')
                    elif rval.startswith('M'):colors.append('m')
                    elif rval.startswith('C'):colors.append('y')
                    elif rval.startswith('B'):colors.append('g')
                    elif rval.startswith('A'):colors.append('b')
                    else: colors.append('k')
                else: #color code by other one
                    if mval.startswith('X'): colors.append('r')
                    elif mval.startswith('M'):colors.append('m')
                    elif mval.startswith('C'):colors.append('y')
                    elif mval.startswith('B'):colors.append('g')
                    elif mval.startswith('A'):colors.append('b')
                    else: colors.append('k')                
            rf=convert_goes2flux(rval)
            mf=convert_goes2flux(mval)

        mvals.append(mf)
        rvals.append(rf)
        ang.append(angle)
        if rf != -1 and rf !=0.:
            labelang.append(angle)
            labelratio.append(mf/rf)
            if labels==1:
                coordlabel.append(ids)
            else: #0 or 2
                coordlabel.append(datetime.strftime(dts,'%D %H:%M'))
        if scatter:
              delta = 50
        elif chisq == '':#notes column is empty
            delta.append(5000*10*2**np.rint(np.log10(np.abs(rf-mf)))) #difference in size between the GOES classes
        else: #notes carries chisq value
            delta.append(50*10*2**np.rint(np.log10(float(chisq))))

    ratio = np.array(mvals)/np.array(rvals)
    full_ratio = ratio #for now ...
    #print list_obj.Flare_properties["RHESSI_GOES"]
    #print list_obj.Notes,delta
    #print colors
    #print sorted(ratio)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(np.array(ang), ratio, s=delta, c=colors,alpha=.75)
    ax1.axhline(y=1,linestyle='dashed',color='k')

    if labels != 0: 
        for x,y,t in zip(np.array(labelang),labelratio,coordlabel):
            #print x,y,t
            ax1.annotate('%s' % t, xy=(x,y), textcoords='data')
        plt.grid()

   
    plt.xlabel('Angle between Earth and Mercury (degrees)')
    plt.ylabel(ylabel)

    if ylog:
        ax1.set_yscale('log')
    ax1.set_ylim([ymin,ymax])
    ax1.set_xlim([0,180])
    plt.title(title)
    X = mpatches.Patch(color='red', label='X')
    M = mpatches.Patch(color='magenta', label='M')
    C = mpatches.Patch(color='yellow', label='C')
    B = mpatches.Patch(color='green', label='B')
    A= mpatches.Patch(color='blue', label='A')
    K= mpatches.Patch(color='black', label='>A')
    ax1.legend(handles=[X,M,C,B,A,K],loc='upper left',fontsize='medium')

    ax1.plot()
    if save:
        plt.savefig(save)
    if show:
        fig.show()
    return ratio,full_ratio

def hist_int_time_ratio(flare_list):
    '''Make histogram to see if ratio depends on integration time'''
    foo=foo
    #do stuff with the datetimes
    return ratio

def hist_ratio(flare_list,title='All flares',gc='all', save=False,show=True):
    '''Make histogram to see distribution of ratios by goes class'''
    ratio=get_ratio(flare_list)
    ratio = ratio[np.where(ratio != -1)]
    #print ratio
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    goes=np.array(convert_goes2flux(flare_list.Flare_properties['GOES_GOES']))
    #print goes
    
    if gc == 'all':
        n, bins, patches = plt.hist(ratio, np.logspace(-3,3,12), facecolor='orange', alpha=0.75)
    if gc == 'A':
        gt,lt =np.where(goes > 10**-8),np.where(goes < 10**-7)
        all=set(gt[0]) & set(lt[0])
        ratio= ratio[list(all)]
        n, bins, patches = plt.hist(ratio, np.logspace(-3,3,12), facecolor='blue', alpha=0.75)
    if gc == 'B':
        gt,lt =np.where(goes > 10**-7),np.where(goes < 10**-6)
        all=set(gt[0]) & set(lt[0])
        ratio= ratio[list(all)]
        n, bins, patches = plt.hist(ratio,np.logspace(-3,3,12) , facecolor='green', alpha=0.75)
    if gc == 'C':
        gt,lt =np.where(goes > 10**-6),np.where(goes < 10**-5)
        all=set(gt[0]) & set(lt[0])
        ratio= ratio[list(all)]
        n, bins, patches = plt.hist(ratio, np.logspace(-3,3,12), facecolor='yellow', alpha=0.75)
    if gc == 'M':
        gt,lt =np.where(goes > 10**-5),np.where(goes < 10**-4)
        all=set(gt[0]) & set(lt[0])
        ratio= ratio[list(all)]
        n, bins, patches = plt.hist(ratio, np.logspace(-3,3,12), facecolor='magenta', alpha=0.75)
    if gc == 'overlay':
        gt,lt =np.where(goes > 10**-8),np.where(goes < 10**-7)
        all=set(gt[0]) & set(lt[0])
        ratioA= ratio[list(all)]
        plt.hist(ratioA, np.logspace(-3,3,12), facecolor='blue', alpha=0.6)
        gt,lt =np.where(goes > 10**-7),np.where(goes < 10**-6)
        all=set(gt[0]) & set(lt[0])
        ratioB= ratio[list(all)]
        plt.hist(ratioB,np.logspace(-3,3,12) , facecolor='green', alpha=0.6)
        gt,lt =np.where(goes > 10**-6),np.where(goes < 10**-5)
        all=set(gt[0]) & set(lt[0])
        ratioC= ratio[list(all)]
        plt.hist(ratioC, np.logspace(-3,3,12), facecolor='yellow', alpha=0.6)
        gt,lt =np.where(goes > 10**-5),np.where(goes < 10**-4)
        all=set(gt[0]) & set(lt[0])
        ratioM= ratio[list(all)]
        plt.hist(ratioM, np.logspace(-3,3,12), facecolor='magenta', alpha=0.6)
        
        M = mpatches.Patch(color='magenta', label='M')
        C = mpatches.Patch(color='yellow', label='C')
        B = mpatches.Patch(color='green', label='B')
        A= mpatches.Patch(color='blue', label='A')
        ax1.legend(handles=[M,C,B,A],loc='upper left',fontsize='medium')
        
    
    plt.xlabel('Ratio between Messenger GOES and actual GOES')
    plt.ylabel('Number of Events')
    plt.gca().set_xscale("log")
    plt.title(title)
    ax1.set_ylim([0,100])
    #ax1.set_xlim([0,150])

    ax1.plot()
    if save:
        plt.savefig(save)
    if show:
        fig.show()
        
    return goes

def select_outliers_lt90(flare_list, ratio,threshold=10.):
    '''Select certain flares from the list to be examined in visually'''
    indices=[]
    newlist=copy.deepcopy(flare_list)
    for i,flare in enumerate(flare_list.ID):
        if np.abs(ratio[i]) > threshold and flare_list.Angle[i] < 90.:
            indices.append(i)
    newlist.slice(indices)    
    return newlist

def select_outliers_gt90(flare_list, ratio,threshold=10.):
    '''Select certain flares from the list to be examined in visually'''
    indices=[]
    newlist=copy.deepcopy(flare_list)
    for i,flare in enumerate(flare_list.ID):
        if np.abs(ratio[i]) > threshold and flare_list.Angle[i] > 90.:
            indices.append(i)
    newlist.slice(indices)    
    return newlist

def select_outliers(flare_list, ratio, angle=90.,threshold=10.): #greater than given angle... need to fix bug/feature with ratio
    '''Select certain flares from the list to be examined in visually'''
    indices=[]
    newlist=copy.deepcopy(flare_list)
    for i,flare in enumerate(flare_list.ID):
        if np.abs(ratio[i]) > threshold and flare_list.Angle[i] < angle:
            indices.append(i)
    newlist.slice(indices)    
    return newlist
    
def filter_source_vis(flare_list, field='Y'):
    '''Select certain flares from the list to be examined in sswidl'''
    vis = flare_list.Flare_properties['source_vis']
    indices=[]
    for i,v in enumerate(vis):
        if v == field:
            indices.append(i)
    flare_list.slice(indices)    
    return flare_list

def two_thermals(flare_list_high,flare_list_low):
    flare_list_2T=copy.deepcopy(flare_list_high)
    fh=flare_list_high.Flare_properties
    fl=flare_list_low.Flare_properties
    for i,t1,t2,em1,em2 in zip(range(0,len(flare_list_low.ID)-1),fh['Messenger_T'],fl['Messenger_T'],fh['Messenger_EM1'],fl['Messenger_EM1']):
        if (type(t2) and type(t1)) == str:
            t2=float(t2)
            t1=float(t1)
        if np.isnan(t2) != 1 and t2 != -1: #if there was no solution, ignore the low-T component (this occurs more because of the EM though, and only in the goes_fluxes function)
            flare_list_2T.Flare_properties['Messenger_T'][i]=t1+t2
        else:
            flare_list_2T.Flare_properties['Messenger_T'][i]=em1
        if (type(em2) and type(em1)) == str:
            em2=float(em2)
            em1=float(em1)
        if np.isnan(em2) != 1 and em2 != -1:        
            flare_list_2T.Flare_properties['Messenger_EM1'][i]=em1+em2
        else:
            flare_list_2T.Flare_properties['Messenger_EM1'][i]=em1

    return flare_list_2T #will need to go to IDL to calculate the goes flux... unless I translate that to Python but I'm lazy

def download_messenger(flare_list):
    '''Downloads Messenger .dat and .lbl files from the database, given event date. Can be used to get missing files too.'''
    import urllib
    dataurl = 'https://hesperia.gsfc.nasa.gov/messenger/'
    #subfolders by year, month,day (except 2001)
    listlen = len(flare_list.ID)
    for i,dt in zip(range(0,listlen -1),flare_list.Datetimes['Messenger_datetimes']):
        datestring=dt.strftime('%Y%j')
        newurl=dataurl+datestring[0:4] #just the year
        filename='/xrs'+datestring
        datfile=filename+'.dat'
        lblfile=filename+'.lbl'
        #first check if the file is already there:
        if not os.path.exists('/Users/wheatley/Documents/Solar/occulted_flares/data/dat_files/'+filename+'.dat'):
            urllib.urlretrieve(newurl+datfile,'/Users/wheatley/Documents/Solar/occulted_flares/data/dat_files/'+filename+'.dat')
        if not os.path.exists('/Users/wheatley/Documents/Solar/occulted_flares/data/dat_files/'+filename+'.lbl'):
            urllib.urlretrieve(newurl+lblfile,'/Users/wheatley/Documents/Solar/occulted_flares/data/dat_files/'+filename+'.lbl')
        #fill in the xrsfilename section
        flare_list.Data_properties['XRS_files'][i]=filename
    return flare_list
        
def do_quicklook(flare_list, download=False, opent=False):
    import webbrowser
    import urllib
    dataurl = 'http://soleil.i4ds.ch/hessidata/metadata/qlook_image_plot'
    #subfolders by year, month,day (except 2001)
    listlen = len(flare_list.ID)
    for i,id,dt in zip(range(0,listlen -1),flare_list.ID,flare_list.Datetimes['RHESSI_datetimes']):
        newurl=dataurl+'/'+dt.strftime('%Y/%m/%d')
        filename='/hsi_qlimg_'+str(id)+'_012025.png' #get the 12-25 keV image
        image=newurl+filename
        if download:
            urllib.urlretrieve(image,'/Users/wheatley/Documents/Solar/occulted_flares/data/round2/'+filename)
        if opent:
            webbrowser.open_new_tab(image)
        else: #just save the data to the object
            flare_list.Data_properties['QL_images'][i] = image 

def open_in_RHESSI_browser(flare_list, opent=False):
    import webbrowser
    browserurl = 'http://sprg.ssl.berkeley.edu/~tohban/browser/?show=grth1+qlpcr+qlpg9+qli02+qli03+qli04+synff&date=' #20120917&time=061500'
    #subfolders by year, month,day (except 2001)
    #warn if you're opening more than 20 tabs
    if opent == True and len(flare_list.ID) > 20:
        ans = raw_input('You are about to open ' + str(len(flare_list.ID)) + ' tabs! Are you sure?')
        if ans != 'Y':
            opent = False        
    for i,dt in enumerate(flare_list.Datetimes['RHESSI_datetimes']):
        try:
            address=browserurl + dt.strftime('%Y%m%d') + '&time=' +dt.strftime('%H%M%S')
        except AttributeError: #datetime for RHESSI is empty
            address=browserurl + flare_list.Datetimes['Messenger_datetimes'][i].strftime('%Y%m%d') + '&time=' +flare_list.Datetimes['Messenger_datetimes'][i].strftime('%H%M%S')            
        flare_list.Data_properties['RHESSI_browser_urls'][i] = address
        if opent:
            webbrowser.open_new_tab(address)

def select_time_intervals(flare_list,before=10,after=10):
    '''Generate time intervals for use in sswidl'''
    #list=pickle.load(open('list_observed.p','rb'))
    Mstarttime,Mendtime,Rstarttime,Rendtime=[],[],[],[]
    #Format for time is YY/MM/DD, HHMM:SS.SSS
    #or alternatively, DD-Mon-YY HH:MM:SS.SSS -> '%d-%b-%Y %H:%M:%S.000'

    for index in range(0,len(flare_list['date'])):
        Mstarttime.append((flare_list['Messenger_datetime'][index] - td(minutes=before)).strftime('%d-%b-%Y %H:%M:%S.000'))
        Mendtime.append((flare_list['Messenger_datetime'][index] + td(minutes=after)).strftime('%d-%b-%Y %H:%M:%S.000'))
        Rstarttime.append((flare_list['Messenger_datetime'][index] - td(minutes=before)).strftime('%d-%b-%Y %H:%M:%S.000'))
        Rendtime.append((flare_list['Messenger_datetime'][index] + td(minutes=after)).strftime('%d-%b-%Y %H:%M:%S.000'))
    print type(Mstarttime),Mendtime
    time_intervals={'Mstart_time':Mstarttime,'Mend_time':Mendtime,'Rstart_time':Rstarttime,'Rend_time':Rendtime}
    flare_list.update(time_intervals) #adds these to the flare list dictionary
    return flare_list

def make_plots(flist='full',wavelength='short',temp='1T',abund=True,ptype='scat',flux=False,pall=False,ymin=0.001,ymax=1000):
    '''Make selected plots cuz I'm lazy'''
    os.chdir('/Users/wheatley/Documents/Solar/occulted_flares/flare_lists')
    if flist == 'full': fliststr='all Messenger'
    if flist == 'joint': fliststr='jointly observed'
    if flist == 'vis': fliststr='visibly occulted'
    
    if not pall:
        if not flux:
            filename=flist+'_'+wavelength+'_'+temp+'_abund.sav'
            savename='../plots/'+flist+'_'+wavelength+'_'+temp+'_abund'
            abundstr = ' with abundances'
        else:
            filename=flist+'_'+wavelength+'_'+temp+'_abund_flux.sav'
            savename='../plots/'+flist+'_'+wavelength+'_'+temp+'_abund_flux'
            abundstr = ' with abundances'            
        if not abund:
            filename=flist+'_'+wavelength+'_'+temp+'.sav'
            abundstr=''
            savename='../plots/'+flist+'_'+wavelength+'_'+temp
        f=da.Data_Study(filename)
        title='Ratio of Messenger GOES to observed GOES for '+ fliststr + ' flares,\n' + wavelength + ' channel flux, ' + temp + ' fit' + abundstr
        if ptype =='scat':
            foo=plot_goes_ratio(f,title=title, mgoes='goes',scatter=True,labels=0,ymin=ymin,ymax=ymax,ylog=True,save=savename+'_scat.png')
        elif ptype == 'scat2': 
            foo=plot_goes_ratio(f,title=title, mgoes='goes',scatter=True,labels=0,ymin=ymin,ymax=ymax,ylog=True,cc='M',save=savename+'_scat2.png')
        else: #type = hist #need these elses to execute when 'all' also...
            foo=hist_ratio(f,title=title, gc='overlay',save=savename+'_hist.png')

    else: #make all the plots
        for fl in ['full','joint','vis']:
            if fl == 'full': fliststr='all Messenger'
            if fl == 'joint': fliststr='jointly observed'
            if fl == 'vis': fliststr='visibly occulted'
            
            for wv in ['short','long']:
                for t in ['1T','2T']:
                   filename1=fl+'_'+wv+'_'+t+'_abund.sav'
                   filename2=fl+'_'+wv+'_'+t+'.sav'
                   savename1='../plots/'+fl+'_'+wv+'_'+t+'_abund'                  
                   savename2='../plots/'+fl+'_'+wv+'_'+t
                   try:
                       f1=da.Data_Study(filename1)
                   except IOError:
                       print filename1 +' was not found!'
                       continue
                   try:
                       f2=da.Data_Study(filename2)
                   except IOError:
                       print filename2 +' was not found!'
                       continue
                       
                   title1='Ratio of Messenger GOES to observed GOES for '+ fliststr + ' flares,\n' + wavelength + ' channel flux, ' + t + ' fit with abundances'
                   title2='Ratio of Messenger GOES to observed GOES for '+ fliststr + ' flares,\n' + wavelength + ' channel flux, ' + t + ' fit'
                   foo=plot_goes_ratio(f1,title=title1, mgoes='goes',scatter=True,labels=0,ymin=ymin,ymax=ymax,ylog=True,save=savename1+'_scat.png',show=False)
                   foo=plot_goes_ratio(f1,title=title1, mgoes='goes',scatter=True,labels=0,ymin=ymin,ymax=ymax,ylog=True,cc='M',save=savename1+'_scat2.png',show=False)
                   foo=hist_ratio(f1,title=title1, gc='overlay',save=savename1+'_hist.png',show=False)
              
                   foo=plot_goes_ratio(f2,title=title2, mgoes='goes',scatter=True,labels=0,ymin=ymin,ymax=ymax,ylog=True,save=savename2+'_scat.png',show=False)
                   foo=plot_goes_ratio(f2,title=title2, mgoes='goes',scatter=True,labels=0,ymin=ymin,ymax=ymax,ylog=True,cc='M',save=savename2+'_scat2.png',show=False)
                   foo=hist_ratio(f2,title=title2, gc='overlay',save=savename2+'_hist.png',show=False)

    os.chdir('/Users/wheatley/Documents/Solar/occulted_flares/code')
        
