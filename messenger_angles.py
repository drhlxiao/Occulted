 #######################################
# messenger_angles.py
# Erica Lastufka 30/11/2016

#Description: Check if Messenger can see an occulted flare
#######################################

#######################################
# Usage:

# for default output: python messenger_angles.py
######################################

#import glob
#import os
import numpy as np
#import grid_eff_area as ga
import pidly #set up pidly
#import csv
import scipy.constants as sc
import pandas as pd
from datetime import datetime
from datetime import timedelta as td
import os
import pickle

def match_flare_list(l1,l2):
    #l1 is OCflare object, l2 is not, it's a dictonary
    l2dt=l2['datetimes']
    for l in l1.list:
        bt=l.Datetimes.Messenger_peak
        mdt=min(l2dt, key=lambda x: abs(x - bt))
        mindex=l2dt.index(mdt)
        l.Properties.dennisT1=l2['T1'][mindex]
        l.Properties.dennisT2=l2['T2'][mindex]
        l.Properties.dennisEM1=l2['EM1'][mindex]
        l.Properties.dennisEM2=l2['EM2'][mindex]
        l.Properties.Fe_abund=l2['Fe'][mindex]

def read_flare_list(filename='test.csv'):
    '''reads in list of occulted flare positions and time intervals'''
    #edit 2020 read the rest of everything as well
    #Number	Date	Start time	End time	Duration	EM1	T1	EM2	T2	Fe	Ca	S	Si	Ar	Mg	drmmod	drmmod	drmmod	chimax			chimax
    #fit_resK = pd.read_csv('Occulted_Results_23.csv', sep=',', index_col=0) - Frederic's data
    #fit_resO = pd.read_csv('Occulted_Results_24.csv', sep=',', index_col=0)
    #fit_res = pd.concat([fit_resK,fit_resO])

    data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
    #print data.keys()
    datelen=data['Date']
    newdates=[]
    datetimes=[]
    starttimes=[]
    endtimes=[]
    strdates=newdates
    for (date,time,i) in zip(data['Date'],data['Start_Time'],range(0,len(datelen)-1)):
        try:
            newdates.append(str(datetime.strptime(date, '%Y %b %d').strftime('%d-%b-%Y')))
        except ValueError:
            newdates=datelen.tolist()
            break

    if data['Start_Time'][0].endswith('M'): #convert to 24 hour time
        for (s,e,i) in zip(data['Start_Time'],data['End_Time'],range(0,len(datelen)-1)):
            if s.endswith('PM'):
                news=datetime.strptime(s,'%I:%M:%S %p')
                new_s = str(news.time())
            else:
                new_s = s[:-2]
            if e.endswith('PM'):
                newe=datetime.strptime(e,'%I:%M:%S %p')
                new_e = str(newe.time())
            else:
                new_e = e[:-2]
            starttimes.append(new_s)
            endtimes.append(new_e)
            datetimes.append(datetime.combine(datetime.strptime(newdates[i],'%d-%b-%Y'),datetime.strptime(new_s.strip(),'%H:%M:%S').time()))
    else:
        starttimes=data['Start_Time'].tolist()
        endtimes=data['End_Time'].tolist()
        for (date,time) in zip(data['Date'],data['Start_Time']):
            ndate = datetime.strptime(date,'%d-%b-%Y')
            ntime = datetime.strptime(time,'%H:%M:%S')
            datetimes.append(datetime.combine(ndate.date(),ntime.time()))

    flare_list = {'date':newdates,'stime':starttimes,'etime':endtimes,'datetimes':datetimes}
    if 'Angle' in data.keys():
        flare_list['Angle']=data['Angle'].tolist()
    if 'ID' in data.keys():
        flare_list['ID']=data['ID'].tolist()

    return flare_list

def is_flare_observed(list1,list2, cutoff=7200):
    '''Checks if Messenger observed anything within the flare time intervals. Use Brian Dennis's flare list. They shouldn't be too different for an occulted flare? Cutoff time should be specified in seconds, default is 2 hours'''

    dt1=list1['datetimes']
    dt2=list2['datetimes']
    id1=list1['ID']
    #test=['03/28/02', '04/04/02', '04/04/02', '04/18/02', '04/22/02', '04/29/02', '04/30/02', '04/30/02', '05/17/02', '05/17/02']
    #test2=['05/28/07', '05/28/07', '04/04/02', '05/29/07', '05/29/07', '05/30/07', '06/01/07', '06/01/07', '06/01/07', '06/01/07']
    #found_list = list(set(aa) & set(bb)) #This misses things that might have happened on different days relative to RHESSI/messenger. Can fix this by just checking diff in datetime? for future test....

    IDs,newfoundlist, nMtime,nRtime,nMdt,nRdt=[],[],[],[],[],[]
    for dt,id in zip(dt1,id1):
        closest = map(lambda d: abs(d-dt), dt2) #this is a list
        for item,index in zip(closest, range(0,len(dt2))):
            if item.total_seconds() < cutoff:
                print dt, '    ',dt2[index]
                newfoundlist.append(str(dt.date()))
                nMtime.append(str(dt2[index].time()))
                nRtime.append(str(dt.time()))
                nMdt.append(dt2[index])
                nRdt.append(dt)
                IDs.append(id)

    list_dict={'ID':IDs,'date':newfoundlist, 'Messenger_time':nMtime, 'RHESSI_time':nRtime, 'Messenger_datetime':nMdt,'RHESSI_datetime':nRdt}
    #list_dict = 'foo'
    return list_dict

def messenger_flare_angle(flare_list):
    '''Executes sswidl function messenger_flare_angle()'''
    angle_list=[]
    idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')

    for flare in flare_list['date'][0:10]:
        #is_observed == False:
        #    angle=-1
        #else:
        #print flare
        #time=flare['date']
        if flare > '8/30/2007' and flare < '4/1/2014':
            #time=flare[0]+','+flare[1]#.tolist()
            print flare
            xloc=-200.
            yloc=200.
            aa=idl('print,"this is idl"')
            print aa
            angle = idl.messenger_flare_angle([-200.,200.], '15-feb-2011') #mess_helio_lonlat=mess_helio_lonlat
            #append angle to flare list
            angle_list.append(angle)
    idl.close()
    print angle_list
    #print time,angle_list
    return angle_list

def save_flare_list(flare_list, fname):
    '''save the final flare list'''
    pickle.dump(flare_list, open(fname, 'wb'))

def save_flare_list_idl(flare_list, fname):
    '''save the final flare list in IDL save file'''
    #start new IDL session
    length=len(flare_list['date'])
    dates=flare_list['date']
    Mtimes=flare_list['Messenger_time']
    Rtimes=flare_list['RHESSI_time']
    if 'Rend_time' in list_selected.keys():
        Mst=flare_list['Mstart_time']
        Met=flare_list['Mend_time']
        Rst=flare_list['Rstart_time']
        Ret=flare_list['Rend_time']

    idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
    #transfer variables there
    idl('len',length)
    idl('dates',dates)
    idl('Mtimes',Mtimes)
    idl('Rtimes',Rtimes)
    idl('fname',fname)
    if 'Rend_time' in list_selected.keys():
        #print 'also here', Mst, Mtimes #type(Mst) -> string!, type(Mtimes) -> list
        idl('Mst',Mst)
        idl('Met',Met)
        idl('Rst',Rst)
        idl('Ret',Ret)
        idl('flare_list={date:dates,Messenger_time:Mtimes,RHESSI_time:Rtimes,Mstart_time:Mst,Mend_time:Met,Rstart_time:Rst,Rend_time:Ret}')

    else:
        idl('flare_list={date:dates,Messenger_time:Mtimes,RHESSI_time:Rtimes}')
    #save them
    idl('save,flare_list, filename=fname')

def plot_angle_distribution(angle_list):
    '''make a histogram of the angle distributions'''

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    n, bins, patches = plt.hist(angle_list['Angle'], 18, facecolor='green', alpha=0.75)

    plt.xlabel('Angle between Earth and Mercury (degrees)')
    plt.ylabel('Number of Events')
    #ax1.set_ylim([0,1])
    #ax1.set_xlim([0,150])
    #plt.title("Effective area of grids "+ substrate_properties['material'] + ' '+ substrate_properties['thickness'] + '$\mu$m, ' + slit_properties['material'] + ' '+ slit_properties['thickness']+'$\mu$m' )
    #ax1.legend(loc='upper right',fontsize='medium')

    ax1.plot()

    fig.show()

def list_selected(indices=[9,10,16,18]):
    '''Select certain flares from the list to be examined in sswidl'''
    #indices=[9,10,16,18]
    list=pickle.load(open('list_observed.p','rb'))
    ldate,lMtime,lRtime,lMdt,lRdt=[],[],[],[],[]
    for index in indices:
        ldate.append(list['date'][index])
        lMtime.append(list['Messenger_time'][index])
        lRtime.append(list['RHESSI_time'][index])
        lMdt.append(list['Messenger_datetime'][index])
        lRdt.append(list['RHESSI_datetime'][index])
    list_selected={'date':ldate,'Messenger_time':lMtime,'RHESSI_time':lRtime, 'Messenger_datetime':lMdt,'RHESSI_datetime':lRdt}
    #print list_selected
    return list_selected

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

# if __name__ == "__main__":
#     os.chdir('/Users/wheatley/Documents/Solar/occulted_flares/')
#     occulted_list = read_flare_list(filename='list_filtered.csv')
#     messenger_list = read_flare_list(filename='messenger_list.csv')
#     list_observed = is_flare_observed(occulted_list, messenger_list,cutoff=3600)
#     save_flare_list(list_observed,'list_observed_1hour.p')
#     pd.DataFrame(list_observed).to_csv('list_observed_1hour.csv')
#     #df.to_csv('list_observed_1hour.csv')
#     #pickle.load(open('list_observed.p','rb'))
#     #list_selected=list_selected()
#     #save_flare_list(list_selected,'list_selected.p')
#     #save_flare_list_idl(list_selected,'list_selected.sav')
#     #time_intervals=select_time_intervals(list_selected)
#     #save_flare_list(time_intervals,'time_intervals.p')
#     #save_flare_list_idl(time_intervals,'time_intervals.sav')

#     #angle_list = messenger_flare_angle(messenger_list)
#     #plot_angle_distribution(angle_list)
#     #save_flare_list(angle_list)
#     #angle_list = read_flare_list(filename='messenger_angles.csv')
#     #plot_angle_distribution(angle_list)


