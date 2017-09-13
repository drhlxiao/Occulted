 #######################################
# OCDatetimes.py
# Erica Lastufka 17/5/17 

#Description: Class of functions to deal with all datetimes for one flare
#######################################

import numpy as np
import os
import glob
import pickle
from datetime import datetime as dt
import re

class OCDatetimes(object):
    ''' Data for this object will have the following attributes:
            Messenger Peak
            Rhessi peak
            Obs_start_time
            Obs_end_time
            Spec_start_time
            Spec_end_time
            lc_start_time
            lc_end_time
            pf_loop_time
            stereo_start_time
            stereo_times
         Methods are:
            format_datetimes() to format for reading from/exporting to csv/sav
            extend_time_int() does the time interval extensions
            '''
    def __init__(self, ID,legacy=True,filename=False, calc_times=False):
        '''Initialize the object, given the flare ID. This requires use of the OCFiles object I think...'''
        if legacy: #get info from legacy OCData object and the associated csv/sav files
            if not filename:
                filename= '/Users/wheatley/Documents/Solar/occulted_flares/flare_lists/list_final.csv'#default file to read
            import pandas as pd
            data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
            i=self.get_index(ID,data) #get the index of the flare if it's in a list
            self.Messenger_peak=data["Messenger_datetimes"][i]
            self.RHESSI_peak=data["RHESSI_datetimes"][i]
            self.Obs_start_time=data["Obs_start_time"][i]
            self.Obs_end_time=data["Obs_end_time"][i]
            try:
                self.Spec_start_time=data["Spec_start_time"][i]
            except KeyError:
                self.Spec_start_time=self.Obs_start_time #actually should modify this to do whatever time extension the code did
            try:
                self.Spec_end_time=data["Spec_end_time"][i]
            except KeyError:
                self.Spec_end_time=self.Obs_end_time
            try:                
                self.lc_start_time=data["lc_start_time"][i]
            except KeyError:
                self.lc_start_time=self.Obs_start_time
            try:                               
                self.lc_end_time=data["lc_end_time"][i]
            except KeyError:
                self.lc_end_time=self.Obs_end_time
            try:                
                self.pf_loop_time=data["pf_loop_time"][i]
            except KeyError:
                self.pf_loop_time=self.Obs_start_time
            try:                
                self.stereo_start_time=data["stereo_start_time"][i]
            except KeyError:
                self.stereo_start_time=self.Obs_start_time
            try:                
                self.stereo_times=data["stereo_times"][i]
            except KeyError:
                self.stereo_times=[self.Obs_start_time,self.Obs_end_time]

        if not legacy:
            #read datetimes csv file (can just restore the pickle file otherwise. Build this into OCFlare class):
            import pandas as pd
            if filename: #it's the big one
                data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
                #i=self.get_index(ID,data) #get the index of the flare if it's in a list
                for key in data.keys(): #only do this if it starts with Observation
                    if key.startswith('Datetimes.'):
                        dat=data[key]
                        key=key[key.find('.')+1:] #trim it
                        setattr(self,key,dat.values[0])
            if not filename:
                filename= '/Users/wheatley/Documents/Solar/occulted_flares/objects/'+str(ID)+'OCDatetimes.csv'#default file to read
                data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
                for key in data.keys(): #only do this if it starts with Observation
                    setattr(self,key,data.values[0])
                 
        self.convert2datetime()
        #update the OCFiles object with the file used to generate this instance of the object
        if calc_times:
            #self.convert2datetime() #in case it doesn't work the first time?
            self.calc_times(i)
            self.read_loop_time()
            self.read_stereo_times()
        
    def __iter__(self):
        '''Returns a generator that iterates over the object - not sure if it's needed here but whatever'''
        for attr, value in self.__dict__.iteritems():
            yield attr, value

    def write(self, picklename=False):
        '''Write object to pickle'''
        import pickle
        if not picklename:
            picklename='/Users/wheatley/Documents/Solar/occulted_flares/objects/'+str(self.ID)+'OCDatetimes.p'
        pickle.dump(self, open(picklename, 'wb'))

    def write_csv(self, csvname=False): #do I need to convert2string first?
        import csv
        if not csvname: csvname='/Users/wheatley/Documents/Solar/occulted_flares/objects/'+ str(ID) + 'OCDatetimes.csv'
        d=self.__dict__
        with open(csvname,'wb') as f:
            w=csv.writer(f)
            w.writerow(d.keys())
            w.writerow(d.values())

    def get_index(self,ID,data):
        '''Write object to pickle'''
        try:
            i= np.where(ID == np.array(data['ID']))
            i=i[0]
        except TypeError:
            i= np.searchsorted(np.array(datadata['ID']), ID)
        return i[0]

    def calc_times(self,i):
        '''Calculate lc and spec times based on parameters I had in the original lightcurves.py and spectrograms.py'''
        import lightcurves as lc #should I read it straight from the data files? perhaps this is safest...
        data,ebins=lc.import_spectrogram('/Users/wheatley/Documents/Solar/occulted_flares/data/rerun_sgrams.sav')
        last=int(data['len'][0]-1)
        self.lc_start_time= dt.strptime(data['UT'][i][0],'%d-%b-%Y %H:%M:%S.000') #let's hope everything has the same index....
        try:
            self.lc_end_time= dt.strptime(data['UT'][i][last],'%d-%b-%Y %H:%M:%S.000') 
        except ValueError:
            import datetime
            self.lc_end_time = start + datetime.timedelta(seconds=time_axis[-1])
        self.spec_start_time=self.lc_start_time
        self.spec_end_time=self.lc_end_time
        self.spec_end_time=self.lc_end_time

    def read_stereo_times(self):
        '''Read map.time from all STEREO fits files '''
        fdir='/Users/wheatley/Documents/Solar/occulted_flares/data/stereo_pfloops/'
        fname= dt.strftime(self.Messenger_peak,'%Y%m%d')
        ffiles=glob.glob(fdir+fname+'*.fts')
        #open it in IDL? or can read a fits file in Python... probably. should be in the header info one would hope
        from sunpy import io as sio
        try:
            for f in ffiles:
                head=sio.fits.get_header(f)
                self.stereo_times.append(dt.strptime(head[0]['DATE-OBS'].replace('T',' ')[:-3]+'000','%Y-%m-%d %H:%M:%S.000'))#convert to datetime
        except IndexError:
            self.stereo_times=[]
        
    def read_loop_time(self,filename=False):
        '''Read map.time from STEREO fits files used to do the loop tracing. filename can be set using the OCFiles object'''
        #find the correct file
        fdir='/Users/wheatley/Documents/Solar/occulted_flares/data/stereo_pfloops/'
        fname= dt.strftime(self.Messenger_peak,'%Y%m%d')#format the datetime correctly
        if not filename: filename=glob.glob(fdir+fname+'*.fts')
        ffile=filename
        #open it in IDL? or can read a fits file in Python... probably. should be in the header info one would hope
        from sunpy import io as sio
        head=sio.fits.get_header(ffile[0])
        self.pf_loop_time = dt.strptime(head[0]['DATE-OBS'].replace('T',' ')[:-3]+'000','%Y-%m-%d %H:%M:%S.000')#convert to datetime
    
    def convert2datetime(self):
        '''If it's a string make it a datetime'''
        for a in dir(self):
            if not a.startswith('__') and not callable(getattr(self,a)) and type(getattr(self,a)) == str:
                val=getattr(self,a)
                #print a,val
                if val == '':
                    fc=''
                elif ('AM' or 'PM') in val:
                    fc = '%m/%d/%y %I:%M:%S %p'
                elif '/' in val:
                    fc = '%m/%d/%y %H:%M:%S.00'
                elif re.search('[a-z]', val) is not None: #not None if there is a letter
                    fc ='%d-%b-%Y %H:%M:%S.000' #current format code - might need to add .000
                else: #it's come straight from idl ....
                    fc ='%Y-%m-%d %H:%M:%S'
                if fc != '':
                    if val.startswith("'"):
                        try:
                            setattr(self,a,dt.strptime(val[1:-1], fc))
                        except ValueError:
                            fc=fc ='%d-%b-%Y %H:%M:%S.000'
                            setattr(self,a,dt.strptime(val[1:-1], fc))
                    else:
                        try:
                            setattr(self,a,dt.strptime(val, fc))
                        except ValueError:
                            pass
                            #fc=fc ='%d-%b-%Y %H:%M:%S.000'
                            #setattr(self,a,dt.strptime(val, fc))

    def convert2string(self):
        '''If it's a datetime make it a string'''
        for a in dir(self):
            if not a.startswith('__') and not callable(getattr(self,a)) and type(getattr(self,a)) == dt:
                val=getattr(self,a)
                fc='%d-%b-%Y %H:%M:%S.000'
                setattr(self,a,dt.strftime(val,fc))
