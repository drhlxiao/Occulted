 #######################################
# OC.Flare.py
# Erica Lastufka 17/5/17

#Description: Class of functions to deal with all data,files,calculations of the study
#######################################

import numpy as np
import glob
import pickle
import sunpy.map
from datetime import datetime as dt
from datetime import timedelta
import os
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd

import OCDatetimes as ocd
import OCFiles as ocf
import OCProperties as ocp
import OCObservation as oco

class OCFlare(object):
    ''' Object will have the following attributes:
            ID: RHESSI flare ID or generated ID
            Datetimes: OCDatetimes object (might inherit from OCData object)
            Files: OCFiles object that contains all information about file names and locations
            Data: OCData object (different from original) that contains methods for loading/accessing the data ->tag this onto Files?
            Properties: OCProperties object that contains all information about the flare itself
            Observation:OCObservation object that contains all information about the observations
            Notes
        Methods will include:
            Method for translating between pickle, .sav and .csv
        A list of OCFlare objects will comprise an OCFlareList object, which will include all the plotting methods
     '''
    def __init__(self, ID, legacy=True, calc_times=False, from_csv=False, filename= False,calc_missing=False,gen=False,freload=False):
        '''Initialize the flare object from an ID, pickle files in given folder, or legacy .csv or .sav'''
        self.ID=ID
        ppath='/Users/wheatley/Documents/Solar/occulted_flares/data/att_data/'+str(ID)+'_att.p'
        ispickle=glob.glob(ppath)
        if ispickle !=[] and freload == False:
            atts=pickle.load(open(ispickle[0],'rb')) #this is the attributes. We need to distribute them now.
            self.Datetimes=ocd.OCDatetimes(ID,legacy=False,calc_times=calc_times,att_dict=atts['Datetimes'][1])
            self.Properties=ocp.OCProperties(ID,legacy=False,calc_missing=calc_missing,att_dict=atts['Properties'][3])
            self.Files=ocf.OCFiles(ID,legacy=False,att_dict=atts['Files'][2],gen=gen)
            self.Observation=oco.OCObservation(ID,legacy=False,att_dict=atts['Observation'][4])
            self.Notes= atts['Notes'][5]
        elif type(from_csv)==pd.core.frame.DataFrame:
             #make the attribute dictionaries and pass them on
            atts={}
            if 'Date' in from_csv.keys():
                atts['Datetimes']={'Date':from_csv['Date'][from_csv[from_csv.keys()[0]].tolist().index(ID)]}
            if 'Time' in from_csv.keys():
                atts['Datetimes']['Time']=from_csv['Time'][from_csv[from_csv.keys()[0]].tolist().index(ID)]

            self.Datetimes=ocd.OCDatetimes(ID,legacy=False,calc_times=calc_times,att_dict=atts['Datetimes'],convert=False)
            self.Properties=ocp.OCProperties(ID,legacy=False,calc_missing=calc_missing)
            self.Files=ocf.OCFiles(ID,legacy=False,gen=gen)
            self.Observation=oco.OCObservation(ID,legacy=False)
            #self.Notes= atts['Notes'][5]

        else: #reload ... ie regenerate the attribute pickle even if there is one. run in //?
            self.Datetimes=ocd.OCDatetimes(ID,legacy=legacy,calc_times=calc_times,filename=filename)
            self.Properties=ocp.OCProperties(ID,legacy=legacy,calc_missing=calc_missing,filename=filename)
            self.Files=ocf.OCFiles(ID,legacy=legacy,filename=filename,gen=gen)
            self.Observation=oco.OCObservation(ID,legacy=legacy,filename=filename)
            self.Notes= [""]
            #self.extract_stereo_times()

    def att2df(self):
        '''Store the flare attributes in the database'''
        #convert to pandas DataFrame
        df=pd.DataFrame.from_dict([{'ID':self.ID},{'Datetimes':self.Datetimes.__dict__},{'Files':self.Files.__dict__},{'Properties':self.Properties.__dict__},{'Observation':self.Observation.__dict__},{'Notes':self.Notes}])
        return df

    def jsonify(self,key):
            '''convert attributes that are numpy.float32 etc to normal floats'''
            obj=getattr(self,key)
            for k in obj.__dict__.keys():
                if type(getattr(obj,k)) == np.float32:
                    setattr(obj,k,float(getattr(obj,k))) #what if it's a list or dict? deal with that later...

    def export2csv(self, filename=False, save=True):
        '''Exports Python dictionary to csv file. First convert to json to enable normalization'''
        import json
        from pandas.io.json import json_normalize

        #prep for json
        self.Datetimes.convert2string()
        self.Datetimes.convert2string()
        #self.jsonify('Properties')
        #self.jsonify('Observation')
        df=self.att2df()

        jsonser={'ID':df['ID'][0],'Datetimes':df['Datetimes'][1],'Files':df['Files'][2],'Properties':df['Properties'][3],'Observation':df['Observation'][4],'Notes':df['Notes'][5]}

        #move to json - guess I don't actually have to do this though? can normalize straight from dataframe
        #json_df=json.dumps(jsonser)

        #normalize
        norm_df=json_normalize(jsonser)

        if save:
            if not filename:
                filename=str(self.ID)+'_att.csv'
            os.chdir(self.Files.dir+self.Files.folders['flare_lists'])

            import csv
            with open(filename,'wb') as f:
                writer=csv.writer(f)
                writer.writerow(norm_df.keys())
                writer.writerow([norm_df[k][0] for k in norm_df.keys()])
            os.chdir(self.Files.dir)
        else: #for writing a whole list
            return jsonser

    def export2idl(self, savname):
        '''because the dictionaries might be to big to send directly, read it in from the csv file'''
        #convert datetimes to IDL-friendly format
        #Format for time is YY/MM/DD, HHMM:SS.SSS
        #or alternatively, DD-Mon-YY HH:MM:SS.SSS -> '%d-%b-%Y %H:%M:%S.000'

        for i,ang in enumerate(self.Angle):
            if np.isnan(ang) == 1:
                self.Angle[i] = 0.0 #to avoid any awkwardness with IDL when it does the calculation

        self.format_datetimes() #this should turn them all into datetimes UNLESS they are empty strings
        data = self.Datetimes
        for key in data.keys():
            #check that it's got quote marks then continue
            newval=[]
            if type(data[key][0]) == str:
                    continue #it's an empty string so let's leave it that way
            else: #make it an IDL-friendly string and put quotes around everything
                for item in data[key]:
                    newval.append("'"+item.strftime('%d-%b-%Y %H:%M:%S.000')+"'")
                self.Datetimes[key]=newval

        #check if the csv file exists, if not, make it
        import glob
        if not savname[:-3]+'csv' in glob.glob('*.csv'):
            print 'no csv file found, converting to csv first'
            self.export2csv(savname[:-3]+'csv')
        else:
            ans=raw_input(savname[:-3]+'csv already exists! Overwrite? (Y/N)')
            if ans == 'Y':
                self.export2csv(savname[:-3]+'csv')

        import pidly
        idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
        idl('savname',savname)
        idl('csvname',savname[:-3]+'csv')
        idl('foo=read_csv(csvname)')
        idl('flare_list={ID:foo.field01,Datetimes:{Messenger_datetimes:foo.field02,RHESSI_datetimes:foo.field03,Obs_start_time:foo.field04,Obs_end_time:foo.field05},Flare_properties:{Messenger_T:foo.field06,Messenger_EM1:foo.field07,Messenger_GOES:foo.field08,Messenger_total_counts:foo.field09,RHESSI_GOES:foo.field10,GOES_GOES:foo.field11,RHESSI_total_counts:foo.field12,Location:foo.field13,source_vis:foo.field14},Data_properties:{Messenger_data_path:foo.field15,RHESSI_data_path:foo.field16,XRS_files:foo.field17,QL_images:foo.field18,RHESSI_browser_urls:foo.field19,csv_name:foo.field20,good_det:foo.field21},Angle:foo.field22,Notes:foo.field23}')
        idl('save,flare_list, filename=savname')

    def export2pickle(self, picklename=False, silent=True):
        '''Replaces save_flare_list. Saves the data frame in a .p file.'''
        import pickle
        if self.Files.att_data == '':
            picklename=str(self.ID)+'_att.p'
        else:
            picklename=self.Files.att_data
        path=self.Files.dir+self.Files.folders['att_data']+'/'
        df=self.att2df()
        pickle.dump(df, open(path+picklename, 'wb'))
        if not silent:
            print picklename + ' was created in '+ path

    #############################################  DATA READING & STORAGE METHODS  #######################################################

    #do I want the data to be an separate object? possibly....it would need to be initiated along with the other classes though to have the datetime/file info
    #it could be an optional attribute? I don't want to clog up memory
    #maybe don't store...

    def get_RHESSI_lightcurve(self, filename=False, save=False):
        '''Get data for this flare's rhessi lightcurve. Save in new file if requested.'''
        os.chdir(self.Files.dir + self.Files.folders['lightcurves'])
        if self.Files.Raw['lightcurves'].endswith('.p'):
            pickle.load(self.Files.Raw['spectrogram']) #restore and be done. should be a mySpectrogram object
            os.chdir(self.Files.dir)
            return #think I probably have to assign a variable name...

        if self.Files.Raw['spectrogram'] == '': #get it from the sav file
            if not filename: filename='rerun_sgrams.sav'
            from scipy.io import readsav
            os.chidr('../')
            d=readsav(filename, python_dict=True)
            Rdata=d['data']['rdata'] #access via Rdata[0]['UT'], etc
            for j in len(data['len']):
                if self.Datetimes.Messenger_peak.date() == dt.strptime(Rdata['UT'][j][0],'%d-%b-%Y %H:%M:%S.000').date():
                    i=j #use this index
                    break
            Rdata=Rdata[i] #how is this organized again???
        if save:
            os.chdir(self.Files.dir + self.Files.folders['spectrograms'])
            newfname=str(self.ID)+'RHESSIlc.p'
            pickle.dump(a,open(newfname,'wb'))
            self.Files.Raw['lightcurves']=newfname
        os.chdir(self.Files.dir)
        return Rdata

    def download_messenger(self): #CHECK THAT THIS IS COMPATIBLE WITH OBJECT CONVENTIONS
        '''Downloads Messenger .dat and .lbl files from the database, given event date. Can be used to get missing files too.'''
        import urllib
        dataurl = 'https://hesperia.gsfc.nasa.gov/messenger/'
        dt =self.Datetimes['Messenger_peak'].date()
        datestring=dt.strftime('%Y%j')
        newurl=dataurl+datestring[0:4] #just the year

        #first check if the file is already there:
        if not self.Files.Raw['xrs_files']:
            filename='/xrs'+datestring
            self.Files.Raw['xrs_files']=filename
        else:
            filename=self.Files.Raw['xrs_files']
        datfile=filename+'.dat'
        lblfile=filename+'.lbl'
        # check if the file is already there:
        if not os.path.exists('/Users/wheatley/Documents/Solar/occulted_flares/data/dat_files/'+filename+'.dat'):
            urllib.urlretrieve(newurl+datfile,'/Users/wheatley/Documents/Solar/occulted_flares/data/dat_files/'+filename+'.dat')
        if not os.path.exists('/Users/wheatley/Documents/Solar/occulted_flares/data/dat_files/'+filename+'.lbl'):
            urllib.urlretrieve(newurl+lblfile,'/Users/wheatley/Documents/Solar/occulted_flares/data/dat_files/'+filename+'.lbl')

    def download_norp(self):
        '''Downloads NORP .xdr files from the database, given event date. Can be used to get missing files too.'''
        import urllib
        dataurl = 'ftp://solar-pub.nao.ac.jp/pub/nsro/norp/xdr/'
        if type(self.Datetimes.Date)==str:
            dtm=dt.strptime(self.Datetimes.Date,'%Y %b %d')
            datestring=dt.strftime(dtm,'%Y%m%d')
        else:
            dtm =self.Datetimes.Date.date()
            datestring=dt.strftime(dtm,'%Y%m%d')
        newurl=dataurl+datestring[0:4]+'/'+datestring[4:6] #year and month

        #first check if the file is already there:
        if not self.Files.Raw['norp_files']:
            filename='/norp'+datestring
            self.Files.Raw['norp_files']=filename
        else:
            filename=self.Files.Raw['norp_files']
        xdrfile=filename+'.xdr'
        # check if the file is already there:
        if not os.path.exists('/Users/wheatley/Documents/Solar/occulted_flares/Nobeyama/NORP/data/'+xdrfile):
            urllib.urlretrieve(newurl+xdrfile,'/Users/wheatley/Documents/Solar/occulted_flares/Nobeyama/NORP/data/'+xdrfile)

    def gen_NORP_lightcurve(self,offsets=[0,.1,.2,.3,.4,.5],mint=10,norm_int=False,write=True):
        '''Make the NORP lightcurves, offset for visibility. Normalise either by given norm_int or the first 20 or so ponits.'''
        import make_timeseries_plot_flux as mtp # plot_norp_lc
        #some basic checks to make life easier
        try:
            getattr(self.Files,'dir')
        except AttributeError:
            self.Files.dir='/Users/wheatley/Documents/Solar/occulted_flares/'
        try:
            getattr(self.Files,'Raw')
        except AttributeError:
            self.Files.Raw={}
        try:
            getattr(self.Files,'Folders')
        except AttributeError:
            self.Files.Folders={}
        if 'norp_files' not in self.Files.Folders.keys():
            self.Files.Folders['norp_files']='Nobeyama/NORP/Data/'
        if 'norp_files' not in self.Files.Raw.keys():
            if type(self.Datetimes.Date) == str:
                datestr0=dt.strptime(self.Datetimes.Date,'%Y %b %d')
                datestr=dt.strftime(datestr0,'%Y%m%d')
            else:
                datestr=dt.strftime(self.Datetimes.Date,'%Y%m%d')
            norpfile='norp'+datestr[2:]+'.xdr'
            self.Files.Raw['norp_files']=[norpfile]
        else:
            norpfile=self.Files.Raw['norp_files'][0]
        if write:
            self.Files.Raw['norp_files'].append(norpfile[:-4]+'_'+str(mint)+'s.p')
        #get the peak time to plot as a vertical line for reference
        pt_time=self.Datetimes.Time
        if type(pt_time) == str:
            peaktime=dt.strptime(self.Datetimes.Date+pt_time,'%Y %b %d%H:%M:%S')
            time_int=[peaktime-timedelta(hours=2), peaktime+timedelta(hours=3)]
        else:
            peaktime=False
            time_int=False

        #defaults
        from_idl=False
        raw=False
        if norpfile[-4:]=='.xdr':
            from_idl=True
        elif 'raw' in norpfile:
            raw=True

        os.chdir(self.Files.dir + self.Files.Folders['norp_files'])
        mtp.plot_norp_lc(norpfile=norpfile,raw=raw,from_idl=from_idl,mint=mint,norm_int=norm_int,plot=False,save=True,write=write,peaktime=peaktime,time_int=time_int)
        os.chdir(self.Files.dir)


    def gen_Messenger_lightcurve(self, countspersec=True, ret=False,save=True,recalc=False):
        '''Get data FROM SCRATCH aka IDL for this flare's Messenger lightcurves, and adjust to counts/s if specified (raw = False). Save in new file if requested. Write peak flux and peak time intervals to object.'''
        import os.path
        import pickle
        try:
            if recalc:
                data=pickle.load(open('notafile.p','rb'))
            data=pickle.load(open(self.Files.folders['lightcurves']+self.Files.Raw['messenger'],'rb'))
            countspersec=False #assume it's already been fixed
            save=False
        except IOError:

            self.Datetimes.convert2string() #to make IDL easier

            import pidly
            idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
            idl(".compile lightcurves.pro")
            idl('ft',[self.Datetimes.lc_start_time,self.Datetimes.lc_end_time])
            idl('ebinlow',[1.5,12.4])
            idl('ebinhigh',[3,25])
            idl('xrs_file',self.Files.Raw['xrs_files'])
            idl('clow=one_curve_M(ft,xrs_file,ebinlow,quiet=1)') #returns data={taxis:taxis,phigh:phigh}
            idl('chigh=one_curve_M(ft,xrs_file,ebinhigh,quiet=1)') #returns data={taxis:taxis,phigh:phigh}
            clow=idl.clow
            chigh=idl.chigh
            idl.close() #could probably speed this up if I get the idl bridge to work, then don't have to open and close a session each time
            Mtiml=[]
            for t in clow['taxis']: Mtiml.append(dt.strptime(t,'%d-%b-%Y %H:%M:%S.%f')) #fix messenger times to datetimes

            #since OSPEX ignores my time interval commands, fix it here:
            startidx=np.searchsorted(Mtiml,dt.strptime(self.Datetimes.lc_start_time,'%d-%b-%Y %H:%M:%S.000'))
            endidx=np.searchsorted(Mtiml,dt.strptime(self.Datetimes.lc_end_time,'%d-%b-%Y %H:%M:%S.000'))

            data={'taxis':Mtiml[startidx-1:endidx+1],'lowEcounts':clow['phigh'][startidx-1:endidx+1],'highEcounts':chigh['phigh'][startidx-1:endidx+1]} # +-1 there because time bins are so big... better to have more than less

        if len(data['taxis']) !=0:
            #adjust to counts per second if requested
            if countspersec:
                #get the time bin for each data point
                tlow=data['taxis']
                mlenl=len(tlow)
                #thigh=chigh['taxis'] #should be the same in theory, but do this just in case
                #mlenh=len(thigh)
                Mtiml,Mtimh,cps1,cps2=[],[],[],[]
                for t in tlow: Mtiml.append(t) #it's already a datetime now
                #for t in thigh: Mtimh.append(dt.strptime(t,'%d-%b-%Y %H:%M:%S.%f')) #fix messenger times to datetimes

                for i in range(mlenl-1):
                    tbin=Mtiml[i+1]-Mtiml[i] #timedelta object in seconds
                    cps1.append(data['lowEcounts'][i]/tbin.total_seconds())
                    cps2.append(data['highEcounts'][i]/tbin.total_seconds())

                data={'taxis':Mtiml,'lowEcounts':cps1,'highEcounts':cps2}

            #now get the peaks and time bins, write to object
            lpeak=np.searchsorted(data['lowEcounts'],np.max(data['lowEcounts'])) #gives index of the peak
            hpeak=np.searchsorted(data['highEcounts'],np.max(data['highEcounts']))
            #get time range around peak...should I take the time before or after the peak? let's do after to be consistent with the plotting method...

            #store stuff
            self.Datetimes.Messenger_peak = data['taxis'][hpeak]
            try:
                self.Datetimes.messenger_peak_bin['long'] = [data['taxis'][lpeak],data['taxis'][lpeak+1]]
                self.Datetimes.messenger_peak_bin['short'] = [data['taxis'][hpeak],data['taxis'][hpeak+1]]#for use with gen_GOES_lightcurve to have a more accurate peak value in GOES
            except IndexError: #in case the last value is the max
                self.Datetimes.messenger_peak_bin['long'] = [data['taxis'][lpeak-1],data['taxis'][lpeak]]
                self.Datetimes.messenger_peak_bin['short'] = [data['taxis'][hpeak-1],data['taxis'][hpeak]]
            #print self.Datetimes.messenger_peak_bin

            self.Properties.Messenger_GOES_long=np.max(data['lowEcounts'])
            #self.Properties.Messenger_GOES=np.max(data['highEcounts']) #for now! since it's counts and we want flux

            #pickle it up if requested, save filename and location to object
            if save:
                import pickle
                fdir=self.Files.dir+self.Files.folders['lightcurves']
                filename= '/'+str(self.ID)+'_Mess_lc.p' #something
                self.Files.Raw['messenger']=filename
                pickle.dump(data,open(fdir+filename,'wb'))
                print 'File '+filename+' saved in '+fdir

        else:
            print 'No data available for flare ' + str(self.ID)

        os.chdir(self.Files.dir)
        self.Datetimes.convert2datetime()
        #print self.Datetimes.messenger_peak_bin

        if ret:
            return data

    def gen_Messenger_spectrum(self, ret=False,save=True):
        '''Get data FROM SCRATCH aka IDL for this flare's Messenger spectrum. Save in new file if requested. '''
        import os.path
        import pickle
#         try:
#             if recalc:
#                 data=pickle.load(open('notafile.p','rb'))
#             data=pickle.load(open(self.Files.folders['lightcurves']+self.Files.Raw['messenger'],'rb'))
#             countspersec=False #assume it's already been fixed
#             save=False
#         except IOError:

        self.Datetimes.convert2string() #to make IDL easier

        import pidly
        idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
        idl(".compile lightcurves.pro")
        idl('ft',[self.Datetimes.Spec_start_time,self.Datetimes.Spec_end_time])
        idl('ebinlow',[1.5,12.4])
        #idl('ebinhigh',[3,25])
        idl('xrs_file',self.Files.Raw['xrs_files'])
        idl('clow=one_spec_M(ft,xrs_file,ebinlow,quiet=1)') #returns data={taxis:taxis,phigh:phigh}
        #idl('chigh=one_curve_M(ft,xrs_file,ebinhigh,quiet=1)') #returns data={taxis:taxis,phigh:phigh}
        clow=idl.clow
        #chigh=idl.chigh
        idl.close() #could probably speed this up if I get the idl bridge to work, then don't have to open and close a session each time
        #    Mtiml=[]
        #    for t in clow['taxis']: Mtiml.append(dt.strptime(t,'%d-%b-%Y %H:%M:%S.%f')) #fix messenger times to datetimes

        #since OSPEX ignores my time interval commands, fix it here:
        #startidx=np.searchsorted(Mtiml,dt.strptime(self.Datetimes.lc_start_time,'%d-%b-%Y %H:%M:%S.000'))
        #endidx=np.searchsorted(Mtiml,dt.strptime(self.Datetimes.lc_end_time,'%d-%b-%Y %H:%M:%S.000'))

        data={'eaxis':clow['eaxis'],'ctot':clow['ctot']} # +-1 there because time bins are so big... better to have more than less

        #store stuff
        #self.Datetimes.Messenger_peak = data['taxis'][hpeak]

        #pickle it up if requested, save filename and location to object
        if save:
            import pickle
            fdir=self.Files.dir+self.Files.folders['lightcurves']
            filename= '/'+str(self.ID)+'_Mess_spec.p' #something
            self.Files.Raw['messenger']=filename
            pickle.dump(data,open(fdir+filename,'wb'))
            print 'File '+filename+' saved in '+fdir

        else:
            print 'No data available for flare ' + str(self.ID)

        os.chdir(self.Files.dir)
        self.Datetimes.convert2datetime()
        #print self.Datetimes.messenger_peak_bin

        if ret:
            return data


    def gen_GOES_lightcurve(self, ret=False,save=True,recalc=True):
        '''Get data FROM SCRATCH aka IDL for this flare's GOES lightcurves. Save in new file if requested. Write peak flux  to object.'''
        import os.path
        try:
            import pickle
            if recalc:
                data=pickle.load(open('notafile.p','rb'))
            data=pickle.load(open(self.Files.folders['lightcurves']+self.Files.Raw['goes'],'rb'))
            save=False
        except IOError:
            self.Datetimes.convert2string() #to make IDL easier

            import pidly
            idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
            idl(".compile lightcurves.pro")
            idl('ft',[self.Datetimes.lc_start_time,self.Datetimes.lc_end_time])
            idl('d=one_curve_G(ft,quiet=1)') #returns data={taxis:taxis,phigh:phigh}
            d=idl.d
            idl.close()
            if d== -1:
                print 'no data in array for flare ' + str(self.ID)
                return
            Mtiml=[]
            for t in d['tarray']: Mtiml.append(dt.strptime(t,'%d-%b-%Y %H:%M:%S.%f')) #fix messenger times to datetimes

            data={'taxis':Mtiml,'lowEflux':d['ydata'][0],'highEflux':d['ydata'][1]}

            self.Datetimes.convert2datetime()

        #now get the peaks and time bins, write to object
        binlong=self.Datetimes.messenger_peak_bin['long']
        lidx=[]
        for i in range(0,len(data['taxis'])-1):
            if data['taxis'][i] > binlong[0] and data['taxis'][i] < binlong[1]:
                lidx.append(i) #the index of the value
        lrange=[lidx[0],lidx[len(lidx)-1]]
        if lrange[0]==lrange[1]: #don't allow this!
            lrange=[lidx[0]-1,lidx[0]]

        lpeak=np.searchsorted(data['lowEflux'][lrange[0]:lrange[1]],np.max(data['lowEflux'][lrange[0]:lrange[1]])) #gives index of peak IN DESIRED RANGE
        lpeaktime=data['taxis'][lidx[0]+lpeak] #gives time at actual peak
        lpeakflux=data['lowEflux'][lidx[0]+lpeak] #gives flux at actual peak

        binshort=self.Datetimes.messenger_peak_bin['short']
        sidx=[]
        for i in range(0,len(data['taxis'])-1):
            if data['taxis'][i] > binshort[0] and data['taxis'][i] < binshort[1]:
                sidx.append(i) #the index of the value
        #sidx=np.where(data['taxis'] > binshort[0] and data['taxis'] < binshort[1])
        srange=[sidx[0],sidx[len(sidx)-1]]
        if srange[0]==srange[1]: #don't allow this!
            srange=[sidx[0]-1,sidx[0]]
        speak=np.searchsorted(data['highEflux'][srange[0]:srange[1]],np.max(data['highEflux'][srange[0]:srange[1]])) #gives index of peak IN DESIRED RANGE
        speaktime=data['taxis'][sidx[0]+speak] #gives time at actual peak
        speakflux=data['highEflux'][sidx[0]+speak] #gives flux at actual peak

        #store stuff
        self.Datetimes.Messenger_peak = speaktime
        self.Datetimes.Messenger_peak_long = lpeaktime
        self.Properties.GOES_GOES_long=lpeakflux
        self.Properties.GOES_GOES=speakflux

        #pickle it up if requested, save filename and location to object
        if save:
            import pickle
            fdir=self.Files.dir+self.Files.folders['lightcurves']
            filename= '/'+str(self.ID)+'_GOES_lc.p' #something
            self.Files.Raw['goes']=filename
            pickle.dump(data,open(fdir+filename,'wb'))
            print 'File '+filename+' saved in '+fdir

        os.chdir(self.Files.dir)
        #self.Datetimes.convert2datetime()

        if ret:
            return data

    def get_Messenger_lightcurve(self, filename=False, raw=False, save=False):
        '''Get data for this flare's Messenger lightcurve, and adjust to counts/s if specified (raw = False). Save in new file if requested.'''
        os.chdir(self.Files.dir + self.Files.folders['lightcurves'])
        if self.Files.Raw['messenger'].endswith('.p'):
            Mdata= pickle.load(self.Files.Raw['messenger']) #restore and be done. should be a mySpectrogram object
            #os.chdir(self.Files.dir)
            #return Mdata#think I probably have to assign a variable name...

        elif self.Files.Raw['spectrogram'] == '': #get it from the sav file
            if not filename: filename='rerun_sgrams.sav'
            from scipy.io import readsav
            os.chidr('../')
            d=readsav(filename, python_dict=True)
            Mdata=d['data']['mdata']
            for j in len(data['len']):
                if self.Datetimes.Messenger_peak.date() == dt.strptime(Mdata['UT'][j][0],'%d-%b-%Y %H:%M:%S.000').date():
                    i=j #use this index
                    break
            Mdata=Mdata[i]

        def counts_ps(Mdata,n):
            '''Adjust messenger data to counts/s'''
            #get the time bin for each data point
            tim=Mdata['taxis'][0][n]
            mlen=Mdata['len'][0][n]
            nflares=np.shape(Mdata['taxis'][0])[0]/2
            Mtim,cps1,cps2=[],[],[]
            for t in tim: Mtim.append(dt.strptime(t,'%d-%b-%Y %H:%M:%S.%f')) #fix messenger times to datetimes
            M1=Mdata['phigh'][0][n][0:mlen-1]
            M2=Mdata['phigh'][0][n+nflares-1][0:mlen-1]

            for i in range(mlen-1):
                tbin=Mtim[i+1]-Mtim[i] #timedelta object in seconds
                cps1.append(M1[i]/tbin.total_seconds())
                cps2.append(M2[i]/tbin.total_seconds())
            return cps1,cps2

        if not Raw:
            #adjust counts
            cps1,cps2=counts_ps(Mdata,i)

        if save:
            os.chdir(self.Files.dir + self.Files.folders['spectrograms'])
            newfname=str(self.ID)+'Messengerlc.p'
            pickle.dump(Mdata,open(newfname,'wb'))
            self.Files.Raw['spectrogram']=newfname

        os.chdir(self.Files.dir)
        return Mdata

    def get_GOES_lightcurve(self, filename=False, save=False):
        '''Get data for this flare's GOES lightcurve. Save in new file if requested.'''
        os.chdir(self.Files.dir + self.Files.folders['lightcurves'])
        if self.Files.Raw['goes'].endswith('.p'):
            Gdata=pickle.load(self.Files.Raw['goes']) #restore and be done. should be a mySpectrogram object
            #os.chdir(self.Files.dir)
            #return Gdata#think I probably have to assign a variable name...

        elif self.Files.Raw['spectrogram'] == '': #get it from the sav file
            if not filename: filename='rerun_sgrams.sav'
            from scipy.io import readsav
            os.chidr('../')
            d=readsav(filename, python_dict=True)
            Gdata=d['data']['gdata']
            for j in len(data['len']):
                if self.Datetimes.Messenger_peak.date() == dt.strptime(Gdata['UT'][j][0],'%d-%b-%Y %H:%M:%S.000').date():
                    i=j #use this index
                    break
            Gdata={}

        if save:
            os.chdir(self.Files.dir + self.Files.folders['spectrograms'])
            newfname=str(self.ID)+'GOESlc.p'
            pickle.dump(a,open(newfname,'wb'))
            self.Files.Raw['spectrogram']=newfname
        os.chdir(self.Files.dir)
        return Gdata

    def get_spectrogram(self, filename=False, save=False):
        '''Get data for this flare's RHESSI spectrogram. Save in new file if requested.'''
        import mySpectrogram as s
        os.chdir(self.Files.dir + self.Files.folders['spectrograms'])
        if self.Files.Raw['spectrogram'].endswith('.p'):
            pickle.load(self.Files.Raw['spectrogram']) #restore and be done. should be a mySpectrogram object
            os.chdir(self.Files.dir)
            return #think I probably have to assign a variable name...
        if self.Files.Raw['spectrogram'] == '': #get it from the sav file
            if not filename: filename='rerun_sgrams.sav'
            from scipy.io import readsav
            os.chdir('../')
            d=readsav(filename, python_dict=True)
            data={'UT':0.,'rate':0.,'erate':0.,'ltime':0.,'len':0.}
            data['UT'] = d['data']['UT']
            data['rate']=d['data']['rate']
            data['erate']=d['data']['erate']
            data['len']=d['data']['len']
            ebins=d['ebins']
            eaxis=ebins[0:-1]
            last=int(data['len'][0]-1)
            time_axis=np.arange(0,last+1)*4.#data[0]['ltime'] #ndarray of time offsets,1D #need the 4 for RHESSI time interval
            for j in len(data['len']):
                if self.Datetimes.Messenger_peak.date() == dt.strptime(data['UT'][j][0],'%d-%b-%Y %H:%M:%S.000').date():
                    i=j #use this index
                    break
            start=dt.strptime(data['UT'][i][0],'%d-%b-%Y %H:%M:%S.000')
            try:
                end= dt.strptime(data['UT'][i][last],'%d-%b-%Y %H:%M:%S.000') #datetime #might want to modify this to be where the data =0
            except ValueError:
                end = start + datetime.timedelta(seconds=time_axis[-1])
                drate=np.transpose(np.log10(data['rate'][i][0:last+1])) #get rid of zeros before taking log10?
                drate[drate == -np.inf] = 0.0
            for n,col in enumerate(drate.T):
                if all([c ==0.0 for c in col]):
                    drate[:,n] = np.nan
            a=s.Spectrogram(data=drate,time_axis=time_axis,freq_axis=eaxis,start=start,end=end,instruments=['RHESSI'],t_label='',f_label='Energy (keV)',c_label='log(counts/cm^2 s)')
            if save:
                os.chdir(self.Files.dir + self.Files.folders['spectrograms'])
                newfname=str(self.ID)+'spectrogram.p'
                pickle.dump(a,open(newfname,'wb'))
                self.Files.Raw['spectrogram']=newfname
        os.chdir(self.Files.dir)
        return a

    def get_ospex_fit(self, filename=False, save=False):
        '''Get data for this flare's OSPEX fit to the Messenger lightcurve. Save in new file if requested. These are IDL sav files..'''

    def _query(self, path=False, AIA=False, STEREO=True, tint= False,wave=False, timedelay = 30,both=False):
        '''Internal method to query VSO database for AIA or STEREO data and download it. Option to set time delay for stereo observations (in minutes). '''
        from sunpy.net import Fido, attrs as a#vso
        from datetime import timedelta as td
        import sunpy.map

        #vc = vso.VSOClient()
        if type(self.Datetimes.Messenger_peak) != dt: self.Datetimes.convert2datetime()
        if AIA:
            instr= a.Instrument('AIA')
            sample = a.Sample(24 * u.hour)
            if wave == '171': wave = a.Wavelength(16.9 * u.nm, 17.2 * u.nm)
            elif wave == '193': wave = a.Wavelength(19.1 * u.nm, 19.45 * u.nm)
            elif wave == '304': wave = a.Wavelength(30 * u.nm, 31 * u.nm)
            elif wave == '131': wave = a.Wavelength(13.0 * u.nm, 13.2 * u.nm)

            if not tint:
                time = a.Time(dt.strftime(self.Datetimes.AIA_start,'%Y-%m-%dT%H:%M:%S'),dt.strftime(self.Datetimes.AIA_end,'%Y-%m-%dT%H:%M:%S')) #should have data to within 1 s #used to be Datetimes.Messenger_peak

        if STEREO:
            if not both:
                source=a.vso.Source('STEREO_'+self.Observation.STEREO) #this should be a string
            instr= a.Instrument('EUVI')
            if wave == '195': wave = a.Wavelength(19.1 * u.nm, 19.45 * u.nm) #193
            elif wave == '171': wave = a.Wavelength(16.9 * u.nm, 17.2 * u.nm) #193
            if not tint:
                time = a.Time(dt.strftime(self.Datetimes.Messenger_peak+td(minutes=timedelay),'%Y-%m-%dT%H:%M:%S'),dt.strftime(self.Datetimes.Messenger_peak +td(minutes=timedelay+15),'%Y-%m-%dT%H:%M:%S')) #have to look 15 minutes after designated timedelay?
        if tint:
             time=a.Time(tint[0],tint[1]) #these are datetimes
        if not both:
            res=Fido.search(wave, source, time, instr)
        else:
            res=Fido.search(wave, time, instr)

        if not path: files = Fido.fetch(res,path='/Users/wheatley/Documents/Solar/occulted_flares/data/stereo-aia/{file}')
        else: files = Fido.fetch(res,path=path+'{file}').wait()
        #put filenames in Files object
        if AIA:
            self.Files.Raw['aia'] = files[0]
        if STEREO:
            self.Files.Raw['stereo'] = files #take all of them in a list

    def get_AIA(self, filename=False, save=False, wave='304'):
        '''Get data for this flare's AIA map of specified wavelength. Query the database if need be. Save in new file if requested.'''
        from astropy.coordinates import SkyCoord

        path=self.Files.dir + self.Files.folders['stereo-aia']+'/'
        os.chdir(path)
        #first check the if there is a local file with filename stored in Files object.
        fname=self.Files.Raw['aia']
        if type(fname) == list:
            try:
                fname=[f for f in fname if wave in f][0]
            except IndexError:
                fname=''
        if fname.endswith('.p'):
            maps= pickle.load(open(fname,'rb')) #if it's a pickle just restore it. (do I need to assign it to maps or not?)
            os.chdir(self.Files.dir)
            return maps
        if fname == '' or not os.path.isfile(path+fname) or wave not in fname: #if not, query database. Now filename should be in the correct place in Files.
            self._query(AIA=True, STEREO=False, wave=wave)

        #convert to maps
        files=self.Files.Raw['aia']
        if type(files) == list:
            files=files[0]
        if files !=[]:
            smap= sunpy.map.Map(files)
            maps = {smap.instrument: smap.submap(SkyCoord((-1100, 1100) * u.arcsec, (-1100, 1100) * u.arcsec,frame=smap.coordinate_frame))} #this could be empty if files is empty
        else: print 'No files found and no maps made!'

        if save: #pickle it?
            newfname=files[files.rfind('/')+1:files.rfind('.')]+'.p'
            pickle.dump(maps,open(newfname,'wb'))
            self.Files.Raw['aia']=newfname
        os.chdir(self.Files.dir)

        self.Observation.Rsun_AIA=maps[maps.keys()[0]].rsun_obs #what are the units?
        return maps

    def get_STEREO(self,tint=False, wave='195',filename=False, save=False):
        '''Get data for this flare's STEREO map. Query the database if need be. Save in new file if requested.'''
        from astropy.coordinates import SkyCoord

        path=self.Files.dir + self.Files.folders['stereo-aia']+'/'
        os.chdir(path)
        #first check the if there is a local file with filename stored in Files object.
        fname=self.Files.Raw['stereo']
        try:
            if fname.endswith('.p'):
                maps = pickle.load(open(fname,'rb')) #if it's a pickle just restore it. (do I need to assign it to maps or not?)
                os.chdir(self.Files.dir)
                return maps
            if fname == '' or not os.path.isfile(path+fname) or fname==[]: #if not, query database. Now filename should be in the correct place in Files.
                self._query(tint=tint,wave=wave)
        except AttributeError: #it's a list of fits filenames
            pass

        #convert to maps
        files=self.Files.Raw['stereo']
        if filename:
            files=[]
            files.append(filename)
        maps=[]
        if files !=[]:
            #smap= sunpy.map.Map(files)
            if len(files) == 1:
                f=sunpy.map.Map(files[0])
                maps.append({f.instrument: f.submap(SkyCoord((-1100, 1100) * u.arcsec, (-1100, 1100) * u.arcsec,frame=f.coordinate_frame))})
            else:
                for f in sunpy.map.Map(files):
                    maps.append({f.instrument: f.submap(SkyCoord((-1100, 1100) * u.arcsec, (-1100, 1100) * u.arcsec,frame=f.coordinate_frame))}) #this could be empty if files is empty
                self.extract_stereo_times()
        else: print 'No files found and no maps made!'

        if save: #pickle it?
            os.chdir(path)
            newfname=files[0][files[0].rfind('/')+1:files[0].rfind('.')]+'.p'
            pickle.dump(maps,open(newfname,'wb'))
            self.Files.Raw['stereo']=newfname
        os.chdir(self.Files.dir)

        return maps

    def do_clean(self,erange):
        '''Do clean using method do_clean from full_disk_image.pro'''
        import pidly
        idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
        idl('search_network,/enable')
        idl('.compile full_disk_image.pro')
        self.Datetimes.convert2string()
        idl.time_int=[self.Datetimes.Messenger_peak_time,self.Datetimes.Messenger_peak_end]
        idl.pix_size=[1,1]
        idl.image_dim=[128,128]
        ldet=list(self.Observation.Det[1:-1])
        idl.det=list(int(l) for l in ldet)
        idl.xyoffset=self.Properties.source_pos
        idl.erange=erange
        idl.outfilename=str(self.ID)+str(erange[0])+'-'+str(erange[1])+'plot.fits'
        #idl.outsavename=str(self.ID)+'clean.sav'
        idl('im=do_clean(time_int,pix_size,det,image_dim,xyoffset,erange,outfilename=outfilename,quiet=quiet)')
        #idl('save,o,filename=outsavname')

    def clean_all(self,eranges): #eranges=[[4,9],[12,18],[18,30],[30,80]]
        for erange in eranges:
            self.do_clean(erange)

    def download_peak_stereo(self,n=3):
        '''Download the n files from the best suited STEREO around the Messeneger peak'''
        from datetime import timedelta as td
        start_time=self.Datetimes.Messenger_peak - td(minutes=6*n/2)
        end_time=self.Datetimes.Messenger_peak + td(minutes=6*n/2)
        tint=[start_time,end_time] #n files around Messenger peak. datetime range.
        #smap=f.get_STEREO(tint=tint)
        self._query(tint=tint,wave='195',both=True)

    def download_preflare_stereo(self,overwrite=True,peak=False):
        '''Download the preflare image - one hour before the peak'''
        from datetime import timedelta as td

        #determine download constraints
        start_time=self.Datetimes.Messenger_peak - td(minutes=65)
        end_time=self.Datetimes.Messenger_peak - td(minutes=55)
        tint=[start_time,end_time] #n files around Messenger peak. datetime range.
        try:
            pffile=self.Files.prepped['stereo-preflare']
            if pffile[-5]!=self.Observation.STEREO:
                overwrite=True
        except AttributeError:
            pass
        if overwrite:
            self._query(tint=tint,wave='195')
            try:
                self.Files.prepped['stereo-preflare']=self.Files.Raw['stereo'][0]
                pfstr=self.Files.prepped['stereo-preflare']
                pftime=dt.strptime(pfstr[pfstr.rfind('/')+1:pfstr.rfind('_')],'%Y%m%d_%H%M%S')
                peaktime=dt.strptime(self.Files.prepped['stereo-peak'][-25:-10],'%Y%m%d_%H%M%S')
                self.Datetimes.peakdiff = (peaktime-pftime).total_seconds()
            except IndexError:
                pass
            if peak:
                try:
                    pfile=self.Files.prepped['stereo-peak']
                    tint=dt.strptime(pfile[:-10],'%Y%m%d_%H%M%S')
                    self._query(tint=[tint-td(minutes=2.5),tint+td(minutes=2.5)],wave='195')
                except AttributeError:
                   pass



    def find_stereo_peak(self,prep=False,minutes=25):
        '''Find the STEREO image with the peak of the flare in it. Could be from A or B.'''
        import sunpy.map
        import pidly
        self.update_stereo_filenames(minutes=minutes)
        ff=self.Files.Raw['stereo']
        if prep:
            preplist=[f0 for f0 in ff] # if f0[-4].islower()]
        else:
            try:
                foo=self.Files.prepped['stereo-l1']
            except AttributeError:
                    self.Files.prepped={'stereo-l1':[],'stereo-peak':'','stereo-preflare':''} #will this transfer back to the original object? I think so

        #now run secchi_prep
        if prep:
            idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
            if ff !=[]: #skip
                idl('flist',ff)
                idl('secchi_prep,flist,/write_fts')
                pl=[l[0:-9]+'1'+l[-8:-5]+l[-5].upper()+l[-4:] for l in ff if not l[-5].isdigit()] #new naming convention
                #should delete all the .0.fts from the directory first
                os.chdir(self.Files.dir + self.Files.folders['stereo-aia'])
                plmaps=sunpy.map.Map(pl,cube=True)
                dav=[plm.meta['dataavg'] for plm in plmaps]
                davmax=dav.index(np.max(dav))
                self.Files.prepped['stereo-l1']=pl
                self.Files.prepped['stereo-peak']=pl[davmax]
                self.Files.prepped['stereo-preflare']=pl[0]
                if pl[davmax][-5] != self.Observation.STEREO:
                    self.Observation.STEREO=pl[davmax][-5]
            idl.close()
        else:
            if ff !=[]: #skip
                #pl=[l[0:-9]+'1'+ll[-8:-5]+l[-5].upper()+l[-4:] for l in ll if not l[-5].isdigit()] #new naming convention
                #should delete all the .0.fts from the directory first
                pl=self.Files.prepped['stereo-l1']
                if pl != [] and pl != '':
                    self.Datetimes.stereo_times=[dt.strptime(sf[:-10],'%Y%m%d_%H%M%S') for sf in pl]
                    difftimes=[np.abs((pftime - self.Datetimes.Messenger_peak)).seconds for pftime in self.Datetimes.stereo_times]
                    darr=np.array(difftimes)
                    try:
                        goodtimes=np.where(darr <=900.) #files within 15 minutes
                        goodpl=[p for i,p in enumerate(pl) if i in goodtimes[0]]
                        os.chdir(self.Files.dir + self.Files.folders['stereo-aia'])
                        plmaps=sunpy.map.Map(goodpl,cube=True)
                        dav=[plm.meta['dataavg'] for plm in plmaps]
                        davmax=dav.index(np.max(dav))
                        #self.list[i+nlow].Files.prepped['stereo-l1']=pl
                        self.Files.prepped['stereo-peak']=pl[davmax]
                        #check the time difference! if more than 15 minutes, just pick the closest one!
                        self.Datetimes.stereo_peak=dt.strptime(self.Files.prepped['stereo-peak'][:-10],'%Y%m%d_%H%M%S')
                        self.Datetimes.peakdiff=np.abs((self.Datetimes.Messenger_peak - self.Datetimes.stereo_peak)).seconds
                        #self.Datetimes.stereo_times=[dt.strptime(sf[:-10],'%Y%m%d_%H%M%S') for sf in self.Files.prepped['stereo-l1']]
                        #self.Files.prepped['stereo-preflare']=pl[0]
                        if pl[davmax][-5] != self.Observation.STEREO:
                            self.Observation.STEREO=pl[davmax][-5]
                    except (TypeError, ValueError):
                        self.Datetimes.peakdiff=100000.


    #############################################  MISC CALCULATION METHODS  #######################################################

    def angle_closest_stereo(self):
        '''Find the angle between Messenger and the closest STEREO spacecraft...might need to re-calc messenger position since the Angle parameter is a bit vague'''
        import stereo_angle as sa
        #if ...
        heexa,heeya,anglea=sa.get_stereo_angle(self.Datetimes.Messenger_peak, stereo='A')
        heexb,heeyb,angleb=sa.get_stereo_angle(self.Datetimes.Messenger_peak, stereo='B')
        mheex,mheey=self.get_messenger_coords()
        dma=np.sqrt((heexa-mheex)**2 + (heeya-mheey)**2)
        dmb=np.sqrt((heexb-mheex)**2 + (heeyb-mheey)**2)
        if dma < dmb:
            closest_stereo='A'
            angles=anglea
        else:
            closest_stereo='B'
            angles=angleb
        angle=angles-self.Observation.Angle
        self.Observation.closest_stereo=closest_stereo #or does it have to be sett_attr?
        self.Observation.msangle=angle
        return closest_stereo, angle

    def dump_messenger_coords(self):
        '''Get messenger coordinates in heliocentric longditude and latitude. For use in angle_closest_stereo'''
        import pidly
        idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
        datestr=dt.strftime(self.Datetimes.Messenger_peak.date(),'%d-%b-%Y')
        if self.Properties.source_pos != '' and self.Properties.source_pos != '[0,0]':
            pos=[float(self.Properties.source_pos[1:self.Properties.source_pos.find(':')]),float(self.Properties.source_pos[self.Properties.source_pos.find(':')+1:-1])] #format?
        else:
            pos=[0,0]
        idl.pos=pos
        idl.datestr=datestr
        angle = idl('messenger_flare_angle(pos, datestr,mess_helio_longlat)')
        idl('dd=anytim(mess_helio_longlat.date,/vms)')
        datet=idl.dd
        idl('ll=mess_helio_longlat.xyz')
        longlat=idl.ll
        idl('rad=mess_helio_longlat.d_km')
        radius=idl.rad
        idl.close()
        coords={'datet':datet,'radius':radius,'longlat':longlat}
        pickle.dump(coords,open('messenger_helio_longlat.p','wb'))
        #return datet,longlat

    def get_messenger_coords(self):
        '''Get messenger coordinates in heliocentric longditude and latitude. For use in angle_closest_stereo'''
        cdict=pickle.load(open('messenger_helio_longlat.p','rb'))
        dates=cdict['datet']
        radius=cdict['radius']
        coords=cdict['longlat']
        fdate=self.Datetimes.Messenger_peak.date()
        #search the dates for the correct index
        for i,d in enumerate(dates):
            idldate=dt.strptime(d,'%d-%b-%Y %H:%M:%S.000').date()
            delt=fdate-idldate
            if delt.days==0: #use timedeltas
                idx=i
                break
        #search the coords for the corresponding x and y
        mheex=coords[idx][0]*radius[idx] #convert to km
        mheey=coords[idx][1]*radius[idx]
        return mheex,mheey

    def update_stereo_filenames(self,minutes=15,clean=False):
        ''' search the stereo-aia folder for all files corresponding to a particular flare'''
        from datetime import timedelta as td
        path=self.Files.dir + self.Files.folders['stereo-aia']
        os.chdir(path)
        self.Datetimes.convert2datetime()
        start_time=self.Datetimes.Messenger_peak - td(minutes=minutes) #should get all the ones from the day then narrow down
        start_string=dt.strftime(start_time,'%Y%m%d_%H%M%S')[:-4]
        #end_time=self.Datetimes.Messenger_peak + td(minutes=minutes)
        end_hour="{:02}".format(int(start_string[-2:])+1)
        end_string=start_string[:-2]+end_hour
        pf_time=self.Datetimes.Messenger_peak - td(hours=1)
        pf_string=dt.strftime(pf_time,'%Y%m%d_%H%M%S')[:-4]

        rawglob= glob.glob(start_string+'*_n*.fts')
        preppedglob= glob.glob(start_string+'*_1*.fts')
        which_stereo=self.Observation.STEREO
        pfglob=glob.glob(pf_string+'*_14eu'+which_stereo+'.fts') #need to make sure this is from the correct STEREO
        if pfglob == []: #need to give it a few more minutes buffer
            pf_timea=pf_time-td(minutes=2.5)
            pf_timeb=pf_time+td(minutes=2.5)
            pfstringa=dt.strftime(pf_timea,'%Y%m%d_%H%M%S')[:-4]
            pfstringb=dt.strftime(pf_timeb,'%Y%m%d_%H%M%S')[:-4]
            pfglob =glob.glob(pfstringa+'*_14eu'+which_stereo+'.fts')
            pfglob.extend(glob.glob(pfstringb+'*_14eu'+which_stereo+'.fts'))
        rawglob.extend(glob.glob(end_string+'*_n*.fts'))
        preppedglob.extend(glob.glob(end_string+'*_1*.fts'))
        #if len(rawglob) > 10:
        #    rawglob=rawglob[:10]
        #if len(preppedglob) > 10:
        #    preppedglob=preppedglob[:10]

        if rawglob !=[]:
            self.Datetimes.stereo_times=[]
            for f in rawglob:
                timestr=f[f.find('_')+1:f.rfind('_')]
                self.Datetimes.stereo_times.append(dt.strptime(f[:f.find('_')+1]+timestr,'%Y%m%d_%H%M%S')) #parse the string and extract dateti

            self.Files.Raw['stereo']=rawglob#put the filenames in Raw
        if preppedglob !=[]:
            try:
                self.Files.prepped['stereo-l1']=preppedglob#put the filenames in Raw
                try:
                    self.Files.prepped['stereo-preflare']=pfglob[0]#put the filenames in Raw
                    self.Datetimes.peakdiff=(dt.strptime(self.Files.prepped['stereo-peak'][-25:-10],'%Y%m%d_%H%M%S')-dt.strptime(self.Files.prepped['stereo-preflare'][-25:-10],'%Y%m%d_%H%M%S')).total_seconds()
                except IndexError:
                    pass
            except AttributeError:
                self.Files.prepped={'stereo-l1':preppedglob,'stereo-peak':'','stereo-preflare':preppedglob[0]}

            if clean: #delete all the files that aren't the peak or the preflare image
                print 'clean not enabled yet'
        else:
             self.Files.prepped={'stereo-l1':'','stereo-peak':'','stereo-preflare':''}


    def extract_stereo_times(self):
        '''Put correct dates and times from downloaded stereo files into the datetimes.stereo_files list (except these are not actually the correct times, just the filenames. Can fix this later though, since right now the filename is more important)'''
        path=self.Files.dir + self.Files.folders['stereo-aia']
        os.chdir(path)
        self.Datetimes.convert2datetime()
        try:
            fname= dt.strftime(self.Datetimes.stereo_start_time,'%Y%m%d')
        except TypeError:
            fname= self.Datetimes.Messenger_peak[0:self.Datetimes.stereo_start_time.find(' ')]
        ffiles=glob.glob(fname+'*eub.fts')
        if glob.glob(fname+'*eua.fts') !=[]:
            ffiles=[]
            ffiles =glob.glob(fname+'*eua.fts')
        self.Datetimes.stereo_times=[]
        for f in ffiles:
            timestr=f[f.find('_')+1:f.rfind('_')]
            self.Datetimes.stereo_times.append(dt.strptime(fname+timestr,'%Y%m%d%H%M%S')) #parse the string and extract datetimes
        self.Files.Raw['stereo']=ffiles#put the filenames in Raw
        self.Files.Raw['aia']=glob.glob('aia_lev1*'+dt.strftime(self.Datetimes.Messenger_peak,'%Y_%m_%d')+'*.fits') #put the filenames in Raw
        os.chdir(self.Files.dir)

    def convert_coords_aia2stereo(self,map_aia=True,map_stereo=True,quiet=True,wave='304'):
        '''convert Properties.source_pos_disk to Properties.source_pos_STEREO for a given AIA and STEREO map'''
        from astropy.coordinates import SkyCoord
        #import astropy.wcs
        import sunpy.coordinates
        import sunpy.coordinates.wcs_utils

        if self.Properties.source_pos_disk == '': #think it should be this if empty
            coords=self.Properties.source_pos #is this in arcsec? how to check?
        else: coords=self.Properties.source_pos_disk
        #coords= [int(coords[:coords.find(',')]),int(coords[coords.find(',')+1:])]
        if not map_aia: #get the map
            if self.Files.Raw['aia'] !='' and glob.glob(self.Files.dir+ self.Files.folders['stereo-aia']+'/'+self.Files.Raw['aia']) !=[]:#first see if the pickle file is there
                os.chdir(self.Files.dir+ self.Files.folders['stereo-aia'])
                map_aia=pickle.load(open(self.Files.Raw['aia'],'rb'))
                os.chdir(self.Files.dir)
            else: #download stuff and pickle it
                map_aia=self.get_AIA(filename=False, save=True, wave=wave)

        if not map_stereo: #get the map. Default is to take the first in the list
            if self.Files.Raw['stereo'] !='' and glob.glob(self.Files.dir+ self.Files.folders['stereo-aia']+'/'+self.Files.Raw['stereo']) !=[]:#first see if the pickle file is there
                os.chdir(self.Files.dir+ self.Files.folders['stereo-aia'])
                map_stereo=pickle.load(open(self.Files.Raw['stereo'],'rb'))
                os.chdir(self.Files.dir)
            else: #download stuff and pickle it
                map_stereo=self.get_STEREO(filename=False, save=True)
                if type(map_stereo) == list:
                    map_stereo=map_stereo[0]

        hpc_aia = SkyCoord(coords[0]*u.arcsec,coords[1]*u.arcsec, frame=map_aia[map_aia.keys()[0]].coordinate_frame)
        if not quiet: print(hpc_aia)

        hgs = hpc_aia.transform_to('heliographic_stonyhurst')
        hgs_val=np.isnan(hgs.lon.value) + np.isnan(hgs.lat.value)
        while hgs_val == 1: #while there are nan's
            coordstr=raw_input('Please enter integer coordinates on the AIA disk. Format: xxx,yyy')
            coords= [int(coordstr[:coordstr.find(',')]),int(coordstr[coordstr.find(',')+1:])]#format the string
            hpc_aia = SkyCoord(coords[0]*u.arcsec,coords[1]*u.arcsec, frame=map_aia[map_aia.keys()[0]].coordinate_frame)
            hgs = hpc_aia.transform_to('heliographic_stonyhurst')
            hgs_val= np.isnan(hgs.lon.value)+ np.isnan(hgs.lat.value)
        if not quiet: print(hgs)

        hgs.D0 = map_stereo['SECCHI'].dsun
        hgs.L0 = map_stereo['SECCHI'].heliographic_longitude
        hgs.B0 = map_stereo['SECCHI'].heliographic_latitude

        hpc_B = hgs.transform_to('helioprojective')
        if not quiet: print(hpc_B)

        self.Properties.source_pos_STEREO = [hpc_B.Tx.value,hpc_B.Ty.value] #is this right?

    def calc_projection_effect(self):
        '''calculate the projection effect using formula in Krucker-St.Hilaire 2015'''
        #get the Sun for that date
        s_ang=sunpy.sun.solar_semidiameter_angular_size(self.Datetimes.current_map_time)
        self.Observation.Rsun_ang=s_ang #but convert this to km...
        if not self.Observation.Rsun_AIA:
            Rsun=s_ang*757. #1arcsec=715km (how exact is this though)?
        else:
            Rsun=self.Observation.Rsun_AIA.value*757.
        arcsec2deg=s_ang.value/90. #since s_ang is the radius not diameter, subtends an angle of 90 degrees not 180
        angd=np.radians(self.Properties.e2['distance2limb'].value/arcsec2deg) #convert distance to degrees
        hkm=Rsun*(1-np.cos(angd))/np.cos(angd)
        angda=np.radians((self.Properties.e2['distance2limb'].value-30.)/arcsec2deg)
        hkma=Rsun*(1-np.cos(angda))/np.cos(angda)
        print 'Projection effect of ',str(hkm), ' km, ', str(self.Properties.e2['distance2limb'].value/arcsec2deg) , ' degrees behind limb'
        print 'Adjusted projection effect of ',str(hkma), ' km', str((self.Properties.e2['distance2limb'].value-30.)/arcsec2deg) , ' degrees behind limb'
        self.Properties.e5={'hkm':hkm,'angd':angd,'hkma':hkma,'angda':angda}

    def get_stereo_flare_loc(self):
        '''Aggressively Gaussian smooth the full disk image, take the center of the brightest contour'''
        from scipy import ndimage as ndi
        from find_hessi_centroids import center_of_mass
        os.chdir(self.Files.dir + self.Files.folders['stereo-aia'])
        peakmap=sunpy.map.Map(self.Files.prepped['stereo-peak'])
        gmap=sunpy.map.Map(ndi.gaussian_filter(peakmap.data,7),peakmap.meta)
        cc=plt.contour(gmap.data,levels=[.85*np.max(gmap.data)],frame=gmap.coordinate_frame)
        try:
            cmass=center_of_mass(cc.allsegs[0][0])
            hpj=gmap.pixel_to_world(cmass[0]*u.pixel,cmass[1]*u.pixel,origin=0)
            self.Properties.source_pos_STEREO=[hpj.Tx.value,hpj.Ty.value]
        except IndexError: #contour level probably too high?
            pass
        os.chdir(self.Files.dir)

    def calc_nitta_flux(self,recalc=False):
        '''calculate Nitta flux via total of background-subtracted flare image
        Based on my scripts, I ran secchi_prep for the flare and pre-flare images and got two fluxes by averaging them.
        The unit of the output from secchi_prep is DN/s.
        Subtract the preflare flux from the flare flux.
        EL: Pretty sure he still means the sum, because this is what gives the correct result for the example of 2010-08-31'''
        os.chdir(self.Files.dir + self.Files.folders['stereo-aia'])
        try:
            peakmap=sunpy.map.Map(self.Files.prepped['stereo-peak'])
            bgmap=sunpy.map.Map(self.Files.prepped['stereo-preflare'])
            #self.Properties.peakav=np.mean(peakmap.data)
            #self.Properties.bgav=np.mean(bgmap.data)
            self.Properties.peaksum=np.sum(peakmap.data)
            self.Properties.bgsum=np.sum(bgmap.data)
            #self.Properties.Nittaflux=self.Properties.peakav-self.Properties.bgav
            self.Properties.Nittaflux=self.Properties.peaksum-self.Properties.bgsum
        except (AttributeError, ValueError):
            if recalc:
               self.Properties.Nittaflux=0.

    def MSobs(self):
        '''is the flare observed by both Messenger and STEREO? (first round: roughly, don't take into account flare position on STEREO disk). also neglect Z-coord'''
        self.Observation.MSvis = ''
        try:
            mlonlat=self.Observation.mlonlat
            if self.Observation.STEREO == 'A':
                slonlat=self.Observation.Alonlat
            elif self.Observation.STEREO == 'B':
                slonlat=self.Observation.Blonlat

            #if type(mlonlat) == np.array and type(slonlat) == np.array: #find the angle
            sang=np.rad2deg(np.arctan(slonlat)[0])
            mang=np.rad2deg(np.arctan(mlonlat)[0])
            abs_ang=np.max([sang,mang])-np.min([sang,mang])
            #print abs_ang
            if abs_ang <=60:
                self.Observation.MSvis=True
            else:
                self.Observation.MSvis = False

        except AttributeError:
            pass

    def GSobs(self):
        '''is the flare observed by both GOES and STEREO? (first round: roughly, don't take into account flare position on STEREO disk)'''
        self.Observation.GSvis = ''
        try:
            if self.Observation.STEREO == 'A':
                slonlat=self.Observation.Alonlat
            elif self.Observation.STEREO == 'B':
                slonlat=self.Observation.Blonlat

            #if type(slonlat == np.array): #find the angle
            sang=np.rad2deg(np.arctan(slonlat)[0])
            if sang <=60:
                self.Observation.GSvis=True
            else:
                self.Observation.GSvis=False
        except AttributeError:
            pass

    def MSobs_with_sloc(self):
        import theta_obs as th
        if self.Observation.STEREO == 'A':
            sx,sy=self.Observation.Alonlat
        else:
            sx,sy=self.Observation.Blonlat
        mx,my=self.Observation.mlonlat
        spx,spy=self.Properties.source_pos_STEREO
        v2=(np.deg2rad(sy),np.deg2rad(sx))
        vec1=(mx,my)
        vec2=(np.cos(v2[0])*np.cos(v2[1]),np.cos(v2[0])*np.sin(v2[1]))
        smap=sunpy.map.Map(self.Files.dir+'/'+self.Files.folders['stereo-aia']+'/'+self.Files.prepped['stereo-peak'])
        stereo_loc=SkyCoord(spx*u.arcsec,spy*u.arcsec,frame=smap.coordinate_frame)
        s1E,s1W=th.senkret(vec1)
        s2E,s2W=th.senkret(vec2)
        thobs,thex,vlimb,vElimb,vfov=th.theta_obs(vec1,s1E,s1W,vec2,s2E,s2W)
        stereo_lon=th.stereo_loc_to_lon(stereo_loc)
        joint_obs=th.is_flare_observed(stereo_lon,vElimb,vfov,vlimb)
        self.Properties.MSobs_with_sloc=joint_obs

    def GSobs_with_sloc(self):
        import theta_obs as th
        if self.Observation.STEREO == 'A':
            sx,sy=self.Observation.Alonlat
        else:
            sx,sy=self.Observation.Blonlat
        mx,my=(1.0,0.000001) # Earth. arctan doesn't do well with zeros
        spx,spy=self.Properties.source_pos_STEREO
        v2=(np.deg2rad(sy),np.deg2rad(sx))
        vec1=(mx,my)
        vec2=(np.cos(v2[0])*np.cos(v2[1]),np.cos(v2[0])*np.sin(v2[1]))
        smap=sunpy.map.Map(self.Files.dir+'/'+self.Files.folders['stereo-aia']+'/'+self.Files.prepped['stereo-peak'])
        stereo_loc=SkyCoord(spx*u.arcsec,spy*u.arcsec,frame=smap.coordinate_frame)
        s1E,s1W=th.senkret(vec1)
        s2E,s2W=th.senkret(vec2)
        thobs,thex,vlimb,vElimb,vfov=th.theta_obs(vec1,s1E,s1W,vec2,s2E,s2W)
        stereo_lon=th.stereo_loc_to_lon(stereo_loc)
        joint_obs=th.is_flare_observed(stereo_lon,vElimb,vfov,vlimb)
        self.Properties.GSobs_with_sloc=joint_obs

    ################################################# INDIVIDUAL PLOT METHODS  ###########################################################
    def plot_GM_otf(self, save=False):
        '''Plot lightcurves on-the-fly'''
        import matplotlib.dates as mdates

        #although no reason not to be smart and look for the files first
        try:
            if self.Files.Raw['messenger'] != '':
                Mdata=pickle.load(open(self.Files.folders['lightcurves']+self.Files.Raw['messenger'],'rb'))
            else:
                Mdata=self.gen_Messenger_lightcurve(ret=True)
        except KeyError:
            Mdata=self.gen_Messenger_lightcurve(ret=True)
        try:
            if self.Files.Raw['goes'] != '':
                Gdata=pickle.load(open(self.Files.folders['lightcurves']+self.Files.Raw['goes'],'rb'))
            else:
                Gdata=self.gen_GOES_lightcurve(ret=True)
        except KeyError:
            Gdata=self.gen_GOES_lightcurve(ret=True)

        fig,ax1=plt.subplots()
        ax2=ax1.twinx()
        l1,=ax1.step(Mdata['taxis'][:-1],Mdata['lowEcounts'],'b',label= '1.5-12.4 keV')
        l2,=ax1.step(Mdata['taxis'][:-1],Mdata['highEcounts'],'g',label= '3-24.8 keV')
        l3,=ax2.plot(Gdata['taxis'],Gdata['lowEflux'],'k',label='GOES 1-8 $\AA$') #goes short - plot with
        l4,=ax2.plot(Gdata['taxis'],Gdata['highEflux'],'m',label='GOES .5-4 $\AA$') #goes long

        myFmt = mdates.DateFormatter('%H:%M')
        ax1.xaxis.set_major_formatter(myFmt)
        plt.gcf().autofmt_xdate()
        #ax1.set_xlabel(dt.strftime(Mtim[0].date(),'%Y-%b-%d'))
        ax1.set_ylabel('Messenger counts $cm^{-2} keV^{-1} s^{-1}$')
        ax1.set_ylim([10**0,10**4])
        ax1.set_yscale('log')
        ax2.set_ylabel('GOES Flux W$m^{-2}$')
        ax2.set_yscale('log')

        plt.title(dt.strftime(Mdata['taxis'][0].date(),'%Y-%b-%d'))
        ax1.set_xlim([Gdata['taxis'][0],Gdata['taxis'][-1]])
        #plt.legend((l1,l2,l3,l4),(l1.get_label(),l2.get_label(),l3.get_label(),l4.get_label()),loc='upper left',prop={'size':12})
        fig.show()
        if save:
            fname='data/lightcurves/'+dt.strftime(Mdata['taxis'][0].date(),'%Y-%b-%d')+'MG.png'
            fig.savefig(fname)
            #return glong,gshort

    def plot_GM(Mdata,Gdata,n): #will probably have to deal with times to make them all the same...
        import matplotlib.dates as mdates
        tim=Mdata['taxis'][0][n]
        mlen=Mdata['len'][0]#[n]
        nflares=np.shape(Mdata['taxis'][0])[0]/2 #assume 2 energy ranges for now
        Mtim=[]
        for t in tim: Mtim.append(dt.strptime(t,'%d-%b-%Y %H:%M:%S.%f')) #fix messenger times to datetimes
        cps1,cps2=counts_ps(Mdata,n)
        print type(Mtim),np.shape(Mtim[:-1]),np.shape(cps1)

        gtim=Gdata['taxis'][0]#[n]
        glen=Gdata['len'][0]#[n]
        Gtim=[]
        for t in gtim: Gtim.append(dt.strptime(t,'%d-%b-%Y %H:%M:%S.%f')) #fix GOES times to datetimes
        glong=Gdata['ydata'][0][1,0:glen-1]#[n,1,0:glen-1] #what's up with these data?
        gshort=Gdata['ydata'][0][0,0:glen-1]#[n,0,0:glen-1]

        fig,ax1=plt.subplots()
        ax2=ax1.twinx()
        l2,=ax1.step(Mtim[:-1],cps2,'g',label= '3-24.8 keV')

        myFmt = mdates.DateFormatter('%H:%M')
        ax1.xaxis.set_major_formatter(myFmt)
        plt.gcf().autofmt_xdate()
        #ax1.set_xlabel(dt.strftime(Mtim[0].date(),'%Y-%b-%d'))
        ax1.set_ylabel('Messenger counts $cm^{-2} keV^{-1} s^{-1}$')
        ax1.set_ylim([10**0,10**4])
        ax1.set_yscale('log')
        ax2.set_ylabel('GOES Flux W$m^{-2}$')
        ax2.set_yscale('log')

        plt.title(dt.strftime(Mtim[0].date(),'%Y-%b-%d'))
        ax1.set_xlim([Gtim[0],Gtim[glen-1]])
        #plt.legend((l1,l2,l3,l4),(l1.get_label(),l2.get_label(),l3.get_label(),l4.get_label()),loc='upper left',prop={'size':12})
        fig.show()
        fname='data/lightcurves/'+dt.strftime(Mtim[0].date(),'%Y-%b-%d')+'MG.png'
        fig.savefig(fname)
        return glong,gshort

    def plotR(Rdata,n):
        import matplotlib.dates as mdates
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        n=n*4
        tim=Rdata['UT'][0][n][:,0]#since there are 4 channels per flare
        #rlen=Rdata['len'][0][n]

        Rtim=[]
        nflares=np.shape(Rdata['rate'][0])[0]/4 #assume 4 energy ranges for now
        for t in tim: Rtim.append(dt.strptime(t,'%d-%b-%Y %H:%M:%S.%f'))

        #get the energy bins - or do I need to do this since they should be the same? check first

        if np.mean(Rdata['rate'][0][n]) != 0.0:
            ax1.plot(Rtim,Rdata['rate'][0][n],'m',label='4-9 keV') #first energy channel
        if np.mean(Rdata['rate'][0][n+1]) != 0.0:
            ax1.plot(Rtim,Rdata['rate'][0][n+1],'g',label='12-18 keV') #second energy channel I think...
        if np.mean(Rdata['rate'][0][n+2]) != 0.0:
            ax1.plot(Rtim,Rdata['rate'][0][n+2],'c',label='18-30 keV') #etc
        if np.mean(Rdata['rate'][0][n+3]) != 0.0:
            ax1.plot(Rtim,Rdata['rate'][0][n+3],'k',label='30-80 keV') #etc
            #ax1.set_xlabel(dt.strftime(Rtim[0].date(),'%Y-%b-%d'))
        ax1.set_yscale('log')
        ax1.set_ylim([0,10**5])
        ax1.legend(loc='upper right')
        myFmt = mdates.DateFormatter('%H:%M')
        ax1.xaxis.set_major_formatter(myFmt)
        plt.gcf().autofmt_xdate()
        plt.title(dt.strftime(Rtim[0].date(),'%Y-%b-%d'))
        #plt.show()
        fname='data/lightcurves/'+dt.strftime(Rtim[0].date(),'%Y-%b-%d')+'R.png'
        fig.savefig(fname)

    def sunpy_spectrogram(filename,i): #since this is for one spectrogram move it to Flare?
        from scipy.io import readsav #import the spectrogram
        #format everything so it's right
        d=readsav(filename, python_dict=True)
        data={'UT':0.,'rate':0.,'erate':0.,'ltime':0.,'len':0.}
        data['UT'] = d['data']['UT']
        data['rate']=d['data']['rate']
        data['erate']=d['data']['erate']
        data['len']=d['data']['len']
        ebins=d['ebins']

        import mySpectrogram as s
        eaxis=ebins[0:-1]
        last=int(data['len'][0]-1)
        time_axis=np.arange(0,last+1)*4.#data[0]['ltime'] #ndarray of time offsets,1D #need the 4 for RHESSI time interval
        start= dt.strptime(data['UT'][i][0],'%d-%b-%Y %H:%M:%S.000') #datetime
        try:
            end= dt.strptime(data['UT'][i][last],'%d-%b-%Y %H:%M:%S.000') #datetime #might want to modify this to be where the data =0
        except ValueError:
            import datetime
            end = start + datetime.timedelta(seconds=time_axis[-1])
            #print dt.strptime(data['UT'][i][last-1],'%d-%b-%Y %H:%M:%S.000')
        drate=np.transpose(np.log10(data['rate'][i][0:last+1])) #get rid of zeros before taking log10?
        #drate=np.nan_to_num(drate)
        drate[drate == -np.inf] = 0.0
        for n,col in enumerate(drate.T):
            if all([c ==0.0 for c in col]):
                drate[:,n] = np.nan
        a=s.Spectrogram(data=drate,time_axis=time_axis,freq_axis=eaxis,start=start,end=end,instruments=['RHESSI'],t_label='',f_label='Energy (keV)',c_label='log(counts/cm^2 s)')
        #create a figure of designated size - does this even work without a subplot?
        f=plt.figure(figsize=(8,5))
        fig=a.plot(figure=f)

        outfilename='data/spectrograms/'+s.get_day(a.start).strftime("%d%b%Y")+'sgram.png'
        fig.figure.savefig(outfilename)
        plt.clf()
        return a

    def plot_aia(self, AIAmap=True,zoom=False,save=False, wave='304',all_wave=False):
        '''Plot the stereo map(s) and overplot the limb and source_pos_STEREO. If all_wave, then subplots'''
        if not AIAmap:
            AIAmap=self.get_AIA(wave=wave)
        m=AIAmap[AIAmap.keys()[0]]

        if all_wave:
            maps = [self.get_AIA(wave=w) for w in ['171','193','304']]
            print maps.keys()
            fig = plt.figure(figsize=(len(maps)*3, 5))
        if not zoom:
            p=m.wcs
            if not all_wave:
                fig = plt.figure(figsize=(6, 5))
                ax = fig.add_subplot(1, 1, 1, projection=p)
                m.plot(axes=ax)
            if all_wave:
                for i in range(len(maps)):
                    ax = fig.add_subplot(1,len(maps),i+1, projection=p) #do the math
                    mm=maps[i][maps.keys()[0]] #I think
                    mm.plot(axes=ax)
                print maps.keys()

        if zoom: #zoom to a rectangle centered on source_pos
            width=400.*u.arcsec
            height=450.*u.arcsec
            coord=self.Properties.source_pos*u.arcsec
            submap = m.submap(u.Quantity((coord[0]-width/2,
                                        coord[0] + width/2)),
                            u.Quantity((coord[1]-height/2,
                                        coord[1] + height/2)))
            p=submap
            if not all_wave:
                fig = plt.figure(figsize=(6, 5))
                ax = fig.add_subplot(1, 1, 1, projection=p.wcs)
                submap.plot(axes=ax)
            if all_wave:
                for i in range(len(maps)):
                    ax = fig.add_subplot(1,len(maps),i+1, projection=p.wcs) #do the math
                    mm=maps[i][maps.keys()[0]] #I think
                    submap= mm.submap(u.Quantity((coord[0]-width/2,
                                        coord[0] + width/2)),
                            u.Quantity((coord[1]-height/2,
                                        coord[1] + height/2)))
                    submap.plot(axes=ax)

        fig.show()
        fname=self.Files.dir+self.Files.folders['plots']+'/'+dt.strftime(self.Datetimes.Messenger_peak.date(),'%Y%m%d')+'AIA_'+str(wave)
        if save:
            tag=raw_input('Identifying tag? ')
            fig.savefig(fname+tag+'.png')

    def plot_stereo(self,AIAmap=True,maps=True,oplotlimb=False,oplotpoint=False,save=False,subplots=False,zoom=False,diff=False,return_limb=True,n=False,xran=False,yran=False):
        '''Plot the stereo map(s) and overplot the limb and source_pos_STEREO. Enable the option to plot the difference of two maps, along with the original (map0)'''
        from astropy.coordinates import SkyCoord
        import copy

        if not AIAmap: AIAmap=self.get_AIA()
        if not maps: maps=self.get_STEREO()
        if subplots: fig = plt.figure(figsize=(len(maps)*3, 5))
        if self.Properties.source_pos_STEREO=='':
            self.convert_coords_aia2stereo(map_aia=AIAmap,map_stereo=maps[0],quiet=True) #because for each map it will be slightly different?
        if n: #limit the number of plots to this number
            maps=maps[0:n]

        for i,m in enumerate(maps):
            m=copy.deepcopy(m['SECCHI'])
            self.Datetimes.current_map_time=m.date
            if diff:
                if i!=0:
                    m0=copy.deepcopy(maps[i-1]['SECCHI'])
                    md=copy.deepcopy(m)
                    m.data=md.data-m0.data #can I add color scaling somehow?
                    print i, i-1

            if not zoom:
                p=m.wcs
                pmap=m
                if not subplots:
                    fig = plt.figure(figsize=(len(maps)*3, 5))
                    ax = fig.add_subplot(1, 1, 1, projection=p)
                if subplots: ax = fig.add_subplot(1,len(maps),i+1, projection=p) #do the math
                m.plot(axes=ax)

            if zoom: #zoom to a rectangle centered on source_pos_STEREO
                width=400.*u.arcsec
                height=450.*u.arcsec
                coord=self.Properties.source_pos_STEREO*u.arcsec
                submap = m.submap(u.Quantity((coord[0]-width/2,
                                        coord[0] + width/2)),
                            u.Quantity((coord[1]-height/2,
                                        coord[1] + height/2)))
                p=submap
                pmap=submap
                if not subplots:
                    fig = plt.figure(figsize=(6, 5))
                    ax = fig.add_subplot(1, 1, 1, projection=p)
                if subplots: ax = fig.add_subplot(1,len(maps),i+1, projection=p) #do the math
                submap.plot(axes=ax)

            ax.set_autoscale_on(False)
            if oplotlimb: #do I need to abjust this for Stereo A? check...
                r = AIAmap[AIAmap.keys()[0]].rsun_obs.to(u.deg)-1*u.arcsec # remove the one arcsec so it's on disk.
                #r = (maps[0]['SECCHI'].rsun_arcseconds*u.arcsec).to(u.deg)-1*u.arcsec # remove the one arcsec so it's on disk.
                import scipy
                #r=astropy.constants.R_sun/scipy.constants.au
                if self.Properties.source_pos[0] >0:
                    th = np.linspace(-360*u.deg, -180*u.deg)
                else:
                    th = np.linspace(-180*u.deg, 0*u.deg) #change from 180

                x = r * np.sin(th)
                y = r * np.cos(th)
                coords = SkyCoord(x, y, frame=AIAmap[AIAmap.keys()[0]].coordinate_frame)
                #print coords._rsun
                coords._rsun=m.rsun_meters
                #print coords._rsun
                hgs = coords.transform_to('heliographic_stonyhurst')
                hgs.D0 = m.dsun
                hgs._rsun=m.rsun_meters
                print hgs._rsun
                if np.mean(m.heliographic_longitude) >0:
                    hgs.L0 = m.heliographic_longitude#+1*u.arcsec #add back the 1 arcsec
                else:
                    hgs.L0 = m.heliographic_longitude#-1*u.arcsec) #*1.01 try this... it's what's done for the euv in pascal's routine
                hgs.B0 = m.heliographic_latitude
                limbcoords = hgs.transform_to(pmap.coordinate_frame)
                ax.plot_coord(limbcoords, color='w') #this is not overplotting in the zoomed version! why not? do I have to do this before the map.plot command?

                #translate Pascal's method and see what it gives...
                #default, roll_angle_range, [0,360]
                #default, npts, 360L
                #default, factor, 1d	;;1.0 for photospheric limb, about 1.01 for EUV limb...?
                #roll_angle_range= [0,360]
                #npts=360
                #factor=1.01

                #Rs=phys_const('R_sun')/phys_const('AU')		;;Rs in [AU]
                #a=psh_binning(roll_angle_range[0],roll_angle_range[1],Npts)/radeg()
				# heeq=[0,cos(a[i]),sin(a[i])]*Rs*factor		;;soft/transparent/translucid limb
                #Rs=astropy.constants.R_sun/astropy.constants.au #Rs in AU
                #a=np.linspace(roll_angle[0]*u.deg,roll_angle[1]*u.deg,npts)

                #;;Now, rotate for L-angle around Z-axis, then B-angle around Y-axis:

				#th=blr[0]/radeg()
				#M = [[cos(th),0,sin(th)],[0,1,0],[-sin(th),0,cos(th)]]
				#heeq = M # heeq

				#th=blr[1]/radeg()
				#M = [[cos(th),sin(th),0],[-sin(th),cos(th),0],[0,0,1]]
				#heeq = M # heeq

				#xxyyd=heeq2hpc(heeq,[map.ROLL_ANGLE,map.b0,map.L0,map.RSUN])


            if oplotpoint: #overplot source_loc_STEREO
                scs=self.Properties.source_pos_STEREO*u.arcsec
                sc=SkyCoord(scs[0],scs[1],frame=pmap.coordinate_frame)
                ax.plot_coord(sc, color='r',marker='o',markersize=10) #might need to make it bigger
                #ax.plot([self.Properties.source_pos_STEREO[0]],[self.Properties.source_pos_STEREO[1]],marker='o',markersize=10,color='red',projection=p)
                print self.Properties.source_pos_STEREO

            #pmap.plot(axes=ax)

            fig.show()

        fname=self.Files.dir+self.Files.folders['plots']+'/'+dt.strftime(self.Datetimes.Messenger_peak.date(),'%Y%m%d')+'EUVI_'
        if save:
            tag=raw_input('Identifying tag? ')
            fig.savefig(fname+tag+'.png')

        if return_limb:
            return limbcoords

    def get_closest_limbcoord(self,limbcoords=True,redef=False,plot=True,smap=False):
        '''Find the closest limb coordinate to source_pos_STEREO'''
        import shapely.geometry as geom
        if not limbcoords:
            limbcoords=self.plot_stereo(AIAmap=False,maps=False,oplotlimb=True,subplots=True,return_limb=True)
        sp=self.Properties.source_pos_STEREO
        if type(sp) != list or sp == [] or redef: #ask for source position
            coordstr=raw_input('Please enter integer coordinates on the disk. Format: xxx,yyy')
            self.Properties.source_pos_STEREO= [int(coordstr[:coordstr.find(',')]),int(coordstr[coordstr.find(',')+1:])]#format the string
        sp=geom.Point(self.Properties.source_pos_STEREO)
        lc=[]
        for x,y in zip(limbcoords.Tx.value,limbcoords.Ty.value):
            lc.append([x,y]) #format limbcoords
        curve=geom.LineString(lc) #these might be in some weird format...lon,lat or something

        distance=curve.distance(sp)*u.arcsec
        if plot:        #draw segment
            fig,ax=plt.subplots()
            ax.set_autoscale_on(False)
            if smap:
                if type(smap) != dict:
                    pass
                else: #plot the map
                    smap['SECCHI'].plot(axes=ax)
            ax.plot(limbcoords.Tx.value,limbcoords.Ty.value)
            ax.axis('equal')
            point_on_line=curve.interpolate(curve.project(sp))
            ax.plot([sp.x, point_on_line.x], [sp.y, point_on_line.y], color='red', marker='o', scalex=False, scaley=False)
            fig.canvas.draw()
            plt.show()
        print 'coordinates of closest point on AIA limb: ', list(point_on_line.coords)[0]
        print 'distance to closest point on AIA limb: ', distance
        self.Properties.e2={'distance2limb':distance,'coordsonlimb': list(point_on_line.coords)[0]}

        return list(point_on_line.coords)[0]

    def plot_RHESSI_AIA(self,AIAmap=False,erange=False,wave='193',title=False):
        '''Overplot RHESSI contours on AIA map'''
        import matplotlib.transforms as transforms
        #from mpl_toolkits.axes_grid.inset_locator import inset_axes
        #default energy range: all
        all_erange=[[4,9],[12,18],[18,30],[30,80]] #deal with choosing energy range later
        if not erange:
            erange=all_erange
        #convert the fits files to maps
        os.chdir(self.Files.dir+self.Files.folders['clean'])
        fitsfiles=glob.glob(str(self.ID)+'*plot.fits')
        #            for f in sunpy.map.Map(files):
        #        maps.append({f.instrument: f.submap((-1100, 1100) * u.arcsec, (-1100, 1100) * u.arcsec)}) #this could be empty if files is empty
        maps=[]
        for i,f in enumerate(fitsfiles):
            mapf=sunpy.map.Map(f)
            maps.append({str(all_erange[i]):mapf}) #try this

        if not AIAmap:
            AIAmap=self.get_AIA(wave=wave)
        #trim the map to RHESSI dimensions
        center=self.Properties.source_pos
        fov= 128 #128 rhessi pixels to arcsec
        AIAmap=AIAmap[AIAmap.keys()[0]].submap((center[0]-fov,center[0]+fov)*u.arcsec,(center[1]-fov,center[1]+fov)*u.arcsec)
        if not title:
            title=''

        strerange=[str(e) for e in erange]
        fig,ax=plt.subplots()
        AIAmap.plot(axes=ax,title=title)
        #maps[0][maps[0].keys()[0]].draw_contours([.75]*u.percent,axes=ax)
        #inset_axes=inset_axes(ax,width="100%",height="100%",)
        a=plt.axes([center[0]-fov,center[1]-fov,center[0]+fov,center[1]+fov])
        a=plt.axes([-1000,-450,-800,-200])
        for m in maps:
            if m.keys()[0] in strerange:
                offset = transforms.ScaledTranslation(center[0], center[1],fig.dpi_scale_trans)
                shadow_transform = ax.transData + offset
                m[m.keys()[0]].draw_contours([.15]*u.percent,axes=a) #contour args ... whatever those are. colors and labels I hope?
                #do I need to make another axis... offset it by center? it's not reading the keywords from fits :(
        plt.show()

    def plot_RHESSI_AIA_IDL(self,AIAmap=False,erange=False,wave='193',title=False,tag=False,fov=100):
        '''Overplot RHESSI contours on AIA map using IDL'''
        #default energy range: all
        import pidly
        all_erange=[[4,9],[12,18],[18,30],[30,80]] #deal with choosing energy range later
        if not erange:
            erange=all_erange
        #convert the fits files to maps
        os.chdir(self.Files.dir+self.Files.folders['clean'])
        rfitsfiles=glob.glob(str(self.ID)+'*plot.fits')

        os.chdir(self.Files.dir)
        afitsfiles=self.Files.Raw['aia']
        for a in afitsfiles:
            if wave in a:
                aa=a
                print aa
        if not aa:
            aa=afitsfiles[0]

        idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
        idl.afile=aa
        if AIAmap:
            idl.afile=AIAmap
        idl.adir=self.Files.dir+self.Files.folders['stereo-aia']+'/'
        idl.rdir=self.Files.dir+self.Files.folders['clean']+'/'
        idl('fits2map,adir+afile,amap')
        idl('loadct,9')
        idl('keywords = PSWINDOW()')
        #DEVICE, _EXTRA=keywords
        idl("set_plot,'ps'")
        if not tag:
            idl.filename=str(self.ID)+'_rhessi_AIA.eps'
        else:
            idl.filename=str(self.ID)+'_rhessi_AIA_'+tag+'.eps'

        idl('device, filename=filename,/encapsulated,/decomposed,/color,bits_per_pixel=8')
        #Idl('loadct,9')
        center=self.Properties.source_pos

        labels=[]
        idl.center=center
        idl.xran=[center[0]-fov,center[0]+fov]
        idl.yran=[center[1]-fov,center[1]+fov]
        idl('plot_map,amap,center=center,xran=xran,yran=yran,/window,/limb,lcolor=cgcolor("White"),lthick=5')
        #idl('loadct,12')
        idl("colors=[cgcolor('White'),cgcolor('Sky Blue'),cgcolor('Lime Green'),cgcolor('Magenta')]")
        ne=len(erange)
        idl.ner=ne
        #idl('ne=fix(ne)')
        idl('lcolor=strarr(ner)')
        for j,ee in enumerate(erange):
            estr=str(ee[0])+"-"+str(ee[1])
            idl.j=j
            for i,r in enumerate(rfitsfiles):
                if estr in r:
                    print r
                    idl.rfiles=r
                    #color=i*75
                    idl.i=i
                    idl("color=colors[i]")
                    idl("lcolor[j]=i")
                    #idl("colors=[cgcolor('White'),cgcolor('Sky Blue'),cgcolor('Lime Green'),cgcolor('Magenta')]")
                    #if i==0:
                    #    color=16
                    #    idl.color=color
                    #    #idl("color=cgcolor(White)")
                    print i,all_erange[i],ee,idl.color
                    #colors.append(int(color))
                    labels.append(estr+' keV')
                    idl('fits2map,rdir+rfiles,rmap')
                    idl('level=mean(rmap.data)+3*stddev(rmap.data)')
                    idl('plot_map,rmap,/cont,/over,levels=[level],color=string(color),thick=10')
        idl.labels=labels
        #idl.colors=colors
        #idl('colors=string(colors)')
        idl('print,lcolor')
        #idl('print,typename(colors),typename(string(colors))')
        idl("al_legend,labels,colors=colors[lcolor],linestyle=0,/bottom,/right,thick=5,textcolors=cgcolor('White')")
        #idl.filename=str(self.ID)+'_rhessi_AIA.png'
        #idl("write_png,filename,tvrd(/true)")
        idl('device,/close_file')



