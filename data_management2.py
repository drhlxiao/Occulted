 #######################################
# data_management.py
# Erica Lastufka 31/1/17 

#Description: Class of functions to deal with data for this study
#######################################

#######################################
# Usage:

# my_flare_list = Data_Study()

######################################
import numpy as np

class Data_Study(object):
    ''' Dota for this study will have the following properties:
            ID: RHESSI flare ID (if applicable)
            Datetimes:
                Messenger datetimes
                RHESSI datetimes
                Observation time intervals
            Angle: Angle between the Earth and Messenger
            Flare properties:
                Messenger T and EM1
                GOES class calculated given messenger T and EM1
                GOES class RHESSI calculated how?
                GOES class from GOES satellites
                source visibility: good or not?
                total counts Messenger
                total counts RHESSI (need to define energy range)
            Data properties:
                XRS .dat and .lbl filenames
                Quicklook image filenames
                Messenger data path
                RHESSI data path
                RHESSI browser URL's
                csv storage filename
                array of good detectors
            Notes'''
    def __init__(self, filename=0, length=0):
        '''Return an empty Data_Study object, unless a filename to a correctly formatted csv is given. Do I need to be able to go from a dictionary to an object too?'''
        self.ID= [0]
        self.Datetimes={"Messenger_datetimes":[''],"RHESSI_datetimes":[''],"Obs_start_time":[''],"Obs_end_time":['']}
        self.Angle = [0.0]
        self.Flare_properties={"Messenger_T":[0.0],"Messenger_EM1":[0.0], "Messenger_GOES":["A.0"],"RHESSI_GOES":["A.0"],"GOES_GOES":["A.0"],"source_vis":["Y"],"Messenger_total_counts":[0],"RHESSI_total_counts":[0],"Location":[0,0]}
        self.Data_properties={"XRS_files":["xrs.dat"],"QL_images":["foo.png"],"Messenger_data_path":["/Users/wheatley/Documents/Solar/occulted_flares/data/dat_files"],"RHESSI_data_path":["/Users/wheatley/Documents/Solar/occulted_flares/data/"],"RHESSI_browser_urls":["http://sprg.ssl.berkeley.edu/~tohban/browser/?show=grth1+qlpcr+qlpg9+qli02+qli03+qli04+synff&date="],"csv_name":["foo.csv"],"good_det":[1,0,1,1,1,1,1,1,1]}
        self.Notes= ["notes"]
        
        if length != 0:
            zerov = np.zeros(length, float)
            zerovs = np.zeros(length, str)
            self.ID = zerov
            self.Datetimes={"Messenger_datetimes":zerov,"RHESSI_datetimes":zerov,"Obs_start_time":zerov,"Obs_end_time":zerov}
            self.Angle = zerov
            self.Flare_properties={"Messenger_T":zerov,"Messenger_EM1":zerov, "Messenger_GOES":zerovs,"RHESSI_GOES":zerovs,"GOES_GOES":zerovs,"source_vis":zerovs,"Messenger_total_counts":zerov,"RHESSI_total_counts":zerov,"Location":zerov}
            self.Data_properties={"XRS_files":zerovs,"QL_images":zerovs,"Messenger_data_path":np.repeat("/Users/wheatley/Documents/Solar/occulted_flares/data/dat_files",length),"RHESSI_data_path":np.repeat("/Users/wheatley/Documents/Solar/occulted_flares/data/",length),"RHESSI_browser_urls":zerovs,"csv_name":zerovs,"good_det":zerovs}
            self.Notes= zerovs       
        
        if filename != 0:
            if filename.endswith('csv'):
                import pandas as pd
                data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
                self.ID = data['ID']
                self.Datetimes={"Messenger_datetimes":data["Messenger_datetimes"],"RHESSI_datetimes":data["RHESSI_datetimes"],"Obs_start_time":data["Obs_start_time"],"Obs_end_time":data["Obs_end_time"]} #might need to convert these all to datetime format. There will be different formats depending on the csv!
                self.format_datetimes()
                self.Angle = data['Angle']
                self.Flare_properties={"Messenger_T":data["Messenger_T"],"Messenger_EM1":data["Messenger_EM1"], "Messenger_GOES":data["Messenger_GOES"],"RHESSI_GOES":data["RHESSI_GOES"],"GOES_GOES":data["GOES_GOES"],"source_vis":data["source_vis"],"Messenger_total_counts":data["Messenger_total_counts"],"RHESSI_total_counts":data["RHESSI_total_counts"],"Location":data["Location"]}
                self.Data_properties={"XRS_files":data["XRS_files"],"QL_images":data["QL_images"],"Messenger_data_path":data["Messenger_data_path"],"RHESSI_data_path":data["RHESSI_data_path"],"RHESSI_browser_urls":data["RHESSI_browser_urls"],"csv_name":data["csv_name"],"good_det":data["good_det"]}
                self.Notes= data["Notes"]
            elif filename.endswith('sav'):
                from scipy.io import readsav
                #format everything so it's right
                data=readsav(filename, python_dict=True)
                data=data['flare_list']
                self.ID = data['ID'][0].tolist()
                self.Datetimes['Messenger_datetimes']=data['datetimes'][0]['messenger_datetimes'][0].tolist()
                self.Datetimes['RHESSI_datetimes']=   data['datetimes'][0]['rhessi_datetimes'][0].tolist()
                self.Datetimes['Obs_start_time']=     data['datetimes'][0]['obs_start_time'][0].tolist()
                self.Datetimes['Obs_end_time']=       data['datetimes'][0]['obs_end_time'][0].tolist()
                self.Angle = data['angle'][0].tolist()
                self.Flare_properties["Messenger_T"]=   data['flare_properties'][0]["messenger_t"][0].tolist()
                self.Flare_properties["Messenger_EM1"]= data['flare_properties'][0]["messenger_em1"][0].tolist()
                self.Flare_properties["Messenger_GOES"]=data['flare_properties'][0]["messenger_goes"][0].tolist()
                self.Flare_properties["RHESSI_GOES"]=   data['flare_properties'][0]["RHESSI_GOES"][0].tolist()
                self.Flare_properties["GOES_GOES"]=   data['flare_properties'][0]["GOES_GOES"][0].tolist()
                self.Flare_properties["source_vis"]=    data['flare_properties'][0]["source_vis"][0].tolist()       
                self.Flare_properties["Messenger_total_counts"]=data['flare_properties'][0]["messenger_total_counts"][0].tolist()
                self.Flare_properties["RHESSI_total_counts"]=data['flare_properties'][0]["rhessi_total_counts"][0].tolist()
                self.Flare_properties["Location"]=data['flare_properties'][0]["location"][0].tolist()
                self.Data_properties["XRS_files"]=data['data_properties'][0]["xrs_files"][0].tolist()
                self.Data_properties["QL_images"]=data['data_properties'][0]["ql_images"][0].tolist()
                self.Data_properties["Messenger_data_path"]=data['data_properties'][0]["messenger_data_path"][0].tolist()
                self.Data_properties["RHESSI_data_path"]=data['data_properties'][0]["rhessi_data_path"][0].tolist()
                self.Data_properties["RHESSI_browser_urls"]=data['data_properties'][0]["rhessi_browser_urls"][0].tolist()
                self.Data_properties["csv_name"]=data['data_properties'][0]["csv_name"][0].tolist()
                self.Data_properties["good_det"]=data['data_properties'][0]["good_det"][0].tolist()
                self.Notes= data["notes"][0].tolist()
                self.format_datetimes()
                
    def __iter__(self):
        '''Returns a generator that iterates over the object'''
        for attr, value in self.__dict__.iteritems():
            yield attr, value

    def slice(self, indices):
        for n,i in enumerate(indices):
            self.ID[n]=self.ID[i]
            self.Angle[n]=self.Angle[i]
            self.Notes[n]=self.Notes[i]
            for key in self.Datetimes.keys():
                #print key,type(key),i
                self.Datetimes[key][n] = self.Datetimes[key][i]
            for key in self.Flare_properties.keys():
                self.Flare_properties[key][n] = self.Flare_properties[key][i]
            for key in self.Data_properties.keys():
                self.Data_properties[key][n] = self.Data_properties[key][i]
        self.ID = self.ID[0:n-1]
        self.Angle=self.Angle[0:n-1]
        self.Notes=self.Notes[0:n-1]
        for key in self.Datetimes.keys():
            self.Datetimes[key]= self.Datetimes[key][0:n-1]
        for key in self.Flare_properties.keys():
            self.Flare_properties[key]= self.Flare_properties[key][0:n-1]
        for key in self.Data_properties.keys():
            self.Data_properties[key]= self.Data_properties[key][0:n-1]
            
    def format_datetimes(self):
        '''Move the datetimes in the self.Datetimes dictionary between formats. This assumes that the datetime is stored in a STRING not a datetime - as it will be when reading from a .csv or .sav file. If it is a datetime, format it for output to csv'''
        from datetime import datetime as dt
        import re
        data = self.Datetimes
        for key in data.keys():
            newdata = []
            if type(data[key][0]) == dt:
                for d in data[key]:
                    newdata.append(d.strftime('%d-%b-%Y %H:%M:%S.000'))
            else:
                try: #it is a string
                    if data[key][0] == '':
                        newdata=data[key]
                        fc=''
                        #for i in range(0,len(data['RHESSI_datetimes'])-1):
                        #   newdata.append('')
                        #continue
                    elif ('AM' or 'PM') in data[key][0]:
                        fc = '%m/%d/%y %I:%M:%S %p'
                    elif '/' in data[key][0]:
                        fc = '%m/%d/%y %H:%M:%S.00'
                    elif re.search('[a-z]', data[key][0]) is not None: #not None if there is a letter
                        fc ='%d-%b-%Y %H:%M:%S' #current format code - might need to add .000
                    else: #it's come straight from idl ....
                        fc ='%Y-%m-%d %H:%M:%S'
                    if fc != '':
                       for string in data[key]:
                           if string.startswith("'"):
                               try:
                                   newdata.append(dt.strptime(string[1:-1], fc))
                               except ValueError:
                                   fc=fc ='%d-%b-%Y %H:%M:%S.000'
                                   newdata.append(dt.strptime(string[1:-1], fc))
                           else:
                               try:
                                   newdata.append(dt.strptime(string, fc))
                               except ValueError:
                                   fc=fc ='%d-%b-%Y %H:%M:%S.000'
                                   newdata.append(dt.strptime(string, fc))
                except TypeError: #if nothing in there - let's change things to empty strings for now so that IDL plays nice
                    for i in range(0,len(data['RHESSI_datetimes'])-1):
                        newdata.append('')
            self.Datetimes[key] = newdata
            #print key, newdata[0:10]
                    
    def export2csv(self, filename):
        '''Exports Python dictionary to csv file'''
        #from pandas.io.json import json_normalize
        #json_normalize(self.__dict__).to_csv(filename) #flattens the dictionaries so it looks ok in a spreadsheet... only it does this wrong!

        #first flatten the dictionaries - I should write a function to do this later:
        self.Data_properties['csv_name']=np.repeat(filename,len(self.ID))
        self.format_datetimes()
        #print self.Datetimes['Messenger_datetimes'][0:10]
        d=self.__dict__
        headers = ['ID','Messenger_datetimes','RHESSI_datetimes','Obs_start_time','Obs_end_time','Messenger_T','Messenger_EM1','Messenger_GOES','Messenger_total_counts','RHESSI_GOES','GOES_GOES','RHESSI_total_counts','Location','source_vis','Messenger_data_path','RHESSI_data_path','XRS_files','QL_images','RHESSI_browser_urls','csv_name','good_det','Angle','Notes']
        #zip them all together
        rows = zip(d['ID'],
                   d['Datetimes']['Messenger_datetimes'],d['Datetimes']['RHESSI_datetimes'],
                   d['Datetimes']['Obs_start_time'],d['Datetimes']['Obs_end_time'],
                   d['Flare_properties']['Messenger_T'],d['Flare_properties']['Messenger_EM1'],
                   d['Flare_properties']['Messenger_GOES'],d['Flare_properties']['Messenger_total_counts'],
                   d['Flare_properties']['RHESSI_GOES'],d['Flare_properties']['GOES_GOES'],d['Flare_properties']['RHESSI_total_counts'],
                   d['Flare_properties']['Location'],d['Flare_properties']['source_vis'],
                   d['Data_properties']['Messenger_data_path'],d['Data_properties']['RHESSI_data_path'],
                   d['Data_properties']['XRS_files'],d['Data_properties']['QL_images'],
                   d['Data_properties']['RHESSI_browser_urls'],d['Data_properties']['csv_name'],d['Data_properties']['good_det'],
                   d['Angle'],d['Notes'])

        import csv
        with open(filename,'wb') as f:
            writer=csv.writer(f)
            writer.writerow(headers)
            for row in rows:
                writer.writerow(row)
        
#    def export2idl_old(self, savname):
#        '''because the dictionaries might be to big to send directly, read it in from the csv file'''
#        #convert datetimes to IDL-friendly format
#        #Format for time is YY/MM/DD, HHMM:SS.SSS
#        #or alternatively, DD-Mon-YY HH:MM:SS.SSS -> '%d-%b-%Y %H:%M:%S.000'
#        nMdt, nRdt,nOst,nOet = [],[],[],[]
#        if type(self.Datetimes['RHESSI_datetimes'][0]) == str:
#            #check that it's got quote marks then continue
#            print 'dates are already strings!'
#            if self.Datetimes['RHESSI_datetimes'][0].startswith("'") or self.Datetimes['RHESSI_datetimes'][0] == '':
#                nMdt =self.Datetimes['Messenger_datetimes']
#                nRdt = self.Datetimes['RHESSI_datetimes']
#            else:
#                #put quotes around everything
#                for Mdt, Rdt in zip(self.Datetimes['Messenger_datetimes'], self.Datetimes['RHESSI_datetimes']): 
#                    nMdt.append("'"+Mdt+"'") 
#                nRdt.append("'"+Rdt+"'")
#        else:
#            for Mdt, Rdt in zip(self.Datetimes['Messenger_datetimes'], self.Datetimes['RHESSI_datetimes']): #will this work even when they're diffe#rent lengths? We'll see
#                nMdt.append("'"+Mdt.strftime('%d-%b-%Y %H:%M:%S.000')+"'") #if you don't have the quotes idl won't read the csv in correctly
#                nRdt.append("'"+Rdt.strftime('%d-%b-%Y %H:%M:%S.000')+"'")

#        if (self.Datetimes['Obs_start_time'][0] != 0.0 and self.Datetimes['Obs_start_time'][0] != ''): #should change the default to be an empty string array...
#            if type(self.Datetimes['Obs_start_time'][0]) == str:
#                if self.Datetimes['Obs_start_time'][0].startswith("'"):
#                    nOst=self.Datetimes['Obs_start_time']
#                    nOet=self.Datetimes['Obs_end_time']
#                else:
#                    for Ost, Oet in zip(self.Datetimes['Obs_start_time'], self.Datetimes['Obs_end_time']): 
#                        nOst.append("'"+Ost+"'") 
#                        nOet.append("'"+Oet+"'")
#            else:
#                for Ost, Oet in zip(self.Datetimes['Obs_start_time'], self.Datetimes['Obs_end_time']): 
#                    nOst.append("'"+Ost.strftime('%d-%b-%Y %H:%M:%S.000')+"'") 
#                    nOet.append("'"+Oet.strftime('%d-%b-%Y %H:%M:%S.000')+"'")
#            self.Datetimes['Obs_start_time']=nOst
#            self.Datetimes['Obs_end_time']=nOet
            
#        self.Datetimes['RHESSI_datetimes']=nRdt
#        self.Datetimes['Messenger_datetimes']=nMdt
#        #print self.Datetimes['Messenger_datetimes'][0:10]
#        #check if the csv file exists, if not, make it
#        import glob
#        if not savname[:-3]+'csv' in glob.glob('*.csv'):
#            print 'no csv file found, converting to csv first'
#            self.export2csv(savname[:-3]+'csv')
#        else:
#            ans=raw_input(savname[:-3]+'csv already exists! Overwrite? (Y/N)')
#            if ans == 'Y':
#                self.export2csv(savname[:-3]+'csv')

#        import pidly
#        idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
#        idl('savname',savname)
#        idl('csvname',savname[:-3]+'csv')
#        idl('foo=read_csv(csvname)')
#        idl('flare_list={ID:foo.field01,Datetimes:{Messenger_datetimes:foo.field02,RHESSI_datetimes:foo.field03,Obs_start_time:foo.field04,Obs_end_time:foo.field05},Flare_properties:{Messenger_T:foo.field06,Messenger_EM1:foo.field07,Messenger_GOES:foo.field08,Messenger_total_counts:foo.field09,RHESSI_GOES:foo.field10,GOES_GOES:foo.field11,RHESSI_total_counts:foo.field12,source_vis:foo.field13},Data_properties:{Messenger_data_path:foo.field14,RHESSI_data_path:foo.field15,XRS_files:foo.field16,QL_images:foo.field17,RHESSI_browser_urls:foo.field18,csv_name:foo.field19},Angle:foo.field20,Notes:foo.field21}') 
#        idl('save,flare_list, filename=savname')

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

    def export2pickle(self, picklename):
        '''Replaces save_flare_list. Saves the data object in a .p file'''
        import pickle
        pickle.dump(self, open(picklename, 'wb'))


def download_messenger(self):
    '''Downloads Messenger .dat and .lbl files from the database, given event date. Can be used to get missing files too.'''
    import urllib
    dataurl = 'https://hesperia.gsfc.nasa.gov/messenger/'
    #subfolders by year, month,day (except 2001)
    listlen = len(self.ID)
    for i,dt in zip(range(0,listlen -1),self.Datetimes['Messenger_datetimes']):
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
        if not self.Data_properties['XRS_files'][i]:
            self.Data_properties['XRS_files'][i]=filename

def open_in_RHESSI_browser(self, opent=False):
    import webbrowser
    browserurl = 'http://sprg.ssl.berkeley.edu/~tohban/browser/?show=grth1+qlpcr+qlpg9+qli02+qli03+qli04+synff&date=' #20120917&time=061500'
    #subfolders by year, month,day (except 2001)
    #warn if you're opening more than 20 tabs
    if opent == True and len(self.ID) > 20:
        ans = raw_input('You are about to open ' + str(len(self.ID)) + ' tabs! Are you sure?')
        if ans != 'Y':
            opent = False        
    for i,dt in enumerate(self.Datetimes['Messenger_datetimes']):
        try:
            address=browserurl + dt.strftime('%Y%m%d') + '&time=' +dt.strftime('%H%M%S')
        except AttributeError: #datetime for RHESSI is empty
            address=browserurl + self.Datetimes['Messenger_datetimes'][i].strftime('%Y%m%d') + '&time=' +self.Datetimes['Messenger_datetimes'][i].strftime('%H%M%S')
        if not self.Data_properties['RHESSI_browser_urls'][i]:
            self.Data_properties['RHESSI_browser_urls'][i] = address
        if opent:
            webbrowser.open_new_tab(address)
            
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

def get_ratio(self):
    '''Get ratio of Messenger goes v actual goes'''
    mvals,rvals=[],[]
    mc = self.Flare_properties["Messenger_GOES"]
    gc=self.Flare_properties["GOES_GOES"]
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
