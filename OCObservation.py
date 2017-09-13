 #######################################
# flare_class.py
# Erica Lastufka 17/5/17 

#Description: Class of functions to deal with all data,files,calculations of the study
#######################################

import numpy as np

class OCObservation(object):
    ''' Object will have the following attributes:
            Angle: Angle between RHESSI and messenger observation
            Det: Good detectors (RHESSI)
            STEREO: Closest satellite to Messenger: A or B
            CLEAN: Dictionary of clean parameters (RHESSI)
            Rcounts: Dictionary of rhessi counts in various energy ranges
            Mcounts: Same for Messenger
            Rcountrate:Rhessi count rate
            Mcountrate:Messenger
            Rpflux: Rhessi photon flux
            Mpflux: Messenger photon flux
            Rebins: Rhessi energy bins
            Mebins: Messenger energy bins
            Swave: Stereo wavelength
            Slos: Stereo line of sight
            AIAwave:AIA wavelength(s)
            Notes
        Methods: 
            Calculate photon flux rates? TBD'''
    def __init__(self, ID, legacy=True,filename=False):
        '''Return an empty Data_Study object, unless a filename to a correctly formatted csv is given. Do I need to be able to go from a dictionary to an object too?'''
        tags=['CLEAN','Rcountrate','Mcountrate','Rpflux','Mpflux','Rebins','Mebins','Swave','Slos','AIAwave']
        if legacy: #get info from legacy OCData object and the associated csv/sav files
            if not filename: filename= '/Users/wheatley/Documents/Solar/occulted_flares/flare_lists/list_final.csv'#default file to read
            import pandas as pd
            data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
            i=self.get_index(ID,data) #get the index of the flare if it's in a list
            self.Angle=data["Angle"][i]
            self.Det=data['good_det'][i]
            self.STEREO=data['RHESSI_datetimes'][i]
            for att,key in zip(tags,tags):
                setattr(self,att,'')
 
        if not legacy:
            #read attributes from csv file (can just restore the pickle file otherwise. Build this into OCFlare class):
            import pandas as pd
            if filename: #it's the big one
                data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
                #i=self.get_index(ID,data) #get the index of the flare if it's in a list
                for key in data.keys(): #only do this if it starts with Observation
                    if key.startswith('Observation.'):
                        dat=data[key]
                        key=key[key.find('.')+1:] #trim it
                        setattr(self,key,dat.values[0])
            if not filename:
                filename= '/Users/wheatley/Documents/Solar/occulted_flares/data/objects/'+str(ID)+'OCObservation.csv'#default file to read - need an except in case it doesn't exist
                data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
                #i=self.get_index(ID,data) #get the index of the flare if it's in a list
                for key in data.keys(): #only do this if it starts with Observation
                    setattr(self,key,data.values[0])

        #set some defaults
        if self.Rebins == '': self.Rebins=[[],[],[],[],[],[]]
        if self.Mebins == '': self.Mebins=[[],[]]
        if self.Swave == '': self.Swave=193
        if self.AIAwave == '': self.AIAwave=[171,193,304]

        self.Notes=''

    def get_index(self,ID,data):
        '''Get index of flare in flare list file'''
        try:
            i= np.where(ID == np.array(data['ID']))
            i=i[0]
        except TypeError:
            i= np.searchsorted(np.array(datadata['ID']), ID)
        return i[0]

    
    def write(self, picklename=False):
        '''Write object to pickle'''
        import pickle
        if not picklename: picklename='/Users/wheatley/Documents/Solar/occulted_flares/data/objects/'+str(self.ID)+'OCObservation.p'
        pickle.dump(self, open(picklename, 'wb'))

    def write_csv(self, csvname=False): #might not flatten the dictionaries correctly
        if not csvname: csvname='/Users/wheatley/Documents/Solar/occulted_flares/data/objects/'+ str(ID) + 'OCObservation.csv'
        d=self.__dict__
        import pandas as pd
        pd.io.json.json_normalize(d).to_csv(csvname)

