 #######################################
# flare_class.py
# Erica Lastufka 17/5/17 

#Description: Class of functions to deal with all data,files,calculations of the study. Maybe merge this with calculations since they are so interrelated
#######################################

import numpy as np

class OCProperties(object):
    ''' Object will have the following attributes:
            T1: 
            EM1:
            T2: 
            EM2: 
            chisq: chisq value from ospex fit
            GOES_GOES: GOES flux 
            Messenger_GOES: Same for Messenger
            source_vis:
            source_pos:
            source_pos_disk: source position projected down onto disk if necessary
            source_pos_STEREO: 
            f1pos: in case this is different than what is used in the loops?
            f2pos:  
            loop: dictionary of loop data
            Rheight: Calculated rhessi height above limb
            occ_percent: occultation percentage
            ratio: GOES ratio
            e1: empty field for addition later (in case I can't figure out how to make this dynamic)
            e2
            e3
            e4
            e5
            e6
            e7
            e8
            e8
            e10
            Notes
        Methods are:
            all the calculation methods...'''
    
    def __init__(self, ID, legacy=True, filename=False,calc_missing=True):
        '''Return an empty Data_Study object, unless a filename to a correctly formatted csv is given. Do I need to be able to go from a dictionary to an object too?'''
        tags=['T2','EM2','chisq','source_pos_STEREO','f1pos','f2pos','loop','Rheight','occ_percent','source_pos_disk','e2','e3','e4','e5','e6','e7','e8','e9','e10']
        if legacy: #get info from legacy OCData object and the associated csv/sav files
            if not filename:
                filename= '/Users/wheatley/Documents/Solar/occulted_flares/flare_lists/list_final.csv'#default file to read
            import pandas as pd
            data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
            i=self.get_index(ID,data) #get the index of the flare if it's in a list
            self.T1=data['Messenger_T'][i]
            self.EM1=data['Messenger_EM1'][i]
            self.GOES_GOES=data['GOES_GOES'][i]
            self.Messenger_GOES=data['Messenger_GOES'][i]
            self.source_vis=data['source_vis'][i]
            self.source_pos=data['Location'][i]
            if not self.source_pos.startswith('['): #trim off all the weird junk
                aa=self.source_pos
                self.source_pos=[float(aa[aa.find('[')+1:aa.find(';')]),float(aa[aa.find(';')+1:aa.find(']')])]
            self.ratio=self.Messenger_GOES/self.GOES_GOES
            for att,key in zip(tags,tags):
                setattr(self,att,'')
 
        if not legacy:
            #read attributes from csv file (can just restore the pickle file otherwise. Build this into OCFlare class):
            import pandas as pd
            if filename: #it's the big one
                data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
                #i=self.get_index(ID,data) #get the index of the flare if it's in a list
                for key in data.keys(): #only do this if it starts with Observation
                    if key.startswith('Properties.'):
                        dat=data[key]
                        key=key[key.find('.')+1:] #trim it
                        setattr(self,key,dat.values[0])
                #self.string_to_dict() #ast doesn't like my data for some reason... the coordinate pairs?
            if not filename:
                filename= '/Users/wheatley/Documents/Solar/occulted_flares/flare_lists/'+str(ID)+'OCProperties.csv'#default file to read - need an except in case it doesn't exist
                data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
                i=self.get_index(ID,data) #get the index of the flare if it's in a list
                for key in data.keys():
                    setattr(self,key,data.values[0])
                #self.string_to_dict()

        if calc_missing: #caluclate missing values if possible
            self.get_T_EM(ID)
            self.convert2stereo()
            self.readloops(ID)
            #etc
            
        self.Notes=''

    def get_index(self,ID,data):
        '''Get index of flare in flare list file'''
        try:
            i= np.where(ID == np.array(data['ID']))
            i=i[0]
        except TypeError:
            i= np.searchsorted(np.array(datadata['ID']), ID)
        return i[0]



    def string_to_dict(self):
        '''convert any string dictionaries to proper dictionaries'''
        import ast
        #check what's a dictionary
        d=dir(self)
        for att in d:
            attv=getattr(self,att)
            if type(attv) == str and attv.startswith('{'):
                setattr(self,att,ast.literal_eval(attv))
    
    def write(self, picklename=False):
        '''Write object to pickle'''
        import pickle
        if not picklename:
            picklename=str(self.ID)+'OCProperties.p'
        pickle.dump(self, open(picklename, 'wb'))

    def write_csv(self, csvname=False): #might need to flatten dictionaries here because of loops
        #replace , in source_pos and source_pos_stereo with ;
        if ',' in self.source_pos:
            self.source_pos.replace(',',';')
        if ',' in self.source_pos_disk:
            self.source_pos_disk.replace(',',';')
        if ',' in self.source_pos_STEREO:
            self.source_pos_STEREO.replace(',',';')
        
        if not csvname: csvname= str(ID) + 'OCProperties.csv'
        d=self.__dict__
        import pandas as pd
        df = pd.io.json.json_normalize(d)
        df.to_csv(csvname)  #dict(orient='records')

    def get_T_EM(self,ID): #read in the fit... use IDL cuz I'm too lazy to sift through the whole .fits file for the right keywords
        #put in error exception
        import pidly
        idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
        idl.filename = str(ID)
        idl("aa=spex_read_fit_results('/Users/wheatley/Documents/Solar/occulted_flares/data/dat_files/'+filename+'aall_a.fits')")
        params = idl.ev('aa.spex_summ_params')
        self.T1=params[1]
        self.T2=params[4]
        self.EM1=params[0]
        self.EM2=params[3]
        self.chisq=idl.ev('aa.spex_summ_chisq')[0]
        
    def convert2stereo(self): #convert RHESSI/AIA source location to stereo location (roughly, since if it's off the limb that will be hard) (How exactly should I go about this?)
        self.source_pos_STEREO = 0.0
        
    def readloops(self,ID,filename=False): #read loops spreadsheet and store values in self. 
        if not filename:
            filename= '/Users/wheatley/Documents/Solar/occulted_flares/data/stereo_pfloops/loops.csv'#default file to read
        import pandas as pd
        data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
        i = np.where(data['ID'].values == ID)[0][0]#the index of the right flare
        #do some formatting of the data before placing it into the dictionary...
        def format_fp(fstr):
            if fstr=='' or type(fstr) != str:
                return fstr
            arr=fstr.split(' ')
            newarr=[]
            for a in arr:
                try:
                    newarr.append(float(a))
                except ValueError:
                    if (a.startswith('x') or a.startswith('y')) and len(a) > 4:
                        newarr.append(float(a[2:]))
                    else:
                        pass
            coords=[[newarr[0],newarr[1]],[newarr[2],newarr[3]]]
            return coords

        fp1=format_fp(data['Footpoints1'][i])
        fp2=format_fp(data['Footpoints2'][i])
        
        def format_coords(cstr):
            if cstr=='' or type(cstr) != str:
                return cstr
            coords=[float(cstr[cstr.find('[')+1:cstr.find(';')]),float(cstr[cstr.find(';')+1:cstr.find(']')])]
            return coords
        
        center=format_coords(data['center'][i])
        e1=format_coords(data['ellipse2'][i])
        e2=format_coords(data['ellipse2'][i])
        self.loop={'footpoints1':fp1,'length1':data['length1'][i],'ellipse1':e1,'eccentricity1':data['eccentricity1'][i],'Footpoints2':fp2,'length2':data['length2'][i],'ellipse2':e2,'eccentricity2':data['eccentricity2'][i],'center':center,'notes':data['notes'][i]}
