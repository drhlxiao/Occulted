 #######################################
# OCFiles.py
# Erica Lastufka 17/5/17 

#Description: Class of functions to deal with all files of the study. Has to be dynamic because I expect there will be a lot of data and plot names yet to come.
#######################################

import numpy as np
import glob

class OCFiles(object):
    ''' object will have the following attributes:
           Naming convention 
           home directory
           subfolders
           special names (xrs, etc)
           Separate images and plots from fits and binary data and spreadsheets/.sav/.p
        Methods:
            generate file names
            test if file exists
            print/returm date file created/modified
            print/return path
    '''

    def __init__(self, ID,legacy=True,filename=False,gen=False, att_dict=False):
        '''Initialize the object, given the flare ID. This requires use of the OCFiles object I think...'''
        if att_dict: #distribute these
            legacy=False
            gen=False
            for k in att_dict.keys():
                setattr(self,k,att_dict[k]) #types?
            
        if legacy: #get info from legacy OCData object and the associated csv/sav files
            self.dir='/Users/wheatley/Documents/Solar/occulted_flares/' #top level directory
            self.folders={'spectrograms':'data/spectrograms','lightcurves':'data/lightcurves','stereo-aia':'data/stereo-aia','stereo_pfloops':'data/stereo_pfloops','ql_images':'data/ql_images','xrs_files':'data/xrs_files','flare_lists':'flare_lists','bproj_vis':'data/bproj_vis','plots':'plots','clean':'data/clean','att_data':'data/att_data'}

            if not filename:
                filename= '/Users/wheatley/Documents/Solar/occulted_flares/flare_lists/list_final.csv'#default file to read
            import pandas as pd
            data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
            i=self.get_index(ID,data) #get the index of the flare if it's in a list
            self.Raw={'xrs_files':data["XRS_files"][i],'aia':'','stereo':'','RHESSI':'','ospex':'','spectrogram':'','messenger':'','goes':'','e4':'','e5':''}
            self.Lists={'flare_list':data["csv_name"][i]}
            self.att_data='' #place where attribute data is stored
            self.Plots={} #all the locations/naming conventions for the plots?

        #if not legacy:
        #    #read datetimes csv file (can just restore the pickle file otherwise. Build this into OCFlare class):
        #    import pandas as pd
        #    if filename: #it's the big one
        #        data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
        #        #i=self.get_index(ID,data) #get the index of the flare if it's in a list
        #        self.Raw=dict()
        #        self.folders=dict()
        #        self.att_data=str(ID)+'_att.p' #place where attribute data is stored
        #        for key in data.keys(): #only do this if it starts with Observation
        #            if key.startswith('Files.') and key !='Files.Raw' and key !='Files.folders':
        #                dat=data[key]
        #                key=key[key.find('.')+1:] #trim it
        #                if key.startswith('Raw'): #do stuff
        #                    key=key[key.find('.')+1:]
        #                    self.Raw[key]=dat.values[0]
        #                if key.startswith('folders'): #do stuff
        #                    key=key[key.find('.')+1:]
        #                    self.folders[key]=dat.values[0]
        #                setattr(self,key,dat.values[0])
        #    if not filename:
        #        filename= '/Users/wheatley/Documents/Solar/occulted_flares/flare_lists/'+str(ID)+'OCFiles.csv' #default file to read
        #        data=pd.read_csv(filename,sep=',', header=0) #column 0 will be NaN because it's text
        #        i=self.get_index(ID,data) #get the index of the flare if it's in a list
        #        self.Raw={'xrs_files':data['Raw']["xrs_files"][i],'aia':data['Raw']["xrs_files"][i],'stereo':data['Raw']["xrs_files"][i],'RHESSI':data['Raw']["xrs_files"][i],'ospex':'','spectrogram':'','messenger':'','e3':'','e4':'','e5':''}
        #        self.Lists={'flare_list':data["csv_name"][i]}
        #        self.att_data=str(ID)+'_att.p' #place where attribute data is stored
        #        self.Plots={} #all the locations/naming conventions for the plots?

        if gen: #generate missing filenames by searching for them
            self.find_ospex(ID)
            self.find_aia(ID,legacy=legacy,filename=filename)
            self.find_stereo(ID,legacy=legacy,filename=filename)
            #self.gen_plotnames(ID) #get a giant dictionary of plot names
            self.gen_lcnames(ID)

        #self.Notes= '' #if I need to write a note about this iteration
 
    def get_index(self,ID,data):
        try:
            i= np.where(ID == np.array(data['ID']))
            i=i[0]
        except TypeError:
            i= np.searchsorted(np.array(datadata['ID']), ID)
        return i[0]


    def find_ospex(self,ID):
        odir=self.dir+'data/dat_files/'
        oname=str(ID)+'aafull_a.fits'
        try:
            res=glob.glob(odir+oname)[0]
        except IndexError:
            res=''
        self.Raw['ospex']=res
        
    def find_aia(self,ID,legacy=False,filename=False):
        odir=self.dir + self.folders['stereo-aia']+'/'
        res=[]
        import OCDatetimes as ocd
        from datetime import datetime as dt
        with ocd.OCDatetimes(ID,legacy=legacy,filename=filename) as dtinst:
            for wave in ['304','171','193']:
                oname= 'aia_lev1_'+wave+'*'+ dt.strftime(dtinst.Messenger_peak,'%Y_%m_%d') +'*.fits' #need to get the datetime and format that...also the wavelength? or make this a list
                try:
                    res.append(glob.glob(odir+oname)[0])
                except IndexError:
                    continue
        self.Raw['aia']=res

    def find_stereo(self,ID,legacy=False,filename=False):
        odir=self.dir + self.folders['stereo-aia']+'/'
        res=[]
        import OCDatetimes as ocd
        from datetime import datetime as dt
        with ocd.OCDatetimes(ID,legacy=legacy,filename=filename) as dtinst:         
            oname= dt.strftime(dtinst.Messenger_peak,'%Y%m%d') +'*.fts'
        try:
            res=glob.glob(odir+oname)[0]
        except IndexError:
            pass
        self.Raw['stereo']=res

    def gen_plotnames(self,ID):
        odir='/Users/wheatley/Documents/Solar/plots/'
        #etc

    def gen_lcnames(self,ID):
        odir='/Users/wheatley/Documents/Solar/occulted_flares/data/'
        mfilename='/'+ str(ID)+'_Mess_lc.p' #something
        self.Raw['messenger']=mfilename
        gfilename= '/'+str(ID)+'_GOES_lc.p' #something
        self.Raw['goes']=gfilename
      
    def rename2convention(self): #fix weird conventions I might have had earlier. 
        #first test if the file exists
        if self.test('Plots','spectrogram',ret=True) == '': #do stuff
            ofile=self.test('Plots','spectrogram',ret=True) 
            date=ofile[ofile.rfind('/'):-9]#parse the filename
            nfile=ofile[:ofile.rfind('/')]+dt.strftime(date,'%Y%m%d') + 'sgram.png' #rearange the date
            os.rename(ofile,nfile)#rename file
            
        if self.test('Plots','spectrogram',ret=True) == '': #next one
            ofile=self.test('Plots','spectrogram',ret=True) 
            date=ofile[ofile.rfind('/'):-9]#parse the filename
            nfile=ofile[:ofile.rfind('/')]+dt.strftime(date,'%Y%m%d') + 'sgram.png' #rearange the date
            os.rename(ofile,nfile)#rename file        

    def test(self,group,key,ret=False):
        fname=getattr(self,group)[key]
        import glob2
        a= glob2.glob(fullpath+'data/*/'+fname)
        if a != '':
            print a + ' exists'
        else:
            print a + ' does not exist!'
        if ret: return a
        
    def write(self, picklename=False):
        '''Write object to pickle'''
        import pickle
        if not picklename:
            picklename=str(self.ID)+'OCFiles.p'
        pickle.dump(self, open(picklename, 'wb'))

    def write_csv(self, csvname=False):
        if not csvname: csvname= str(ID) + 'OCFiles.csv'
        d=self.__dict__
        import pandas as pd
        df = pd.io.json.json_normalize(d)
        df.to_csv(csvname) #use pandas to write csv
