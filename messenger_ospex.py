 #######################################
# messenger_ospex.py
# Erica Lastufka 13/2/17  

#Description: Find peak time intervals in messenger data. Maybe even do the fit in ospex?
#######################################

#######################################
# Usage:

# for default output: python imager_widget.py
######################################
import os
import Datetime.Datetimes as dt
import pidly

class messenger_data(object): #do I need to put the id arg here?

    def __init__(self,id):
        self.datfile=id+'.dat'
        self.lblfile=id+'.lbl'
        self.read_lbl()
        #self.read_dat()
        #self.find_messenger_peak()
 
    #first find the maximum of a given spectrum
    def find_messenger_peak(self,erange = erange):    
        
        #data=o->getdata(class='spex_data')
        #eaxis=o->getaxis(/ct_energy)
        #taxis=o->getaxis(/mean)
        #elist=where( (eaxis ge min(er_high)) AND (eaxis le max(er_high)) )
        #phigh=total(data.data(elist,*),1)
        #tlist=where( (taxis ge min(this_tr)) AND (taxis le max(this_tr)) )
        #phigh=phigh(tlist)
        #this_max=max(phigh,this_max_ind)
        #this_tr0=this_tr
        #this_tr1=taxis(tlist(this_max_ind))+[-8.,8.]
        #this_tr=[max([this_tr0(0),this_tr1(0)]),min([this_tr0(1),this_tr1(1)])]

        self.START_TIME=dt.strftime(self.START_TIME,'%Y-%M-%DT%H:%M:%S')
        self.STOP_TIME=dt.strftime(self.STOP_TIME,'%Y-%M-%DT%H:%M:%S')

        #now find the maxima
        
    def read_lbl(self):
        with open(self.lblfile) as lbl:
            for line in lbl:
                if line.startswith('/***'):
                    continue
                else:
                    key=line[:line.find('=')].rstrip()
                    value = line[line.find('=')+2:-2].rstrip()#take away the quotation marks?
                    #strip a
                    setattr(self,key,value)

    def read_dat(self): #this is binary...
        #do this in pidly, get the structure back and convert to attributes?
        idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
        os.chdir('/Users/wheatley/Documents/Solar/occulted_flares/data/dat_files')
        idl('specfile',self.datfile)
        idl('.compile read_messenger_pds')
        idl('read_messenger_pds, files=specfile, data_str=data_str')
        idl('x=data_str')
        foo=idl.x
        idl.close()

filenames=['xrs2012201']
os.chdir('/Users/wheatley/Documents/Solar/occulted_flares/data/dat_files')
for f in filenames:
    obj = messenger_data(f)
   #find the peak
os.chdir('/Users/wheatley/Documents/Solar/MiSolFA/code/Occulted')          
   #save the results for inspection

