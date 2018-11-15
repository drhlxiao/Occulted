 #######################################
# messenger_goes_class.py
# Erica Lastufka 27/07/2018 

#Description: Determine GOES class of messenger flare from T and EM
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
import scipy.constants as sc

def messenger_goes_class(tlist,emlist):
    #in idl:
    idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
    idl('tlist',tlist)
    idl('emlist',emlist)
    idl('run_goes_fluxes, tlist,emlist,flong,fshort')
    flong=idl.flong
    fshort=idl.fshort
    idl.close()
    return flong,fshort
