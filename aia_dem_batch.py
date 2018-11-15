 #######################################
#aia_dem_batch.py
# Erica Lastufka 29/03/2018  

#Description: Run Sparse DEM for AIA data in given folders, store, clean up
#######################################

#######################################
# Usage:

######################################

import numpy as np
import numpy.ma as ma

import os
from scipy.ndimage.filters import generic_filter as gf
import matplotlib.cm as cm
import matplotlib.colors as colors
import glob
from sunpy.net import Fido, attrs as a
import sunpy.map
from scipy.io import readsav
import astropy.units as u
from astropy.coordinates import SkyCoord
#import display_aia_dem as disp
import zipfile
import pidly
from datetime import datetime as dt

def dem_from_sav(filename):
    dem_dict=readsav(filename,python_dict=True)
    return dem_dict

def aia_prep(files,zip_old=True):
    '''run AIA prep on given files, clean up'''
    idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
    idl('files',files)
    idl('aia_prep,files,indgen(size(files,/dim)),/do_write_fits,/normalize,/use_ref') #what is input 2? input2 - List of indices of FITS files to read 
    #rename files
    prepped_files=glob.glob('AIA2*.fits')
    for pf in prepped_files:
        wave=pf[-8:-5]
        newf=pf[:3]+'_'+wave+'_'+pf[3:-10]+'.fits'
        #print newf
        os.rename(pf,newf)
    if zip_old:
        #archive level 0 data
        zipf = zipfile.ZipFile('AIALevel0.zip', 'w', zipfile.ZIP_DEFLATED)
        zipdir(files, zipf)
        zipf.close()

def zipdir(files, ziph):
    # ziph is zipfile handle
    root=os.getcwd()
    for f in files:
        ziph.write(os.path.join(root, f))

def group6(files):
    '''group AIA data into blocks of 6: '''
    #get time tags for each file, put in a dict: {filename:x,time:y,closest_files:[(file1,td1),(file2,td2)]}
    fdict={'094':[],'131':[],'171':[],'193':[],'211':[],'335':[]}
    wavelengths=['094','131','171','193','211','335']
    for f in files:
        ttag=dt.strptime(f[8:-5],'%Y%m%d_%H%M%S')
        wave=f[4:7]
        tdict={'filename':f,'time':ttag,'closest_files':[],'tds':[],'group_id':0}
        fdict[wave].append(tdict)
        #find the closest files in other wavelengths
        #for w in fdict.keys():
        #    wavelengths=['094','131','171','193','211','335']
        #    wavelengths.remove(wave) #all the others remain
    for f in fdict['094']: #get the closest file in other wavelengths
        #print f,ttag
        ttag=f['time']
        f['closest_files']=[find_closest_file(ttag,fdict,wave=ww) for ww in wavelengths]
        #don't forget to sort
    groups=[f['closest_files'] for f in fdict['094']]
    #return list of groups ...
    return groups

def find_closest_file(ttag,fdict,wave='094'):
    '''what it sounds like'''
    tvec=[fdict[wave][idx]['time'] for idx in range(0,len(fdict[wave]))]
    fvec=[fdict[wave][idx]['filename'] for idx in range(0,len(fdict[wave]))] #hopefully this doesn't mess up the sorting
    #find closest t in tvec to given ttag
    closest_t= min(tvec, key=lambda x: abs(x - ttag))
    ct_idx=tvec.index(closest_t)
    closest_file=fvec[ct_idx]
    #print closest_t,ct_idx,closest_file
    #return corresponding filename... figure out a way to return timedelta as well
    return closest_file
        
def run_sparse_dem(group,submap):
    idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
    idl('.compile /Users/wheatley/Documents/Solar/DEM/tutorial_dem_webinar/aia_sparse_em_init.pro')
    idl('.compile /Users/wheatley/Documents/Solar/occulted_flares/code/run_sparse_dem.pro')
    idl('group6', group)
    idl('submap',submap)
    idl('result=run_sparse_dem(group6, submap)')
    result=idl.result
    #get out the result
    return result 

if __name__ == '__main__':
    os.chdir('low_cadence_cutout')
    #os.chdir('longbefore')
    files=glob.glob('ssw_cutout*.fts')
    #files=glob.glob('aia_lev1*.fits')
    aia_prep(files,zip_old=False)
    preppedfiles=glob.glob('AIA_*.fits')    
    groups=group6(preppedfiles)
    #print groups #should be sorted already
    do_over=[]
    for g in groups[2:]:
        res=run_sparse_dem(g,[-1220,-800,0,400])
        if res !=1:
            do_over.append(g)
    pickle.dump(do_over,open('do_over.p','wb'))
