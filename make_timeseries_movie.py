 #######################################
#make_timeseries_movie.py
# Erica Lastufka 5/03/2018  

#Description: Make edge-enhanced timeseries of AIA data following example here: http://docs.sunpy.org/en/stable/generated/gallery/loop_edge_enhance.html
#######################################

#######################################
# Usage:

######################################
from __future__ import print_function, division

import numpy as np
import numpy.ma as ma
import scipy.constants as sc
import os
from scipy import ndimage
from scipy.ndimage.filters import generic_filter as gf
import matplotlib.cm as cm
import make_aligned_movie as mov
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal


import glob
import sunpy.map
import astropy.units as u
from astropy.coordinates import SkyCoord
from ffmpy import FFmpeg #ffmpeg python wrapper

def edge_filter(smap):
    sx = ndimage.sobel(smap.data, axis=0, mode='constant')
    sy = ndimage.sobel(smap.data, axis=1, mode='constant')
    edge_enhanced_im = np.hypot(sx, sy)
    edge_map = sunpy.map.Map(edge_enhanced_im, smap.meta)
    return edge_map

def get_submap(smap,bl,tr):
    bl = SkyCoord(bl[0] * u.arcsec, bl[1] * u.arcsec, frame=smap.coordinate_frame)
    tr = SkyCoord(tr[0] * u.arcsec, tr[1] * u.arcsec, frame=smap.coordinate_frame)
    ssmap = smap.submap(bl, tr)
    return ssmap

def difference_map2(m1,m2):
    data=m1.data-m2.data
    diff_map=sunpy.map.Map((data,m1.meta))
    return diff_map

def plot_log_cmap(maps):
    '''generate a logarithmic colormap based on a given minimum and maximum data value'''
    #N = 100
    #X, Y = np.mgrid[-3:3:complex(0, N), -2:2:complex(0, N)]
    fig,ax=plt.subplots()
    palette_name=maps.plot_settings['cmap']
    im = ndimage.gaussian_filter(maps.data, 3)
    levels=[np.min(im),.6*np.min(im)+1,.6*np.max(im)-1,np.max(im)]
    pcm = ax.contour(im, levels=levels, cmap=palette_name)#,vmin=maps.min(), vmax=maps.max(), cmap=palette_name)
    #pcm = ax.pcolormesh(maps.data,
    #                   norm=colors.SymLogNorm(linthresh=0.0001, linscale=0.0001,vmin=maps.min(), vmax=maps.max()),
    #                   cmap=palette_name)
    #maps.plot(pcm)
    fig.colorbar(pcm, ax=ax, extend='max')

    #pcm = ax[1].pcolor(X, Y, Z1, cmap='PuBu_r')
    #fig.colorbar(pcm, ax=ax[1], extend='max')
    fig.show()
    #return pcm

#def get_extra_maps(files):
#        maps=sunpy.map.Map(files)
#        return maps

def fname2png(fname, imname,submap=False,overwrite=False,edges=True):
    ''' process one file->map->png+pickle at a time'''
    #check if any of the map images already exist
    maps=sunpy.map.Map(fname)
    if submap:
        bl=submap[0]
        tr=submap[1]
        mapss=get_submap(maps,bl,tr) 
    else:
        mapss=maps

    #pickle
    if edges:
        mapss=edge_filter(mapss)
        pname=fname[:-4]+'_edges.p'
    else:
        pname=fname[:-4]+'.p'
    pickle.dump(mapss,open(pname,'wb'))
    plt.clf()
    mapss.plot()
    plt.savefig(imname)

        
def make_timeseries_movie(ddir,moviename,framerate=12,tags=False,n=False,imtags=False,imdir=False,ext='png',submap=False,overwrite=False,edges=True):
    '''Make the timeseries movie. input for kwarg submap is [(bottom_left_tuple), (top_righ_tuple)]. should probably mpi this '''
    #check if any of the map images already exist
    cwd=os.getcwd()
    if imtags:
        imnames='map_'+imtags
    else:
        imnames='map_'
    if not tags:
        tags='*'
    if not imdir:
        imdir='.'
    os.chdir(ddir)
    all_fits=glob.glob('*'+tags+'*.fts') #or .fits
    fitsnames=[af[:-4] for af in all_fits]
    all_ims=glob.glob(imnames+'*.'+ext)
    all_pickles=glob.glob('*'+tags+'*.p')
    len_fits=len(all_fits)
    len_ims=len(all_ims)
    if len_fits != len_ims: #get the indexes of where the png images end? assuming only later images will be added
        overwrite=True
        remaining_fits=all_fits[len_ims:]
    else:
        remaining_fits=all_fits
    if overwrite:
        #get the maps
        for i,f in enumerate(remaining_fits):
            if edges:
                pname=f[:-4]+'_edges.p'
            else:
                pname=f[:-4]+'.p'                
            imname=imdir+'/'+imnames+str(i+len_ims)+'.'+ext
            if pname in all_pickles:
                #restore and make the image from here
                pmap=pickle.load(open(pname,'rb'))
                plt.clf()
                pmap.plot()
                plt.savefig(imname)
            else:
                fname2png(f,imname,submap=submap,overwrite=overwrite,edges=edges)
            
    os.chdir(imdir)
    #make movie. add try and except if images not found... 
    mov.run_ffmpeg(imnames+'%d.png',moviename,framerate=framerate)
    #print 'movie ', moviename, ' is stored in ', mdir
        
