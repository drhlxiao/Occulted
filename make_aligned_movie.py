 #######################################
#make_aligned_movie.py
# Erica Lastufka 28/02/2018  

#Description: Because jhelioviewer was down
#######################################

#######################################
# Usage:

######################################

import numpy as np
import numpy.ma as ma
import scipy.constants as sc
from datetime import datetime as dt
from datetime import timedelta as td
import os
from scipy.ndimage.filters import generic_filter as gf
from plot_composite import composite_plot
import matplotlib.cm as cm

import glob
from sunpy.net import Fido, attrs as a
import sunpy.map
import astropy.units as u
from astropy.coordinates import SkyCoord
from ffmpy import FFmpeg #ffmpeg python wrapper

def get_maps_and_headers(path,tags=False,number=False):
    '''Convert fits files into maps'''
    os.chdir(path)
    if tags:
        files=glob.glob('*'+tags+'*.fts')
    else:
        files=glob.glob('*.fts')
    files.sort()
    if not number:
        maps=sunpy.map.Map(files)
        headers=[sunpy.io.fits.get_header(f) for f in files]
    else:
        maps=sunpy.map.Map(files[:number])
        headers=[sunpy.io.fits.get_header(f) for f in files[:number]]
    return maps,headers
    
def align_maps(m1,m2,a1=1,a2=1): #should m1 and m2 be the masked arrays?
    comp_map=sunpy.map.Map(m1,m2,composite=True)
    comp_map.set_alpha(0,a1)
    comp_map.set_alpha(1,a2)
    return comp_map

def difference_map(m1,m2,h1,h2,logscale=True,lowmask=False,highmask=False):
    data=m1.data-m2.data
    if lowmask:
        if logscale:
            lowmask=10**lowmask
        if lowmask < 1:
            data=np.ma.masked_less(data,lowmask*np.max(data))
        else: #it's just a pixel value
            data=np.ma.masked_less(data,lowmask)
    if highmask:
        if logscale:
            highmask=10**highmask
        if highmask < 1:
            data=np.ma.masked_greater(data,highmask*np.max(data))
        else: #it's just a pixel value
            data=np.ma.masked_greater(data,highmask)
    if logscale:
        data=np.log10(data)
    header=h1
    diff_map=sunpy.map.Map((data,header))
    return diff_map

def running_diff(maps,headers,logscale=True,lowmask=False,highmask=False):
    '''make base difference map. basically a wrapper for difference_map when applied to a bunch of maps'''
    outmaps=[]
    for j,m,h in zip(enumerate(maps),headers):
        while j < len(m[:-1]):
            diff_map=difference_map(maps[j+1],m,headers[j+1],h,logscale=logscale,lowmask=lowmask,highmask=highmask)
            outmaps.append(diff_map)
    return outmaps
    
def base_diff(maps,headers,logscale=True,lowmask=False,highmask=False):
    '''make base difference map. basically a wrapper for difference_map when applied to a bunch of maps'''
    bm=maps[0]
    bh=headers[0]#[0] #need to remove this for not-COR (also for h below)
    outmaps=[]
    for m,h in zip(maps,headers):
        diff_map=difference_map(m,bm,h,bh,logscale=logscale,lowmask=lowmask,highmask=highmask)
        outmaps.append(diff_map)
    return outmaps

###############Image Processing###################

def mask_circular(imap, radius,less=True,blue=False):
    '''Make a circular mask of given radius centered at (a,b)'''
    x, y = np.meshgrid(*[np.arange(v.value) for v in imap.dimensions]) * u.pixel
    hpc_coords = imap.pixel_to_world(x, y)
    r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / (radius*imap.rsun_obs)
    if less:
        mask = ma.masked_less_equal(r, 1)
    else:
        mask = ma.masked_greater_equal(r, 1)        
    palette = imap.plot_settings['cmap']
    if blue:
        palette=cm.Blues
    palette.set_bad(alpha=0.0)
    return mask,palette

def mask_solar_disk(imap,radius=1.5,plot=False):
    '''mask the solar disk in the coronograph images. Units of Rsun.'''
    mask,palette=mask_circular(imap,radius)
    if plot:
        plot_with_mask(imap,mask,palette,show=True,save=False)

def mask_stereo_circular(imap,radius=1.5,plot=False):
    '''Mask the STEREO image OUTSIDE a certain radius. Units of Rsun.'''
    mask,palette=mask_circular(imap,radius,less=False)
    if plot:
        plot_with_mask(imap,mask,palette,show=True,save=False)

#def scale_palette(map,palette,plow=2,phigh=98,log=False):
#    '''Scale the color palette to best display the data'''
#    pl = np.percentile(map.data, plow)
#    ph = np.percentile(map.data, phigh)
#    im2 = exposure.rescale_intensity(masked_array, in_range=(pl, ph)) #also a masked array
    
#    return scaled_palette #returnes matplotlib.colors.**map object

def process_map(imap,mask,type=False,plot=False,blue=False):
    mmask=[]
    for i,m in enumerate(mask):
        #print i,m
        if m > 0:
            if i==1:
                less=False
            else:
                less=True
            foo=mask_circular(imap,m,less=less,blue=blue)
            ppalette=foo[1]
            mmask.append(foo[0])
    #if more than 2 masks, combine
    #print 'here',np.shape(mmask)
    if np.shape(mmask)[0] > 1:
        tmask=mmask[0]+mmask[1]
    else:
        tmask=mmask
    #scale the color palette
    #spalette=scale_palette(ppalette)
    if plot:
        plot_with_mask([imap],[tmask],ppalette,show=True,save=False)            
    return tmask,ppalette #scaling should be done in the Axes object?
#################################################

def plot_with_mask(maplist,mask,palette,save=True,show=False):
    scaled_maps=[]
    for m,mmask in zip(maplist,mask):
        try:
            #scaled_maps.append(sunpy.map.Map(m.data, m.meta, mask=mmask.mask))
            #masked_map=ma.fix_invalid(m.data,mmask.mask)
            scaled_maps.append(sunpy.map.Map(m.data, m.meta, mask=mmask.mask))
        except AttributeError:
            #masked_map=ma.fix_invalid(m.data,mmask[0].mask)
            scaled_maps.append(sunpy.map.Map(m.data, m.meta, mask=mmask[0].mask))
            
    if len(scaled_maps) > 1:
        scaled_map=align_maps(scaled_maps[0],scaled_maps[1])#[0]#align_maps(scaled_maps[0],scaled_maps[1])
    else:
        scaled_map = scaled_maps[0]
    fig = plt.figure()
    #plt.subplot(projection=scaled_map)
    #can I add in the axes scaling here?
    composite_plot(scaled_map)
    #scaled_map[1].plot(cmap=palette)
    ax=plt.gca()
    #scaled_map[1].plot(ax,cmap=palette)
    #scaled_map.draw_limb()
    #fig,ax=plt.subplots()
    #ax.imshow(np.ma.masked_where(masked_map>0,masked_map))
    #ax.imshow(scaled_maps[0].data)

    if show:
        fig.show()
    if save:
        return fig
    
def sunpy2png(maplist,outname,mask=False,outpath=False,ext='png',colorbar=True,lownorm=False,highnorm=False,big=True):
    '''input to mask kward is [[list of masks], palette]'''
    if big:
        plt.figure(figsize=(16,12))
    plt.clf()        
    if lownorm or highnorm:
        if not lownorm:
            norm=plt.Normalize(0,highnorm)
        else:
            norm=plt.Normalize(lownorm,highnorm)
    if mask:
        fig=plot_with_mask(maplist,mask[0],mask[1])
    else:
        try:
            fig=maplist[0].plot(cmap=cm.rainbow,norm=norm)
        except UnboundLocalError:
            fig=maplist[0].plot(cmap=cm.rainbow)
    if not outpath:
        outpath='.'
    if colorbar:
        plt.colorbar()            
    plt.savefig(outpath+'/'+outname+ext)
    
def run_ffmpeg(imwild,mname,mdir='.',framerate=24):
    if mdir != '.':
        os.chdir(mdir)
    ff=FFmpeg(inputs={imwild:None},outputs={mname:'-framerate '+str(framerate)})
    #ff.cmd()
    ff.run()

def make_difference_movie(ddir,mtype,moviename,mdir=False,tags=False,n=False,imtags=False,imdir=False,ext='png',submap=False,logscale=True,lownorm=False,highnorm=False,colorbar=True,big=True):
    #get the maps. submap is in form [(xrange_tuple),(yrange_tuple)]. low and highmask are in percentage of max pixel or absolute pixel cutoffs (>1)
    cwd=os.getcwd()
    maps1,h1=get_maps_and_headers(ddir,tags=tags,number=n)
    if submap:
        submaps,h1=[],[]
        for m in maps1:
            bl=submap[0]*u.arcsec#(submap[0][0]*u.arcsec,submap[0][1]*u.arcsec)
            tr=submap[1]*u.arcsec#(submap[1][0]*u.arcsec,submap[1][1]*u.arcsec)
            sm=m.submap(SkyCoord((bl),(tr),frame=m.coordinate_frame))
            submaps.append(sm)
            h1.append(sm.meta)
        maps1=None
        maps1=submaps
    if mtype=='base':
        outmaps=base_diff(maps1,h1,logscale=logscale,lowmask=lownorm,highmask=highnorm)
    elif mtype == 'running':
        outmaps=running_diff(maps1,h1,logscale=logscale,lowmask=lownorm,highmask=highnorm)
    if not imdir:
        imdir='.'
    #convert aligned maps to png for use in ffmpeg
    for i,om in enumerate(outmaps):
        if imtags:
            outname='map'+imtags+'_'+'{0:03d}'.format(i)+'.'
        else:
            outname='map'+'{0:03d}'.format(i)+'.'
        #sunpy2png([om],outname,**kwargs)            
        sunpy2png([om],outname,outpath=imdir,ext=ext,colorbar=colorbar,lownorm=lownorm,highnorm=highnorm,big=big)
    os.chdir(imdir)
    #if imtags:
    #    imlist=glob.glob('*'+imtags+'*.'+ext)
    #else:
    #    imlist=glob.glob('*.'+ext)
    #make movie
    run_ffmpeg('map'+imtags+'_%d.png',moviename)
    if not mdir:
        mdir='.'
    print 'movie ', moviename, ' is stored in ', mdir
    
def make_aligned_movie(dir1,dir2,moviename,tags1=False,tags2=False,n1=False,n2=False,submap=False,imdir=False,ext='png',mdir='.',framerate=24,mask1=[0,1.5],mask2=[1.5,0]):
    #get the maps
    cwd=os.getcwd()
    maps1,h1=get_maps_and_headers(dir1,tags=tags1,number=n1)
    os.chdir(cwd)
    maps2,h2=get_maps_and_headers(dir2,tags=tags2,number=n2)
    os.chdir(cwd)
    #process maps - mask out unnecessary bits ##, contrast stretch, fix color palette
    masks1,masks2=[],[]
    for m in maps1:
        masks1.append(process_map(m,mask1)[0])
    for m in maps2:
        masks2.append(process_map(m,mask2,blue=True)[0])
        palette=process_map(m,mask2)[1]
    #align maps
    aligned_maps=[]
    #if submap requested, treat that here
    if submap: #pretty useless for coronograph so leave it out
        submaps=0
    if not imdir:
        imdir='.'
    i=0
    for m1,m2,mm1,mm2 in zip(maps1,maps2,masks1,masks2): #maybe do this later.... 
        #comp_map=align_maps(m1,m2)
        #aligned_maps.append(comp_map)
        outname='map'+str(i)+'.'
        #fig=plot_with_mask([m1,m2],[mm1[0],mm2],palette,show=False,save=False)
        
        sunpy2png([m1,m2],outname, mask=[[mm1[0],mm2],palette],outpath=imdir,ext=ext)
        i+=1
    os.chdir(imdir)
    #imlist=glob.glob('*.'+ext)
    #make movie
    run_ffmpeg('map%d.'+ext,moviename,framerate=framerate)
    print 'movie ', moviename, ' is stored in ' , mdir

    
    
