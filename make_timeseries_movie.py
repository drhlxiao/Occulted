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
import make_timeseries_plot_flux as ts
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
import matplotlib.ticker as tk
import matplotlib.pyplot as plt

from datetime import datetime as dt
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

def check_masks(check_flux=True):
    '''these are the images following convol()'''
    adata=readsav('norh_convol.sav',python_dict=True)
    masklist=pickle.load(open('masks.p','rb'))
    f34=adata['f34_c']
    f17=adata['f17_c']

    for i,m17,m34,mask in zip(range(0,len(f17)),f17,f34,masklist):
        imname='norhcheck'+'{0:03d}'.format(i)+'.png'
        fig=plt.figure(figsize=[12,5])
        ax1=plt.subplot(121)
        ax2=plt.subplot(122)
        norm1=colors.Normalize(vmin=np.min(m17),vmax=np.max(m17))
        norm2=colors.Normalize(vmin=np.min(m34),vmax=np.max(m34))
        if check_flux:
            ax1.imshow(m17,cmap=sunpy.cm.cmlist['sdoaia304'])
            ax2.imshow(m34,cmap=sunpy.cm.cmlist['sdoaia304'])
        else:
            ax1.imshow(np.ma.masked_array(m17[80:185,0:125],mask=mask),cmap=sunpy.cm.cmlist['sdoaia304'],norm=norm1)
            ax2.imshow(np.ma.masked_array(m34[80:185,0:125],mask=mask),cmap=sunpy.cm.cmlist['sdoaia304'],norm=norm2)
            
        fig.savefig(imname)
        plt.clf()
    mov.run_ffmpeg('norhcheck'+'%03d.png','check_flux.mp4',framerate=12)
 
    

def norh_together(maplist17,maplist34,moviename='norh_together.mp4',cmap='sdoaia304',framerate=12,bsub=0,norm1=[2500,25000],norm2=[0,10000],mask_disk=True,creverse=True,rhessi='soft'):
    #maplist17.sort()
    #maplist34.sort()
    #build the pre-flare image for subtraction
    bl17=maplist17[:bsub]
    bl34=maplist34[:bsub]
    bav17=np.sum([m.data for m in bl17],axis=0)/float(bsub)
    bav34=np.sum([m.data for m in bl34],axis=0)/float(bsub)
    print( np.shape(bav17),np.mean(bav17))
    norhtimes17,norhtimes34=pickle.load(open('lcdates.p','rb'))
    
    if rhessi == 'soft':
        rhessi_files=glob.glob('../Xray/rhsi_low*.fits')
    elif rhessi == 'hard':
        rhessi_files=glob.glob('../Xray/rhsi_high*.fits')
        
    try:    
        rhessi_maps=[sunpy.map.Map(f) for f in rhessi_files]
        rtimes,srhsi=[],[]
        for rhsi_map in rhessi_maps:
            rhsi_map.meta['CTYPE1'] = 'HPLN-TAN'
            rhsi_map.meta['CTYPE2'] = 'HPLT-TAN'
            rhsi_map.meta['cunit1']='arcsec'
            rhsi_map.meta['cunit2']='arcsec'
            rhsi_map.meta['date-obs']=rhsi_map.meta['date-obs']
            rhsi_map.plot_settings['cmap']=cm.Greens
            rhsi_map.rotate()
            rtimes.append(dt.strptime(rhsi_map.meta['date-obs'][:-4],'%Y-%m-%dT%H:%M:%S'))
            srhsi.append(rhsi_map.submap(SkyCoord([-1410,-750]*u.arcsec,[-50,500]*u.arcsec,frame=rhsi_map.coordinate_frame)))

        #make a of -1 (no rhessi) or index (rhessi) corresponding to NORH times
        rindex17,rindex34=[],[]
        current_idx=0
        next_idx=1
        for t17 in norhtimes17[bsub:]:
            if next_idx >= len(rtimes):
                rindex17.append(-1)
            elif t17 < rtimes[current_idx]:
                rindex17.append(-1)
            elif t17 < rtimes[next_idx]: #t17 >= rtimes[current_idx] and 
                rindex17.append(current_idx)
                #print(t17,rtimes[current_idx],current_idx)
            elif t17 >= rtimes[next_idx]:
                current_idx+=1
                next_idx+=1
                rindex17.append(current_idx)
        current_idx=0
        next_idx=1
        for t34 in norhtimes34[bsub:]:
            if next_idx >= len(rtimes):
                rindex34.append(-1)
            elif t34 < rtimes[current_idx]:
                rindex34.append(-1)
            elif t34 < rtimes[next_idx]:
                rindex34.append(current_idx)
            elif t34 >= rtimes[next_idx]:
                current_idx+=1
                next_idx+=1
    except NameError:
        pass
                
    for i,m17,m34 in zip(range(0,len(maplist17[bsub:])),maplist17[bsub:],maplist34[bsub:]):
        imname='norh2map'+'{0:03d}'.format(i)+'.png'
        fig=plt.figure(figsize=[12,5])
        ax1=plt.subplot(121)
        ax2=plt.subplot(122)
        #subtract background
        if bsub !=0:
            new17=m17.data-bav17
            meta17=m17.meta
            new34=m34.data-bav34
            meta34=m34.meta
            m17=sunpy.map.Map(new17,meta17)
            m17.meta['CTYPE1'] = 'HPLN-TAN'
            m17.meta['CTYPE2'] = 'HPLT-TAN'
            #m17.plot_settings['cmap']=cmap
            m17.rotate()
            m34=sunpy.map.Map(new34,meta34)
            m34.meta['CTYPE1'] = 'HPLN-TAN'
            m34.meta['CTYPE2'] = 'HPLT-TAN'
            #m17.plot_settings['cmap']=cmap
            m34.rotate()
            m17s=m17.submap(SkyCoord([-1410,-750]*u.arcsec,[-50,500]*u.arcsec,frame=m17.coordinate_frame))
            m34s=m34.submap(SkyCoord([-1410,-750]*u.arcsec,[-50,500]*u.arcsec,frame=m34.coordinate_frame))
            #print(np.min(new17),np.max(new17))

        palette_name=sunpy.cm.cmlist[cmap]#m17.plot_settings['cmap']
        if creverse and not palette_name.name.endswith('_r'): #true means reverse the color table
            #palette_name=mapss.plot_settings['cmap']
            new_cdata=cm.revcmap(palette_name._segmentdata)
            new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name+'_r',new_cdata)
        elif not creverse and palette_name.name.endswith('_r'):
            new_cdata=cm.revcmap(palette_name._segmentdata)
            new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name[:-2],new_cdata)

        if mask_disk:
            masked_map17=ts.mask_disk(m17s,fac=10.,plot=False)
            masked_map34=ts.mask_disk(m34s,fac=10.,plot=False)
            m17s=masked_map17#.filled(0)
            m34s=masked_map34#.filled(0)
    
        m17s.plot_settings['cmap']=new_cmap
        m17s.plot_settings['norm']=colors.Normalize(norm1[0],norm1[1])
        m34s.plot_settings['cmap']=new_cmap
        m34s.plot_settings['norm']=colors.Normalize(norm2[0],norm2[1])

        try:
            if rindex17[i] != -1: #make composite map with RHESSI
                #print(rindex17[i])
                comp_map17 = sunpy.map.Map(m17s,srhsi[rindex17[i]], composite=True)
                comp_map17.set_alpha(0, 1)
                comp_map17.set_alpha(1, 1)
                comp_map17.set_levels(1, [50,51,69,70,89,90],percent=True)
                plt.sca(ax1)
                #ax1.set_xticklabels(default_xticks)
                #ax1.set_yticklabels(default_yticks)
                ax1.get_xaxis().set_visible(False)
                ax1.get_yaxis().set_visible(False)
                rplot_settings=comp_map17.get_plot_settings(1)
                rplot_settings['cmap']=cm.Greens
                comp_map17.set_plot_settings(1, rplot_settings)
                comp_map17.plot(title='17 GHz 2013-05-01 '+m17s.meta['time-obs'])
            else:
                m17s.plot(axes=ax1,title='17 GHz 2013-05-01 '+m17s.meta['time-obs'])
            
            if rindex34[i] != -1:
                comp_map34 = sunpy.map.Map(m34s,srhsi[rindex34[i]], composite=True)
                comp_map34.set_alpha(0, 1)
                comp_map34.set_alpha(1, 1)
                comp_map34.set_levels(1, [50,51,69,70,89,90],percent=True)
                plt.sca(ax2)
                #ax2.set_xticklabels(default_xticks)
                #ax2.set_yticklabels(default_yticks)
                ax2.get_xaxis().set_visible(False)
                ax2.get_yaxis().set_visible(False)
                rplot_settings=comp_map34.get_plot_settings(1)
                rplot_settings['cmap']=cm.Greens
                comp_map34.set_plot_settings(1, rplot_settings)

                comp_map34.plot(title='34 GHz 2013-05-01 '+m34s.meta['time-obs'])
            else:
                m34s.plot(axes=ax2,title='34 GHz 2013-05-01 '+m34s.meta['time-obs'])
        except NameError:
            m17s.plot(axes=ax1,title='17 GHz 2013-05-01 '+m17s.meta['time-obs'])
            m34s.plot(axes=ax2,title='34 GHz 2013-05-01 '+m34s.meta['time-obs'])

        #fig.show()
        plt.savefig(imname)
        plt.clf()
    #mov.run_ffmpeg('norhmap'+'%03d.png',moviename,framerate=framerate)
    #return rindex17,srhsi

def norh_on_304(maplist17,file304='AIA/304/mapcube304_peak.p',moviename='norh_on_304.mp4',framerate=12,bsub=0,levels=[19,20,29,30,39,40],creverse=True):
    #maplist17.sort()
    #maplist34.sort()
    #build the pre-flare image for subtraction
    maplist304=pickle.load(open(file304,'rb'))
    times304=[dt.strptime(m.meta['date-obs'][:-3],'%Y-%m-%dT%H:%M:%S') for m in maplist304.maps]
    bl17=maplist17[:6]
    bav17=np.sum([m.data for m in bl17],axis=0)/float(6.)
    #norhtimes17,norhtimes34=pickle.load(open('NORH/lcdates.p','rb'))
    
    submaps17=[]
    norhtimes17=[]
    for m17 in maplist17[bsub:]:
        if bsub !=0:
            new17=m17.data-bav17
            meta17=m17.meta
            m17=sunpy.map.Map(new17,meta17)
            m17.meta['CTYPE1'] = 'HPLN-TAN'
            m17.meta['CTYPE2'] = 'HPLT-TAN'
            #m17.plot_settings['cmap']=cmap
            m17.rotate()
            m17s=m17.submap(SkyCoord([-1210,-750]*u.arcsec,[-50,500]*u.arcsec,frame=m17.coordinate_frame))
            submaps17.append(m17s)
            norhtimes17.append(dt.strptime(m17s.meta['date-obs']+'T'+m17s.meta['time-obs'][:-4],'%Y-%m-%dT%H:%M:%S'))
    
    #make a of -1 (no rhessi) or index (rhessi) corresponding to NORH times
    rindex17=range(0,len(norhtimes17))
    current_idx=0
    next_idx=1
    for i,t17 in enumerate(norhtimes17):
        closest_t= min(times304, key=lambda x: abs(x - t17))
        ct_idx=times304.index(closest_t)
        rindex[ct_idx]=i
    print(rindex17)
    #return rindex17,times304,norhtimes17,maplist304,submaps17

        
    for i,m304 in enumerate(maplist304[:len(submaps17)]):
        s304=m304.submap(SkyCoord([-1210,-750]*u.arcsec,[-50,500]*u.arcsec,frame=m304.coordinate_frame))
        palette_name=m304.plot_settings['cmap']#m17.plot_settings['cmap']
        if creverse and not palette_name.name.endswith('_r'): #true means reverse the color table
            #palette_name=mapss.plot_settings['cmap']
            new_cdata=cm.revcmap(palette_name._segmentdata)
            new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name+'_r',new_cdata)
        elif not creverse and palette_name.name.endswith('_r'):
            new_cdata=cm.revcmap(palette_name._segmentdata)
            new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name[:-2],new_cdata)

        s304.plot_settings['cmap']=new_cmap    
        imname='norh304_'+'{0:03d}'.format(i)+'.png'
        fig,ax=plt.subplots()
    
        if rindex17[i] != -1: #make composite map with AIA
            #print(rindex17[i])
            comp_map17 = sunpy.map.Map(s304,submaps17[rindex17[i]], composite=True)
            comp_map17.set_alpha(0, 1)
            comp_map17.set_alpha(1, 1)
            comp_map17.set_levels(1, levels,percent=True)
            plt.sca(ax)
            #ax1.set_xticklabels(default_xticks)
            #ax1.set_yticklabels(default_yticks)
            #ax1.get_xaxis().set_visible(False)
            #ax1.get_yaxis().set_visible(False)
            rplot_settings=comp_map17.get_plot_settings(1)
            rplot_settings['cmap']=cm.Greens
            comp_map17.set_plot_settings(1, rplot_settings)
            comp_map17.plot(title='AIA 304 '+s304.meta['date-obs']+' NORH '+submaps17[rindex17[i]].meta['time-obs'])
        else:
            #s304.plot(axes=ax,title='AIA 304'+s304.meta['date-obs'])
            break
        plt.savefig(imname)
        plt.clf()
    #mov.run_ffmpeg('norhmap'+'%03d.png',moviename,framerate=framerate)
    return rindex17,maplist304,submaps17


def rhessi_on_304(file304='AIA/304/mapcube304_peak.p',moviename='rhessi_on_304.mp4',framerate=12,bsub=0,norm=[2500,25000],mask_disk=True,creverse=True,rhessi='soft'):
    #maplist17.sort()
    #maplist34.sort()
    maplist304=pickle.load(open(file304,'rb'))
    times304=[dt.strptime(m.meta['date-obs'][:-3],'%Y-%m-%dT%H:%M:%S') for m in maplist304.maps]
    #build the pre-flare image for subtraction
    bl304=maplist304[:bsub]
    bav304=np.sum([m.data for m in bl304],axis=0)/float(bsub)
   
    if rhessi == 'soft':
        rhessi_files=glob.glob('Xray/rhsi_low*.fits')
    elif rhessi == 'hard':
        rhessi_files=glob.glob('Xray/rhsi_hi*.fits')
        
    rhessi_maps=[sunpy.map.Map(f) for f in rhessi_files]
    rtimes,srhsi=[],[]
    for rhsi_map in rhessi_maps:
        rhsi_map.meta['CTYPE1'] = 'HPLN-TAN'
        rhsi_map.meta['CTYPE2'] = 'HPLT-TAN'
        rhsi_map.meta['cunit1']='arcsec'
        rhsi_map.meta['cunit2']='arcsec'
        rhsi_map.meta['date-obs']=dt.strftime(dt.strptime('0'+rhsi_map.meta['date-obs'][1:-4],'%d-%b-%Y %H:%M:%S'),'%Y-%m-%dT%H:%M:%S')
        rhsi_map.plot_settings['cmap']=cm.Greens
        rhsi_map.rotate()
        #rtimes.append(dt.strptime(rhsi_map.meta['date-obs'][:-4],'%Y-%m-%dT%H:%M:%S'))
        rtimes.append(dt.strptime(rhsi_map.meta['date-obs'],'%Y-%m-%dT%H:%M:%S'))
        srhsi.append(rhsi_map.submap(SkyCoord([-1210,-750]*u.arcsec,[-50,500]*u.arcsec,frame=rhsi_map.coordinate_frame)))

    #make a of -1 (no rhessi) or index (rhessi) corresponding to AIA times
    rindex304=[]
    current_idx=0
    next_idx=1
    for t304 in times304[bsub:]:
        if next_idx >= len(rtimes):
            rindex304.append(-1)
        elif t304 < rtimes[current_idx]:
            rindex304.append(-1)
        elif t304 < rtimes[next_idx]:
            rindex304.append(current_idx)
        elif t304 >= rtimes[next_idx]:
            current_idx+=1
            next_idx+=1
    for i in range(0,5):
        rindex304.append(-1)
        
    print(len(rindex304),rindex304)
        
    for i,m304 in zip(range(0,len(maplist304[bsub:])),maplist304[bsub:]):
        imname='rhsi_aia_hi'+'{0:03d}'.format(i)+'.png'
        fig,ax=plt.subplots()
        #subtract background
        if bsub !=0:
            new304=m304.data-bav304
            meta304=m304.meta
        m304.rotate()
        m304s=m304.submap(SkyCoord([-1210,-750]*u.arcsec,[-50,500]*u.arcsec,frame=m304.coordinate_frame))
            #print(np.min(new17),np.max(new17))

        palette_name=m304.plot_settings['cmap']#sunpy.cm.cmlist[cmap]#m17.plot_settings['cmap']
        if creverse and not palette_name.name.endswith('_r'): #true means reverse the color table
            #palette_name=mapss.plot_settings['cmap']
            new_cdata=cm.revcmap(palette_name._segmentdata)
            new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name+'_r',new_cdata)
        elif not creverse and palette_name.name.endswith('_r'):
            new_cdata=cm.revcmap(palette_name._segmentdata)
            new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name[:-2],new_cdata)

        #if mask_disk:
            #masked_map304=ts.mask_disk(m304s,fac=0.,plot=False)
            #m304s=masked_map304#.filled(0)
        m304s.plot_settings['cmap']=new_cmap
        #m304s.plot_settings['norm']=colors.Normalize(norm[0],norm[1])
        
        if rindex304[i] != -1:
            comp_map304 = sunpy.map.Map(m304s,srhsi[rindex304[i]], composite=True)
            comp_map304.set_alpha(0, 1)
            comp_map304.set_alpha(1, 1)
            comp_map304.set_levels(1, [50,51,69,70,89,90],percent=True)
            plt.sca(ax)
            #ax2.set_xticklabels(default_xticks)
            #ax2.set_yticklabels(default_yticks)
            #ax.get_xaxis().set_visible(False)
            #ax.get_yaxis().set_visible(False)
            rplot_settings=comp_map304.get_plot_settings(1)
            rplot_settings['cmap']=cm.Greens
            comp_map304.set_plot_settings(1, rplot_settings)

            comp_map304.plot(title='AIA 304 '+m304s.meta['date-obs'])
        else:
            m304s.plot(axes=ax,title='AIA 304 '+m304s.meta['date-obs'])
        #fig.show()
        plt.savefig(imname)
        plt.clf()
    #mov.run_ffmpeg('norhmap'+'%03d.png',moviename,framerate=framerate)
    return rindex304,srhsi

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

def fname2png(fname, imname,cmap=False,norm=False,submap=False,overwrite=False,edges=True,retitle=False):
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
    #pickle.dump(mapss,open(pname,'wb'))
    plt.clf()
    if retitle:
        rtitle=mapss.meta['date-obs']+ ' ' + mapss.meta['time-obs']
        if cmap and norm:
            mapss.plot(cmap=cmap,norm=norm,title=rtitle)
        elif cmap and not norm:
            mapss.plot(cmap=cmap,title=rtitle)
        elif norm and not cmap:
            mapss.plot(norm=norm,title=rtitle)       
        else:
            mapss.plot(title=rtitle)
    else:
        if cmap and norm:
            mapss.plot(cmap=cmap,norm=norm)
        elif cmap and not norm:
            mapss.plot(cmap=cmap)
        elif norm and not cmap:
            mapss.plot(norm=norm)       
        else:
            mapss.plot()
        
    plt.savefig(imname)

        
def make_timeseries_movie(ddir,moviename,framerate=12,tags=False,n=False,imtags=False,imdir=False,ext='png',submap=False,overwrite=False,edges=False,cmap=False, norm=False,from_fits=False,retitle=True):
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
    if from_fits:
        all_fits=glob.glob(tags+'*')#+'*.fits') #or .fits
        fitsnames=[af[:-5] for af in all_fits]
        all_ims=glob.glob(imnames+'*.'+ext)
        len_fits=len(all_fits)
        len_ims=len(all_ims)
        if len_fits != len_ims: #get the indexes of where the png images end? assuming only later images will be added
            overwrite=True
            remaining_fits=all_fits[len_ims:]
    else:
        all_pickles=glob.glob('*'+tags+'*.p')
        remaining_fits=all_pickles
        len_ims=len(all_pickles)
    if overwrite:
        #get the maps
        if norm:
            norm=colors.Normalize(vmin=norm[0],vmax=norm[1])
        for i,f in enumerate(remaining_fits):
            if edges:
                pname=f[:-4]+'_edges.p'
            if from_fits:
                pname=f[:-4]+'.p'                  
            else:
                pname=f                
            imname=imdir+'/'+imnames+'{0:03d}'.format(i)+'.'+ext
            if not from_fits and pname in all_pickles:
                #restore and make the image from here
                pmap=pickle.load(open(pname,'rb'))
                plt.clf()
                pmap.plot(cmap=cmap,norm=norm)
                plt.savefig(imname)
            else:
                fname2png(f,imname,submap=submap,cmap=cmap,norm=norm,overwrite=overwrite,edges=edges,retitle=retitle)
            
    os.chdir(imdir)
    #make movie. add try and except if images not found... 
    mov.run_ffmpeg(imnames+'%03d.png',moviename,framerate=framerate)
    #print 'movie ', moviename, ' is stored in ', mdir
        
def make_timeseries_movie_mapcubes(ddir,moviename,framerate=12,tags=False,n=False,imtags=False,imdir=False,ext='png',submap=False,overwrite=False,edges=False):
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
    all_maps=glob.glob('mapcube'+tags+'*.p')
    all_maps.sort()
    #fitsnames=[af[:-5] for af in all_fits]
    #all_ims=glob.glob(imnames+'*.'+ext)
    #all_pickles=glob.glob('*'+tags+'*.p')
    len_cubes=len(all_maps)
    #len_ims=len(all_ims)
    idx=0
    #get the maps
    for mc in all_maps:
        #if edges:
        #    pname=f[:-4]+'_edges.p'
        #else:
        #    pname=f[:-4]+'.p'
        mapcube= pickle.load(open(mc,'rb'))
        for m in mapcube:
            imname=imdir+'/'+imnames+'{0:03d}'.format(idx)+'.'+ext
            plt.clf()
            m.plot()
            plt.savefig(imname)
            idx+=1
            #else:
            #    fname2png(f,imname,submap=submap,overwrite=overwrite,edges=edges)
            
    os.chdir(imdir)
    #make movie. add try and except if images not found... 
    mov.run_ffmpeg(imnames+'%03d.png',moviename,framerate=framerate)
    #print 'movie ', moviename, ' is stored in ', mdir

def make_joint_movies_mapcubes(ddir,moviename,framerate=12,tags=False,n=-1,imtags=False,imdir=False,ext='png',submap=False,overwrite=False,creverse=False,edges=False,show=False,extrema2=False,normalize=True):
    '''Make joint movies from pickles or mapcubes and time syncronize them as best as possible.'''
    #check if any of the map images already exist
    cwd=os.getcwd()
    all_maps=[]
    if imtags:
        imnames='map_'+imtags
    else:
        imnames='map_'
    if not tags:
        tags=['*' for d in ddir]
    if not imdir:
        imdir='.'
    for dd,tt in zip(ddir,tags):
        #os.chdir(ddir)
        lmaps=glob.glob(dd+'/mapcube'+tt+'*.p')
        if lmaps == []:
            lmaps=glob.glob(dd+'/*'+tt+'.p')
        else:
            lmaps.append(dd+'/AIA_lev1-5_names.p')
            
        lmaps.sort()
        all_maps.append(lmaps)
        #os.chdir('../')

    idx=0
    #get the maps
    synced_maps,file_groups,extrema=sync_maps(all_maps)
    movie_from_synced_maps(synced_maps,moviename,framerate=framerate,n=n,imdir=imdir,ext=ext,submap=submap,overwrite=overwrite,creverse=creverse,edges=edges,show=show,extrema=extrema2,normalize=normalize)
    return synced_maps,file_groups,extrema

def movie_from_synced_maps(synced_maps,moviename,framerate=12,imdir='.',n=-1,ext='png',submap=False,overwrite=False,creverse=False,edges=False,show=False,extrema=False,normalize=True):
    #synced_maps,file_groups=sync_maps(all_maps)
    idx=0
    def my_func(x, pos):
        return str(x/1000)

    fmt1 = tk.FuncFormatter(my_func)

    for sm in synced_maps[:n]:
        imname=imdir+'/map_'+'{0:03d}'.format(idx)+'.'+ext
        if submap:
            pms=[process_map(m,submap=submap[i],creverse=creverse,edges=edges,pickle=False,normalize=normalize) for i,m in enumerate(sm)]
        else:
            pms=[process_map(m,creverse=creverse,edges=edges,pickle=False,normalize=normalize) for i,m in enumerate(sm)]

        pms.reverse()#for when the AIA is on top ... put it at the bottom
        pms[0],pms[1]=pms[1],pms[0]
        #plt.clf()
        fig=plt.figure(figsize=(7,15))
        for a,pm in enumerate(pms):
            if extrema:
                norm=colors.Normalize(vmin=extrema[a][0],vmax=extrema[a][1])
            else:
                norm=colors.Normalize(vmin=np.min(pm.data),vmax=np.max(pm.data))
                
            ax=fig.add_subplot(len(pms),1,a+1)
            cf=pm.plot(norm=norm) #later add option for colormap normalization
            ax.xaxis.label.set_visible(False)
            ax.yaxis.label.set_visible(False)
            #cbar_ax = fig.add_axes([0.85, .7-(float(a)/3.), 0.02, 0.25])
            #cbar=fig.colorbar(cf, cax=cbar_ax,format=fmt1)
            #cbar.set_label('x 1000')
            if a == 0:
                box1=ax.get_position()
                #print(box1)
                #ax.set_position([box1.x0+.05,box1.y0,box1.width,box1.height])
                cbar_ax = fig.add_axes([0.85, .675, 0.03, 0.2])
                cbar=fig.colorbar(cf, cax=cbar_ax,format=fmt1)
                cbar.set_label('x1000')
                ax.set_xlim([250,750])
                ax.set_ylim([0,400])
                #fig.subplots_adjust(top=.9)
            elif a == 1:
                #box1=ax.get_position()
                ax.set_position([box1.x0,box1.y0-.3,box1.width,box1.height])
                cbar_ax = fig.add_axes([0.85, box1.y0-.3, 0.03, 0.23])
                cbar=fig.colorbar(cf, cax=cbar_ax,format=fmt1)
                cbar.set_label('x1000')
                ax.set_xlim([750,1250])
                ax.set_ylim([140,540])
                #fig.subplots_adjust(left=.005, top=.9)
            else:
                #box1=ax.get_position()
                ax.set_position([box1.x0,box1.y0-.6,1.01*box1.width,1.01*box1.height])
                cbar_ax = fig.add_axes([0.85, .07, 0.03, 0.23])
                cbar=fig.colorbar(cf, cax=cbar_ax,format=fmt1)
                cbar.set_label('x1000')
                ax.set_xlim([-1223,-723])
                ax.set_ylim([40,440])
            #fig.subplots_adjust(bottom=.7-(float(a)/3.))
        #fig.tight_layout()
        #fig.subplots_adjust(hspace=.1)
        if show:
            fig.show()
        plt.savefig(imname)
        idx+=1
    
    #os.chdir(imdir)
    #make movie. add try and except if images not found... 
    mov.run_ffmpeg('map_'+'%03d.png',moviename,framerate=framerate)
    #print 'movie ', moviename, ' is stored in ', mdir
    #return synced_maps,file_groups

def sync_maps(all_maps,norh=False):
    '''time syncronize the maps. assume AIA is first... since it is higher cadence'''
    synced_maps=[]
    nmaps=np.shape(all_maps)[0]
    namearr=[]
    for am in all_maps:
        if am[0].endswith('AIA_lev1-5_names.p'): #it's a mapcube
            all_aia=pickle.load(open(am[0],'rb'))
            wavestr=am[1][-6:-3]
            if wavestr.startswith('e'):
                wave='AIA_094'
            else:
                wave='AIA_'+wavestr
            namearr.append([aia for aia in all_aia if wave in aia])
        else:
            namearr.append(am) #list of pickle files

    #now parse the names like in aia_dem_batch
    basenvec=namearr[0] #vector to match the names to.
    basetimevec=time_vecs(basenvec,norh=norh)
    other_tvec=[time_vecs(nv,norh=norh) for nv in namearr[1:]]

    file_groups=[]
    for bt in basetimevec:
        closest_file=[basenvec[basetimevec.index(bt)]]
        for j,tvec in enumerate(other_tvec):
            closest_t= min(tvec, key=lambda x: abs(x - bt))
            ct_idx=tvec.index(closest_t)
            closest_file.append(namearr[j+1][ct_idx])
        file_groups.append(closest_file)

    #get the actual maps now
    mapcube=pickle.load(open(all_maps[0][1],'rb')) #for now hard code that AIA map is last
    ncube=1
    extrema=[[50,51],[50,51],[50,51]]
    for i, cf in enumerate(file_groups):
        gmaps=[]
        mapdatelist=[m.name[-19:] for m in mapcube]
 
        #now restore the maps
        for ci,c in enumerate(cf):
            if c.endswith('.p'): #make sure the ddir info is carried over with c
                newmap=pickle.load(open(c,'rb'))
                    
            else:
                mdate=dt.strftime(dt.strptime(c[8:-5],'%Y%m%d_%H%M%S'),'%Y-%m-%d %H:%M:%S')
                test=True
                while test:
                    try:
                        aia_idx= mapdatelist.index(mdate)
                        test=False
                    except ValueError: #load the next mapcube
                        ncube+=1
                        mapcube=pickle.load(open(all_maps[0][ncube],'rb'))
                        mapdatelist=[m.name[-19:] for m in mapcube]
                newmap=mapcube[aia_idx]
            gmaps.append(newmap)
            if np.max(newmap.data) > extrema[ci][1]:
                extrema[ci][1]=np.max(newmap.data)
            if np.min(newmap.data) < extrema[ci][0]:
                extrema[ci][0]=np.min(newmap.data)
        synced_maps.append(gmaps)
        
    return synced_maps, file_groups,extrema

def time_vecs(namevec,norh=False):
    timevec=[]
    for n in namevec:
        if n.endswith('.p'):
            timevec.append(dt.strptime(n[-21:-8],'%y%m%d_%H%M%S'))
        elif norh:
            timevec.append(dt.strptime(n[3:],'%y%m%d_%H%M%S'))
        else:
            timevec.append(dt.strptime(n[8:-5],'%Y%m%d_%H%M%S'))
    return timevec   
    
def process_map(smap,submap=False,creverse=False,edges=False,pickle=False,normalize=True):
    '''do the individual map processing to make life easier'''
    if submap:
        bl=submap[0]
        tr=submap[1]
        mapss=get_submap(smap,bl,tr) 
    else:
        mapss=smap

    if normalize:
        if mapss.meta['exptime'] > 8.0: #this really should only be the case for STEREO
            newdata=mapss.data/2.
            newmap=sunpy.map.Map(newdata,mapss.meta)
            del(mapss)
            mapss=newmap

    if edges:
        mapss=edge_filter(mapss)
        #pname=fname[:-4]+'_edges.p'
    #else:
        #pname=fname[:-4]+'.p'

    palette_name=mapss.plot_settings['cmap']
    if creverse and not palette_name.name.endswith('_r'): #true means reverse the color table
        #palette_name=mapss.plot_settings['cmap']
        new_cdata=cm.revcmap(palette_name._segmentdata)
        new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name+'_r',new_cdata)
        mapss.plot_settings['cmap']=new_cmap
    elif not creverse and palette_name.name.endswith('_r'):
        new_cdata=cm.revcmap(palette_name._segmentdata)
        new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name[:-2],new_cdata)
        mapss.plot_settings['cmap']=new_cmap
        
        
    if pickle:
        pickle.dump(mapss,open(pname,'wb'))
    return mapss

def map_resize(filelist,bl,tr):
    for f in filelist:
        newname=f[:-4]+'.p'
        fmap=sunpy.map.Map(f)
        fmaps=get_submap(fmap, bl,tr)
        pickle.dump(fmaps,open(newname,'wb'))
    
