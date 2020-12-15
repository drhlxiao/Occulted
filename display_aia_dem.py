 #######################################
#display_aia_dem.py
# Erica Lastufka 15/03/2018  

#Description: Because OSX doesn't play well with XQuartz and IDL sucks
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
from scipy import ndimage as ndi
import astropy.units as u
from astropy.coordinates import SkyCoord
import make_aligned_movie as mov

#global emcube
#global lgtaxis
#should store all this stuff in common block equivalents...how to deal with the map meta?

def dem_from_sav(filename):
    dem_dict=readsav(filename,python_dict=True)
    return dem_dict

def plot_whole_cube(dem_dict=True,picklename=False,meta=False,lgtaxis=False,extent=False,ccmap=False):
    '''what it sounds like. update so that there's only one colobar extending the entire y axis...'''
    try:
        lgtaxis=dem_dict['lgtaxis']
    except KeyError:
        pass
    #if not lgtaxis:
    #    lgt=readsav('lgtaxis.sav')
    #    #global lgtaxis
    #    lgtaxis=lgt['lgtaxis']
    xylabels=['Log T='+str(np.round(lg,2)) for lg in lgtaxis]
    if not picklename:
        emcube=dem_dict['emcube']
        mapcube,mins,maxs=[],[],[]
        if not meta:
            try:
                meta=pickle.load(open('map_meta.p','rb'))
            except IOError:
                #get it from a fitsfile itself
                fitsf=glob.glob('AIA_094*.fits')
                metamap=sunpy.map.Map(fitsf[0])
                meta=metamap.meta
        #convert to sunpy maps
        for n in range(0,21):
            data=emcube[n,:,:]
            data=np.ma.masked_less(data,0.00001)
            mmap=sunpy.map.Map(data,meta)
            if type(extent) == list:
                mmap=mmap.submap(SkyCoord(extent[0]*u.arcsec,extent[2]*u.arcsec,frame=mmap.coordinate_frame),SkyCoord(extent[1]*u.arcsec,extent[3]*u.arcsec,frame=mmap.coordinate_frame))
            mapcube.append(mmap)
            mins.append(mmap.min())
            maxs.append(mmap.max())
    else:
        mapcube=pickle.load(open(picklename,'rb'))
        mins,maxs=[],[]
        for mmap in mapcube:
            mins.append(mmap.min())
            maxs.append(mmap.max())
        meta=mapcube[0].meta
    if not ccmap:
        ccmap=cm.rainbow

    norm=colors.Normalize(vmin=np.min(mins),vmax=np.max(maxs))
        
    fig=plt.figure(figsize=(6,10))
    for i,im in enumerate(mapcube):
        #ii=i%3
        #jj=i%2
        ax=fig.add_subplot(7,3,i+1)
        cf=im.plot(cmap=ccmap,axes=ax,norm=norm)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.title.set_visible(False)
        #cf=ax[ii][jj].imshow(immasked,cmap=cm.rainbow,norm=norm,origin='lower',extent=extent)
        ax.annotate(xylabels[i],xy=(im.bottom_left_coord.Tx.value+5,im.bottom_left_coord.Ty.value+5),xytext=(im.bottom_left_coord.Tx.value+5,im.bottom_left_coord.Ty.value+5))#,textcoords='axes points')
        #print i,i%3
        #if i%3==2:
        #

    fig.tight_layout()
    fig.subplots_adjust(top=0.95)
    fig.suptitle('AIA DEM Analysis Results '+meta['date_obs'])
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.875, 0.15, 0.05, 0.7])
    fig.colorbar(cf, cax=cbar_ax)
    fig.show()
    return mapcube

def get_max_T_EM(vec,lgtaxis):
    maxem=np.max(vec)
    maxx=vec.tolist().index(maxem)
    maxt=lgtaxis[maxx]
    return maxt,maxem
    
def plot_T_EM_simple(dem_dict,tm=False,emm=False,test=False,meta=False,lgtaxis=False,extent=False,tmin=False,EMrange=False,Trange=False,EMonly=False,plotname=False,gauss=False,plot=True):
    '''Make maps showing the max T and EM for each pixel in emcube. No bells and whistles'''
    #first test for a single pixel
    #global emcube
    emcube=dem_dict['emcube'] # [ndem, y,x] #need it global for plot_DEM_curve? not sure
    status=dem_dict['status'] #[y,x]
    try:
        lgtaxis=dem_dict['lgtaxis']
    except KeyError:
        pass
    if not meta:
        try:
            meta=pickle.load(open('map_meta.p','rb'))
        except IOError:
            pass
    if type(lgtaxis) != np.ndarray:
        lgt=readsav('lgtaxis.sav')
        #global lgtaxis
        lgtaxis=lgt['lgtaxis']
    if not tm and not emm:    
        #make whole map
        maxTmap=np.zeros(np.shape(np.transpose(status)))
        maxEMmap=np.zeros(np.shape(np.transpose(status)))
        emx=range(0,np.shape(status)[1])
        emy=range(0,np.shape(status)[0])
        for coordpair in np.transpose([np.tile(emx,len(emy)),np.repeat(emy,len(emx))]):#zip(range(np.shape(status)[1]),range(np.shape(status)[0])):
            mx=coordpair[0]
            my=coordpair[1]
            vec=emcube[:,my,mx]
            if np.max(vec)!=0.0:
                maxt,maxem=get_max_T_EM(vec,lgtaxis)
                maxTmap[mx,my]=maxt
                maxEMmap[mx,my]=maxem
        #fig,ax=plt.subplots(1,2)
        #do some masking
        plotTmap=np.ma.masked_less(np.transpose(maxTmap),lgtaxis[0])
        plotEMmap=np.ma.masked_less(np.transpose(maxEMmap),0.001)
        if not tmin:
            norm1=colors.Normalize(vmin=np.min(plotTmap),vmax=np.max(plotTmap))
        else:
            plotTmap.data=np.ma.masked_less(plotTmap.data,tmin)
            norm1=colors.Normalize(vmin=tmin,vmax=np.max(plotTmap))        
        norm2=colors.Normalize(vmin=np.min(plotEMmap),vmax=np.max(plotEMmap))

        if not meta:
            try:
                meta=pickle.load(open('map_meta.p','rb'))
            except IOError:
                #get it from a fitsfile itself
                fitsf=glob.glob('AIA_094*.fits')
                metamap=sunpy.map.Map(fitsf[0])
                meta=metamap.meta

        #convert into SunPy maps and return
        sTmap=sunpy.map.Map(plotTmap,meta)
        if gauss:
            plotEMmap_smoothed=ndi.gaussian_filter(plotEMmap,gauss)
            #re-mask
            plotEMmap_smoothed=np.ma.masked_less(plotEMmap_smoothed,EMrange[0])
            sEMmap=sunpy.map.Map(plotEMmap_smoothed,meta)
            plotTmap_smoothed=ndi.gaussian_filter(plotTmap,gauss)
            #re-mask
            plotTmap_smoothed=np.ma.masked_less(plotTmap_smoothed,Trange[0])
            sTmap=sunpy.map.Map(plotTmap_smoothed,meta)
            sEMmap=sunpy.map.Map(plotEMmap_smoothed,meta)
        else:
            sEMmap=sunpy.map.Map(plotEMmap,meta)
        if extent:
            #print(sEMmap.bottom_left_coord)
            #print('here')
            sTmap=sTmap.submap(SkyCoord(extent[0]*u.arcsec,extent[2]*u.arcsec),SkyCoord(extent[1]*u.arcsec,extent[3]*u.arcsec,frame=sEMmap.coordinate_frame))
            sEMmap=sEMmap.submap(SkyCoord(extent[0]*u.arcsec,extent[2]*u.arcsec),SkyCoord(extent[1]*u.arcsec,extent[3]*u.arcsec,frame=sEMmap.coordinate_frame))
            #print(sEMmap.bottom_left_coord)
            #if extent:
    else: #sTmap and sEM map are inputs
        sTmap=tm
        sEMmap=emm
    norm1=colors.Normalize(vmin=np.min(lgtaxis[0]),vmax=np.max(sTmap.data))
    if EMrange:
        #mask below
        mdata=np.ma.masked_less(sEMmap.data,EMrange[0])
        sEMmap=sunpy.map.Map(mdata,sEMmap.meta)
        norm2=colors.Normalize(vmin=EMrange[0],vmax=EMrange[1])
    else:
        norm2=colors.Normalize(vmin=np.min(sEMmap.data),vmax=np.max(sEMmap.data))
        
    if EMonly:
        fig,ax=plt.subplots(figsize=(6,8))
        #print(sEMmap.bottom_left_coord)
        cf=sEMmap.plot(cmap=cm.rainbow,norm=norm2,axes=ax)
        ax.xaxis.label.set_visible(False)
        ax.yaxis.label.set_visible(False)
        ax.set_title('Max Emission Measure (Log $cm^{-5}$) '+meta['date-obs'][11:-3])
        fig.colorbar(cf,ax=ax,fraction=.03,pad=.07)      

    elif plot:
        fig,ax=plt.subplots(1,2,sharex=True,sharey=True,figsize=(12,6))
        cf1=sTmap.plot(cmap=cm.rainbow,norm=norm1,axes=ax[0])#,extent=extent)
        cf2=sEMmap.plot(cmap=cm.rainbow,norm=norm2,axes=ax[1])#,extent=extent)
        for a in ax:
            a.xaxis.label.set_visible(False)
            a.yaxis.label.set_visible(False)
        #ax.title.set_visible(False)
        ax[0].set_title('Max log T')
        ax[1].set_title('Max Emission Measure (Log $cm^{-5}$)')
        fig.colorbar(cf1,ax=ax[0],fraction=.03,pad=.07)
        fig.colorbar(cf2,ax=ax[1],fraction=.03,pad=.07)
    fig.tight_layout()
    if plotname:
        plt.savefig(plotname)
    elif plot:
        fig.show()
    return sTmap,sEMmap

def plot_T_EM(dem_dict,tm=False,emm=False,test=False,meta=False,lgtaxis=False,extent=False):
    '''Make maps showing the max T and EM for each pixel in emcube'''
    #first test for a single pixel
    #global emcube
    emcube=dem_dict['emcube'] # [ndem, y,x] #need it global for plot_DEM_curve? not sure
    status=dem_dict['status'] #[y,x]
    try:
        lgtaxis=dem_dict['lgtaxis']
    except KeyError:
        pass
    if not meta:
        try:
            meta=pickle.load(open('map_meta.p','rb'))
        except IOError:
            pass
    if not lgtaxis:
        lgt=readsav('lgtaxis.sav')
        #global lgtaxis
        lgtaxis=lgt['lgtaxis']
    if test:
        #find a place where there is a solution
        for sx,sy in zip(range(0,np.shape(status)[0]),range(0,np.shape(status)[1])):
            if np.max(emcube[:,sy,sx])!=0.0: #currently something must be messed up with the status array...
                testpc=sx,sy
                break
        vec=emcube[:,sy,sx]
        fig,ax=plt.subplots()
        ax.plot(vec)
        fig.show()
        maxv=np.max(vec)
        maxx=vec.tolist().index(maxv)
        maxt=lgtaxis[maxx]
        print( sx,sy,maxv,maxx,maxt)
    if not tm and not emm:    
        #make whole map
        maxTmap=np.zeros(np.shape(np.transpose(status)))
        maxEMmap=np.zeros(np.shape(np.transpose(status)))
        emx=range(0,np.shape(status)[1])
        emy=range(0,np.shape(status)[0])
        for coordpair in np.transpose([np.tile(emx,len(emy)),np.repeat(emy,len(emx))]):#zip(range(np.shape(status)[1]),range(np.shape(status)[0])):
            mx=coordpair[0]
            my=coordpair[1]
            vec=emcube[:,my,mx]
            if np.max(vec)!=0.0:
                maxt,maxem=get_max_T_EM(vec,lgtaxis)
                maxTmap[mx,my]=maxt
                maxEMmap[mx,my]=maxem
        #fig,ax=plt.subplots(1,2)
        #do some masking
        plotTmap=np.ma.masked_less(np.transpose(maxTmap),lgtaxis[0])
        plotEMmap=np.ma.masked_less(np.transpose(maxEMmap),0.001)    
        norm1=colors.Normalize(vmin=np.min(plotTmap),vmax=np.max(plotTmap))
        norm2=colors.Normalize(vmin=np.min(plotEMmap),vmax=np.max(plotEMmap))

        if not meta:
            try:
                meta=pickle.load(open('map_meta.p','rb'))
            except IOError:
                #get it from a fitsfile itself
                fitsf=glob.glob('AIA_094*.fits')
                metamap=sunpy.map.Map(fitsf[0])
                meta=metamap.meta
        #convert into SunPy maps and return
        sTmap=sunpy.map.Map(plotTmap,meta)
        sEMmap=sunpy.map.Map(plotEMmap,meta)
        if extent:
            sTmap=sTmap.submap(SkyCoord(extent[0]*u.arcsec,extent[1]*u.arcsec),SkyCoord(extent[2]*u.arcsec,extent[3]*u.arcsec))
            sEMmap=sEMmap.submap(SkyCoord(extent[0]*u.arcsec,extent[1]*u.arcsec),SkyCoord(extent[2]*u.arcsec,extent[3]*u.arcsec))
    else: #sTmap and sEM map are inputs
        sTmap=tm
        sEMmap=emm
        norm1=colors.Normalize(vmin=np.min(lgtaxis[0]),vmax=np.max(sTmap.data))
        norm2=colors.Normalize(vmin=np.min(sEMmap.data),vmax=np.max(sEMmap.data))

    def plot_dem_curve(startx,starty,ymax):
        ax_list = fig.axes
        #print ax_list
        ax=ax_list[2]
        #print ax
        #ax.cla()
        xy = plt.ginput(1)
        #print xy, type(xy[0][0])
        if type(xy[0][0]) != int:
            x = xy[0][0] #have to convert back to image coords that are not in arcsec...
            y = xy[0][1]
            #have to convert back to image coords that are not in arcsec...
            xpix,ypix=sTmap.world_to_pixel(SkyCoord(x*u.arcsec,y*u.arcsec,frame=sTmap.coordinate_frame))
            xpix=int(xpix.value)
            ypix=int(ypix.value)
        else:
            xpix=startx
            ypix=starty
        vec=emcube[:,ypix,xpix]
        #fig,ax=plt.subplots()
        ax.plot(lgtaxis,vec,'b.-')
        ax.set_ylim([0,ymax])
        fig.show()
        #fig.canvas.draw()
        maxv=np.max(vec)
        maxx=vec.tolist().index(maxv)
        maxt=lgtaxis[maxx]
        print( xpix,ypix,maxv,maxt,vec)

    def onclick(event):
        if event.dblclick:
            #if event.button == 1:
            # Draw line
            plot_dem_curve(event.xdata,event.ydata,ymax) # here you click on the plot
            #fig.canvas.draw()
        else: #if event.button == 3:
            ax_list = fig.axes
            #print ax_list
            ax=ax_list[2]
            ax.cla()
            
            #    fig.canvas.mpl_disconnect(cid)
            #    #pass # Do nothing 
        
    fig,ax=plt.subplots(figsize=(25,6))
    ax0=plt.subplot(131)
    ax1=plt.subplot(132,sharex=ax0,sharey=ax0)
    ax2=plt.subplot(133)
    cf1=sTmap.plot(cmap=cm.rainbow,norm=norm1,axes=ax0)
    cf2=sEMmap.plot(cmap=cm.rainbow,norm=norm2,axes=ax1)
    ax0.xaxis.label.set_visible(False)
    ax0.yaxis.label.set_visible(False)
    ax1.xaxis.label.set_visible(False)
    ax1.yaxis.label.set_visible(False)
    ax2.set_xlabel('Log T')
    ax2.set_ylabel('EM ($cm^{-5}$)')
    #ax.title.set_visible(False)
    ax0.set_title('Max log T')
    ax1.set_title('Max Emission Measure (Log $cm^{-5}$)')
    ymax=np.max(sEMmap.data)
    #plot_dem_curve(415,100)
    fig.colorbar(cf1,ax=ax0,fraction=.03,pad=.07)
    fig.colorbar(cf2,ax=ax1,fraction=.03,pad=.07)
    fig.tight_layout()
    #fig.show()
    #interactive plotting of DEM curve for selected pixel
    plt.axes(ax2)
    #print ax[2]
    ax_list= fig.axes
    cid=fig.canvas.mpl_connect('button_press_event', onclick)
    #cid = fig.canvas.mpl_connect('key_press_event', on_key)
    fig.show()
    #fig.canvas.draw()
    #fig.canvas.mpl_disconnect(cid)
    return sTmap,sEMmap
    
def plot6(dem_dict, extent=False,maxval=False,minval=False,logscale=False,lgtaxis=False,meta=False):
    imarr=dem_dict['image']
    try:
        lgtaxis=dem_dict['lgtaxis']
    except KeyError:
        pass
    if not meta:
        try:
            meta=pickle.load(open('map_meta.p','rb'))
        except IOError:
            pass
    if type(lgtaxis) != np.ndarray:
        lgt=readsav('lgtaxis.sav')
        lgtaxis=lgt['lgtaxis']
    xylabels=['EM in lgt ['+str(np.round(bb,3))+','+str(np.round(lgtaxis[i+3],3))+']' for i,bb in enumerate(lgtaxis[:-3]) if i % 3 == 0]
    #print xylabels

    if minval:
        vmin=minval
    else:
        vmin=26
    if maxval:
        vmax=maxval
    else:
        vmax=29.5

    if not logscale:
        norm=colors.Normalize(vmin=vmin,vmax=vmax)
        cblabel='EM [$cm^{-5}$]'
    else:
        norm=colors.Normalize(vmin=np.log10(vmin),vmax=np.log10(vmax))
        cblabel='Log EM [$cm^{-5}$]'    
    fig=plt.figure()
    for i,im in enumerate(imarr):
        #ii=i%3
        #jj=i%2
        ax=fig.add_subplot(3,2,i+1)
        immasked=ma.masked_less(im,.1)
        if minval:
            immasked=ma.masked_less(immasked,minval)
        if maxval:
            immasked=ma.masked_greater(immasked,maxval)
        if logscale:
            immasked=np.log10(immasked)

        if not meta:
            try:
                meta=pickle.load(open('map_meta.p','rb'))
            except IOError:
                #get it from a fitsfile itself
                fitsf=glob.glob('AIA_094*.fits')
                metamap=sunpy.map.Map(fitsf[0])
                meta=metamap.meta
        #convert to map
        smap=sunpy.map.Map(immasked,meta)
        if extent:
            print( extent)
            smap=smap.submap(SkyCoord(extent[0]*u.arcsec,extent[1]*u.arcsec),SkyCoord(extent[2]*u.arcsec,extent[3]*u.arcsec))
        cf=smap.plot(cmap=cm.rainbow,norm=norm)
        ax.xaxis.label.set_visible(False)
        ax.yaxis.label.set_visible(False)
        ax.title.set_visible(False)
        #cf=ax[ii][jj].imshow(immasked,cmap=cm.rainbow,norm=norm,origin='lower',extent=extent)
        ax.annotate(xylabels[i],xy=(smap.bottom_left_coord.Tx.value+5,smap.bottom_left_coord.Ty.value+5),xytext=(smap.bottom_left_coord.Tx.value+5,smap.bottom_left_coord.Ty.value+5))#,textcoords='axes points')

        if i%2==1:
            fig.colorbar(cf,ax=ax,label=cblabel)

    #fig.colorbar(immasked,ax=ax[0][0])
    fig.tight_layout()
    fig.show()

def maxEM_movie(savfiles,moviename,framerate=False,EMrange=False,extent=False,gauss=4):
    '''make a movie showing where the max EM is'''
    for i,s in enumerate(savfiles):
        dem_dict=dem_from_sav(s)
        #get map metadata from fits
        rootname=s[4:-4]
        fitsname=glob.glob('AIA'+rootname+'*.fits')[0]
        fmap=sunpy.map.Map(fitsname)
        if gauss:
            names='maxEM_gauss'+str(gauss)+'_'
            plotname=names+'{0:03d}'.format(i)+'.png'
        else:
            names='maxEM'
            plotname=names+'{0:03d}'.format(i)+'.png'
        plt.clf()
        foo=plot_T_EM_simple(dem_dict,meta=fmap.meta,lgtaxis=dem_dict['lgtaxis'],extent=extent,EMrange=EMrange,plotname=plotname,EMonly=True,gauss=gauss)
    #run ffmpeg
    mov.run_ffmpeg(names+'%03d.png',moviename,framerate=framerate)

def maxEM_Fe18_movie(savfiles,sidx,eidx,moviename,framerate=12,EMrange=False,extent=False,gauss=4):
    '''make a movie showing where the max EM is'''
    import FeXVIII as Fe
    preppedfiles=glob.glob('AIA_*.fits')    
    groups=aia.group6(preppedfiles)

    for i,(s,g) in enumerate(zip(savfiles,groups[sidx:eidx])):
        map94=sunpy.map.Map(g[0])
        map171=sunpy.map.Map(g[2])
        map211=sunpy.map.Map(g[4])
        map18=Fe.get_Fe18(map94,map171,map211,submap=[(extent[0],extent[1]),(extent[2],extent[3])]) #save to listl/mapcube? do something
        
        dem_dict=dem_from_sav(s)
        #get map metadata from fits
        rootname=s[:-4]
        fitsname=glob.glob('AIA_*'+rootname+'.fits')[0]
        fmap=sunpy.map.Map(fitsname)
        if gauss:
            names='maxT_gauss'+str(gauss)+'_'
            plotname=names+'{0:03d}'.format(i)+'.png'
        else:
            names='maxT'
            plotname=names+'{0:03d}'.format(i)+'.png'
        plt.clf()
        Tplot,EMplot=plot_T_EM_simple(dem_dict,meta=fmap.meta,lgtaxis=dem_dict['lgtaxis'],extent=extent,EMrange=EMrange,plotname=False,EMonly=True,gauss=gauss,plot=False)
        #plot stuff
        fig,[ax1,ax2]=plt.subplots(ncols=2,figsize=(10,6))
        #ax1=fig.add_subplot(2,1,1)
        #map18.plot(axes=ax1)
        ax1.xaxis.label.set_visible(False)
        ax1.yaxis.label.set_visible(False)
        ax1.title.set_visible(False)
        #ax1.set(aspect=.55)
        map18.plot(axes=ax1)
        norm2=colors.Normalize(vmin=EMrange[0],vmax=EMrange[1])
        #ax2=fig.add_subplot(2,1,2)
        cf=Tplot.plot(cmap=cm.rainbow,norm=norm2,axes=ax2)
        ax2.xaxis.label.set_visible(False)
        ax2.yaxis.label.set_visible(False)
        ax2.title.set_visible(False)
        #ax2.set_ylim([0,400])
        ax2.set(aspect=1.15)
        #ax2.set_title('Max Emission Measure (Log $cm^{-5}$) '+fmap.meta['date-obs'][11:-3])
        plt.subplots_adjust(bottom=0.1, right=0.9, top=0.9)
        #cax = plt.axes([0.85, 0.1, 0.075, 0.8])
        fig.colorbar(cf,ax=ax2,fraction=.05,pad=.09)
        #fig.tight_layout()
        plt.suptitle(fmap.meta['date-obs'][11:-3])
        #fig.subplots_adjust(wspace=0.4)
        plt.savefig(names+'{0:03d}'.format(i)+'.png')
        #save
        
    #run ffmpeg
    mov.run_ffmpeg(names+'%03d.png',moviename,framerate=framerate)

    def maxEM_maxT_movie(savfiles,moviename,framerate=12,EMrange=False,Trange=False,extent=False,gauss=4):
        '''make a movie showing where the max EM is'''
        for i,s in enumerate(savfiles):
            dem_dict=dem_from_sav(s)
            #get map metadata from fits
            rootname=s[:-4]
            fitsname=glob.glob('AIA_*'+rootname+'.fits')[0]
            fmap=sunpy.map.Map(fitsname)
            if gauss:
                names='maxEM_maxT_gauss'+str(gauss)+'_'
                plotname=names+'{0:03d}'.format(i)+'.png'
            else:
                names='maxEM_maxT_'
                plotname=names+'{0:03d}'.format(i)+'.png'
            plt.clf()
            Tplot,EMplot=plot_T_EM_simple(dem_dict,meta=fmap.meta,lgtaxis=dem_dict['lgtaxis'],extent=extent,EMrange=EMrange,Trange=Trange,plotname=False,EMonly=True,gauss=gauss,plot=False)
            #plot stuff
            fig,[ax1,ax2]=plt.subplots(ncols=2,figsize=(10,6))
            #ax1=fig.add_subplot(2,1,1)
            #map18.plot(axes=ax1)
            norm1=colors.Normalize(vmin=Trange[0],vmax=Trange[1])
            ax1.xaxis.label.set_visible(False)
            ax1.yaxis.label.set_visible(False)
            ax1.title.set_visible(False)
            #ax1.set(aspect=.55)
            cf1=Tplot.plot(cmap=cm.rainbow,norm=norm1,axes=ax1)
            norm2=colors.Normalize(vmin=EMrange[0],vmax=EMrange[1])
            #ax2=fig.add_subplot(2,1,2)
            cf2=EMplot.plot(cmap=cm.rainbow,norm=norm2,axes=ax2)
            ax2.xaxis.label.set_visible(False)
            ax2.yaxis.label.set_visible(False)
            ax2.title.set_visible(False)
            #ax2.set_ylim([0,400])
            #ax2.set(aspect=1.15)
            #ax2.set_title('Max Emission Measure (Log $cm^{-5}$) '+fmap.meta['date-obs'][11:-3])
            plt.subplots_adjust(bottom=0.1, right=0.9, top=0.9)
            #cax = plt.axes([0.85, 0.1, 0.075, 0.8])
            fig.colorbar(cf2,ax=ax2,fraction=.05,pad=.09)
            fig.colorbar(cf1,ax=ax1,fraction=.05,pad=.09)        #fig.tight_layout()
            plt.suptitle(fmap.meta['date-obs'][11:-3])
            #fig.subplots_adjust(wspace=0.4)
            plt.savefig(names+'{0:03d}'.format(i)+'.png')
            #save
            
        #run ffmpeg
        mov.run_ffmpeg(names+'%03d.png',moviename,framerate=framerate)
