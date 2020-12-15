import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord
#import wcsaxes
from astropy.wcs import WCS

import sunpy.map
import sunpy.coordinates
import sunpy.coordinates.wcs_utils
from sunpy.net import vso
import numpy as np
import numpy.ma as ma
import matplotlib.dates as mdates
import pandas as pd

from datetime import datetime as dt
import glob
import plotly.graph_objects as go
import matplotlib
from matplotlib import cm
import pidly
from sunpy.physics.differential_rotation import solar_rotate_coordinate, diffrot_map

def get_limbcoords(aia_map):
    r = aia_map.rsun_obs - 1 * u.arcsec  # remove one arcsec so it's on disk.
    # Adjust the following range if you only want to plot on STEREO_A
    th = np.linspace(-180 * u.deg, 0 * u.deg)
    x = r * np.sin(th)
    y = r * np.cos(th)
    limbcoords = SkyCoord(x, y, frame=aia_map.coordinate_frame)
    return limbcoords

def draw_circle(sf,center,r):
    r = r* u.arcsec  # remove one arcsec so it's on disk.
    # Adjust the following range if you only want to plot on STEREO_A
    th = np.linspace(-180 * u.deg, 180 * u.deg)
    x = (r * np.sin(th))+center.Tx
    y = (r * np.cos(th))+center.Ty
    ccoords = SkyCoord(x, y, frame=sf.coordinate_frame)
    return ccoords

# from make_timeseries_plot_flux.py
def mask_disk(submap,limbcoords=False,plot=False,fac=50.,filled=False,greater=False):
    x, y = np.meshgrid(*[np.arange(v.value) for v in submap.dimensions]) * u.pixel
    if limbcoords:
        r=limbcoords.to_pixel(submap.wcs)
    else:
        hpc_coords = submap.pixel_to_world(x, y)
        r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / (submap.rsun_obs+fac*u.arcsec)
    if not greater:
        mask = ma.masked_less_equal(r, 1)
    else:
        mask = ma.masked_greater_equal(r, 1)

    palette = submap.plot_settings['cmap']
    palette.set_bad('black')
    if filled:
        data_arr=np.ma.masked_array(submap.data,mask=mask.mask)
        filled_data=ma.filled(data_arr,0)
        scaled_map = sunpy.map.Map(filled_data, submap.meta)
    else:
        scaled_map = sunpy.map.Map(submap.data, submap.meta, mask=mask.mask)
    if plot:
        fig = plt.figure()
        ax=fig.add_subplot(projection=scaled_map)
        scaled_map.plot(cmap=palette,axes=ax)
        if limbcoords:
            ax.plot_coord(limbcoords, color='w')
        else:
            scaled_map.draw_limb(axes=ax)
        fig.show()
    return scaled_map

def get_corresponding_maps(mtimes,fhead='AIA/AIA20200912_2*_0335.fits'):
    ffiles=glob.glob(fhead)
    ftimes=[dt.strptime(f[-25:-10],'%Y%m%d_%H%M%S')  for f in ffiles]
    matchfiles=[]
    for m in mtimes:
        #find closest t in tvec to given ttag
        closest_t= min(ftimes, key=lambda x: abs(x - m))
        ct_idx=ftimes.index(closest_t)
        closest_file=ffiles[ct_idx]
        matchfiles.append(closest_file)
    matchmaps=[sunpy.map.Map(aa) for aa in matchfiles]
    return matchfiles,matchmaps


def fits2mapIDL(files,coreg=True):
    '''run fits2map in idl. List of files should be: AIA094.fits, AIA171.fits,AIA211.fits'''
    idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
    idl('files',files)
    idl('fits2map,files,maps')
    if coreg:
        idl('refmap=maps[0]')
        idl('c171=coreg_map(maps[1],refmap)')
        idl('c211=coreg_map(maps[2],refmap)')
        coreg_map171= 'AIA171_coreg94_' + files[1][8:-5]#assuming AIA files are in the usual level 1.5 naming convention
        coreg_map211= 'AIA211_coreg94_' +files[2][8:-5]
        idl('coreg_map171',coreg_map171)
        idl('coreg_map211',coreg_map211)
        idl('map2fits,c171,coreg_map171')
        idl('map2fits,c211,coreg_map211')

def coreg_maps_IDL(fits1,fits2,coreg=True,coreg_mapname='coreg_map.fits'):
    '''run fits2map in idl. List of files should be: AIA094.fits, AIA171.fits,AIA211.fits'''
    idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
    idl('fits1',fits1)
    idl('fits2',fits2)
    idl('fits2map,fits1,map1')
    idl('fits2map,fits2,map2')
    if coreg:
        idl('refmap=map1')
        idl('coreg_map=coreg_map(map2,refmap)')
        #idl('help,/st,coreg_map')
        #print(coreg_mapname)
        idl('coreg_mapname',coreg_mapname)
        #idl('print,coreg_mapname')
        idl('map2fits,coreg_map,coreg_mapname')


def get_Fe18(map94,map171,map211,submap=False,save2fits=False):
    if type(map94) == str:
        map94=sunpy.map.Map(map94)
    if type(map171) == str:
        map171=sunpy.map.Map(map171) #the coregistered map
    if type(map211) == str:
        map211=sunpy.map.Map(map211) #the coregistered map
    if submap:
        map94=map94.submap(SkyCoord(submap[0]*u.arcsec,submap[1]*u.arcsec,frame=map94.coordinate_frame))
        map171=map171.submap(SkyCoord(submap[0]*u.arcsec,submap[1]*u.arcsec,frame=map171.coordinate_frame))
        map211=map211.submap(SkyCoord(submap[0]*u.arcsec,submap[1]*u.arcsec,frame=map211.coordinate_frame))
    try:
        map18data=map94.data-map211.data/120.-map171.data/450.
    except ValueError: #shapes are different- is 94 always the culprit?
        d94=map94.data[1:,:]
        map18data=d94-map211.data/120.-map171.data/450.
    map18=sunpy.map.Map(map18data,map94.meta)
    #if submap:
    #    map18=map18.submap(SkyCoord(submap[0]*u.arcsec,submap[1]*u.arcsec,frame=map18.coordinate_frame))
    if save2fits:
        filename='AIA_Fe18_'+map18.meta['date-obs']+'.fits'
        map18.save(filename)
    return map18

def int_image(filelist, nint=15,outname='AIA_Fe18_2020-09-12_'):
    for i,f in enumerate(filelist):
        st=i*nint
        nd=(i+1)*nint
        try:
            print(filelist[st],filelist[nd])
        except IndexError:
            print(nd,' out of index!')
            break #exits the loop?
        mdata=[]
        mapname=outname+str(i).zfill(2)+'.fits'
        for f in filelist[st:nd]:
            m=sunpy.map.Map(f)
            if m.meta['exptime'] != 0.0:
                mdata.append(m.data)
        intdata=np.mean(mdata,axis=0)
        map_int=sunpy.map.Map(intdata,m.meta)
        map_int.save(mapname)

#difference maps
def diff_maps(mlist):
    mdifflist=[]
    m0=mlist[0].data
    for m in mlist[1:]:
        try:
            diff=m.data-m0
        except ValueError: #size mismatch!
            m1s=m.data.shape
            m0s=m0.shape
            if m1s[0] > m0s[0]:
                mdata=m.data[1:][:]
            m.data.reshape(m0.shape)
            diff=m.data-m0
        map_diff=sunpy.map.Map(diff,m.meta)
        mdifflist.append(map_diff)
    return mdifflist

def fix_units(map_in):
    map_in.meta['cunit1']= 'arcsec'
    map_in.meta['cunit2']= 'arcsec'
    map_in.meta['date_obs']= dt.strftime(dt.strptime(map_in.meta['date_obs'],'%d-%b-%Y %H:%M:%S.%f'),'%Y-%m-%dT%H:%M:%S')#fix the time string
    map_in.meta['date-obs']= dt.strftime(dt.strptime(map_in.meta['date-obs'],'%d-%b-%Y %H:%M:%S.%f'),'%Y-%m-%dT%H:%M:%S')#fix the time string
    return map_in
    
def hand_coalign(map_in,crpix1_off,crpix2_off):
    map_in.meta['crpix1']=map_in.meta['crpix1']+ crpix1_off
    map_in.meta['crpix2']=map_in.meta['crpix2']+ crpix2_off
    return map_in

#get actual contour centroids...from find_hessi_centroids.py
def center_of_mass(X):
    # calculate center of mass of a closed polygon
    x = X[:,0]
    y = X[:,1]
    g = (x[:-1]*y[1:] - x[1:]*y[:-1])
    A = 0.5*g.sum()
    cx = ((x[:-1] + x[1:])*g).sum()
    cy = ((y[:-1] + y[1:])*g).sum()
    return 1./(6*A)*np.array([cx,cy])

def find_centroid_from_map(m,levels=[90],idx=0,show=False):
    cs,hpj_cs=[],[]
    ll=np.max(m.data)*np.array(levels)

    fig,ax=plt.subplots()
    ax.imshow(m.data,alpha=.75,cmap=m.plot_settings['cmap'])
    contour=m.draw_contours(levels=levels*u.percent,axes=ax,frame=m.coordinate_frame)
    print(len(contour.allsegs[-1]))
    c =  center_of_mass(contour.allsegs[-1][idx])
    cs.append(c)
    hpj=m.pixel_to_world(c[0]*u.pixel,c[1]*u.pixel,origin=0)
    hpj_cs.append(hpj)
    ax.plot(c[0],c[1], marker="o", markersize=12, color="red")
    if show:
        fig.show()
    return cs,hpj_cs,contour
    
def centroids_from_XRT(show=False):
    xrt8=fix_units(sunpy.map.Map('/Users/wheatley/Documents/Solar/NuStar/Xray/XRT_Be_orbit8.fits'))
    xbl1=SkyCoord(-850*u.arcsec,300*u.arcsec,frame=xrt8.coordinate_frame)
    xtr1=SkyCoord(-750*u.arcsec,400*u.arcsec,frame=xrt8.coordinate_frame)
    xrts8=xrt8.submap(xbl1,xtr1)
    _,c1f,_=find_centroid_from_map(xrts8,show=show)
    xbl2=SkyCoord(-900*u.arcsec,200*u.arcsec,frame=xrt8.coordinate_frame)
    xtr2=SkyCoord(-850*u.arcsec,300*u.arcsec,frame=xrt8.coordinate_frame)
    xrts2=xrt8.submap(xbl2,xtr2)
    _,c2f,_=find_centroid_from_map(xrts2,show=show)
    return c1f,c2f
    
def str_to_SkyCoord(in_str): #for when restoring from json etc
    clist=in_str[in_str.rfind('(')+1:-2].split(',')
    #print(clist)
    Tx=float(clist[0])*u.arcsec
    Ty=float(clist[1])*u.arcsec
    framestr=in_str[in_str.find('(')+1:in_str.find('observer=')]
    fr_name=framestr[:framestr.find(':')].lower().replace(' ','_')
    obstime=framestr[framestr.find('=')+1:framestr.find(',')]
    rsun=float(framestr[framestr.find('rsun=')+5:framestr.find(', obs')-3])*u.km
    obstr=in_str[in_str.find('observer=')+10:in_str.find(')>)')]
    observer=obstr[:obstr.find('(')]
    obscoords=obstr[obstr.rfind('(')+2:].split(',')
    #print(obscoords)
    lon=float(obscoords[0])*u.deg
    lat=float(obscoords[1])*u.deg
    radius=float(obscoords[2])*u.m
    obsloc=SkyCoord(lon,lat,radius,frame=fr_name,obstime=obstime)
    sc_out=SkyCoord(Tx,Ty,frame=obsloc.frame)
    #sc_out.observer.lat=lat
    #sc_out.observer.radius=radius

return sc_out

def flux_from_box(smap,bottom_left,top_right,mask=False, maps=False, expnorm=True):
    ssmap=smap.submap(bottom_left,top_right)
    data=ssmap.data
    flux=np.ma.sum(data)
    if type(mask) == np.ndarray:
        print('flux_from_box ', np.shape(data),np.shape(mask))
        data=np.ma.masked_array(data,mask=mask) #data.T?
        if expnorm:
            data=data/ssmap.meta['exptime']
        flux=np.ma.sum(data)
    if maps:
        newmap=sunpy.map.Map(data,ssmap.meta)
        newmap.peek()
        return flux,data,newmap
    return flux,data

def derot_box(bl,tr,start_time, mlist):
    bl_rot_all=[bl]
    tr_rot_all=[tr]
    for i,s in enumerate(mlist[1:]):
        end_time = s.meta['date-obs']
        bl_rot=solar_rotate_coordinate(bl, new_observer_time=end_time,new_observer_location=s.observer_coordinate)
        tr_rot=solar_rotate_coordinate(tr, new_observer_time=end_time,new_observer_location=s.observer_coordinate)
        bl_rot_all.append(bl_rot)
        tr_rot_all.append(tr_rot)
    return bl_rot_all,tr_rot_all

def fit_to_mask(s,bl_rot,tr_rot,mask):
    nx,ny=np.shape(mask)

    #print(bl_rot)
    #check shape
    ssmap=s.submap(bl_rot,tr_rot)
    data=ssmap.data
    meta=ssmap.meta
    sx,sy=np.shape(data)
    #print('mask shape: ',nx,ny)
    print('before: ',nx,ny, sx,sy)
    countx,county=0,0
    if sx == nx and sy == ny:
        maskt=mask
        print('no changes to mask or data')
        
    else:
    
        while sx != nx and countx < 30: #x in array changes Y in SkyCoords
            if sx > nx: #should probably do in a while
                #1 px ~ .6" for AIA
                bl_rot=SkyCoord(bl_rot.Tx,bl_rot.Ty+0.01*u.arcsec,frame=bl_rot.frame)
                tr_rot=SkyCoord(tr_rot.Tx,tr_rot.Ty-0.01*u.arcsec,frame=tr_rot.frame)
            elif sx < nx:
                bl_rot=SkyCoord(bl_rot.Tx,bl_rot.Ty-0.01*u.arcsec,frame=bl_rot.frame)
                tr_rot=SkyCoord(tr_rot.Tx,tr_rot.Ty+0.01*u.arcsec,frame=tr_rot.frame)
            newdata=s.submap(bl_rot,tr_rot).data
            sx,sy=np.shape(newdata)
            countx+=1
            #print(np.shape(s.submap(bl_rot,tr_rot).data), tr_rot.Tx - bl_rot.Tx)
        while sy !=ny and county < 30: #y in array changes X in SkyCoords
            if sy > ny:
                bl_rot=SkyCoord(bl_rot.Tx+0.01*u.arcsec,bl_rot.Ty,frame=bl_rot.frame)
                tr_rot=SkyCoord(tr_rot.Tx-0.01*u.arcsec,tr_rot.Ty,frame=tr_rot.frame)
            elif sy < ny:
                bl_rot=SkyCoord(bl_rot.Tx-0.01*u.arcsec,bl_rot.Ty,frame=bl_rot.frame)
                tr_rot=SkyCoord(tr_rot.Tx+0.01*u.arcsec,tr_rot.Ty,frame=tr_rot.frame)
            newdata=s.submap(bl_rot,tr_rot).data
            sx,sy=np.shape(newdata)
            county+=1
            #print('after: ',sx,sy,tr_rot.Ty - bl_rot.Ty)
        print('new data shape: ', sx,sy,np.shape(newdata))
        data=newdata
        maskt=mask
        #what if you need to pad both? grrr
        while sx !=nx: #should be within at least 1 row/col now, trim mask to data size
            try:
                mask[sx-1]
                mask2=mask[:sx,:]
                print('Cropping mask X')
            except IndexError: #sx too big, so pad with True
                pad=np.ones(ny,dtype=bool)
                #pad.resize(1,sy)
                mask2=np.vstack((mask,pad))
                print('Padding mask X')
                #print('sx',sx,sy,np.shape(mask2))
            nx,ny=np.shape(mask2)
            maskt=mask2
        while sy !=ny:
            try:
                mask[:,sy-1]
                mask2=maskt[:,:sy]
                #print('sy',sx,sy,np.shape(mask2))
                print('Cropping mask Y')
            except IndexError: #sx too big, so pad with True
                pad=np.ones(nx,dtype=bool)
                pad.resize(1,nx)
                mask2=np.hstack((maskt,pad.T))
                print('Padding mask Y')
                #print('sy',sx,sy,np.shape(mask2))
            nx,ny=np.shape(mask2)
            maskt=mask2
                
#        if sx == nx and sy == ny:
#            maskt=mask
#            print('current mask shape: ',np.shape(maskt))
#        else:
#            maskt=mask2
#            print('new mask shape: ', np.shape(mask2))
#            fig,ax=plt.subplots()
#            ax.imshow(maskt,origin='bottom left')
#            fig.show() #masks look right at least, so problem occurs elsewhere
    #do flux from box operations here, see if that makes things work
    data=np.ma.masked_array(data,mask=maskt) #data.T?
    flux=np.ma.sum(data)
    newmap=sunpy.map.Map(data,meta)
    #newmap.peek()
    return flux,data,newmap
    #return bl_rot,tr_rot,maskt

def track_region_box(submaps, circle, xoff=0.,yoff=0.,start_idx=0,bin=False, mask=False, plot=False, expnorm=True):
    '''add binning option, return timestamps too for greater ease'''
    map0=submaps[0]
    start_time = map0.meta['date-obs']
    timestamps=[]#[start_time]
    bl=SkyCoord(np.min(circle.Tx)+xoff*u.arcsec,np.min(circle.Ty)+yoff*u.arcsec,frame=map0.coordinate_frame) #note the +10!
    tr=SkyCoord(np.max(circle.Tx)+xoff*u.arcsec,np.max(circle.Ty)+yoff*u.arcsec,frame=map0.coordinate_frame)
    bl.set_obstime=start_time
    tr.set_obstime=start_time
    bwidth=tr.Tx-bl.Tx
    bheight=tr.Ty-bl.Ty
    #flux=[flux_from_box(map0,bl,tr,mask=mask)]
    #rotcoords1=[(bl,tr)]
    #data=[map0.submap(bl,tr).data]
    #nx,ny=np.shape(map0.submap(bl,tr).data)
    flux,data,rotcoords1,maps=[],[],[],[]
    for s in submaps:
        if s.meta['exptime'] !=0.0:
            end_time = s.meta['date-obs']
            timestamps.append(end_time)
            bl_rot=solar_rotate_coordinate(bl, new_observer_time=end_time,new_observer_location=s.observer_coordinate)
            tr_rot=solar_rotate_coordinate(tr, new_observer_time=end_time,new_observer_location=s.observer_coordinate)
            if np.isnan(bl_rot.Tx.value): #need to re-define in terms of box width and height
                bl_rot=SkyCoord(tr_rot.Tx - bwidth,tr_rot.Ty - bheight, frame=s.coordinate_frame)

            if type(mask) == np.ndarray:
                fl,dat,current_map=fit_to_mask(s,bl_rot, tr_rot,mask)
            else:
                fl,dat,current_map=flux_from_box(s,bl_rot,tr_rot,mask=False,maps=True, expnorm=expnorm)
            #print(nx,ny,sx,sy,np.shape(maskt),np.shape(dat))
            flux.append(fl)
            data.append(dat)
            maps.append(current_map)
            rotcoords1.append((bl_rot,tr_rot))

    dout={'timestamps':timestamps,'fluxes':flux,'rotated_coords':rotcoords1,'data':data,'maps':maps}
    if bin: #integer number of maps to bin
        newts,newflux,newcoords=[],[],[]
        for j in range(int(len(submaps)/bin)): #+1?
            newts.append(timestamps[j*bin])
            newcoords.append(rotcoords1[j*bin])
            sumflux=np.sum([flux[j*bin+k] for k in range(bin)])
            newflux.append(sumflux)
        dout={'timestamps':newts,'fluxes':newflux,'rotated_coords':newcoords,'data':data}
    return dout

def plot_and_save(mlist,subcoords=False,outname='STEREO_orbit8_',creverse=True,vrange=False):
    '''assume these are de-rotated difference images. bcoords is circle coords'''
    map0=mlist[0]

    palette_name=map0.plot_settings['cmap']
    #if creverse and not palette_name.name.endswith('_r'):
    #    new_cdata=cm.revcmap(palette_name._segmentdata)
    #else:
    #    new_cdata=palette_name._segmentdata
    #new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name+'_r',new_cdata)
    #new_cmap='Greys_r'
    for i,s in enumerate(mlist):
        fig,ax=plt.subplots(figsize=[6,6])
        #ax=fig.add_subplot(111,projection=s.wcs)
        #if creverse and not palette_name.name.endswith('_r'):
        #    s.plot_settings['cmap']=new_cmap
        if subcoords != False:
            s=s.submap(subcoords[i][0],subcoords[i][1])
            #print(subcoords[i])
        if vrange:
            s.plot_settings['norm']=matplotlib.colors.Normalize(vmin=vrange[0],vmax=vrange[1])
        s.plot(axes=ax,cmap='Greys_r')
        plt.colorbar()
        fig.savefig(outname+str(i).zfill(2)+'.png')
    #return submaps


def save_boxfig(submaps,rotcoords1,rotcoords2,outname='STEREO_orbit8_',creverse=False,vrange=False):
    palette_name=submaps[0].plot_settings['cmap']
    if creverse and not palette_name.name.endswith('_r'):
        new_cdata=cm.revcmap(palette_name._segmentdata)
        new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name+'_r',new_cdata)
    for i,s in enumerate(submaps):
        fig=plt.figure(figsize=[8,8])
        ax=fig.add_subplot(111,projection=s.wcs)
        if creverse and not palette_name.name.endswith('_r'):
            s.plot_settings['cmap']=new_cmap
        if vrange:
            s.plot_settings['norm']=matplotlib.colors.Normalize(vmin=vrange[0],vmax=vrange[1])
        s.plot(axes=ax,alpha=0.7)
        b1l=rotcoords1[i][0]
        w1=rotcoords1[i][1].Tx-rotcoords1[i][0].Tx
        h1=rotcoords1[i][1].Ty-rotcoords1[i][0].Ty
        b2l=rotcoords2[i][0]
        w2=rotcoords2[i][1].Tx-rotcoords2[i][0].Tx
        h2=rotcoords2[i][1].Ty-rotcoords2[i][0].Ty

        s.draw_rectangle(b1l, w1,h1,axes=ax,color='b')
        s.draw_rectangle(b2l, w2,h2,axes=ax,color='m')
        fig.savefig(outname+str(i).zfill(2)+'.png')
        
