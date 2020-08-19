 #######################################
#find_hessi_centroids.py
# Erica Lastufka 14/06/2018

#Description: find centroids of RHESSI image contours
#######################################

#######################################
# Usage:

######################################

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
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
from datetime import datetime as dt
from datetime import timedelta
import matplotlib.dates as mdates
import sunpy.physics.differential_rotation as diffrot
import make_timeseries_plot_flux as mtp#import mask_disk
from matplotlib.path import Path
from matplotlib.patches import Wedge
from matplotlib.collections import PatchCollection


#In [24]: f.get_closest_limbcoord(limbcoords=coords)
#coordinates of closest point on AIA limb:  (-927.3644490660265, 216.14800663328145)
#distance to closest point on AIA limb:  11.1450513653 arcsec
#Out[24]: (-927.3644490660265, 216.14800663328145)

def plot_instr_loc():
    '''duplicate plot from stereo science center, only better resolution'''
    #v2=(np.deg2rad(sy),np.deg2rad(sx))
    #vec2=(np.cos(v2[0])*np.cos(v2[1]),np.cos(v2[0])*np.sin(v2[1]))
    heeA=[np.cos(np.deg2rad(.077))*np.cos(np.deg2rad(135.236)),np.cos(np.deg2rad(.077))*np.sin(np.deg2rad(135.236))]
    heeB=[np.cos(np.deg2rad(.288))*np.cos(np.deg2rad(-141.839)),np.cos(np.deg2rad(.288))*np.sin(np.deg2rad(-141.839))]
    heeM=[.1,-1.5]
    #heeE=

    fig,ax=plt.subplots(figsize=[6,6])
    ax.plot([0,0],[0,1],'--k')
    ax.plot([0,heeA[1]],[0,heeA[0]],'--k')
    ax.plot([0,heeB[1]],[0,heeB[0]],'--k')
    ax.scatter([0],[0],color='darkorange',s=400)
    ax.scatter(heeA[1],heeA[0],color='g',s=200)
    ax.scatter(heeB[1],heeB[0],color='b',s=200)
    ax.scatter(heeM[0],heeM[1],color='darkred',s=200)
    ax.scatter([0],[1],color='c',s=300)
    ax.text(0.1,0,'Sun',color='darkorange')
    ax.text(0.1,1,'Earth',color='c')
    ax.text(0.2,-1.5,'Mars',color='darkred')
    ax.text(.5,-1,'STEREO-A',color='g')
    ax.text(-0.75,-1,'STEREO-B',color='b')
    ax.set_xlabel('Y (HEE)')
    ax.set_ylabel('X (HEE)')

    ax.set_ylim([1.5,-2])
    ax.set_xlim([-1.5,1.5])
    fig.show()

def fix_fits_map(fitsfile):
    '''fix units so it can be read as a sunpy map'''
    rhsi_map=sunpy.map.Map(fitsfile)#.rotate()
    rhsi_map.meta['CTYPE1'] = 'HPLN-TAN'
    rhsi_map.meta['CTYPE2'] = 'HPLT-TAN'
    rhsi_map.meta['cunit1']='arcsec'
    rhsi_map.meta['cunit2']='arcsec'
    #dates=dt.strptime('0'+rhsi_map.meta['date-obs'][1:-4],'%d-%b-%Y %H:%M:%S')
    dates=dt.strptime(rhsi_map.meta['date-obs'][:-4],'%Y-%m-%dT%H:%M:%S')
    rhsi_map.meta['date-obs']=dt.strftime(dates,'%Y-%m-%d %H:%M:%S')
    return rhsi_map


def find_centroid_from_map(m,levels=[90],idx=0,show=False):
    '''more generic'''
    cs,hpj_cs=[],[]
    ll=np.max(m.data)*np.array(levels)

    fig,ax=plt.subplots()
    #m.plot(axes=ax,alpha=.5)
    ax.imshow(m.data,alpha=.75,cmap=m.plot_settings['cmap'])
    contour=m.draw_contours(levels=levels*u.percent,axes=ax,frame=m.coordinate_frame)
    #contour=ax.contour(m.data, levels=ll,frame=m.coordinate_frame)
    print len(contour.allsegs[-1])
    c =  center_of_mass(contour.allsegs[-1][idx])
    cs.append(c)
    hpj=m.pixel_to_world(c[0]*u.pixel,c[1]*u.pixel,origin=0)
    hpj_cs.append(hpj)
    ax.plot(c[0],c[1], marker="o", markersize=12, color="red")
    if show:
        fig.show()
    #fig=plt.figure()
    #ax=plt.subplot(projection=m)
    #m.plot(axes=ax,alpha=.5)
    #ax.plot_coord(hpj, marker="o", markersize=12, color="blue")
    ##ax.plot(c[0],c[1], marker="o", markersize=12, color="red")
    #fig.show()
    return cs,hpj_cs,contour

def AIA131_from_rhessi(rhsifits,AIAmap,bl=False,tr=False):
    '''Get the AIA 131 emission inside given RHESSI HXR/SXR contour. Compare with emission at similar heights that is not in the contour '''
    meanin,meanout=[],[]
    rdates,aiadates=[],[]
    rmarea,csum=[],[]

    for rh,aia in zip(rhsifits,AIAmap):
        rmap=fix_fits_map(rh)
        rdates.append(rmap.meta['date-obs'])
        if bl and tr:
            bls = SkyCoord(bl[0] * u.arcsec, bl[1] * u.arcsec, frame=aia.coordinate_frame)
            trs = SkyCoord(tr[0] * u.arcsec, tr[1] * u.arcsec, frame=aia.coordinate_frame)
        aiamap=aia.submap(bls,trs)#sunpy.map.Map(aia)
        aiadates.append(dt.strptime(aiamap.meta['date-obs'][:-3],'%Y-%m-%dT%H:%M:%S'))
        _,hpj_cs,contour=find_centroid_from_map(rmap,levels=[80],show=False)
        #now get the contour coordinates in pixelcoords.
        ccoords=contour.allsegs[-1][0]
        #convert to wcs
        wcs_con=[rmap.pixel_to_world(xy[0]*u.pixel,xy[1]*u.pixel) for xy in ccoords]
        #make into composite map now to line things up
        #cmap=sunpy.map.Map(rmap,aiamap,comp=True)

        #now get corresponding AIA pixels
        #change date info on the skycoords
        aiask=[SkyCoord(wcs.Tx,wcs.Ty,frame=aiamap.coordinate_frame) for wcs in wcs_con]
        aia_datacoords=[aiamap.world_to_pixel(sk) for sk in aiask]
        #change to list of tuples
        tuples=[]
        for pp in aia_datacoords:
            tuples.append((int(pp.x.value),int(pp.y.value)))

        codes=[Path.MOVETO]
        for a in range(1,len(aia_datacoords)-1):
            codes.append(Path.LINETO)
        codes.append(Path.CLOSEPOLY)
        path_aia=Path(tuples,codes)

        #now mask the rest of the AIA image based on outermost and innermost radial boundaries of the rhessi contours
        x=np.linspace(0,np.shape(aiamap.data)[0]-1,np.shape(aiamap.data)[0])
        y=np.linspace(0,np.shape(aiamap.data)[1]-1,np.shape(aiamap.data)[1])
        X,Y=np.meshgrid(x,y)
        points=np.array((X.flatten(),Y.flatten())).T
        boolmask=path_aia.contains_points(points).reshape(X.shape)
        csum.append(np.sum(boolmask))
        maia=np.ma.masked_array(aiamap.data,mask=boolmask)

        mmaia=sunpy.map.Map(maia,aiamap.meta)
        #mask from r_sun to somewhere high in the corona
        #mask_radius(submap,plot=False,radius=1100.,greater=True):
        rmasked=mtp.mask_disk(aiamap,fac=150,greater=True,filled=True) #mask out to +150"
        #rmasked.peek()
        rmasked_all=mtp.mask_disk(rmasked,fac=1.1,greater=False,filled=True)
        #rmasked_all.peek()
        #flatten
        print type(rmasked_all.data),np.mean(rmasked_all.data)
        rmasked_filled=rmasked_all.data#filled(0)
        rmasked_filled[rmasked_filled == 0.0] = np.nan
        rmeans = np.nanmean(rmasked_filled)
        rmsum=np.sum(np.isfinite(rmasked_filled))
        rmarea.append(rmsum)
        meanout.append(rmeans)
        #finally average and store
        aiaflat=maia.filled(0)
        meanin.append(np.mean(aiaflat))
        #mmaia.peek()

    #plot ratio
    fig,ax=plt.subplots()
    ax.plot_date(aiadates, np.array(meanin)/np.array(meanout),'b-')
    dax=ax.twinx()
    dax.plot_date(aiadates, np.array(csum)/np.array(rmarea),color='g')
    ax.set_ylabel("Ratio of mean flux in contour vs. outside contour")
    dax.set_ylabel("Ratio of contour area vs. box area")
    fig.show()

    rmasked_all.peek()

    #return data
    return rdates,aiadates,ccoords,wcs_con, meanin, meanout,csum,rmarea#,ratio

def plot_ratio_better():
    from datetime import timedelta as td
    params=pickle.load(open('rhessi_hi_131_ratio.p','rb'))
    meanin=params['meanin']
    meanout=params['meanout']
    aiadt=params['aia_dt']
    csum=params['contour_area'][:-1]
    rmarea=params['box_area'][:-1]
    aratio=np.array(csum)/np.array(rmarea)
    mratio=np.array(meanin)/np.array(meanout)
    fig,ax=plt.subplots()
    ax.plot_date(aiadt[:-1], mratio[:-1],'b-')
    ax.plot_date(aiadt[:-1], mratio[:-1])
    #ax.plot_date(aiadt[:-1], mratio[:-1]*aratio*100.,'r-')
    #ax.plot_date(aiadt[:-1], mratio[:-1]*aratio*100.,color='r')
    dax=ax.twinx()
    dax.plot_date(aiadt[:-1], aratio,'g-')
    dax.plot_date(aiadt[:-1], aratio,color='g')
    dax.set_ylabel("Ratio of contour area vs. box area")
    ax.set_ylabel("Ratio of mean flux in contour vs. outside contour")
    myFmt = mdates.DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_formatter(myFmt)
    ax.set_xlabel('Time on 2013 May 01')
    ax.set_xlim([aiadt[0]-td(seconds=30),aiadt[-1]])
    #ax.set_ylim([.5,1.0])
    dax.set_ylim([0,0.1])
    plt.gcf().autofmt_xdate()

    fig.show()


def find_hessi_centroids(savfile='rhsi_1min_may1.sav',aiafile='AIA193test37.fits',offset=(-927.364,216.148),plot=False):
    hsi_dict=readsav(savfile,python_dict=True)
    aiamap=sunpy.map.Map(aiafile)
    cs=[]
    hpj_cs=[]
    dates=[]
    dist,speed=[],[]
    c0=offset
    karr=hsi_dict.keys()
    #karr.remove('n234')
    karr.sort()
    for i,k in enumerate(karr):
        im=hsi_dict[k].field(0)[0]
        dates.append(dt.strptime('0'+hsi_dict[k].field(5)[0][1:-4],'%d-%b-%Y %H:%M:%S'))
        levels=np.max(im)*np.array([.3,.5,.7,.9])
        #get new x-axis
        xc=hsi_dict[k].field(1)[0]
        yc=hsi_dict[k].field(2)[0]
        dx=hsi_dict[k].field(3)[0]
        dy=hsi_dict[k].field(4)[0]
        fig=plt.figure()
        ax=plt.subplot(projection=aiamap)
        contour=ax.contour(im, levels=levels)
        c =  center_of_mass(contour.allsegs[-1][0])
        cs.append(c)
        hpj=transform_to_wcs(c,xc,yc,dx,dy,np.shape(im)) #convert to hpc coords ... hopefully this works!
        hpj_cs.append(hpj)
        if k == 'c233':
            timed=60.
        else:
            timedt=dates[i]-dates[i-1]
            timed=timedt.seconds
        atts=calc_distance(hpj,c0,tds=timed)
        dist.append(atts[0])
        speed.append(atts[1])
        #ax.plot(c[0],c[1], marker="o", markersize=12, color="red")
        if plot:
            fig.show()
    #plot in data coords. do the conversion to actual coords here or in IDL?
    print len(dates),len(cs)
    #fig,ax=plt.subplots()
    #ax.plot_date(dates,hpj_cs)
    #fig.show()
    fig,ax=plt.subplots()
    hx=[h[0] for h in hpj_cs]
    hy=[h[1] for h in hpj_cs]
    ax.scatter(hx,hy)
    #for i,l in enumerate(dates):
    #    ax.annotate(dt.strftime(l,'%d-%b-%Y %H:%M:%S'),(hx[i],hy[i]))
    ax.set_xlim([-1200,-700])
    ax.set_ylim([-70,450])
    fig.show()

    return dates,hpj_cs,dist,speed

def hessi_on_norh(savfile='rhsi_hi_45s.sav',norhfile='../NORH/130501_17g10s/ifa130501_023205',offset=(-927.364,216.148),plot=False):
    #hsi_dict=readsav(savfile,python_dict=True)
    norhmap=sunpy.map.Map(norhfile)#.rotate()  #probably need to transform the coordinates
    norhmap.meta['CTYPE1'] = 'HPLN-TAN'
    norhmap.meta['CTYPE2'] = 'HPLT-TAN'
    norhmap.rotate()
    norhmap.plot_settings['cmap']=sunpy.cm.cmlist['sdoaia304']
    snorh=norhmap.submap(SkyCoord([-1410,-750]*u.arcsec,[-50,500]*u.arcsec,frame=norhmap.coordinate_frame))

    rhsi_map=sunpy.map.Map('c233.fits')#.rotate()
    rhsi_map.meta['CTYPE1'] = 'HPLN-TAN'
    rhsi_map.meta['CTYPE2'] = 'HPLT-TAN'
    rhsi_map.meta['cunit1']='arcsec'
    rhsi_map.meta['cunit2']='arcsec'
    rhsi_map.meta['date-obs']=rhsi_map.meta['date-obs'][1:]
    rhsi_map.plot_settings['cmap']=cm.winter
    rhsi_map.rotate()
    srhsi=rhsi_map.submap(SkyCoord([-1410,-750]*u.arcsec,[-50,500]*u.arcsec,frame=norhmap.coordinate_frame))

    comp_map = sunpy.map.Map(snorh,srhsi, composite=True)
    comp_map.set_alpha(0, 1)
    #comp_map.set_alpha(1, 0.5)
    #comp_map.set_colors(1,'winter')
    comp_map.set_levels(1, [50,70,90],percent=True)

    #plot settings
    #nps=comp_map.get_plot_settings(0)
    #nps['norm']=colors.Normalize(vmin=5000,vmax=15000)
    #comp_map.set_plot_settings(0,nps)

    fig=plt.figure()
    comp_map.plot()
    if plot:
        fig.show()
    return comp_map


def hessi_on_stereo(savfile='rhsi_hi_45s.sav',norhfile='../STEREO/A/20130501_023445_15euA.fts',offset=(-927.364,216.148),plot=False):
    #hsi_dict=readsav(savfile,python_dict=True)
    norhmap=sunpy.map.Map(norhfile).rotate()  #probably need to transform the coordinates
    norhmap.meta['CTYPE1'] = 'HPLN-TAN'
    norhmap.meta['CTYPE2'] = 'HPLT-TAN'
    norhmap.rotate()
    norhmap.plot_settings['cmap']=sunpy.cm.cmlist['sdoaia304']
    #snorh=norhmap.submap(SkyCoord([-1410,-750]*u.arcsec,[-50,500]*u.arcsec,frame=norhmap.coordinate_frame))
    #flip the x-axis ... how can I do this? What about we mirror the data array? try that...
    data=norhmap.data
    #print np.shape(data)
    fmap= sunpy.map.Map(np.fliplr(data),norhmap.meta)
    snorh=fmap
    #now get the submap
    snorh=fmap.submap(SkyCoord([-1410,-750]*u.arcsec,[-50,500]*u.arcsec,frame=fmap.coordinate_frame))
    rhsi_map=sunpy.map.Map('rhsi_hi6.fits')#.rotate()
    rhsi_map.meta['CTYPE1'] = 'HPLN-TAN'
    rhsi_map.meta['CTYPE2'] = 'HPLT-TAN'
    rhsi_map.meta['cunit1']='arcsec'
    rhsi_map.meta['cunit2']='arcsec'
    rhsi_map.meta['date-obs']='2013/05/01T'+rhsi_map.meta['date-obs'][12:]
    rhsi_map.plot_settings['cmap']=cm.autumn
    #rhsi_map.rotate()
    map304=sunpy.map.Map('AIA193test35.fits')
    srhsi=rhsi_map.submap(SkyCoord([-1410,-750]*u.arcsec,[-50,500]*u.arcsec,frame=map304.coordinate_frame))
    #srhsi=rhsi_map
    comp_map = sunpy.map.Map(snorh,srhsi, composite=True)
    comp_map.set_alpha(0, 1)
    #comp_map.set_alpha(1, 0.5)
    #comp_map.set_colors(1,'winter')
    comp_map.set_levels(1, [79,80,89,90,91],percent=True)

    #plot settings
    nps=comp_map.get_plot_settings(0)
    nps['norm']=colors.Normalize(vmin=-1000,vmax=3000)
    comp_map.set_plot_settings(0,nps)

    fig=plt.figure()
    comp_map.plot()
    if plot:
        fig.show()
    return fmap,rhsi_map,comp_map

def calc_collision_time(E, delta, n):
    '''E in keV, delta is spectral index of electrons, n in units of cm^-3. From Lin 1974, eq. 3.4 '''
    tc=(E**1.5)/(4.9*10**-9 *(delta + 0.5)*n)
    print tc, 'seconds'
    return tc

def both_hessi(levels=[49,50,69,70,89,90],save=False,legend=True):
    comp_maps=glob.glob('comp_map*.p')
    newhxrmaps=glob.glob('nine_hi_clean*.fits')
    marr,dmaps,centroids=[],[],[]
    for i,c in enumerate(comp_maps):
        imname='comp_map'+str(i+2)+'.png'
        m=pickle.load(open(c,'rb'))
        m.remove_map(1)
        hxrmap=fix_fits_map(newhxrmaps[i])
        _,cs,_=find_centroid_from_map(hxrmap)
        centroids.append(cs)
        m.add_map(hxrmap)
        m.set_alpha(0, .75)
        m.set_alpha(1, 1)
        datemap=m.get_map(0)
        #dmaps.append(datemap)
        palette_name=datemap.plot_settings['cmap']#sunpy.cm.cmlist[cmap]#m17.plot_settings['cmap']
        new_cdata=cm.revcmap(palette_name._segmentdata)
        new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name+'_r',new_cdata)
        newsettings=datemap.plot_settings
        newsettings['cmap']=new_cmap
        newsettings['norm']=colors.Normalize(0,1000)
        m.set_plot_settings(0,newsettings)

        mycolors = [(1, 1, 1), (0, 0, 0), (0, 0,0)]
        #n_bins = [1, 1, 1, 2]  # Discretizes the interpolation into bins
        onecm = colors.LinearSegmentedColormap.from_list('mycmap', mycolors, N=1)
        newsettings=m.get_plot_settings(1)
        newsettings['cmap']=onecm
        m.set_plot_settings(1,newsettings)
        m.set_levels(1, levels[2:],percent=True)

        try:
            m.set_alpha(2,1)
            m.set_levels(2,levels ,percent=True)
            mycolors = [(1, 0, 1), (1, 1, 0), (1, 0, 1)]  # R -> G -> B

             # R -> G -> B
            #n_bins = [1, 1, 1, 2]  # Discretizes the interpolation into bins
            onecm2 = colors.LinearSegmentedColormap.from_list('mycmap', mycolors, N=1)
            newsettings2=m.get_plot_settings(2)
            newsettings2['cmap']=onecm
            m.set_plot_settings(2,newsettings2)
            m.set_levels(1, levels,percent=True)
            newsettings=m.get_plot_settings(1)
            newsettings['cmap']=onecm2
            m.set_plot_settings(1,newsettings)
            marr.append(m)
        except IndexError:
            pass
    #print centroids
    fig=plt.figure(figsize=[15,6])
    codes=[131,132,133]
    ecen=[[-991.86,266.05],[-997.50,267.33],[-957.29,221.86]]
    ex=[17.43,30.51,14.41]
    ey=[12.21,18.38,9.06]
    size_err=[27.5,1.5,37.93]
    sizes=[136.44,200.,44.49]
    for i in range(0,3):
        m=marr[i]
        d=m.get_map(0)
        code=codes[i]
        cs=centroids[i]
        ax=fig.add_subplot(code)
        m.plot(axes=ax,title=d.meta['date-obs'])
        #plot error bars from forward fit
        ax.errorbar(ecen[i][0]/3600.,ecen[i][1]/3600.,xerr=ex[i]/3600.,yerr=ey[i]/3600.,color='k')
        #plot error bars for source size - or plot two circles indicating upper and lower bounds?
        #ax.errorbar(cs[0].Tx.value/3600.,cs[0].Ty.value/3600.,xerr=size_err[i]/3600.,yerr=size_err[i]/3600.,color='k')
        print cs[0].Tx.value/3600.,cs[0].Ty.value/3600.
        p=Wedge((cs[0].Tx.value/3600.,cs[0].Ty.value/3600.), sizes[i]/7200., 0, 360, width=size_err[i]/3600.)  # need to see how it calculates the width
        pc = PatchCollection([p], alpha=0.2,color='k')
        ax.set_ylim([0.02778,0.097223]) #100 - 350 arcsec
        ax.set_xlim([-.3055,-.23])
        #convert labels to arcseconds instead of degrees, every other one is empty
        ax.set_xticks([-0.2917,-0.2639,-0.2362])
        ax.set_yticks([0.0278,0.041667,0.0556,0.069445,0.08334,0.097223])
        xlab=ax.get_xticks()
        ylab=ax.get_yticks()
        newxlab=[str(int(xl*3600.)) for xl in xlab]
        newylab=[str(int(xl*3600.)) for xl in ylab]
        ax.set_xticklabels(newxlab)
        ax.set_yticklabels(newylab)
        ax.set_xlabel('X (arcsec)')
        ax.set_ylabel('Y (arcsec')
        ax.set_aspect('equal')
        ax.add_collection(pc)

        if legend and i == 0:
            labels=['4-8 keV','13-30 keV']
            proxy = [plt.Rectangle((0,0),1,1,color='m'),plt.Rectangle((0,0),1,1,color='k')] #how to make it green or blue?
            ax.legend(proxy, labels,loc="lower left")

    #ax.xaxis.set_visible(False)
    #ax.yaxis.set_visible(False)
    if save:
        fig.savefig(imname)
    else:
        fig.show()
    return marr

def both_hessi2(levels=[49,50,69,70,89,90],save=False,legend=True):
    comp_maps=glob.glob('comp_map_final*.p')
    #newhxrmaps=['nine_hi_clean1.fits','nine_hi_clean1.fits','nine_hi_clean2.fits','nine_hi_clean2.fits']
    highfile='clean_20130501_nt_19ttt_13-30keV_u89_bw14.fits'

    hrhsi_map=sunpy.map.Map(highfile)#.rotate()
    hrhsi_map.meta['CTYPE1'] = 'HPLN-TAN'
    hrhsi_map.meta['CTYPE2'] = 'HPLT-TAN'
    hrhsi_map.meta['cunit1']='arcsec'
    hrhsi_map.meta['cunit2']='arcsec'
    hrhsi_map.meta['date-obs']=hrhsi_map.meta['date-obs']
    hrhsi_map.plot_settings['cmap']=cm.Greens
    hrhsi_map.rotate()
    #hsrhsi=hrhsi_map.submap(SkyCoord([-1210,-850]*u.arcsec,[50,450]*u.arcsec,frame=saia.coordinate_frame))

    marr,dmaps,centroids=[],[],[]
    for i,c in enumerate(comp_maps):
        imname='comp_map'+str(i+2)+'.png'
        m=pickle.load(open(c,'rb'))
        #m.remove_map(1)
        if i < 3:
            #m.remove_map(1)
            hxrmap=hrhsi_map#fix_fits_map(newhxrmaps[i])
            _,cs,_=find_centroid_from_map(hxrmap)
            centroids.append(cs)
            m.add_map(hxrmap)

        m.set_alpha(0, .75)
        m.set_alpha(1, 1)
        datemap=m.get_map(0)
        #dmaps.append(datemap)
        palette_name=datemap.plot_settings['cmap']#sunpy.cm.cmlist[cmap]#m17.plot_settings['cmap']
        new_cdata=cm.revcmap(palette_name._segmentdata)
        new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name+'_r',new_cdata)
        newsettings=datemap.plot_settings
        newsettings['cmap']=new_cmap
        newsettings['norm']=colors.Normalize(0,1000)
        m.set_plot_settings(0,newsettings)

        blackcolors = [(1, 1, 1), (0, 0, 0), (0, 0,0)]
        #n_bins = [1, 1, 1, 2]  # Discretizes the interpolation into bins
        blackcm = colors.LinearSegmentedColormap.from_list('mycmap', blackcolors, N=1)
        pinkcolors = [(0,.5,0),(1,.5,0),(0,.5,0)]#[(1, 0, 1), (1, 1, 0), (1, 0, 1)]
        pinkcm = colors.LinearSegmentedColormap.from_list('mycmap', pinkcolors, N=1)
        newsettings=m.get_plot_settings(1)
        newsettings['cmap']=pinkcm
        m.set_plot_settings(1,newsettings)
        m.set_levels(1, levels,percent=True)
        marr.append(m)

        try:
            m.set_alpha(2,1)
            m.set_levels(2,levels ,percent=True)
            newsettings2=m.get_plot_settings(2)
            newsettings2['cmap']=blackcm
            m.set_plot_settings(2,newsettings2)
            if i != 0:
                m.set_levels(1, levels,percent=True)
            else:
                m.set_levels(1, levels[1:],percent=True)
            newsettings=m.get_plot_settings(1)
            newsettings['cmap']=pinkcm
            m.set_plot_settings(1,newsettings)
            #marr.append(m)
        except IndexError:
            pass
    #print centroids
    fig=plt.figure(figsize=[15,8])
    codes=[241,242,243,244,245,246,247,248]
    ecen=[[-991.86,266.05],[-997.50,267.33],[-957.29,221.86]]
    ex=[17.43,30.51,14.41]
    ey=[12.21,18.38,9.06]
    size_err=[27.5,1.5,37.93]
    sizes=[136.44,200.,44.49]
    for i in range(0,8):
        m=marr[i]
        d=m.get_map(0)
        code=codes[i]
        #cs=centroids[i]
        ax=fig.add_subplot(code)
        m.plot(axes=ax,title=d.meta['date-obs'])
        #plot error bars from forward fit
        #ax.errorbar(ecen[i][0]/3600.,ecen[i][1]/3600.,xerr=ex[i]/3600.,yerr=ey[i]/3600.,color='k')
        #plot error bars for source size - or plot two circles indicating upper and lower bounds?
        #ax.errorbar(cs[0].Tx.value/3600.,cs[0].Ty.value/3600.,xerr=size_err[i]/3600.,yerr=size_err[i]/3600.,color='k')
        #print cs[0].Tx.value/3600.,cs[0].Ty.value/3600.
        #p=Wedge((cs[0].Tx.value/3600.,cs[0].Ty.value/3600.), sizes[i]/7200., 0, 360, width=size_err[i]/3600.)  # need to see how it calculates the width
        #pc = PatchCollection([p], alpha=0.2,color='k')
        ax.set_ylim([0.02778,0.097223]) #100 - 350 arcsec
        ax.set_xlim([-.3055,-.23])
        #convert labels to arcseconds instead of degrees, every other one is empty
        ax.set_xticks([-0.2917,-0.2639,-0.2362])
        ax.set_yticks([0.0278,0.041667,0.0556,0.069445,0.08334,0.097223])
        xlab=ax.get_xticks()
        ylab=ax.get_yticks()
        newxlab=[str(int(xl*3600.)) for xl in xlab]
        newylab=[str(int(xl*3600.)) for xl in ylab]
        ax.set_xticklabels(newxlab)
        ax.set_yticklabels(newylab)
        ax.set_xlabel('X (arcsec)')
        ax.set_ylabel('Y (arcsec)')
        if i not in [0,4]: #turn off y axis
            ax.yaxis.set_visible(False)
        if i not in [4,5,6,7]:
            ax.xaxis.set_visible(False)
        ax.set_aspect('equal')
        #ax.add_collection(pc)

        if legend and i == 4:
            labels=['4-8 keV','13-30 keV']
            proxy = [plt.Rectangle((0,0),1,1,color='g'),plt.Rectangle((0,0),1,1,color='k')] #how to make it green or blue?
            ax.legend(proxy, labels,loc="lower left")

    plt.tight_layout()
    #ax.xaxis.set_visible(False)
    #ax.yaxis.set_visible(False)
    if save:
        fig.savefig(imname)
    else:
        fig.show()
    return marr


def norh_fig1():
    from calc_norh_alpha import cmap_plot
    cmap=pickle.load(open('norh_304_fig1.p','rb'))
    fig,ax=plt.subplots()
    palette_name=sunpy.cm.cmlist['sdoaia304']#cmap.get_plot_settings(0)['cmap']#sunpy.cm.cmlist[cmap]#m17.plot_settings['cmap']
    new_cdata=cm.revcmap(palette_name._segmentdata)
    new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name+'_r',new_cdata)
    #newsettings=cmap.get_plot_settings(0)
    #newsettings['cmap']=new_cmap
    #newsettings['norm']=colors.Normalize(0,1000)
    newsettings={'cmap': new_cmap,
 'interpolation': 'nearest',
 'norm': colors.Normalize(-10,550),
 'origin': 'lower'}
    cmap.set_plot_settings(0,newsettings)
    #norhset=cmap.get_plot_settings(1)
    #norhset['linewidths']=2
    #norhset['norm']=colors.Normalize(0,100)
    #cmap.set_plot_settings(1,norhset)
    #cmap_plot(cmap,axes=ax)
    cmap.set_levels(1,[5.45,5.5,6.9,7,9.9,10],percent=True)
    cmap.plot(axes=ax)
    ax.set_ylim([0.02778,0.097223]) #100 - 350 arcsec
    ax.set_xlim([-.3055,-.23])
    ax.set_xticks([-0.2917,-0.2639,-0.2362])
    ax.set_yticks([0.0278,0.041667,0.0556,0.069445,0.08334,0.097223])
    xlab=ax.get_xticks()
    ylab=ax.get_yticks()
    newxlab=[str(int(xl*3600.)) for xl in xlab]
    newylab=[str(int(xl*3600.)) for xl in ylab]
    ax.set_xticklabels(newxlab)
    ax.set_yticklabels(newylab)
    ax.set_xlabel('X (arcsec)')
    ax.set_ylabel('Y (arcsec')
    ax.set_aspect('equal')
    ax.set_title('2013-05-01 02:37:15')
    labels=['NORH 17 GHz']
    proxy = [plt.Rectangle((0,0),.5,.5,facecolor='none')] #how to make it green or blue?
    ax.legend(proxy, labels,loc="lower right")
    fig.show()


def hessi_on_304(savfile='c233.fits',aiamap='../AIA/304/mapcube304_peak.p',midx=0,offset=(-927.364,216.148),plot=False):
    #hsi_dict=readsav(savfile,python_dict=True)
    aiamaps=pickle.load(open(aiamap,'rb'))
    norhmap=aiamaps[midx]
    saia=norhmap.submap(SkyCoord([-1410,-750]*u.arcsec,[-50,500]*u.arcsec,frame=norhmap.coordinate_frame))

    rhsi_map=sunpy.map.Map(savfile)#.rotate()
    rhsi_map.meta['CTYPE1'] = 'HPLN-TAN'
    rhsi_map.meta['CTYPE2'] = 'HPLT-TAN'
    rhsi_map.meta['cunit1']='arcsec'
    rhsi_map.meta['cunit2']='arcsec'
    rhsi_map.meta['date-obs']=rhsi_map.meta['date-obs']#[1:]
    rhsi_map.plot_settings['cmap']=cm.winter
    rhsi_map.rotate()
    srhsi=rhsi_map.submap(SkyCoord([-1410,-750]*u.arcsec,[-50,500]*u.arcsec,frame=norhmap.coordinate_frame))
    #snorh=norhmap.submap(SkyCoord([-1410,-750]*u.arcsec,[-50,500]*u.arcsec,frame=norhmap.coordinate_frame))

    comp_map = sunpy.map.Map(saia,srhsi, composite=True)
    comp_map.set_alpha(0, 1)
    #comp_map.set_alpha(1, 0.5)
    #comp_map.set_colors(1,'winter')
    comp_map.set_levels(1, [50,70,90],percent=True)

    #plot settings
    #nps=comp_map.get_plot_settings(0)
    #nps['norm']=colors.Normalize(vmin=5000,vmax=15000)
    #comp_map.set_plot_settings(0,nps)

    print 'AIA time ',saia.meta['date-obs']
    print 'RHESSI time ',srhsi.meta['date-obs']
    fig=plt.figure()
    comp_map.plot()
    if plot:
        fig.show()
    return comp_map

def vff_contours(gkey='vff_t0*.fits'):
    t0=glob.glob(gkey)
    t0maps=[]
    keys,cvals=[],[]
    cmaps=[cm.Greys, cm.Blues, cm.Greens,cm.Greens, cm.OrRd, cm.Oranges,cm.Reds,cm.Purples]
    for i,t in enumerate(t0):
        rhsi_map=sunpy.map.Map(t)#.rotate()
        rhsi_map.meta['CTYPE1'] = 'HPLN-TAN'
        rhsi_map.meta['CTYPE2'] = 'HPLT-TAN'
        rhsi_map.meta['cunit1']='arcsec'
        rhsi_map.meta['cunit2']='arcsec'
        rhsi_map.meta['date-obs']=rhsi_map.meta['date-obs']
        rhsi_map.meta['rsun_obs']=952.32
        rhsi_map.plot_settings['cmap']=cmaps[i]
        #a=rhsi_map.plot_settings.pop('norm')
        rhsi_map.plot_settings['label']=t[t.find('n_')+2:t.rfind('_')]
        rhsi_map.rotate()
        if np.sum(rhsi_map.data) != 0.0:
            t0maps.append(rhsi_map)
            keys.append(t[t.find('n_')+2:t.rfind('_')])
            cvals.append(cmaps[i])

    comp_map0 = sunpy.map.Map(t0maps, composite=True)
    comp_map0.set_alpha(0, 1)
    #comp_map.set_alpha(1, 0.5)
    #comp_map.set_colors(1,'winter')
    #ps=comp_map0.get_plot_settings(0)
    #ps.pop('norm')
    for i in range(0,len(t0maps)):
        if i !=0:
            comp_map0.set_levels(i, [69,70,89,90],percent=True)
        ps=comp_map0.get_plot_settings(i)
        lm=comp_map0.get_map(i)
        ps['norm']=matplotlib.colors.Normalize(np.min(lm.data),np.max(lm.data))
        #print ps
        comp_map0.set_plot_settings(i,ps)
    #should put a legend?
    fig,ax=plt.subplots()
    comp_map0.plot(axes=ax)
    comp_map0.draw_limb(axes=ax,**{'color':'k'})
    ax.legend()
    fig.show()
    print ps
    #print keys,cvals



def fig4hessi_on_AIA(lowfile='clean_20130501_nt_19ttt_6-8keV_u6789.fits',highfile='clean_20130501_nt_19ttt_13-30keV_u89_bw14.fits',wave=304,offset=(-927.364,216.148),plot=False,creverse=True,legend=False):
    #hsi_dict=readsav(savfile,python_dict=True)
    if wave == 304:
        aiafile='peak304submap.p'
    else:
        aiafile='peak94submap.p'

    saia=pickle.load(open(aiafile,'rb'))

    #peakmap=aiamaps.maps[2]
    #saia=peakmap.submap(SkyCoord([-1210,-850]*u.arcsec,[50,450]*u.arcsec,frame=peakmap.coordinate_frame))

    lrhsi_map=sunpy.map.Map(lowfile)#.rotate()
    lrhsi_map.meta['CTYPE1'] = 'HPLN-TAN'
    lrhsi_map.meta['CTYPE2'] = 'HPLT-TAN'
    lrhsi_map.meta['cunit1']='arcsec'
    lrhsi_map.meta['cunit2']='arcsec'
    lrhsi_map.meta['date-obs']=lrhsi_map.meta['date-obs']
    lrhsi_map.plot_settings['cmap']=cm.Blues
    lrhsi_map.rotate()
    lsrhsi=lrhsi_map.submap(SkyCoord([-1210,-850]*u.arcsec,[50,450]*u.arcsec,frame=saia.coordinate_frame))

    hrhsi_map=sunpy.map.Map(highfile)#.rotate()
    hrhsi_map.meta['CTYPE1'] = 'HPLN-TAN'
    hrhsi_map.meta['CTYPE2'] = 'HPLT-TAN'
    hrhsi_map.meta['cunit1']='arcsec'
    hrhsi_map.meta['cunit2']='arcsec'
    hrhsi_map.meta['date-obs']=hrhsi_map.meta['date-obs']
    hrhsi_map.plot_settings['cmap']=cm.Greens
    hrhsi_map.rotate()
    hsrhsi=hrhsi_map.submap(SkyCoord([-1210,-850]*u.arcsec,[50,450]*u.arcsec,frame=saia.coordinate_frame))

    comp_map = sunpy.map.Map(saia,lsrhsi,hsrhsi,composite=True)

    if creverse:
        aiaps=comp_map.get_plot_settings(0)
        if wave != 304:
            palette_name = sunpy.cm.cmlist['sdoaia193']
        else:
            palette_name=aiaps['cmap']#.name#m17.plot_settings['cmap']

        new_cdata=cm.revcmap(palette_name._segmentdata)
        new_cmap=colors.LinearSegmentedColormap(palette_name.name+'_r',new_cdata)
        newnorm=colors.Normalize(vmin=.1*np.min(saia.data),vmax=.9*np.max(saia.data))
        aiaps['cmap']=new_cmap
        aiaps['norm']=newnorm
        comp_map.set_plot_settings(0,aiaps)

    lps=comp_map.get_plot_settings(1)
    lps['cmap']=cm.Blues
    lps['norm']=colors.Normalize(vmin=-1500,vmax=500)
    comp_map.set_plot_settings(1,lps)
    hps=comp_map.get_plot_settings(2)
    hps['norm']=colors.Normalize(vmin=-1000,vmax=500)
    hps['cmap']=cm.Greens
    comp_map.set_plot_settings(2,hps)

    comp_map.set_alpha(0, 1)
    #comp_map.set_alpha(1, 0)
    #comp_map.set_colors(1,'Blues')
    #comp_map.set_colors(1,'Greens')
    comp_map.set_levels(1, [29,29.5,30,49,50,69,70,89,90],percent=True)
    comp_map.set_levels(2, [49,50,69,70,89,90],percent=True)

    #plot settings
    #nps=comp_map.get_plot_settings(0)
    #nps['norm']=colors.Normalize(vmin=5000,vmax=15000)
    #comp_map.set_plot_settings(0,nps)

    fig=plt.figure()
    comp_map.plot()
    if legend:
        labels=['6-8 keV','13-30 keV']
        proxy = [plt.Rectangle((0,0),1,1),plt.Rectangle((0,0),1,1)] #how to make it green or blue?
        axes.legend(proxy, labels,loc="upper left")

    if plot:
        fig.show()
    return comp_map


def plot_coords_on_aia(hpj_cs,hpj2=False,aiafile='AIA193test37.fits'):
    aiamap=sunpy.map.Map(aiafile)
    norm=colors.Normalize(vmin=np.min(aiamap.data),vmax=.9*np.max(aiamap.data))
    #cmap=sunpy.cm.cmlist['sdoaia193']
    palette_name=aiamap.plot_settings['cmap']
    new_cdata=cm.revcmap(palette_name._segmentdata)
    new_cmap=mpl.colors.LinearSegmentedColormap(palette_name.name+'_r',new_cdata)
    aiamap.plot_settings['cmap']=new_cmap

    sk_cs,sk_cs2=[],[]
    #convert coords to SkyCoords
    for cs in hpj_cs:
        sk_cs.append(SkyCoord(cs[0]*u.arcsec,cs[1]*u.arcsec, frame=aiamap.coordinate_frame))
    if hpj2:
        for h in hpj2:
            sk_cs2.append(SkyCoord(h[0]*u.arcsec,h[1]*u.arcsec, frame=aiamap.coordinate_frame))

    fig=plt.figure()
    ax=plt.subplot(projection=aiamap)
    aiamap.plot(axes=ax,norm=norm)
    aiamap.draw_limb(color='k')
    for i,cs in enumerate(sk_cs):
        if i==0:
            ax.plot_coord(cs,'o',color='k',markersize=8,label='RHESSI')
        else:
            ax.plot_coord(cs,'o',color='k',markersize=8)

    if hpj2:
        for j,cs in enumerate(sk_cs2):
            if j == 0:
                ax.plot_coord(cs,'o',color='g',markersize=8,label='AIA blob')
            else:
                ax.plot_coord(cs,'o',color='g',markersize=8)

        ax.legend()
    fig.show()
    return sk_cs

def plot_45s_atts(low=True):
    if low:
        fitsf=glob.glob('low_clean*.fits')
    else:
        #fitsf=glob.glob('clean_hi*.fits')
        #fitsf.pop(2)
        fitsf = pickle.load(open('cmaps_90s.p','rb'))

    d,s,t=[],[],[]
    for f in fitsf:
        #amap=sunpy.map.Map(f)
        if  low:
            amap=sunpy.map.Map(f)
            amap=fix_fits_map(amap)
        else:
            amap=f.get_map(2)
        _,cs,_=find_centroid_from_map(amap)
        dd,ss=calc_distance([cs[0].Tx.value,cs[0].Ty.value],[-923,214],tds=45.)
        d.append(dd)
        s.append(ss)
        if not low:
            t.append(dt.strptime(amap.meta['date-obs'],'%Y-%m-%dT%H:%M:%S.000'))
        else:
            t.append(dt.strptime(amap.meta['date-obs'],'%Y-%m-%d %H:%M:%S'))


    plot_atts(t,d,s,fit_accel=False,fit_speed=True)

def plot_atts(dates,dist,speed,params2=False, fit_accel=True,recalc_speed=True, fit_speed=True):
    if recalc_speed:
        speed=[(dist[i+1]-d)*757./(dates[i+1]-dates[i]).seconds for i,d in enumerate(dist[:-1])]
    if fit_accel==1:
        xvec=mdates.date2num(dates[1:])
        fit=np.polyfit(xvec,speed,deg=fit_accel)
        yvec=fit[0]*xvec+fit[1]
        print fit[0]*1.1574*10**-5,fit[1] #need to convert back to datetimes...

        #fit = [ 2.53306701e+02, -1.86177666e+08] what are the units here? number of days on the time axis
        #rough calculation gives slope = 2.966 km/s^2. So maybe slope = 2.533 km/s^2?
    elif fit_accel==2:
        xvec=mdates.date2num(dates)
        fit=np.polyfit(xvec,speed,deg=fit_accel)
        yvec=fit[0]*xvec**2.+fit[1]*xvec + fit[2]
    elif fit_accel==3:
        xvec=mdates.date2num(dates)
        fit=np.polyfit(xvec,speed,deg=fit_accel)
        yvec=fit[0]*xvec**3.+fit[1]*xvec**2. + fit[2]*xvec+fit[3]

    if fit_speed:
        xvec=mdates.date2num(dates[1:])
        fit=np.polyfit(xvec,np.array(dist[1:])*.757,deg=1)
        yvec2=fit[0]*xvec+fit[1]
        print fit[0]*1.1574*10**-5,fit[1] #need to convert back to datetimes...

    if params2:
        dates2,dist2,speed2=params2
        if recalc_speed:
            speed2=[(dist2[i+1]-d)*757./(dates2[i+1]-dates2[i]).seconds for i,d in enumerate(dist2[:-1])]

    fig,ax=plt.subplots()
    ax.plot_date(dates, np.array(dist)*.757,'v',color='m',markersize=10,label='RHESSI source height')
    ax.plot_date(dates, np.array(dist)-100.,'o',markersize=10,label='RHESSI source speed')
    if fit_speed:
        ax.plot_date(dates[1:],yvec2,'m--')
    dax=ax.twinx()
    dax.plot_date(dates[1:],speed,'o',markersize=10,label='Speed')
    if fit_accel:
        dax.plot_date(dates[1:],yvec,'b--')
    if params2:
        #ax.plot_date(dates2, np.array(dist2)*.757,'v',color='g',markersize=10,label='EUV blob height')
        ax.plot_date(dates2, np.array(dist2)-150,'v',color='g',markersize=10,label='EUV blob speed')
        dax.plot_date(dates2[1:],speed2,'o',color='k',markersize=10,label='EUV blob speed')

    myFmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    ax.set_xlabel('Time on 2013 May 01')
    ax.set_ylabel('Height above limb (Mm)')
    dax.set_ylabel('Speed (km/s)')
    ax.set_xlim([dates[0]-timedelta(minutes=1),dates[-1]+timedelta(minutes=1)])
    ax.set_ylim([0,100])
    plt.gcf().autofmt_xdate()
    ax.legend(loc='upper left')
    #dax.legend(loc=(xvec[0],70))
    fig.show()

def calc_distance(cpair,c0,speed=True,tds=60.):
    distance=np.sqrt((cpair[0]-c0[0])**2.+(cpair[1]-c0[1])**2.)
    if speed: # units are arcseconds. Return Mm/s
        speed=distance*757./tds
        return distance,speed
    else:
        return distance

def transform_to_wcs(cpair,xc,yc,dx,dy,imshape):
    xcd=imshape[0]/2
    ycd=imshape[1]/2
    dxd=cpair[0]-xcd
    dyd=cpair[1]-ycd
    wcsx=xc+dxd*dx
    wcsy=yc+dyd*dy
    wcs_cpair=[wcsx,wcsy]
    return wcs_cpair

def center_of_mass(X):
    # calculate center of mass of a closed polygon
    x = X[:,0]
    y = X[:,1]
    g = (x[:-1]*y[1:] - x[1:]*y[:-1])
    A = 0.5*g.sum()
    cx = ((x[:-1] + x[1:])*g).sum()
    cy = ((y[:-1] + y[1:])*g).sum()
    return 1./(6*A)*np.array([cx,cy])

def plot_stereo_blobs():
    print foo
    #convert STEREO to AIA ... lol . Maybe try a direct correction for occultation height?
    #or just plot height above limb, not the actual positions

#def estimate_aia_coords(stereoA_coords,):

def stereo_diffcube(mc,returnm='running'):
    base_diffmap = []
    running_diffmap = []
    for i, map_i in enumerate(mc[1:]):
        mc_rot = diffrot.diffrot_map(map_i, time=mc[0].date)
        mc_rot = diffrot.diffrot_map(mc[i+1], time=mc[i].date)
        diffdata = map_i.data - mc_rot.data
        smap_base = sunpy.map.Map(diffdata, map_i.meta)
        diffdata = mc_rot.data - map_i.data
        smap_run = sunpy.map.Map(diffdata, map_i.meta)
        #smap_base.plot_settings['cmap'] = cm.get_cmap('Greys_r')
        #smap_base.plot_settings['norm'] = colors.LogNorm(100, smap_base.max())
        #smap_run.plot_settings['cmap'] = cm.get_cmap('Greys_r')
        #smap_run.plot_settings['norm'] = colors.LogNorm(100, smap_run.max())
        base_diffmap.append(smap_base)
        running_diffmap.append(smap_run)
    if returnm == 'running':
        return sunpy.map.MapCube(running_diffmap)
    else:
        return sunpy.map.MapCube(base_diffmap)

#earlier ... none?
#better to find something in 193 though, then don't have to convert coords
#find it in the stack - have to convert it to r, theta though. shouldn't be a problem. just use height
#34:45  (964,319)
#36: (994,333)
#37: (1056,371) (appears to split in top and bottom though so not sure)
#38: (,)
#38:30 - new blob is starting to show up over the limb - or is this the same one on a return journey?

#say the bright thing in the middle is sunward-moving. Then its positions are:
#34:45  (964,319)
#36: (994,333)
#37 (1007,345)
#38 (960, 315)

#and the general top of the filament:
#lc0 (956,312)
#34:45  (964,319)
#36: (1026,355)
#37: (1056,371) (the extra bright bit)
#38: (1086,385) (assuming it got spit out)
