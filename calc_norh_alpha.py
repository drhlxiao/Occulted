from scipy.io import readsav
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pidly
from matplotlib.dates import DateFormatter
from PIL import Image
from skimage import feature, exposure
import pickle
from make_timeseries_plot_flux import mask_disk
import glob
import sunpy.map
import sunpy.cm
from matplotlib import cm
from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
from scipy import ndimage as ndi
from sunpy.visualization import axis_labels_from_ctype


def get_flat(flatdir='darksandflats',key='flp00'):
    '''Average all the flats. assume flats are after darks in sequential numbering'''
    files=glob.glob(flatdir+'/SFTH_*'+key+'*.fts')
    flatmap=sunpy.map.Map(files[0])
    return flatmap.data

def get_dark(darkdir='darksandflats',key='dlxxx'):
    '''Average all the darks'''
    files=glob.glob(darkdir+'/SFTH_*'+key+'*.fts')
    darkmap=sunpy.map.Map(files[0])
    return darkmap.data

def correct_image(fitsfile,dark,flat,saveim=False, writefits=True,show=False,fac=32.):
    fmap=sunpy.map.Map(fitsfile)
    fmap.meta['wavelnth']=fmap.meta['wavelen']
    dark=dark/fac
    flat=flat/fac
    im=fmap.data
    imc=im-dark
    flatc=flat-dark
    corrected_image=imc/flatc
    corrected_map=sunpy.map.Map(corrected_image,fmap.meta)
    sm=mask_radius(corrected_map)
    sm2=mask_radius(sm,radius=880.,greater=False)
    #other corrections
    sm2=fix_mapcoords(sm2)

    if show:
        sm2.peek()
    if saveim:
        savec.save(saveim)
    if writefits:
        print('saving to: '+fitsfile[:-4]+'corrected.fits')
        try:
            sm2.save(fitsfile[:-4]+'corrected.fits')
        except IOError:
            pass
    return sm2

def correct_images():
    cmlist,xlabels=[],[]
    l1=glob.glob('minus08/*.fts')
    l2=glob.glob('minus05/*.fts')
    l3=glob.glob('center/*.fts')
    l4=glob.glob('plus05/*.fts')
    l5=glob.glob('plus08/*.fts')
    l6=glob.glob('plus35/*.fts')

    d1=get_dark()
    d2=get_dark(key='dsxxx')

    f1=get_flat(key='fsm08')
    f2=get_flat(key='fsm05')
    f3=get_flat(key='flp00')
    f4=get_flat(key='fsp05')
    f5=get_flat(key='fsp08')
    f6=get_flat(key='fsp35')

    #ln are filename lists
    for i,a,b,c,d,e,f in zip(range(0,len(l1)),l1,l2,l3,l4,l5,l6):
        if 'corrected' not in a:
            m1=correct_image(a,d2,f1)
        if 'corrected' not in b:
            m2=correct_image(b,d2,f2)
        if 'corrected' not in c:
            m3=correct_image(c,d1,f3)
        if 'corrected' not in d:
            m4=correct_image(d,d2,f4)
        if 'corrected' not in e:
            m5=correct_image(e,d2,f5)
        if 'corrected' not in f:
            m6=correct_image(f,d2,f6)

def fix_mapcoords(m):
    m.meta['ctype1']='HPLN-TAN'
    m.meta['ctype2']='HPLT-TAN'
    m.meta['cunit1']='arcsec'
    #m.meta['wavelnth']=m.meta['wavelen']
    m.meta['cunit2']='arcsec'
    if m.meta['date-obs'].startswith(' '):
        m.meta['date-obs']='2013-05-01T'+ m.meta['date-obs'][-12:-4]
    return m

def mask_radius(submap,plot=False,radius=1100.,greater=True):
    x, y = np.meshgrid(*[np.arange(v.value) for v in submap.dimensions])
    xoff,yoff=submap.dimensions[0]/2.,submap.dimensions[1]/2.
    r = np.sqrt((x-xoff.value) ** 2 + (y-yoff.value) ** 2) / (radius)
    if greater:
        mask = ma.masked_greater_equal(r, 1)
    else:
        mask = ma.masked_less_equal(r, 1)
    palette = submap.plot_settings['cmap']
    palette.set_bad('black')
    #scaled_map = sunpy.map.Map(submap.data, submap.meta, mask=mask.mask)
    mdata=ma.masked_array(submap.data,mask=mask.mask)
    scaled_map = sunpy.map.Map(ma.filled(mdata,0), submap.meta)
    if plot:
        fig = plt.figure()
        fig.add_subplot(projection=scaled_map)
        scaled_map.plot(cmap=palette)
        scaled_map.draw_limb()
        fig.show()
    return scaled_map

def prep_ha(box=[(-1200,0),(-900,450)],aia=True):
    cmlist,xlabels=[],[]
    l1=glob.glob('minus08/*rot.fits')
    l2=glob.glob('minus05/*rot.fits')
    l3=glob.glob('center/*rot.fits')
    l4=glob.glob('plus05/*rot.fits')
    l5=glob.glob('plus08/*rot.fits')
    l6=glob.glob('plus35/*rot.fits')

    #d1=get_dark()
    #d2=get_dark(key='dsxxx')

    #f1=get_flat(key='fsm08')
    #f2=get_flat(key='fsm05')
    #f3=get_flat(key='flp00')
    #f4=get_flat(key='fsp05')
    #f5=get_flat(key='fsp08')
    #if aia:
    l7=glob.glob('../AIA/304/mtk*.fts')

    #ln are filename lists
    for i,a,b,c,d,e,f in zip(range(0,len(l1)),l1,l2,l3,l4,l5,l6):
        #m1=correct_image(a,d2,f1)
        #m2=correct_image(b,d2,f2)
        #m3=correct_image(c,d1,f3)
        #m4=correct_image(d,d2,f4)
        #m5=correct_image(e,d2,f5)
        maia0=sunpy.map.Map(l7[i])
        bl=SkyCoord(box[0][0]*u.arcsec,box[0][1]*u.arcsec,frame=maia0.coordinate_frame)
        ur=SkyCoord(box[1][0]*u.arcsec,box[1][1]*u.arcsec,frame=maia0.coordinate_frame)
        maia=mask_disk(maia0,fac=10.,filled=True).submap(bl,ur)
        maia.plot_settings['cmap']=cm.binary
        #print maia.bottom_left_coord, maia.top_right_coord
        m1a=mask_disk(fix_mapcoords(sunpy.map.Map(a)).rotate(),fac=10.,filled=True)
        m1=m1a.submap(bl,ur)
        m2a=mask_disk(fix_mapcoords(sunpy.map.Map(b)).rotate(),fac=10.,filled=True)
        m2=m2a.submap(bl,ur)
        m3a=mask_disk(fix_mapcoords(sunpy.map.Map(c)).rotate(),fac=10.,filled=True)
        m3=m3a.submap(bl,ur)
        m4a=mask_disk(fix_mapcoords(sunpy.map.Map(d)).rotate(),fac=10.,filled=True)
        m4=m4a.submap(bl,ur)
        m5a=mask_disk(fix_mapcoords(sunpy.map.Map(e)).rotate(),fac=10.,filled=True)
        m5=m5a.submap(bl,ur)
        m6a=mask_disk(fix_mapcoords(sunpy.map.Map(f)).rotate(),fac=10.,filled=True)
        m6=m6a.submap(bl,ur)
        #maia.peek()
        #m1.peek()
        #make composite maps or submaps
        if aia:
            cm1=sunpy.map.Map(maia,m1,composite=True)
            #cm1=cm1a.submap(maia.bottom_left_coord,maia.top_right_coord)
            cm2=sunpy.map.Map(maia,m2,composite=True)
            cm3=sunpy.map.Map(maia,m3,composite=True)
            #cm3=cm3a.submap(maia.bottom_left_coord,maia.top_right_coord)
            cm4=sunpy.map.Map(maia,m4,composite=True)
            cm5=sunpy.map.Map(maia,m5,composite=True)
            cm6=sunpy.map.Map(maia,m6,composite=True)
            for comp_map in [cm1,cm2,cm3,cm4,cm5,cm6]:
                comp_map.set_levels(0,[9,10,25,26],percent=True)
                comp_map.set_alpha(0,1)
                #aia_set=comp_map.get_plot_settings(1)
                #aia_set['cmap']=cm.binary
                #comp_map.set_plot_settings(0,aia_set)
                mtk_set=comp_map.get_plot_settings(1)
                mtk_set['cmap']=cm.RdBu#sunpy.cm.cmlist['sdoaia304']
                comp_map.set_plot_settings(1,mtk_set)
                comp_map.set_alpha(1,.95)

            cmlist.append([cm1,cm2,cm3,cm4,cm5,cm6])
        else:
            cmlist.append([maia,m1,m2,m3,m4,m5,m6])


        #fig,ax=plt.subplots(1,5,figsize=[15,5])
        ##ax[0].imshow(m1.data[1400:1830,150:600],cmap=sunpy.cm.cmlist['sdoaia304'],origin='bottom left')
        ##ax[0].set_title('Halpha -0.8 A')
        #plt.sca(ax[0])
        #cm1.plot(title='Halpha -0.8 A')
        #ax[0].set_xlabel(a[-24:-18])
        #ax[1].imshow(m2.data[1400:1830,150:600],cmap=sunpy.cm.cmlist['sdoaia304'],origin='bottom left')
        #plt.sca(ax[1])
        #cm2.plot(title='Halpha -0.5 A')
        #ax[1].set_xlabel(b[-24:-18])
        ##ax[2].imshow(m3.data[1400:1830,150:600],cmap=sunpy.cm.cmlist['sdoaia304'],origin='bottom left')
        #plt.sca(ax[2])
        #cm3.plot(title='Halpha center')
        #ax[2].set_xlabel(c[-24:-18])
        ##ax[3].imshow(m4.data[1400:1830,150:600],cmap=sunpy.cm.cmlist['sdoaia304'],origin='bottom left')
        #plt.sca(ax[3])
        #cm4.plot(title='Halpha +0.5 A')
        #ax[3].set_xlabel(d[-24:-18])
        ##ax[4].imshow(m5.data[1400:1830,150:600],cmap=sunpy.cm.cmlist['sdoaia304'],origin='bottom left')
        #plt.sca(ax[4])
        #cm5.plot(title='Halpha +0.8 A')
        #ax[4].set_xlabel(e[-24:-18])
        #for aa in ax:
        #    #aa.xaxis.set_visible(False)
        #    aa.yaxis.set_visible(False)
        #fig.show()
        #fig.savefig()
        xlabels.append([a[-24:-18],b[-24:-18],c[-24:-18],d[-24:-18],e[-24:-18],f[-24:-18]])
    return cmlist,xlabels

def plot_contours(cmlist, xlabels,levels=[[19,20,49,50],[19,20,49,50],[39,40,49,50],[39,40,49,50],[39,40,49,50],[39,40,49,50],[19,20,39,40]],smooth=3,cmap_check=False, linewidth=2.,single_map=False):
    '''stack all the contours, color code as red- or blue-shifted, overlaid on AIA 304'''
    #ckeys=[3,4,5,7,8,9]
    colormaps=[cm.RdBu]
    mkeys=['RdBu','-.8','-.5','0.0','+.5','+.8','+3.5']
    for j,c in enumerate(cmlist):
        if smooth: #gaussian smooth the Halpha maps
            newc=[sunpy.map.Map(ndi.gaussian_filter(cc.data,smooth),cc.meta) for cc in c]
            cmap=sunpy.map.Map(newc,composite=True)
        else:
            cmap=sunpy.map.Map(c,composite=True)
        #do all the stuff with the plot settings
        csettings=cmap.get_plot_settings(0)
        #csettings['cmap']= sunpy.cm.cmlist['sdoaia304']
        new_cdata=cm.revcmap(sunpy.cm.cmlist['sdoaia171']._segmentdata)
        new_cmap=LinearSegmentedColormap('sdoaia171_r',new_cdata)
        csettings['cmap']=new_cmap

        cmap.set_plot_settings(0,csettings)
        cmap.set_alpha(0,.75)
        if single_map:
            r=single_map
        else:
            r=[1,7]
        for i in range(r[0],r[1]):
            cmap.set_levels(i,levels[i],percent=True)
            cmap.set_alpha(i,1)
            csettings=cmap.get_plot_settings(i)
            #new_cmap=single_cmap(ckeys[i-1])
            new_cmap=single_cmap(i-1)
            #print(new_cmap._segmentdata)
            csettings['cmap']= new_cmap #cm.RdBu#sunpy.cm.cmlist['sdoaia304']
            csettings['linewidths']=linewidth
            cmap.set_plot_settings(i,csettings)
            colormaps.append(new_cmap)
        fig,ax=plt.subplots()
        #cmap.plot(axes=ax,title=c[1].meta['date-obs'])
        foo=cmap_plot(cmap,axes=ax,title=c[1].meta['date-obs'],labels=['-.8 $\AA$','-.5 $\AA$','6535 $\AA$','+.5 $\AA$','+.8 $\AA$','+3.5 $\AA$'],legend=True)
        #ax.legend(loc="upper left", labels=['red','blue'],color=['r','b'])
        #ax.xaxis.set_visible(False)
        #ax.yaxis.set_visible(False)
        fig.show()

        if cmap_check: #plot the contour color tables
            lcols=[]
            fig,ax=plt.subplots()
            x = np.linspace(0.0, 256., 256.)
            for j in range(0,7):
                ax.scatter(x,x+20.*j, c=x, s=100, cmap=colormaps[j],linewidths=0.0,label=mkeys[j])
                #lcols.append(colormaps[j]._segmentdata[0])
            ax.legend(loc='lower right')
            #ax.set_xlim([0,20])
            #ax.set_ylim([0,20])
            fig.show()
        #return cmap,colormaps

def plot_fig4(levels=[[50],[50],[50],[50],[50],[50],[50]],smooth=3, linewidth=1.,single_map=False):
    '''stack all the contours, color code as red- or blue-shifted, overlaid on AIA 304'''
    #ckeys=[3,4,5,7,8,9]
    cmlist=pickle.load(open('compound_maps.p','rb'))
    cmlist=cmlist[2]
    colormaps=[cm.RdBu]
    mkeys=['RdBu','-.8','-.5','0.0','+.5','+.8','+3.5']
    newc=[pickle.load(open('Ha_peak304submap.p','rb'))]#[cmlist[-1]._maps[0]]
    if smooth: #gaussian smooth the Halpha maps
        for cc in cmlist:
            newc.append(sunpy.map.Map(ndi.gaussian_filter(cc._maps[1].data,smooth),cc._maps[1].meta).submap(SkyCoord([-1050,-900]*u.arcsec,[150,300]*u.arcsec,frame=newc[0].coordinate_frame)))
        print( newc[0].meta['date-obs'],newc[1].meta['date-obs'])
    else:
        for c in cmlist:
            newc.append(c._maps[1])
    cmap=sunpy.map.Map(newc,composite=True)
    #do all the stuff with the plot settings
    csettings=cmap.get_plot_settings(0)
    #csettings['cmap']= sunpy.cm.cmlist['sdoaia304']
    new_cdata=cm.revcmap(sunpy.cm.cmlist['sdoaia304']._segmentdata)
    new_cmap=LinearSegmentedColormap('sdoaia304_r',new_cdata)
    newnorm=colors.Normalize(vmin=.1*np.min(cmap._maps[0].data),vmax=.9*np.max(cmap._maps[0].data))
    csettings['norm']=newnorm
    csettings['cmap']=new_cmap

    cmap.set_plot_settings(0,csettings)
    #cmap.set_alpha(0,.75)

    print( len(cmap._maps))
    for i in range(1,7):
        cmap.set_levels(i,levels[i],percent=True)
        cmap.set_alpha(i,1)
        csettings=cmap.get_plot_settings(i)
        ##new_cmap=single_cmap(ckeys[i-1])
        new_cmap=single_cmap(i-1)
        ##print(new_cmap._segmentdata)
        csettings['cmap']= new_cmap #cm.RdBu#sunpy.cm.cmlist['sdoaia304']
        csettings['linewidths']=linewidth
        cmap.set_plot_settings(i,csettings)
        #colormaps.append(new_cmap)
    fig,ax=plt.subplots()
    #cmap.plot(axes=ax,title=c[1].meta['date-obs'])
    foo=cmap_plot(cmap,axes=ax,labels=['-.8 $\AA$','-.5 $\AA$','6535 $\AA$','+.5 $\AA$','+.8 $\AA$','+3.5 $\AA$'],legend=True)
    #cmap.plot(axes=ax)
      #ax.legend(loc="upper left", labels=['red','blue'],color=['r','b'])
      #ax.xaxis.set_visible(False)
      #ax.yaxis.set_visible(False)
    fig.show()
    return cmap

def cmap_plot(cmap,axes=None, annotate=True,title='',labels='',legend=False,**matplot_args):
    '''mod of sunpy.compositemap.plot() to allow line thickness modifications in contours'''
    if not axes:
        axes = plt.gca()

    if annotate:
        axes.set_xlabel(axis_labels_from_ctype(cmap._maps[0].coordinate_system[0],
                                                   cmap._maps[0].spatial_units[0]))
        axes.set_ylabel(axis_labels_from_ctype(cmap._maps[0].coordinate_system[1],
                                               cmap._maps[0].spatial_units[1]))
        axes.set_title(title)

    # Define a list of plotted objects
    ret = []
    # Plot layers of composite map
    for i,m in enumerate(cmap._maps):
        # Parameters for plotting
        bl = m._get_lon_lat(m.bottom_left_coord)
        tr = m._get_lon_lat(m.top_right_coord)
        x_range = list(u.Quantity([bl[0], tr[0]]).to(m.spatial_units[0]).value)
        y_range = list(u.Quantity([bl[1], tr[1]]).to(m.spatial_units[1]).value)
        params = {
                "origin": "lower",
                "extent": x_range + y_range,
                "cmap": m.plot_settings['cmap'],
                "norm": m.plot_settings['norm'],
                "alpha": m.alpha,
                "zorder": m.zorder,
            }
        params.update(matplot_args)

        # The request to show a map layer rendered as a contour is indicated by a
        # non False levels property.  If levels is False, then the layer is
        # rendered using imshow.
        if m.levels is False:
            # Check for the presence of masked map data
            if m.mask is None:
                ret.append(axes.imshow(m.data, **params))
            else:
                ret.append(axes.imshow(np.ma.array(np.asarray(m.data), mask=m.mask), **params))
        else:
            try:
                aa=m.plot_settings['linewidths']
                params = {
                "origin": "lower",
                "extent": x_range + y_range,
                "cmap": m.plot_settings['cmap'],
                "norm": m.plot_settings['norm'],
                "linewidths":m.plot_settings['linewidths'],
                "alpha": m.alpha,
                "zorder": m.zorder#,"label":labels[i-1]
                          }

            except KeyError:
                params = {
                "origin": "lower",
                "extent": x_range + y_range,
                "cmap": m.plot_settings['cmap'],
                "norm": m.plot_settings['norm'],
                #"linewidths":m.plot_settings['linewidths'],
                "alpha": m.alpha,
                "zorder": m.zorder#,"label":labels[i-1]
                          }
            params.update(matplot_args)

            # Check for the presence of masked map data
            if m.mask is None:
                ret.append(axes.contour(m.data, m.levels, **params))
            else:
                ret.append(axes.contour(np.ma.array(np.asarray(m.data), mask=m.mask), m.levels, **params))

            # Set the label of the first line so a legend can be created
            ret[-1].collections[0].set_label(m.name)
    # Adjust axes extents to include all data
    axes.axis('image')
    if legend:
        proxy = [plt.Rectangle((0,0),1,1,fc = (float(j)/6.,0,1.-float(j)/6.)) for j in range(0,6)]
        axes.legend(proxy, labels,loc="upper left")

    # Set current image (makes colorbar work)
    #plt.sci(ret[0])
    return ret

def single_cmap(n):
    '''n is scale on red-blue cmap'''
    #rb=cm.RdBu
    #rval=rb._segmentdata['red'][n]
    #gval=rb._segmentdata['green'][n]#[n]
    #bval=rb._segmentdata['blue'][0]#[n]
    #colors=[rval,gval,bval]
    #nbins=[3,6,10,11]
    #single_cmap=LinearSegmentedColormap.from_list(str(n),colors[n],N=1)
    single_cmap=ListedColormap([(float(n)/6.,0.,1.-float(n)/6.)])
    return single_cmap

def just_plot(cmlist,levels,xlabels):
    newc=[]
    for c in cmlist:
        aiamap=c.get_map(0)
        #smooth the aia
        map1=sunpy.map.Map(ndi.gaussian_filter(aiamap.data,3),aiamap.meta)
        map2=c.get_map(1)
        newc.append(sunpy.map.Map(map2,map1,composite=True))

    cm1,cm2,cm3,cm4,cm5,cm6=newc
    for comp_map in [cm1,cm2,cm3,cm4,cm5,cm6]:
        comp_map.set_levels(1,levels,percent=True)
        aia_set=comp_map.get_plot_settings(0)
        aia_set['cmap']=cm.Greys
        aia_set['linewidths']=2.
        comp_map.set_plot_settings(1,aia_set)
        mtk_set=comp_map.get_plot_settings(0)
        mtk_set['cmap']=cm.RdBu#sunpy.cm.cmlist['sdoaia304']
        mtk_set['norm']=mpl.colors.Normalize(vmin=np.min(comp_map.get_map(0).data),vmax=np.max(comp_map.get_map(0).data))
        comp_map.set_plot_settings(0,mtk_set)
        comp_map.set_alpha(1,1)

    fig,ax=plt.subplots(1,6,figsize=[18,5])
    #ax[0].imshow(m1.data[1400:1830,150:600],cmap=sunpy.cm.cmlist['sdoaia304'],origin='bottom left')
    #ax[0].set_title('Halpha -0.8 A')
    plt.sca(ax[0])
    cmap_plot(cm1,title='H$\\alpha$ -0.8 $\AA$')
    ax[0].set_xlabel(xlabels[0])
    #ax[1].imshow(m2.data[1400:1830,150:600],cmap=sunpy.cm.cmlist['sdoaia304'],origin='bottom left')
    plt.sca(ax[1])
    cmap_plot(cm2,title='H$\\alpha$ -0.5 $\AA$')
    ax[1].set_xlabel(xlabels[1])
    #ax[2].imshow(m3.data[1400:1830,150:600],cmap=sunpy.cm.cmlist['sdoaia304'],origin='bottom left')
    plt.sca(ax[2])
    cmap_plot(cm3,title='H$\\alpha$ 6535 $\AA$')
    ax[2].set_xlabel(xlabels[2])
    #ax[3].imshow(m4.data[1400:1830,150:600],cmap=sunpy.cm.cmlist['sdoaia304'],origin='bottom left')
    plt.sca(ax[3])
    cmap_plot(cm4,title='H$\\alpha$ +0.5 $\AA$')
    ax[3].set_xlabel(xlabels[3])
    #ax[4].imshow(m5.data[1400:1830,150:600],cmap=sunpy.cm.cmlist['sdoaia304'],origin='bottom left')
    plt.sca(ax[4])
    cmap_plot(cm5,title='H$\\alpha$ +0.8 $\AA$')
    ax[4].set_xlabel(xlabels[4])
    plt.sca(ax[5])
    cmap_plot(cm6,title='H$\\alpha$ +3.5 $\AA$')
    ax[5].set_xlabel(xlabels[5])
    for aa in ax:
        #aa.xaxis.set_visible(False)
        aa.yaxis.set_visible(False)
        aa.tick_params(labelbottom='off')
    fig.show()

def convert_to_sfu(nmap,freq=17):
    c=1222
    bt=np.array(nmap.data) #map in brightness temperature units
    intensity=bt*(freq**2 * theta_maj * theta_min)/c
    newmap=sunpy.map.Map(intensity, nmap.meta) #should I put on a mask to screen out negative values?
    return newmap

def calc_alpha(ret_mask=True):
    adata=readsav('norh_convol.sav',python_dict=True)
    f34=adata['f34_c']
    f17=adata['f17_c']

    #first get the subarrays of the region of interest
    f34c=[f[80:185,0:125] for f in f34]
    f17c=[f[80:185,0:125] for f in f17]
    masklist=[]
    #determine thresholds
    f34m,f17m=np.zeros([174,105,125]),np.zeros([174,105,125])
    for k,i,j in zip(range(0,174),f34c,f17c):
        imin=np.min(i)
        imean=np.mean(i)
        idev=np.std(i)
        #if idev < np.abs(imean):
        #    if imean < 0:
        #        ithresh=(imean-imin)/2.
        #    else:
        #        ithresh=imin+idev
        #else:
        #ithresh=(imean-idev)
        #if k in range(0,10):
        #    print(ithresh,imin,np.max(i),imean,idev)
        f34a=np.ma.masked_where(i,i !=0)

        #jmin=np.min(j)
        #jmean=np.mean(j)
        #jdev=np.std(j)
        #if jdev < np.abs(jmean):
        #    jthresh=jmin+jdev
        #else:
        #    jthresh=(jmean-jmin)/2.
        #f17a=np.ma.masked_less(j,jthresh)

        #find out which pixels above threshold are common to both arrays
        #do this by multiplying the masks
        m34=np.ma.getmask(f34a)
        #m17=np.ma.getmask(f17a)
        #newmask=np.multiply(m34,m17) #but this needs to be elementwise not array multiplication
        #if k == 110:
        #    fig=plt.figure(figsize=[16,5])
        #    ax1=plt.subplot(141)
        #    ax2=plt.subplot(142)
        #    ax3=plt.subplot(143)
        #    ax4=plt.subplot(144)
        #    ax1.imshow(m34)
        #    ax2.imshow(m17)
        #    ax3.imshow(newmask)
        #    ax4.imshow(i,alpha=.5)
        #    ax4.contour(i,levels=np.array([1e6,1e5,1e4,1e3,100,10])*np.min(i))
        #    fig.show()
        #nullmask=np.ones([256,256])
        #nullmask[80:185,0:125]=newmask
        masklist.append(m34)

        f34b=np.ma.masked_array(i,mask=m34)
        f17b=np.ma.masked_array(j,mask=m34)
        f34m[k]=f34b.filled(0.)
        f17m[k]=f17b.filled(0.)
    if ret_mask:
        return masklist

    #start idl and run the calculation given these inputs
    idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
    idl("restore, 'norh_convol.sav'")
    idl('f17_c',f17m)
    idl('f34_c',f34m)
    #idl('index17_c', adata['index17_c'])
    #idl('index34_c', adata['index34_c'])
    idl('norh_alpha,index17_c,f17_c,index34_c,f34_c,indexal,alpha,mvda')
    #indexal=idl.indexal
    alpha=idl.alpha
    mvda=idl.mvda
    idl.close()
    print(np.shape(alpha))

    #get mean and plot
    mean_alpha=np.mean(alpha,axis=1)
    fig,ax=plt.subplots()
    ax.plot(range(0,len(alpha)),mean_alpha)
    fig.show()

    #maybe show some visualizations too

    return alpha

def plot_alpha(with_lc=True):
    dates=pickle.load(open('lcdates.p','rb'))
    alpha=pickle.load(open('alpha.p','rb'))
    m1=np.mean(alpha,axis=1)
    mean_alpha=np.mean(m1,axis=1)
    if with_lc: #plot these on the right hand axis
        try:
            lcdata=pickle.load(open('lcdata.p','rb'))
            lc17=lcdata['lc17']
            lc34=lcdata['lc34']
        except IOError:
            norh17file='130501_17g10s/sfu17GHz.sav'
            norh34file='130501_34g10s/sfu34GHz.sav'
            lc17=readsav(norh17file)['fi2']
            lc34=readsav(norh34file)['fi2']

    fig,ax=plt.subplots()
    ax.plot(dates[1],mean_alpha)
    ax.set_xlim([dates[1][0],dates[1][-1]])
    if with_lc:
        dax=ax.twinx()
        dax.plot(dates[1],lc34,'b',label='$\\alpha$')
        dax.plot(dates[0],lc17,'r',label='17 GHz')
        dax.plot(dates[1],lc34,'g',label='34 GHz')
        dax.set_ylabel('Flux (sfu/steradian)')
        dax.legend(loc='upper left')
        dax.set_xlim([dates[1][0],dates[1][-1]])
    ax.set_ylabel('Mean Spectral Index $\\alpha$')
    ax.set_xlabel('Time on 01 May 2013')
    myFmt = DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    plt.gcf().autofmt_xdate()
    if with_lc:
        fig.savefig('alpha_and_lc.png')
    else:
        fig.savefig('alpha.png')
    fig.show()

