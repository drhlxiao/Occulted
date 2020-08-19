def fig1_stereoB():
    peakf = pickle.load(open('20130501_023030_14euB.p','rb'))
    smap=get_submap(peakf,[250,95],[600,355])
    palette_name=smap.plot_settings['cmap']#sunpy.cm.cmlist[cmap]#m17.plot_settings['cmap']
    new_cdata=cm.revcmap(palette_name._segmentdata)
    new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name+'_r',new_cdata)
    newsettings=smap.plot_settings
    newsettings['cmap']=new_cmap
    newsettings['norm']=colors.Normalize(20,2000)

    fig,ax=plt.subplots()
    ax.set_aspect('equal')
    smap.plot()
    ax.set_xticks([250,300,350,400,450,500,550,600])
    ax.set_yticks([100,150,200,250,300,350])
    ax.set_xlabel("X (arcsec)")
    ax.set_ylabel("Y (arcsec)")
    ax.set_xlim([250,600])
    ax.set_ylim([100,350])
    fig.show()

def fig1_stereoA():
    peakf = pickle.load(open('20130501_023715_15euA.p','rb'))
    smap=get_submap(peakf,[850,195],[1200,455])
    palette_name=smap.plot_settings['cmap']#sunpy.cm.cmlist[cmap]#m17.plot_settings['cmap']
    new_cdata=cm.revcmap(palette_name._segmentdata)
    new_cmap=matplotlib.colors.LinearSegmentedColormap(palette_name.name+'_r',new_cdata)
    newsettings=smap.plot_settings
    newsettings['cmap']=new_cmap
    newsettings['norm']=colors.Normalize(30,3000)

    fig,ax=plt.subplots()
    #ax.set_aspect('equal')
    smap.plot()
    ax.set_xticks([850,900,950,1000,1050,1100,1150,1200])
    ax.set_yticks([200,250,300,350,400,450])
    ax.set_xlabel("X (arcsec)")
    ax.set_ylabel("Y (arcsec)")
    ax.set_xlim([850,1200])
    ax.set_ylim([200,450])
    fig.show()

def both_hessi2(levels=[49,50,69,70,89,90],save=False,legend=True):
    comp_maps=glob.glob('comp_map_final*.p')
    newhxrmaps=['nine_hi_clean1.fits','nine_hi_clean1.fits','nine_hi_clean2.fits','nine_hi_clean2.fits']
    marr,dmaps,centroids=[],[],[]
    for i,c in enumerate(comp_maps):
        imname='comp_map'+str(i+2)+'.png'
        m=pickle.load(open(c,'rb'))
        #m.remove_map(1)
        if i < 4:
            #m.remove_map(1)
            hxrmap=fix_fits_map(newhxrmaps[i])
            _,cs,_=find_centroid_from_map(hxrmap)
            centroids.append(cs)
            #m.add_map(hxrmap)

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

        if i == 0:
            levels=[52,53,69,70,89,90]


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


def norh_fig1(levels=[5.45,5.5,6.9,7,9.9,10]):
    from calc_norh_alpha import cmap_plot
    cmap=pickle.load(open('norh_304_fig1b.p','rb'))
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
    cmap.set_levels(1,levels,percent=True)
    cmap.plot(axes=ax)
    ax.set_ylim([0.02778,0.097223]) #100 - 350 arcsec
    ax.set_xlim([-.3055,-.236]) #-1000 - -800 arcsec? Want -1100 to -850
    ax.set_xticks([-0.2917,-.277, -0.2639,-.25,-0.2362])
    ax.set_yticks([0.0278,0.041667,0.0556,0.069445,0.08334,0.097223])
    xlab=ax.get_xticks()
    ylab=ax.get_yticks()
    newxlab=[str(int(xl*3600.)) for xl in xlab]
    newylab=[str(int(xl*3600.)) for xl in ylab]
    ax.set_xticklabels(newxlab)
    ax.set_yticklabels(newylab)
    ax.set_xlabel('X (arcsec)')
    ax.set_ylabel('Y (arcsec)')
    ax.set_aspect('equal')
    ax.set_title('AIA 304$\AA$ 2013-05-01 02:37:15')
    labels=['NORH 17 GHz']
    proxy = [plt.Rectangle((0,0),.5,.5,facecolor='none')] #how to make it green or blue?
    ax.legend(proxy, labels,loc="lower right")
    fig.show()


