def fig7():
    nitta,goes=pickle.load(open('../flare_lists/nitta_fig7.p','rb'))
    nitta=np.array(nitta)
    goes=np.array(goes)
    nzeros=np.where(nitta != 0.0)[0]
    nnonzero=nitta[nzeros]
    gnonzero=goes[nzeros]
    nittaflux=np.array(nnonzero)*1.39e-11
    #error bars
    ylow,yhigh=[],[]
    for n,nf in zip(nnonzero,nittaflux):
        if n > 4e6:
            ylow.append(np.abs(nf*.3-nf))
            yhigh.append(nf*3.)
        else:
            ylow.append(np.abs(nf*.5-nf))
            yhigh.append(nf*1.5)
    fig,ax=plt.subplots()
    ax.errorbar(nittaflux,gnonzero,xerr=[ylow,yhigh],fmt='o')
    #ax.scatter(nittaflux,ylow,c='m')
    #ax.scatter(nittaflux,yhigh,c='g')
    ax.plot(np.linspace(1e-8,1e-3,10),np.linspace(1e-8,1e-3,10))
    ax.set_xlim([1e-8,1e-3])
    ax.set_ylim([1e-8,1e-3])
    ax.set_xlabel('Nitta GOES EUVI proxy (W/m^2)')
    ax.set_ylabel('GOES 1-8 Flux (W/m^2)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.show()
    return nittaflux,ylow,yhigh

def fix_pt():
    for f in fl.list:
        if f.Files.prepped['stereo-peak'] !='':
            f.Datetimes.peakdiff=np.abs((f.Datetimes.Messenger_peak - f.Datetimes.stereo_peak).seconds)

def fix_pt2():
    for f in fl.list:
        if f.Files.prepped['stereo-peak'] !='':
            f.Datetimes.stereo_peak=dt.strptime(f.Files.prepped['stereo-peak'][:-10],'%Y%m%d_%H%M%S')

def Mproxy():
    for f in fl.list:
        if f.Properties.Messenger_GOES_long != '':
            f.Properties.Messenger_GOES_selfproxy=(1.13081*10**-8)*float(f.Properties.Messenger_GOES_long)**.96135

def fig7_all(fl,sat='M',fitline=False,tol=500.):
    nitta=fl.isolate_att('Nittaflux')
    nobs,gobs=[],[]
    if sat == 'G':
        goes=fl.isolate_att('GOES_GOES_long')
        ylabel='GOES 1-8 A Flux ($\\rm{Wm^{-2}}$)'
    if sat == 'M':
        goes=fl.isolate_att('Messenger_GOES_selfproxy')
        ylabel='Messenger Flux ($\\rm{Wm^{-2}}$)'

    #if from_T_EM:
    #    import messenger_goes_class as mgc
    #    #if from_T_EM == '1':
    #    tt=fl.isolate_att('T1')
    #    em=fl.isolate_att('EM1')
    #    #else:
    #    #    tt=self.isolate_att('T2')
    #    #    em=self.isolate_att('EM2')

    #    mcl,mcs=mgc.messenger_goes_class(tt,em)
    #    goes=mcl

    for f in fl.list:
        if sat == 'G':
            try:
                #if f.Properties.GSobs_with_sloc == True and f.Datetimes.peakdiff <=tol:
                #if f.Datetimes.peakdiff <=tol:
                nobs.append(float(f.Properties.Nittaflux))
                #if not from_T_EM:
                gobs.append(float(f.Properties.GOES_GOES_long))
            except AttributeError:
                pass

        if sat == 'M':
            #if f.Observation.MSvis==True:
            try:
                if f.Properties.MSobs_with_sloc == True and f.Datetimes.peakdiff <=tol:
                #if f.Datetimes.peakdiff <=tol:
                    nobs.append(float(f.Properties.Nittaflux))
                    #if not from_T_EM:
                    gobs.append(float(f.Properties.Messenger_GOES_proxy))
            except AttributeError:
                pass

    nitta=np.array(nitta)
    goes=np.array(goes)
    nzeros=np.where(nitta != 0.0)[0]
    nnonzero=nitta[nzeros]
    gnonzero=goes[nzeros]
    if type(gnonzero[0])==str:
        gnonzero=[float(g) for g in gnonzero]
    nittaflux=np.array(nnonzero)*1.39e-11
    nfobs=np.array(nobs)*1.39e-11
    #error bars
    ylow,yhigh=[],[]
    for n,nf in zip(nobs,nfobs):
        if n > 4e6:
            ylow.append(np.abs(nf*.3-nf))
            yhigh.append(nf*3.)
        else:
            ylow.append(np.abs(nf*.5-nf))
            yhigh.append(nf*1.5)
    fig,ax=plt.subplots()
    #ax.scatter(nittaflux,gnonzero,color='gray')
    ax.scatter(nfobs,gobs)
    #ax.errorbar(nfobs,gobs,xerr=[ylow,yhigh],markersize=10.5,fmt='*',c='r')
    #ax.scatter(nittaflux,ylow,c='m')
    #ax.scatter(nittaflux,yhigh,c='g')
    xvec=np.linspace(1e-8,1e-3,10)
    if fitline:
        yvec=10**(fitline[0]*np.log10(xvec)+fitline[1])
        ax.plot(xvec,yvec,'r')
    ax.plot(xvec,xvec,'--k')
    ax.set_xlim([1e-8,1e-3])
    ax.set_ylim([1e-8,1e-3])
    ax.set_xlabel('Nitta EUVI GOES proxy ($\\rm{Wm^{-2}}$)')
    ax.set_ylabel(ylabel)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_aspect('equal')
    fig.show()
    return nfobs,gobs

def overlay(set1,set2,fitline=False,title=''):
    fig,ax=plt.subplots()
    #ax.scatter(nittaflux,gnonzero,color='gray')
    ax.scatter(set1[0],set1[1],color='gray')
    #ax.errorbar(nfobs,gobs,xerr=[ylow,yhigh],markersize=10.5,fmt='*',c='r')
    ax.scatter(set2[0],set2[1],c='r',marker='v',s=60)
    #ax.scatter(nittaflux,yhigh,c='g')
    xvec=np.linspace(1e-8,1e-3,10)
    if fitline:
        yvec=10**(fitline[0]*np.log10(xvec)+fitline[1])
        ax.plot(xvec,yvec,'b')
    ax.plot(xvec,xvec,'--k')
    ax.set_xlim([1e-8,1e-3])
    ax.set_ylim([1e-8,1e-3])
    ax.set_xlabel('GOES Flux ($\\rm{Wm^{-2}}$)')
    ax.set_ylabel('GOES Flux Calculated from Messenger T, EM ($\\rm{Wm^{-2}}$)')
    ax.set_title(title)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_aspect('equal')
    fig.show()
    #return nfobs,gobs

def overlay2(set1,set2,fitline=False,title='',yran=[1e-8,1e-3],ylab='GOES Flux ($\\rm{Wm^{-2}}$)',xran=[1e4,1e9],xlab='Nitta Flux (photons/s)',dashed=True,nfline=False,save=False):
    fig,ax=plt.subplots()
    #ax.scatter(nittaflux,gnonzero,color='gray')
    ax.scatter(set1[0],set1[1],color='gray')
    #ax.errorbar(nfobs,gobs,xerr=[ylow,yhigh],markersize=10.5,fmt='*',c='r')
    ax.scatter(set2[0],set2[1],c='r',marker='v',s=60)
    #ax.scatter(nittaflux,yhigh,c='g')
    xvec=np.linspace(xran[0],xran[1],10)

    #xvec=np.linspace(1e-7,1e-2,10)
    if fitline:
        #yvec=1.39*10**-11*xvec
        yvec=10**(fitline[0]*np.log10(xvec)+fitline[1])
        ax.plot(xvec,yvec,'b')
    if nfline:
        ax.plot(xvec,1.39*10**-11*xvec,'--k')
    if dashed:
        ax.plot(xvec,xvec,'--k')
    ax.set_xlim(xran)
    #ax.set_xlim([1e-7,1e-2])
    ax.set_ylim(yran)
    #ax.set_xlabel('GOES Flux ($\\rm{Wm^{-2}}$)')
    #ax.set_xlabel('Nitta Proxy ($\\rm{Wm^{-2}}$)')
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_aspect('equal')
    fig.show()
    if save:
        plt.savefig(save)
    #return nfobs,gobs

def overlayN(set1,setN,labels=False,fitline=False,title='',yran=[1e-8,1e-3],ylab='GOES Flux ($\\rm{Wm^{-2}}$)',xran=[1e4,1e9],xlab='Nitta Flux (photons/s)',dashed=True,nfline=False,save=False):
    if not labels:
        labels=[str(i) for i in range(0,len(setN))]
    cols=['c','r','orange']
    markers=['d','v','s']
    fig,ax=plt.subplots()
    #ax.scatter(nittaflux,gnonzero,color='gray')
    ax.scatter(set1[0],set1[1],color='gray')
    #ax.errorbar(nfobs,gobs,xerr=[ylow,yhigh],markersize=10.5,fmt='*',c='r')
    for i in range(0,len(setN)):
        ax.scatter(setN[i][0],setN[i][1],c=cols[i],marker=markers[i],s=65,label=labels[i])
    #ax.scatter(nittaflux,yhigh,c='g')
    xvec=np.linspace(xran[0],xran[1],10)

    #xvec=np.linspace(1e-7,1e-2,10)
    if fitline:
        #yvec=1.39*10**-11*xvec
        yvec=10**(fitline[0]*np.log10(xvec)+fitline[1])
        ax.plot(xvec,yvec,'b')
    if nfline:
        ax.plot(xvec,1.39*10**-11*xvec,'--k')
    if dashed:
        ax.plot(xvec,xvec,'--k')
    ax.set_xlim(xran)
    #ax.set_xlim([1e-7,1e-2])
    ax.set_ylim(yran)
    #ax.set_xlabel('GOES Flux ($\\rm{Wm^{-2}}$)')
    #ax.set_xlabel('Nitta Proxy ($\\rm{Wm^{-2}}$)')
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_aspect('equal')
    ax.legend(loc='upper left')
    fig.show()
    if save:
        plt.savefig(save)
    #return nfobs,gobs



def make_goes_plots():
    from scipy import stats
    ddict=pickle.load(open('plot_params.p','rb'))
    kk=ddict.keys()[0]
    ggf,gnf,gvis,mvis,mcounts,msp,mt1=ddict[kk]

    #nittaflux v. GOES
    set1=[gnf,ggf]

    #all visible flares
    s2x,s2y=[],[]
    for g,n,v in zip(ggf,gnf,gvis):
        if v==True:
            s2x.append(n)
            s2y.append(g)
    overlay2(set1,[s2x,s2y],dashed=False, nfline=True,save='GN_vis_recalc.png')
    print stats.spearmanr(s2x,s2y)

    #all visible flares, with Nitta counts condition
    s2x,s2y=[],[]
    for g,n,v in zip(ggf,gnf,gvis):
        if v==True and n >=4e6:
            s2x.append(n)
            s2y.append(g)
    fdatax=np.log10(np.array(s2x))
    fdatay=np.log10(np.array(s2y))
    line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
    print line
    overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=False, nfline=True,save='GN_vis_Ncounts.png')
    print stats.spearmanr(s2x,s2y)

    #all visible above C5.5, with Nitta counts condition > 4e6
    s2x,s2y=[],[]
    for g,n,v in zip(ggf,gnf,gvis):
        if v==True and g >=5.5e-6 and n >=1e5:
            s2x.append(n)
            s2y.append(g)
    #fit line
    fdatax=np.log10(np.array(s2x))
    fdatay=np.log10(np.array(s2y))
    line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
    print line
    overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=False, nfline=True,save='GN_vis_gtC5p5_Ncounts.png')
    print stats.spearmanr(s2x,s2y)

    #all visible above C5.5
    s2x,s2y=[],[]
    for g,n,v in zip(ggf,gnf,gvis):
        if v==True and g >=5.5e-6:
            s2x.append(n)
            s2y.append(g)
    #fit line
    fdatax=np.log10(np.array(s2x))
    fdatay=np.log10(np.array(s2y))
    line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
    print line
    overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=False, nfline=True,save='GN_vis_gtC5p5.png')
    print stats.spearmanr(s2x,s2y)

    #all visible above M, with Nitta counts condition > 4e6
    s2x,s2y=[],[]
    for g,n,v in zip(ggf,gnf,gvis):
        if v==True and g >=1e-5 and n >=1e5:
            s2x.append(n)
            s2y.append(g)
    #fit line
    fdatax=np.log10(np.array(s2x))
    fdatay=np.log10(np.array(s2y))
    line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
    print line
    overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=False, nfline=True,save='GN_vis_gtM_Ncounts.png')
    print stats.spearmanr(s2x,s2y)

    #all visible above M
    s2x,s2y=[],[]
    for g,n,v in zip(ggf,gnf,gvis):
        if v==True and g >=1e-5:
            s2x.append(n)
            s2y.append(g)
    #fit line
    fdatax=np.log10(np.array(s2x))
    fdatay=np.log10(np.array(s2y))
    line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
    print line
    overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=False, nfline=True,save='GN_vis_gtM.png')
    print stats.spearmanr(s2x,s2y)


def make_messenger_plots():
    from scipy import stats
    ddict=pickle.load(open('plot_params.p','rb'))
    kk=ddict.keys()[0]
    ggf,gnf,gvis,mvis,mcounts,msp,mt1=ddict[kk]

    for i,m in enumerate(mcounts):
        if m == '':
          mcounts[i]=0.0
    #nittaflux v. Messenger counts
    set1=[gnf,mcounts]

    #all visible flares
    s2x,s2y=[],[]
    for m,n,v in zip(mcounts,gnf,mvis):
        if v==True and m !='':
            s2x.append(n)
            s2y.append(m)
    overlay2(set1,[s2x,s2y],dashed=False,yran=[1,1e5],ylab='Messenger Peak (Counts/s)',save='MN_vis_recalc.png')
    print stats.spearmanr(s2x,s2y)

    #all visible flares, with Nitta counts condition
    s2x,s2y=[],[]
    for m,n,v in zip(mcounts,gnf,gvis):
        if v==True and m != '' and n >=4e6:
            s2x.append(n)
            s2y.append(m)
    fdatax=np.log10(np.array(s2x))
    fdatay=np.log10(np.array(s2y))
    line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
    print line
    overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=False,yran=[1,1e5],ylab='Messenger Peak (Counts/s)',save='MN_vis_recalc_Ncounts.png')
    print stats.spearmanr(s2x,s2y)


    #Nitta proxy vs. Messenger proxy
    set1=[list(1.39*10**-11*np.array(gnf)),msp]

    #all visible above C5.5, with Nitta counts condition > 4e6
    s2x,s2y=[],[]
    for m,n,v in zip(msp,gnf,gvis):
        if v==True and m!= '' and m >=5.5e-6 and n >=4e6:
            s2x.append(1.39*10**-11*n)
            s2y.append(m)
    #fit line
    fdatax=np.log10(np.array(s2x))
    fdatay=np.log10(np.array(s2y))
    line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
    print line
    overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=True, xran=[1e-7,1e-2],xlab='Nitta Proxy ($\\rm{Wm^{-2}}$)',yran=[1e-7,1e-2],ylab='Messenger Proxy ($\\rm{Wm^{-2}}$)',save='MN_vis_gtC5p5_Ncounts.png')
    print stats.spearmanr(s2x,s2y)

    #all visible above C5.5
    s2x,s2y=[],[]
    for m,n,v in zip(msp,gnf,gvis):
        if v==True and m!='' and m >=5.5e-6:
            s2x.append(1.39*10**-11*n)
            s2y.append(m)
    #fit line
    fdatax=np.log10(np.array(s2x))
    fdatay=np.log10(np.array(s2y))
    line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
    print line
    overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=True,xran=[1e-7,1e-2],xlab='Nitta Proxy ($\\rm{Wm^{-2}}$)',yran=[1e-7,1e-2],ylab='Messenger Proxy ($\\rm{Wm^{-2}}$)',save='MN_vis_gtC5p5.png')
    print stats.spearmanr(s2x,s2y)

    #all visible above M, with Nitta counts condition > 4e6
    s2x,s2y=[],[]
    for m,n,v in zip(msp,gnf,gvis):
        if v==True and m!='' and m >=1e-5 and n >=4e6:
            s2x.append(1.39*10**-11*n)
            s2y.append(m)
    #fit line
    fdatax=np.log10(np.array(s2x))
    fdatay=np.log10(np.array(s2y))
    line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
    print line
    overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=True,xran=[1e-7,1e-2],xlab='Nitta Proxy ($\\rm{Wm^{-2}}$)', yran=[1e-7,1e-2],ylab='Messenger Proxy ($\\rm{Wm^{-2}}$)',save='MN_vis_gtM_Ncounts.png')
    print stats.spearmanr(s2x,s2y)

    #all visible above M
    s2x,s2y=[],[]
    for m,n,v in zip(msp,gnf,gvis):
        if v==True and m!='' and m >=1e-5:
            s2x.append(1.39*10**-11*n)
            s2y.append(m)
    #fit line
    fdatax=np.log10(np.array(s2x))
    fdatay=np.log10(np.array(s2y))
    line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
    print line
    overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=True,xran=[1e-7,1e-2],xlab='Nitta Proxy ($\\rm{Wm^{-2}}$)',yran=[1e-7,1e-2],ylab='Messenger Proxy ($\\rm{Wm^{-2}}$)',save='MN_vis_gtM.png')
    print stats.spearmanr(s2x,s2y)


def make_messenger_Tplots():
    from scipy import stats
    ddict=pickle.load(open('plot_params.p','rb'))
    kk=ddict.keys()[0]
    ggf,gnf,gvis,mvis,mcounts,msp,mt1=ddict[kk]
    tt=11.6045*np.array(mt1)
    #Nitta proxy vs. Messenger proxy
    set1=[list(1.39*10**-11*np.array(gnf)),msp]

    #all visible with T>20MK, Nitta > M
    s2x,s2y=[],[]
    for m,n,v,t in zip(msp,gnf,gvis,tt):
        if v==True and m!= '' and t >=20.:# and 1.39e-11*n>1e-5: #and 13.9e-11*n <1e-3:
            s2x.append(1.39*10**-11*n)
            s2y.append(m)
    #fit line
    fdatax=np.log10(np.array(s2x))
    fdatay=np.log10(np.array(s2y))
    line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
    print line
    overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=True, xran=[1e-7,1e-2],xlab='Nitta Proxy ($\\rm{Wm^{-2}}$)',yran=[1e-7,1e-2],ylab='Messenger Proxy ($\\rm{Wm^{-2}}$)',save='MN_vis_gt20MK.png')
    print stats.spearmanr(s2x,s2y)
    s2x20,s2y20=s2x,s2y

    #all visible T>25MK
    s2x,s2y=[],[]
    for m,n,v,t in zip(msp,gnf,gvis,tt):
        if v==True and m!='' and t >=25. and 1.39e-11*n>=2e-6:# and 13.9e-11*n <1e-3:
            s2x.append(1.39*10**-11*n)
            s2y.append(m)
    #fit line
    fdatax=np.log10(np.array(s2x))
    fdatay=np.log10(np.array(s2y))
    line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
    print line
    overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=True,xran=[1e-7,1e-2],xlab='Nitta Proxy ($\\rm{Wm^{-2}}$)',yran=[1e-7,1e-2],ylab='Messenger Proxy ($\\rm{Wm^{-2}}$)',save='MN_vis_gt25MK.png')
    print stats.spearmanr(s2x,s2y)
    s2x25,s2y25=s2x,s2y

    #all visible T>30MK
    s2x,s2y=[],[]
    for m,n,v,t in zip(msp,gnf,gvis,tt):
        if v==True and m!='' and t >=30. and 1.39e-11*n>=2e-6:
            s2x.append(1.39*10**-11*n)
            s2y.append(m)
    #fit line
    fdatax=np.log10(np.array(s2x))
    fdatay=np.log10(np.array(s2y))
    line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
    print line
    overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=True,xran=[1e-7,1e-2],xlab='Nitta Proxy ($\\rm{Wm^{-2}}$)', yran=[1e-7,1e-2],ylab='Messenger Proxy ($\\rm{Wm^{-2}}$)',save='MN_vis_gt30MK.png')
    print stats.spearmanr(s2x,s2y)
    s2x30,s2y30=s2x,s2y

    #all visible, color coded T
    #fit line
    fdatax=np.log10(np.array(s2x))
    fdatay=np.log10(np.array(s2y))
    line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
    print line
    overlayN(set1,[[s2x20,s2y20],[s2x25,s2y25],[s2x30,s2y30]],labels=['20MK','25MK','30MK'],dashed=True,xran=[1e-7,1e-2],xlab='Nitta Proxy ($\\rm{Wm^{-2}}$)',yran=[1e-7,1e-2],ylab='Messenger Peak (Counts/s)',save='MN_vis_allT.png')
    print stats.spearmanr(s2x,s2y)


def where_GSobs(fl,sat='G'):
    #get the flares where GSobs == TRUE
    nitta,goes=[],[]
    for f in fl.list:
        if f.Observation.GSvis==True:
            nitta.append(f.Properties.Nittaflux)
            if sat == 'G':
                goes.append(f.Properties.GOES_GOES_long)
                ylabel='GOES 1-8 A Flux (W/m^2)'
            if sat == 'M':
                goes.append(f.Properties.Messenger_GOES_proxy)
                ylabel='Messenger Flux (W/m^2)'


    nitta=np.array(nitta)
    goes=np.array(goes)
    nzeros=np.where(nitta != 0.0)[0]
    nnonzero=nitta[nzeros]
    gnonzero=goes[nzeros]
    nittaflux=np.array(nnonzero)*1.39e-11
    #error bars
    ylow,yhigh=[],[]
    for n,nf in zip(nnonzero,nittaflux):
        if n > 4e6:
            ylow.append(np.abs(nf*.3-nf))
            yhigh.append(nf*3.)
        else:
            ylow.append(np.abs(nf*.5-nf))
            yhigh.append(nf*1.5)
    fig,ax=plt.subplots()
    ax.errorbar(nittaflux,gnonzero,xerr=[ylow,yhigh],fmt='o')
    #ax.scatter(nittaflux,ylow,c='m')
    #ax.scatter(nittaflux,yhigh,c='g')
    ax.plot(np.linspace(1e-8,1e-3,10),np.linspace(1e-8,1e-3,10))
    ax.set_title('Flares observed by GOES and STEREO')
    ax.set_xlim([1e-8,1e-3])
    ax.set_ylim([1e-8,1e-3])
    ax.set_xlabel('Nitta GOES EUVI proxy (W/m^2)')
    ax.set_ylabel(ylabel)
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.show()
    return nitta,goes

def where_MSobs(fl,sat='M'):
    nitta,goes=[],[]
    for f in fl.list:
        if f.Observation.MSvis==True:
            nitta.append(f.Properties.Nittaflux)
            if sat == 'M':
                goes.append(f.Properties.Messenger_GOES_long)#proxy)
                ylabel='Messenger Flux (W/m^2)'
            elif sat == 'G':
                goes.append(f.Properties.GOES_GOES_long)
                ylabel='GOES 1-8 A Flux (W/m^2)'

    nitta=np.array(nitta)
    goes=np.array(goes)
    nzeros=np.where(nitta != 0.0)[0]
    nnonzero=nitta[nzeros]
    gnonzero=goes[nzeros]
    nittaflux=np.array(nnonzero)*1.39e-11
    #error bars
    ylow,yhigh=[],[]
    for n,nf in zip(nnonzero,nittaflux):
        if n > 4e6:
            ylow.append(np.abs(nf*.3-nf))
            yhigh.append(nf*3.)
        else:
            ylow.append(np.abs(nf*.5-nf))
            yhigh.append(nf*1.5)
    #fig,ax=plt.subplots()
    #ax.errorbar(nittaflux,gnonzero,xerr=[ylow,yhigh],fmt='o')
    ##ax.scatter(nittaflux,ylow,c='m')
    ##ax.scatter(nittaflux,yhigh,c='g')
    #ax.plot(np.linspace(1e-8,1e-3,10),np.linspace(1e-8,1e-3,10))
    #ax.set_title('Flares observed by Messenger and STEREO')
    #ax.set_xlim([1e-8,1e-3])
    #ax.set_ylim([1e-8,1e-3])
    #ax.set_xlabel('Nitta GOES EUVI proxy (W/m^2)')
    #ax.set_ylabel(ylabel)
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    #fig.show()
    return nitta,goes

def line_of_best_fit(datax,datay,loglog= True,minx=False,maxx=False,minratio=False,maxratio=False):        #if list convert to np.arra
    fdatax,fdatay=[],[]

    if type(datax) == list: datax=np.array(datax)
    if type(datay) == list: datay=np.array(datay)
    if not minx and not maxx: fdatax=datax
    #if not miny and not maxy:
    fdatay=datay
    #print np.shape(fdatax),np.shape(fdatay)

    for x,y in zip(datax,datay):
        if minx and x > minx and maxx and x < maxx:
            fdatax.append(x)
            fdatay.append(y)

    if type(fdatax) == list: fdatax=np.array(fdatax)
    if type(fdatay) == list: fdatay=np.array(fdatay)
    #print np.shape(fdatax),np.shape(fdatay)

    if loglog:
        #filter out values of 0 and nan
        nz=np.where(fdatax !=0.0)
        fdatax=fdatax[nz]
        fdatay=fdatay[nz]

        fdatax=np.log10(fdatax)
        fdatay=np.log10(fdatay)

        fnan=np.where(~np.isfinite(fdatax))
        fdatax=np.delete(fdatax,fnan[0])
        fdatay=np.delete(fdatay,fnan[0])

        fnan=np.where(~np.isfinite(fdatay))
        fdatax=np.delete(fdatax,fnan[0])
        fdatay=np.delete(fdatay,fnan[0])

    #ratio=datay/datax deal with this later
    print np.shape(fdatax),np.shape(fdatay)
    line, res, _, _, _ =np.polyfit(fdatax,fdatay,1, full=True)
    print 'line of best fit: 10**('+str(line[0])+'*log10(x) + ' +str(line[1]) +')'
    print 'error: ',res
    x=10**(fdatax) #assume log
    yfit = lambda x: 10**(line[0]*np.log10(x)+line[1])
    return x,yfit(x),res,fdatax,fdatay
