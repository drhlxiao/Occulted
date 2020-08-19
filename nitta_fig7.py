def fig7(nitta=False, goes=False):
    if not nitta and not goes:
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

def Mproxy(fl):
    for f in fl.list:
        if f.Properties.Messenger_counts_hi != '':
            #f.Properties.Messenger_GOES_selfproxy=(1.13081*10**-8)*float(f.Properties.Messenger_GOES_long)**.96135
            #f.Properties.Messenger_GOES_selfproxy=(0.272*10**-9)*float(f.Properties.Messenger_counts_hi)**1.176
            #f.Properties.Messenger_GOES_selfproxy=(2.36*10**-8)*float(f.Properties.Messenger_counts_hi)**1.015
            f.Properties.Messenger_GOES_selfproxy=(1.4*10**-9)*float(f.Properties.Messenger_counts_hi)**1.22
            #f.Properties.Messenger_GOES_selfproxy=(2.8944*10**-9)*float(f.Properties.Messenger_counts_hi)**1.16


def Mproxy_list(xx):
    #return (2.36*10**-8)*np.array(xx)**1.015
    #return (2.59*10**-9)*np.array(xx)**1.18
    return (1.4*10**-9)*np.array(xx)**1.22

def fitline2eq(fitline):
    '''assume fitline is list'''
    const=10.**(-1*fitline[1]/fitline[0])
    exp=1./fitline[0]
    print 'F_goes = '+str(const)+'x 10^'+str(exp)

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

def fig_occ(fl=False):
    #leave out entries without occultation factors
    #sort by occultation factor, label by ... date?
    import matplotlib.dates as mdates
    if not fl:
        fl=pickle.load(open('occulted_list_final.p','rb'))
    oc=fl.isolate_att('e2',key='distance2limb')
    ol=[o.value for o in oc]
    occ=np.array(ol)/10.208 #rsun (arcsec)*2 /180 formerly/5.33
    #glong=fl.isolate_att('GOES_GOES_long')
    #proxy=fl.isolate_att('Messenger_GOES_selfproxy') #check that this is the right proxy...
    glong=fl.isolate_att('glong_table')#[4.9e-6,2.3e-6,9.6e-7,1.4e-6,8.7e-7,7.2e-6,5.4e-5,7.8e-4,9.1e-5,7.1e-5,5.6e-5,7.4e-6,8.8e-4,5.6e-6,3.12e-7,6.5e-6] #from table - where did 2012-05-05 go?
    proxy=Mproxy_list(fl.isolate_att('Messenger_counts_hi')) #also from table
    dates=fl.isolate_att('Messenger_peak')#[f.Datetimes.Messenger_peak for f in fl.list if 'e2' in f.Properties.__dict__.keys()]
    datestr=[d.strftime("%Y/%m/%d") for d in dates]
    #zip and sort? or not, scatter takes care of it
    yy=[(i,j) for i,j in zip(glong,proxy)]
    #print yy
    fig,ax=plt.subplots()
    ax2=ax.twiny()
    ax.plot((occ,occ),([i for (i,j) in yy], [j for (i,j) in yy]),c='black')
    ax.scatter(occ,glong,c='r', label='GOES') #dt, goes_goes_long
    ax.scatter(occ,proxy,c='b', label='Messenger') #dt, Mproxy
    #twin x axis and put date labels there...
    ax.set_xlabel('Occultation Factor (degrees)')
    ax.set_ylabel('SXR Flux')
    ax.set_yscale('log')
    ax.set_ylim([1e-7,1e-3])
    ax.set_xlim([0,180])
    ax2.set_xlim([0,180])
    ax.set_xticks([0,20,40,60,80,100,120,140,160,180])
    ax.set_xticklabels([0,20,40,60,80,100,120,140,160,180])
    ax2.set_xticks(occ)
    ax2.set_xticklabels(datestr, rotation=70)
    ax.legend(loc='lower right')
    print occ, datestr
    #myFmt = mdates.DateFormatter('%Y/%m/%D')
    #ax2.xaxis.set_major_formatter(myFmt)
    #plt.gcf().autofmt_xdate()
    fig.show()

def goes_temp_from_ratio(fl,int=False):
    import pidly
    idl = pidly.IDL('/Users/wheatley/Documents/Solar/sswidl_py.sh')
    tarr,emarr=[],[]
    for f in fl.list:
        f.Datetimes.convert2string()
        tt=f.Datetimes.messenger_peak_bin['long'][0]
        f.Datetimes.convert2datetime()
#         if int:
#             gs=f.Properties.GOES_short_int
#             gl=f.Properties.GOES_long_int
#             ydata=[gl,gs]
#         else:
        gs=f.Properties.GOES_GOES
        gl=f.Properties.GOES_GOES_long
        ydata=[gl,gs]
        idl('tarray',[tt])
        idl('ydata',ydata)
        idl('goes_tem,tarray=tarray,yclean=ydata,tempr=t,emis=em') #returns data={taxis:taxis,phigh:phigh}
        t=idl.t
        em=idl.em
        tarr.append(t)
        emarr.append(em)
        f.Properties.GOES_T_peak=t
        f.Properties.GOES_EM_peak=em
        gsi=f.Properties.GOES_short_int
        gli=f.Properties.GOES_long_int
        ydatai=[gli,gsi]
        idl('ydatai',ydatai)
        idl('goes_tem,tarray=tarray,yclean=ydatai,tempr=t1,emis=em1') #returns data={taxis:taxis,phigh:phigh}
        t1=idl.t1
        em1=idl.em1
        f.Properties.GOES_T_int=t1
        f.Properties.GOES_EM_int=em1

    idl.close()
    return tarr,emarr


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

def noverlay(set1,fitline=False,title='',yran=[10,1e4],ylab='Messenger Peak (Counts/s)',xran=[1e-7,1e-4],xlab='GOES Flux ($\\rm{Wm^{-2}}$)',dashed=True,nfline=False):
    fig,ax=plt.subplots()
    #ax.scatter(nittaflux,gnonzero,color='gray')
    ax.scatter(set1[0],set1[1],color='gray')
    #ax.errorbar(nfobs,gobs,xerr=[ylow,yhigh],markersize=10.5,fmt='*',c='r')
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


def overlay2(set1,set2,fitline=False,title='',yran=[1e-7,1e-3],ylab='GOES Flux ($\\rm{Wm^{-2}}$)',xran=[1e-7,1e-3],xlab='Nitta Flux (photons/s)',dashed=True,nfline=False,save=False,legend=False,lab1='set 1',lab2='set 2',log=True):
    fig,ax=plt.subplots()
    #ax.scatter(nittaflux,gnonzero,color='gray')
    ax.scatter(set1[0],set1[1],color='gray',label=lab1)
    #ax.errorbar(nfobs,gobs,xerr=[ylow,yhigh],markersize=10.5,fmt='*',c='r')
    ax.scatter(set2[0],set2[1],c='r',marker='v',s=60,label=lab2)
    #ax.scatter(nittaflux,yhigh,c='g')
    xvec=np.linspace(xran[0],xran[1],10)

    #xvec=np.linspace(1e-7,1e-2,10)
    if fitline:
        #yvec=1.39*10**-11*xvec
        yvec=10**(fitline[0]*np.log10(xvec)+fitline[1])
        ax.plot(xvec,yvec,'b')
    if nfline:
        xm5=np.linspace(4e6,xran[1],10) #>4e6
        ax.plot(xm5,1.39*10**-11*xm5,'--k')
    if type(dashed) !=bool:
        #print dashed
        ax.plot(xvec,dashed[0]*(xvec**dashed[1]),'--k')
    elif dashed==True:
        ax.plot(xvec,xvec,'--k')
    ax.set_xlim(xran)
    #ax.set_xlim([1e-7,1e-2])
    ax.set_ylim(yran)
    #ax.set_xlabel('GOES Flux ($\\rm{Wm^{-2}}$)')
    #ax.set_xlabel('Nitta Proxy ($\\rm{Wm^{-2}}$)')
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    if log:
        ax.set_xscale('log')
        ax.set_yscale('log')
    ax.set_aspect('equal')
    if legend:
        ax.legend()
    fig.show()
    if save:
        plt.savefig(save)
    #return nfobs,gobs

def plot_diff_from_fit(set1, fitline,ylab='Difference from fitline',xlab='foo',xran=[1e4,1e9],yran=[1e-8,1e-3],title='foo',log=False, percent=True, hist=True, stacked=True):
    '''plot difference of individual points with fitline '''
    st=False
    xx=set1[0]
    yy=set1[1]
    yvec=10**(fitline[0]*np.log10(np.array(xx))+fitline[1])
    ydiff=np.array(yy)-yvec
    if percent:
        ydiff=(ydiff/xx)*100.
    bins=[-400,-350,-300,-250,-200,-175,-150,-125,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200,250,300,350,400]
    fig,ax=plt.subplots()
    if hist:
        #ydiff=np.abs(ydiff)
        if stacked: #divide up ydiff based on goes class?
            yA,yB,yC,yM=[],[],[],[]
            for x,y in zip(set1[0],set1[1]):
                if x <=1e-7 and x!=0.:
                    yA.append(((y-x)/x)*100)
                elif x <=1e-6 and x >=1e-7:
                    yB.append(((y-x)/x)*100)
                elif x <=1e-5 and x >=1e-6:
                    yC.append(((y-x)/x)*100)
                elif x >=1e-5:
                    yM.append(((y-x)/x)*100)
            #newx=[x for x in [yA,yB,yC,yM] if len(x) >0]
            #length = max(map(len, newx))
            #ydiff=np.array([xi+0.*(length-len(xi)) for xi in newx])
            #ydiff=np.array([np.array(yA),np.array(yB),np.array(yC),np.array(yM)])
            #st=True
            print np.mean(yA),np.std(yA)
            print np.mean(yB),np.std(yB)
            print np.mean(yC),np.std(yC)
            print np.mean(yM),np.std(yM)
            ax.hist([yA,yB,yC,yM],color=["red","blue","green","violet"],bins=bins,stacked=True,histtype='bar', label=["A","B","C",">M"])
            ax.legend(loc='upper right')
            #ax.legend({label1: "A", label2: "B", label3: "C",label4:">M"})
        else:
            ydiff = [((y-x)/x)*100 for x,y in zip(set1[0],set1[1])]
            print "mean",np.mean(ydiff),"dev",np.std(ydiff)
            #do plus-minus
            yplus=[((y-x)/x)*100 for x,y in zip(set1[0],set1[1]) if y>=x]
            yminus=[((x-y)/x)*100 for x,y in zip(set1[0],set1[1]) if y<x]
            print "+ mean",np.mean(yplus),"+ dev",np.std(yplus)
            print "- mean",np.mean(-1*np.array(yminus)),"- dev",np.std(-1*np.array(yminus))
            print "|| mean",np.mean(yplus+yminus),"|| dev",np.std(yplus+yminus)
            ax.hist(ydiff,bins=bins,color='g',histtype='step', label='overall')
            ax.hist(yplus,bins=bins,color='m',alpha=0.5,label='plus')
            ax.hist(yminus,bins=bins,color='c',alpha=0.5,label='minus')
            ax.legend(loc='upper right')
        ax.set_xlabel('percent difference from fitline')
    else:
        ax.scatter(xx,ydiff)
        ax.set_xlim(xran)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.set_xscale('log')
        if log:
            ax.set_yscale('log')
        elif percent:
            ax.set_ylim([-1000,1000])
    #ax.set_aspect('equal')
    fig.show()

def plot_all_hist(flmg=False,flmgs=False,fl_counts_long=[.816,7.225],fl_counts_short=[.65,6.876],fl_flux_long=[.931,-.748],fl_flux_short=[.908,-.879],proxy=[.845,-1.56],fluxone=True):
    if not flmg:
        flmg=pickle.load(open('TEM_recalc_newerange.p','rb'))
    if not flmgs:
        flmgs=pickle.load(open('fig5list.p','rb'))
    if fluxone:
        fl_flux_long=[1.,1.]
        fl_flux_short=[1.,1.] #fl-counts_short[.769,7.126]
        proxy=[1.,1.]

    clong=flmg.isolate_att('Messenger_counts_hi')
    cshort=flmg.isolate_att('Messenger_counts_lo')
    flong=flmg.isolate_att('Messenger_GOES_TEM_long')
    fshort=flmg.isolate_att('Messenger_GOES_TEM_short')
    glong=flmg.isolate_att('GOES_long_int')
    gshort=flmg.isolate_att('GOES_short_int')

    Mproxy(flmgs)
    mproxy=flmgs.isolate_att('Messenger_GOES_selfproxy')
    nproxy=(1.39*10**-11)*np.array(flmgs.isolate_att('Nittaflux'))

    def _get_diff(xx,yy,line):
        ydiffplus,ydiffminus=[],[]
        yfit = lambda x: 10**(line[0]*np.log10(x)+line[1])
        ypred=yfit(xx)
        ydiff=np.array(yy)-ypred
        ydiff=(ydiff/ypred) #y_actual - y_predicted / y_predicted
        ydiff=ydiff[~np.isnan(ydiff)]
        ydiff=ydiff[~np.isinf(ydiff)]
        ydiffminus=[y for y in ydiff if y < 0.0]
        ydiffplus=[y for y in ydiff if y > 0.0]
        print np.shape((ydiff)),np.min(ydiff),np.max(ydiff)
        return ydiffplus, ydiffminus

    bns=np.logspace(-3,2,20,endpoint=True)#range(-10,10,1)
    bins=list(-1.*(bns))
    bins.reverse()
    bins.extend(list(bns))
    bins=np.linspace(-200,200,50)
    normed=True

    #take absolute value of difference
    #color code thos that were negative...or flip so that it's pred-val instead of the other way?
    gcldp,gcldm=_get_diff(glong,clong,fl_counts_long)
    gcsdp,gcsdm=_get_diff(gshort,cshort,fl_counts_short)
    gfldp,gfldm=_get_diff(glong,flong,fl_flux_long)
    gfsdp,gfsdm=_get_diff(gshort,fshort,fl_flux_short)
    npdp,npdm=_get_diff(nproxy,mproxy,proxy)

    fig,ax=plt.subplots()
    ax.hist(gcldp,bins=bins,label='counts vs. flux (1-8$\AA$),',linewidth=2,histtype='step',normed=normed)
    ax.hist(gcldm,bins=bins,linewidth=2,linestyle=('dashed'),histtype='step',normed=normed)
    ax.hist(gcsdp,bins=bins,label='counts vs. flux (.5-4$\AA$)',facecolor=None,linewidth=2,histtype='step',normed=normed)
    ax.hist(gcsdm,bins=bins,linestyle=('dashed'),facecolor=None,linewidth=2,histtype='step',normed=normed)
    ax.hist(gfldp,bins=bins,label='flux vs. flux (1-8$\AA$)',facecolor=None,linewidth=2,histtype='step',normed=normed)
    ax.hist(gfldm,bins=bins,linewidth=2,linestyle=('dashed'),histtype='step',normed=normed)
    ax.hist(gfsdp,bins=bins,label='flux vs. flux (.5-4$\AA$)',facecolor=None,linewidth=2,histtype='step',normed=normed)
    ax.hist(gfsdm,bins=bins,linewidth=2,linestyle=('dashed'),histtype='step',normed=normed)
    ax.hist(npdp,bins=bins,label='proxy vs. proxy',facecolor=None,linewidth=2,histtype='step',normed=normed)
    ax.hist(npdm,bins=bins,linewidth=2,linestyle=('dashed'),histtype='step',normed=normed)

    ax.set_xlabel('Percent Difference from Fitline')
    print np.min(bins),np.max(bins)
    ax.set_xlim(np.min(bins),np.max(bins))
    #ax.set_xscale('log')
    #ax.set_aspect('equal')
    ax.legend()
    fig.show()



def counts_v_counts(fl1,fl2,glim=[10,1e6],channel='long',diff=False,hist=True,stacked=True):
    if channel =='long':
        title='1-8 $\AA$'
    else:
        title='0.5-4 $\AA$'
    s1x,s1y,_=fl1.plot_goes_messenger(channel=channel,gcounts=True)
    s2x,s2y,lparams=fl2.plot_goes_messenger(channel=channel,line=True,linexlim=glim,gcounts=True) #need to return line too
    print stats.spearmanr(s2x,s2y)
    overlay2([s1x,s1y],[s2x,s2y],xran=[10,1e6],yran=[10,1e6],fitline=list(lparams),ylab='Messenger Peak (Counts/s)',xlab='GOES Counts/s',dashed=False,title=title)
    fitline2eq(list(lparams))
    if diff:
        plot_diff_from_fit([s2x,Mproxy_list(s2y)],fitline=list(lparams),xran=[1e-8,1e-3],yran=[10,1e6],xlab='GOES Flux ($\\rm{Wm^{-2}}$)',hist=hist,stacked=stacked)


def flux_v_counts(fl1,fl2,glim=[2e-6,1e-3],yran=[10,1e5],xran=[1e-7,1e-3],channel='long',diff=True,hist=True,stacked=True,swpc=False,gint=False,nofit=False):
    if channel =='long':
        title='1-8 $\AA$'
    else:
        title='0.5-4 $\AA$'
    s1x,s1y,_=fl1.plot_goes_messenger(channel=channel,swpc=swpc,gint=gint)
    s2x,s2y,lparams=fl2.plot_goes_messenger(channel=channel,line=True,linexlim=glim,swpc=swpc,gint=gint) #need to return line too
    print stats.spearmanr(s2x,s2y)
    if nofit:
        fitline=False
    else:
        fitline=list(lparams)
    overlay2([s1x,s1y],[s2x,s2y],xran=xran,yran=yran,fitline=fitline,ylab='Messenger Peak (Counts/s)',xlab='GOES Flux ($\\rm{Wm^{-2}}$)',dashed=False,title=title)
    fitline2eq(list(lparams))
    if diff:
        plot_diff_from_fit([s2x,Mproxy_list(s2y)],fitline=list(lparams),xran=xran,yran=yran,xlab='GOES Flux ($\\rm{Wm^{-2}}$)',hist=hist,stacked=stacked)


def flux_v_flux(fl1,fl2,glim=[2e-6,1e-3],channel='long', useHEK=False, diff= True, one=True, log=True,xran=[1e-8,1e-3],yran=[1e-8,1e-3],swpc=False, gint=False):
    if channel =='long':
        title='1-8 $\AA$'
    else:
        title='0.5-4 $\AA$'
    s1x,s1y,_=fl1.plot_goes_messenger(channel=channel,from_T_EM='1',swpc=swpc,gint=gint)
    s2x,s2y,lparams=fl2.plot_goes_messenger(channel=channel,line=True,linexlim=glim,from_T_EM='1',useHEK=useHEK,swpc=swpc,gint=gint) #need to return line too
    overlay2([s1x,s1y],[s2x,s2y],xran=xran,yran=yran,fitline=list(lparams),ylab='GOES Flux calculated from Messenger T,EM ($\\rm{Wm^{-2}}$)',xlab='GOES Flux ($\\rm{Wm^{-2}}$)',dashed=True,title=title)
    fitline2eq(list(lparams))
    if diff:
        if one:
            lparams=(1.0,0.)
        plot_diff_from_fit([s2x,s2y],fitline=list(lparams),xran=[1e-8,1e-3],yran=[1e-8,1e-3],xlab='GOES Flux ($\\rm{Wm^{-2}}$)',title='flux_v_flux',log=log,stacked=False)

def nf_v_counts(fl1,fl2,glim=[1e3,1e7],channel='long',line2=False,fitline=True):
    s1x=fl1.isolate_att('Nittaflux')
    s1y=fl1.isolate_att('Messenger_counts_hi')
    s1xr,s1yr=[],[]

    for x,y in zip(s1x,s1y):
        if not np.isnan(x):
            s1xr.append(x)
            s1yr.append(y)

    s2x=fl2.isolate_att('Nittaflux')
    s2y=fl2.isolate_att('Messenger_counts_hi')
    s2xr,s2yr=[],[]

    for x,y in zip(s2x,s2y):
        if not np.isnan(x) and x > 0.0:
            s2xr.append(x)
            s2yr.append(y)

    if fitline:
        lparams=line_of_best_fit(s2xr,s2yr,minx=glim[0],maxx=glim[1])
        print stats.spearmanr(s2xr,s2yr)
        flist=list(lparams)
        fitline2eq(flist)
    else:
        flist=False
    if type(line2) !=bool:
        dashed=line2
    overlay2([s1xr,s1yr],[s2xr,s2yr],xran=[1e4,1e8],yran=[10,1e5],fitline=flist,ylab='Messenger Peak Counts/s',xlab='EUV Flux (photons/s)',dashed=line2)
    return s2xr,s2yr

def nf_v_proxy(fl1,fl2,glim=[1e-7,1e-4],channel='long', diff=True,stacked=False):
    s1x=fl1.isolate_att('Nittaflux')
    s1y=fl1.isolate_att('Messenger_GOES_selfproxy')
    s1xr,s1yr=[],[]

    for x,y in zip(s1x,s1y):
        if not np.isnan(x):
            s1xr.append(1.39*10**-11*x)
            s1yr.append(y)

    s2x=fl2.isolate_att('Nittaflux')
    s2y=fl2.isolate_att('Messenger_GOES_selfproxy')
    s2xr,s2yr=[],[]

    for x,y in zip(s2x,s2y):
        if not np.isnan(x) and x > 0.0:
            s2xr.append(1.39*10**-11*x)
            s2yr.append(y)

    lparams=line_of_best_fit(s2xr,s2yr,minx=glim[0],maxx=glim[1])
    print stats.spearmanr(s2xr,s2yr)
    print stats.linregress(np.log10(s2xr),np.log10(s2yr))
    fitline2eq(list(lparams))
    overlay2([s1xr,s1yr],[s2xr,s2yr],xran=[1e-7,1e-3],yran=[1e-7,1e-3],fitline=list(lparams),ylab='Messenger GOES Proxy ($\\rm{Wm^{-2}}$)',xlab='EUV GOES Proxy ($\\rm{Wm^{-2}}$)',dashed=True)
    if diff:
        plot_diff_from_fit([s2xr,s2yr],fitline=[1.,0.],xran=[1e-7,1e-2],yran=[1e-7,1e-2],xlab='EUV GOES Proxy ($\\rm{Wm^{-2}}$)',stacked=stacked)

def nf_combine(fl,glim=[1e-6,1e-3]):
    ''' combine plots 5a and 5b. need 2 axes for this...alignment might be tricky'''
    from matplotlib.transforms import Transform
    from mpl_toolkits.axes_grid1 import host_subplot
    s1x=fl.isolate_att('Nittaflux')
    s1y=fl.isolate_att('Messenger_GOES_selfproxy')
    s1xr,s1yr=[],[]


    def Mproxy(x):
        return (1.4*10**-9)*x**1.22
    def Mproxyr(x):
        return (x/(1.4*10**-9))**(1./1.22)

    def Nproxy(x):
        return 1.39*10**-11*x
    def Nproxyr(x):
        return x/(1.39*10**-11)

    for x,y in zip(s1x,s1y):
        if not np.isnan(x) and x >0.0:
            s1xr.append(1.39*10**-11*x)
            s1yr.append(y)

    lparams=line_of_best_fit(s1xr,s1yr,minx=glim[0],maxx=glim[1])
    print stats.spearmanr(s1xr,s1yr)
    print stats.linregress(np.log10(s1xr),np.log10(s1yr))
    fitline2eq(list(lparams))
    #overlay2([s1xr,s1yr],[s2xr,s2yr],xran=[1e-7,1e-3],yran=[1e-7,1e-3],fitline=list(lparams),ylab='Messenger SXR Proxy ($\\rm{Wm^{-2}}$)',xlab='Nitta SXR Proxy ($\\rm{Wm^{-2}}$)',dashed=True)

    s2x=fl.isolate_att('Nittaflux')
    s2y=fl.isolate_att('Messenger_counts_hi')
    s2xr,s2yr=[],[]

    for x,y in zip(s2x,s2y):
        if not np.isnan(x) and x > 0.0:
            s2xr.append(x)
            s2yr.append(y)

    lparams=line_of_best_fit(s2xr,s2yr)
    print stats.spearmanr(s2xr,s2yr)
    fitline2eq(list(lparams))
    #overlay2([s1xr,s1yr],[s2xr,s2yr],xran=[1e4,1e8],yran=[10,1e5],fitline=list(lparams),ylab='Messenger Peak Counts/s',xlab='Nitta EUV Flux (photons/s)',dashed=False)

    #plot
    fig, ax = plt.subplots()
    #ax = host_subplot(111, adjustable='box-forced', aspect='equal')
    secxax = ax.secondary_xaxis('top', functions=(Nproxy, Nproxyr))
    secyax = ax.secondary_yaxis('right', functions=(Mproxy, Mproxyr))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([1e4,1e8])
    #ax.set_xlim([1e-7,1e-2])
    ax.set_ylim([10,1e5])

    #ax.set_aspect('equal')
    #ax2= ax.twinx().twiny()
    #ax2.set_xscale('log')
    #ax2.set_yscale('log')
    #ax2.set_aspect('equal')
    x2lim=[1.39e-7,1.39e-3]
    #y2lim=1.39*10**-11*np.array(x2lim)
    #ax2.set_xlim(x2lim)
    #ax2.set_ylim(x2lim)
    ax.scatter(s2xr,s2yr,c='r',marker='v',s=60)
    #ax2.scatter(s1xr,s1yr,c='b',marker='v',s=60)
    #ax.scatter(nittaflux,yhigh,c='g')
    #xvec=np.linspace(xran[0],xran[1],10)

    xvec=np.linspace(1e-7,1e-2,10)
        #yvec=1.39*10**-11*xvec
    #yvec=10**(fitline[0]*np.log10(xvec)+fitline[1])
    #ax.plot(xvec,yvec,'b')
    #if nfline:
    #    ax.plot(xvec,1.39*10**-11*xvec,'--k')
    #ax2.plot(xvec,xvec,'--k')
    secyax.set_ylabel('Messenger SXR Proxy ($\\rm{Wm^{-2}}$)')
    secxax.set_xlabel('Nitta SXR Proxy ($\\rm{Wm^{-2}}$)')
    ax.set_ylabel('Messenger Peak counts/s')
    ax.set_xlabel('Nitta EUV Flux (photons/s)')
    #ax.set_title(title)


    fig.show()
    #return ax2.get_xticks(),ax2.get_yticks()

def proxy_over_under(fl):
    '''plot messenger-to-nitta proxy vs nitta proxy'''
    s1x=fl.isolate_att('Nittaflux')
    s1y=fl.isolate_att('Messenger_GOES_selfproxy')
    s1xr,s1yr=[],[]

    for x,y in zip(s1x,s1y):
        if not np.isnan(x) and x >0.0:
            s1xr.append(1.39*10**-11*x)
            s1yr.append(y)

    fig,ax=plt.subplots()
    ratio=np.array(s1yr)/np.array(s1xr)
    print np.min(ratio),np.max(ratio)
    ax.scatter(s1xr,ratio,c='r',marker='v',s=60)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([1e-2,1e2])
    ax.set_xlim([1e-7,1e-3])
    plt.axhline(1.0)
    ax.set_ylabel('Messenger Proxy/Nitta Proxy')
    ax.set_xlabel('Nitta SXR Proxy ($\\rm{Wm^{-2}}$)')
    fig.show()

def together_with_em(fl):
    tem=fl.isolate_att('Messenger_GOES_TEM_long')
    gint=fl.isolate_att('GOES_long_int')
    mem=fl.isolate_att('EM2')
    gem=fl.isolate_att('GOES_EM_int')

    fig,ax=plt.subplots()
    ax1=ax.twiny()
    ax.scatter(tem,gint,marker='v',c='r',s=60)
    ax1.scatter(mem,gem,c='grey',s=60)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax.set_xlim([1e-7,1e-3])
    ax.set_ylim([1e-7,1e-3])
    ax1.set_xlim([1e-3,10])
    ax1.set_ylim([1e-3,10])
    fig.show()


def nf_v_gf(fl1,fl2,glim=[1e-6,1e-3],nlim=[1e3,1e7],channel='long',tmax=False):
    s1x=fl1.isolate_att('Nittaflux')
    s1y=fl1.isolate_att('GOES_long_int')
    s1xr,s1yr=s1x,s1y#[],[]

#     for x,y in zip(s1x,s1y):
#         if not np.isnan(x):
#             s1xr.append(x)
#             s1yr.append(y)

    s2x=fl2.isolate_att('Nittaflux')
    s2y=fl2.isolate_att('GOES_long_int')
    s2xr,s2yr=[],[]

    for x,y in zip(s2x,s2y):
        if not np.isnan(x) and x>0.0:
            s2xr.append(x)
            s2yr.append(y)

    if glim:
        s3xr=[]
        s3yr=[]
        for x,y in zip(s2xr,s2yr):
            if y >= glim[0] and y <=glim[1]:
                s3xr.append(x)
                s3yr.append(y)
        s2xr=s3xr
        s2yr=s3yr

    if tmax:
        tt=11.6045*np.array(fl2.isolate_att('T1'))
        s2xr,s2yr=[],[]
        for x,y,t in zip(s2x,s2y,tt):
           if t > tmax and not np.isnan(x):
               s2xr.append(x)
               s2yr.append(y)

    lparams=line_of_best_fit(s2xr,s2yr,minx=nlim[0],maxx=nlim[1])
    print stats.spearmanr(s2xr,s2yr)
    print stats.linregress(np.log10(s2xr),np.log10(s2yr))
    fitline2eq(list(lparams))
    overlay2([s1xr,s1yr],[s1xr,s1yr],xran=[1e4,1e8],yran=[1e-7,1e-3],ylab='GOES Flux ($\\rm{Wm^{-2}}$)',xlab='EUV Flux (photons/s)',dashed=False,nfline=True,fitline=list(lparams))

def occulted_flares_sum(fl):
    ''' plot summary figure for occulted flares'''
    gc0=fl.isolate_att('GOES_GOES_long') #what about the catalog values? want to use those...edit the flare list itself then
    gcp=fl.isolate_att('Messenger_GOES_selfproxy')
    #occultation factor
    ocf=fl.isolate_att('occultation_factor')
    fig,ax=plt.subplots()
    ax.scatter(ocf,gc0, label='') #occultation factor vs gc0
    ax.scatter(ocf,gcp,label='') #occultation factor vs gc1 -- how to make vertical bars between them?
    ax.set_xlim([0,90])
    ax.set_ylim([])
    ax.set_yscale('log')
    ax.set_xlabel('Occultation Factor (degrees)')
    ax.set_ylabel('')
    ax.legend()
    fig.show()


def euvi_txt2dict(txt='data/stereo-aia/euvi_events.txt'):
    dts,cadence,ab,enum,wave,pfbg,fenhance,bgdelta,fluence,rt,decay,tdelay,gc,l0,b0,l_SC,b_SC,les,b,r=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    with open(txt,'r') as tfile:
        elist=tfile.readlines()
    for e in elist[4:]:
        lstr=e.split(' ')
        ll=[l for l in lstr if l !='']
        dts.append(dt.strptime(ll[0]+ll[1]+ll[2]+' '+ll[3]+ll[4],'%Y%m%d %H%M'))
        cadence.append(float(ll[5]))
        ab.append(ll[6])
        enum.append(int(ll[7]))
        wave.append(int(ll[8]))
        pfbg.append(float(ll[9]))
        fenhance.append(float(ll[10]))
        bgdelta.append(float(ll[11]))
        fluence.append(float(ll[12]))
        rt.append(int(ll[13]))
        decay.append(int(ll[14]))
        tdelay.append(int(ll[15]))
        gc.append(ll[16])
        l0.append(int(ll[17]))
        b0.append(int(ll[18]))
        l_SC.append(int(ll[19]))
        b_SC.append(int(ll[20]))
        les.append(int(ll[21]))
        b.append(int(ll[22]))
        r.append(float(ll[23][:-1]))

    edict={'event_datetime':dts,'cadence':cadence,'spacecraft':ab,'event_num':enum,'lambda':wave,'preflare_bg':pfbg,'flux_enhance':fenhance,'bgdelta':bgdelta,'fluence':fluence,'rise_time':rt,'decay_time':decay,'time_delay':tdelay,'GOES_class':gc,'l0':l0,'b0':b0,'l_SC':l_SC,'b_SC':b_SC,'l':les,'b':b,'r':r}
    return edict


def gc_str2flux(gcstr):
    gcflt=float(gcstr[1:])
    if gcstr.startswith('X'):
        return gcflt*1e-4
    elif gcstr.startswith('M'):
        return gcflt*1e-5
    elif gcstr.startswith('C'):
        return gcflt*1e-6
    elif gcstr.startswith('B'):
        return gcflt*1e-7
    else:
        return gcflt*1e-8

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

def hek_or_stereo(fl):
    '''determine if stereo has better view than goes for flares on earth limb'''
    clist=[f.Observation.HEK_coords[0] for f in fl.list]
    svis=[f.Observation.MSvis for f in fl.list]
    bvis=[c for i,c in enumerate(clist) if svis[i] == True ] #coords
    fvis=[f for i,f in enumerate(fl.list) if svis[i] == True ] #flare objects
    xx,yy=[],[]

    for b in bvis:
        xx.append(float(b[b.find('(')+1:b.find(' ')]))
        yy.append(float(b[b.find(' ')+1:b.find(')')]))

    llist=[]
    for x,y,f in zip(xx,yy,fvis):
        if np.abs(x)>800.:
            llist.append(f)

    return llist


def download_stereo_ex(f,overwrite=True,peak=True):
  '''Download the preflare image - one hour before the peak'''
  from datetime import timedelta as td

  #determine download constraints
  start_time=f.Datetimes.Messenger_peak - td(minutes=65)
  end_time=f.Datetimes.Messenger_peak - td(minutes=55)
  tint=[start_time,end_time] #n files around Messenger peak. datetime range.

  try:
    pffile=f.Files.prepped['stereo-preflare']
    if pffile[-5]!=f.Observation.STEREO:
      overwrite=True
  except AttributeError:
      pass
  if overwrite and type(f.Observation.STEREO) == str:
      #first check ext. hard drive for a file that meets the requirements
      lacie='/Volumes/LaCie/occulted_flares/data/stereo-aia/'
      pfc_sameday=glob.glob(lacie+dt.strftime(start_time,'%Y%m%d_')+'*'+f.Observation.STEREO +'.fts') #let's hope nothing is around midnight
      #get times from all of these filenames and see if they are within the right time range
      try:
          pfc_times=[dt.strptime(p[47:-10],'%Y%m%d_%H%M%S') for p in pfc_sameday]
          in_int=[p for p in pfc_times if p > start_time and p < end_time]
          #print in_int, start_time,end_time
          if len(in_int) > 0:
            f.Files.prepped['stereo-preflare']=pfc_sameday[pfc_times.index(in_int[0])]
            return
      except TypeError: #there's no such files
          pass
          print 'running query, these fies will need to be secchi_prepped'
          f._query(tint=tint,wave='195')
          #need to put these back into Files.prepped!
      try:
          f.Files.prepped['stereo-preflare']=f.Files.Raw['stereo'][0]
          pfstr=f.Files.prepped['stereo-preflare']
          pftime=dt.strptime(pfstr[pfstr.rfind('/')+1:pfstr.rfind('_')],'%Y%m%d_%H%M%S')
          peaktime=dt.strptime(f.Files.prepped['stereo-peak'][-25:-10],'%Y%m%d_%H%M%S')
          f.Datetimes.peakdiff = (peaktime-pftime).total_seconds()
      except IndexError:
          print "nothing in Files.Raw['stereo']!"
          pass
      if peak:
          try:
              pfile=f.Files.prepped['stereo-peak']
              if pfile[-5] !=f.Observation.STEREO:
                  opp=pfile[:-9]+'14eu'+f.Observation.STEREO+'.fts'
                  if os.path.isfile('/Volumes/LaCie/occulted_flares/data/stereo-aia/'+opp):
                      f.Files.prepped['stereo-peak']=opp
                  else:
                      tint=dt.strptime(pfile[:-10],'%Y%m%d_%H%M%S')
                      print 'running query'
                      f._query(tint=[tint-td(minutes=2.5),tint+td(minutes=2.5)],wave='195')
          except AttributeError:
              pass

def rewrite_peak(f,infiles):
    #find the right file
    st=f.Observation.STEREO
    rightfiles=[ll for ll in infiles if ll.endswith(st+'.fts')]
    pftimes=[]
    for pfile in rightfiles:
        pftimes.append(dt.strptime(pfile[pfile.rfind('/')+1:pfile.rfind('_')],'%Y%m%d_%H%M%S'))

    #get the nearest time
    mt=f.Datetimes.Messenger_peak
    nearest= min(pftimes, key=lambda x: abs(x - mt))
    idx=pftimes.index(nearest)
    f.Files.prepped['stereo-peak']=rightfiles[idx]

def rewrite_preflare(f,infiles):
    #find the right file
    st=f.Observation.STEREO
    rightfiles=[ll for ll in infiles if ll.endswith(st+'.fts')]
    pftimes=[]
    for pfile in rightfiles:
        pftimes.append(dt.strptime(pfile[pfile.rfind('/')+1:pfile.rfind('_')],'%Y%m%d_%H%M%S'))

    #get the nearest time
    mt=f.Datetimes.Messenger_peak- td(minutes=60)
    nearest= min(pftimes, key=lambda x: abs(x - mt))
    idx=pftimes.index(nearest)
    f.Files.prepped['stereo-preflare']=rightfiles[idx]


def make_goes_plots(fl=False): #separate by temperature now
    from scipy import stats
    if not fl:
        ddict=pickle.load(open('plot_params.p','rb'))
        kk=ddict.keys()[0]
        ggf,gnf,gvis,mvis,mcounts,msp,mt1=ddict[kk]
    #otherwise it's a flare list object
    else:
        ggf=fl.isolate_att('GOES_GOES_long')
        gnf=fl.isolate_att('Nittaflux')
        gAvis=fl.isolate_att('GAobs_with_gloc')
        gBvis=fl.isolate_att('GBobs_with_gloc')
        mvis=fl.isolate_att('MSvis')
        mcounts=fl.isolate_att('Messenger_GOES_long')
        msp=fl.isolate_att('Messenger_GOES_selfproxy')
        mt1=fl.isolate_att('Messenger_peak')
        tt=fl.isolate_att('T1')

    #nittaflux v. GOES
    set1=[gnf,ggf]

    #all visible flares
    s2x,s2y=[],[]
    for g,n,a,b in zip(ggf,gnf,gAvis,gBvis):
        if a==True or b == True:
            s2x.append(n)
            s2y.append(g)
    overlay2(set1,[s2x,s2y],dashed=False, nfline=True,save='GN_vis_recalc.png')
    print stats.spearmanr(s2x,s2y)

#     #all visible flares, with Nitta counts condition
#     s2x,s2y=[],[]
#     for g,n,v in zip(ggf,gnf,gvis):
#         if v==True and n >=4e6:
#             s2x.append(n)
#             s2y.append(g)
#     fdatax=np.log10(np.array(s2x))
#     fdatay=np.log10(np.array(s2y))
#     line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
#     print line
#     overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=False, nfline=True,save='GN_vis_Ncounts.png')
#     print stats.spearmanr(s2x,s2y)

#     #all visible above C5.5, with Nitta counts condition > 4e6
#     s2x,s2y=[],[]
#     for g,n,v in zip(ggf,gnf,gvis):
#         if v==True and g >=5.5e-6 and n >=1e5:
#             s2x.append(n)
#             s2y.append(g)
#     #fit line
#     fdatax=np.log10(np.array(s2x))
#     fdatay=np.log10(np.array(s2y))
#     line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
#     print line
#     overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=False, nfline=True,save='GN_vis_gtC5p5_Ncounts.png')
#     print stats.spearmanr(s2x,s2y)

#     #all visible above C5.5
#     s2x,s2y=[],[]
#     for g,n,v in zip(ggf,gnf,gvis):
#         if v==True and g >=5.5e-6:
#             s2x.append(n)
#             s2y.append(g)
#     #fit line
#     fdatax=np.log10(np.array(s2x))
#     fdatay=np.log10(np.array(s2y))
#     #line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
#     #print line
#     overlay2(set1,[s2x,s2y],dashed=False, nfline=True,save='GN_vis_gtC5p5.png')
#     print stats.spearmanr(s2x,s2y)

#     #all visible above M, with Nitta counts condition > 4e6
#     s2x,s2y=[],[]
#     for g,n,v in zip(ggf,gnf,gvis):
#         if v==True and g >=1e-5 and n >=1e5:
#             s2x.append(n)
#             s2y.append(g)
#     #fit line
#     fdatax=np.log10(np.array(s2x))
#     fdatay=np.log10(np.array(s2y))
#     line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
#     print line
#     overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=False, nfline=True,save='GN_vis_gtM_Ncounts.png')
#     print stats.spearmanr(s2x,s2y)

#     #all visible above M
#     s2x,s2y=[],[]
#     for g,n,v in zip(ggf,gnf,gvis):
#         if v==True and g >=1e-5:
#             s2x.append(n)
#             s2y.append(g)
#     #fit line
#     fdatax=np.log10(np.array(s2x))
#     fdatay=np.log10(np.array(s2y))
#     #line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
#     #print line
#     overlay2(set1,[s2x,s2y],dashed=False, nfline=True,save='GN_vis_gtM.png')
#     print stats.spearmanr(s2x,s2y)


def make_messenger_plots(fl=False):
    from scipy import stats
    if not fl:
        ddict=pickle.load(open('plot_params.p','rb'))
        kk=ddict.keys()[0]
        ggf,gnf,gvis,mvis,mcounts,msp,mt1=ddict[kk]
    #otherwise it's a flare list object
    else:
        ggf=fl.isolate_att('GOES_GOES_long')
        gnf=fl.isolate_att('Nittaflux')
        gvis=fl.isolate_att('MAobs_with_gloc')
        mvis=fl.isolate_att('MBobs_with_gloc')
        mcounts=fl.isolate_att('Messenger_GOES_long')
        msp=fl.isolate_att('Messenger_GOES_selfproxy')
        mt1=fl.isolate_att('Messenger_peak')


    for i,m in enumerate(mcounts):
        if m == '':
          mcounts[i]=0.0
    #nittaflux v. Messenger counts
    set1=[gnf,mcounts]

    #all visible flares
    s2x,s2y=[],[]
    for m,n,v,v2 in zip(mcounts,gnf,mvis,gvis):
        if v==True or v2 == True and m !='' and not np.isnan(n):
            s2x.append(n)
            s2y.append(m)
    overlay2(set1,[s2x,s2y],dashed=False,yran=[1,1e5],ylab='Messenger Peak (Counts/s)',save='MN_vis_recalc.png',title='no restrictions')
    print stats.spearmanr(s2x,s2y)

#     #all visible flares, with Nitta counts condition
#     s2x,s2y=[],[]
#     for m,n,v in zip(mcounts,gnf,gvis):
#         if v==True and m != '' and n >=4e6:
#             s2x.append(n)
#             s2y.append(m)
#     fdatax=np.log10(np.array(s2x))
#     fdatay=np.log10(np.array(s2y))
#     line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
#     print line
#     overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=False,yran=[1,1e5],ylab='Messenger Peak (Counts/s)',save='MN_vis_recalc_Ncounts.png',title='nitta counts > 4e6')
#     print stats.spearmanr(s2x,s2y)


#     #Nitta proxy vs. Messenger proxy
#     set1=[list(1.39*10**-11*np.array(gnf)),msp]

#     #all visible above C5.5, with Nitta counts condition > 4e6
#     s2x,s2y=[],[]
#     for m,n,v in zip(msp,gnf,gvis):
#         if v==True and m!= '' and m >=5.5e-6 and n >=4e6:
#             s2x.append(1.39*10**-11*n)
#             s2y.append(m)
#     #fit line
#     fdatax=np.log10(np.array(s2x))
#     fdatay=np.log10(np.array(s2y))
#     line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
#     print line
#     overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=True, xran=[1e-7,1e-2],xlab='Nitta Proxy ($\\rm{Wm^{-2}}$)',yran=[1e-7,1e-2],ylab='Messenger Proxy ($\\rm{Wm^{-2}}$)',save='MN_vis_gtC5p5_Ncounts.png',title='>C5.5, Nitta counts > 4e6')
#     print stats.spearmanr(s2x,s2y)

#     #all visible above C5.5
#     s2x,s2y=[],[]
#     for m,n,v in zip(msp,gnf,gvis):
#         if v==True and m!='' and m >=5.5e-6 and not np.isnan(n):
#             s2x.append(1.39*10**-11*n)
#             s2y.append(m)
#     #fit line
#     fdatax=np.log10(np.array(s2x))
#     fdatay=np.log10(np.array(s2y))
#     #line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
#     #print line
#     overlay2(set1,[s2x,s2y],dashed=True,xran=[1e-7,1e-2],xlab='Nitta Proxy ($\\rm{Wm^{-2}}$)',yran=[1e-7,1e-2],ylab='Messenger Proxy ($\\rm{Wm^{-2}}$)',save='MN_vis_gtC5p5.png',title='>C5.5, jointly visible STEREO/GOES')
#     print stats.spearmanr(s2x,s2y)

#     #all visible above M, with Nitta counts condition > 4e6
#     s2x,s2y=[],[]
#     for m,n,v in zip(msp,gnf,gvis):
#         if v==True and m!='' and m >=1e-5 and n >=4e6:
#             s2x.append(1.39*10**-11*n)
#             s2y.append(m)
#     #fit line
#     fdatax=np.log10(np.array(s2x))
#     fdatay=np.log10(np.array(s2y))
#     line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
#     print line
#     overlay2(set1,[s2x,s2y],fitline=[line[0],line[1]],dashed=True,xran=[1e-7,1e-2],xlab='Nitta Proxy ($\\rm{Wm^{-2}}$)', yran=[1e-7,1e-2],ylab='Messenger Proxy ($\\rm{Wm^{-2}}$)',save='MN_vis_gtM_Ncounts.png',title='GSvis, >M, counts >4e6')
#     print stats.spearmanr(s2x,s2y)

#     #all visible above M
#     s2x,s2y=[],[]
#     for m,n,v in zip(msp,gnf,gvis):
#         if v==True and m!='' and m >=1e-5 and not np.isnan(n):
#             s2x.append(1.39*10**-11*n)
#             s2y.append(m)
#     #fit line
#     fdatax=np.log10(np.array(s2x))
#     fdatay=np.log10(np.array(s2y))
#     #line,res,_,_,_=np.polyfit(fdatax,fdatay,1,full=True)
#     #print line
#     overlay2(set1,[s2x,s2y],dashed=True,xran=[1e-7,1e-2],xlab='Nitta Proxy ($\\rm{Wm^{-2}}$)',yran=[1e-7,1e-2],ylab='Messenger Proxy ($\\rm{Wm^{-2}}$)',save='MN_vis_gtM.png',title='GSvis, >M')
#     print stats.spearmanr(s2x,s2y)


def make_messenger_Tplots(fl=False):
    if not fl:
        ddict=pickle.load(open('plot_params.p','rb'))
        kk=ddict.keys()[0]
        ggf,gnf,gvis,mvis,mcounts,msp,mt1=ddict[kk]
    #otherwise it's a flare list object
    else:
        ggf=fl.isolate_att('GOES_GOES_long')
        gnf=fl.isolate_att('Nittaflux')
        gvis=fl.isolate_att('GSobs_with_gloc')
        mvis=fl.isolate_att('MSvis')
        mcounts=fl.isolate_att('Messenger_GOES_long')
        msp=fl.isolate_att('Messenger_GOES_selfproxy')
        mt1=fl.isolate_att('T1')
    tt=11.6045*np.array(mt1)
    #Nitta proxy vs. Messenger proxy
    set1=[list(1.39*10**-11*np.array(gnf)),msp]

    #all visible with T>20MK, Nitta > M
    s2x,s2y=[],[]
    for m,n,v,t in zip(msp,gnf,gvis,tt):
        if v==True and m!= '' and t >=20. and np.isfinite(n) and n !=0.0:# and 1.39e-11*n>1e-5: #and 13.9e-11*n <1e-3:
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
    if not minx and not maxx:
          fdatax=datax
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
    return line#x,yfit(x),res,fdatax,fdatay
