 #######################################
# OCFlareList.py
# Erica Lastufka 17/5/17 

#Description: Class of functions to deal with all data,files,calculations of the study
#######################################

import numpy as np
import os
import pickle
import glob
import OCFlare
import pandas as pd


class OCFlareList(object):
    ''' This object will have the following attributes:
            List/Series of OCFlare objects
            Notes
        Methods will include:
            All selection/filtering/sorting methods
            All the plotting methods
            All statistics methods
     '''
    def __init__(self, IDs=True, folder=False, freload=False,calc_times=False, calc_missing=False,gen=False,filename=False):
        '''Options to initialize from list of flare IDs, file with column of flare IDs, or directory with object pickle files in them?'''
        self.list=[]
        if IDs and type(IDs) == list:
            #initialize from list of flare IDs
            for i in IDs:
                if os.path.isfile(str(i)+'OCFlare.p') and not freload:
                    self.list.append(pickle.load(open(str(i)+'OCFlare.p','rb')))
                else:
                    #initialize all the other objects from the default legacy files?
                    self.list.append(OCFlare.OCFlare(i,calc_times=calc_times, calc_missing=calc_missing,gen=gen,filename=filename))
        elif IDs and type(IDs) == str:
            #assume it's a csv or sav file that has a column of flare ids in it
            import data_management2 as d
            #with d.OCData(IDs) as o: #sadly doesn't work, let's do manual cleanup then:
            o = d.OCData(IDs)
            for i in o.ID:
                self.list.append(OCFlare.OCFlare(i,calc_times=calc_times, calc_missing=calc_missing,gen=gen, legacy=True, filename=IDs))
            del o

        if folder and not freload:
            #restore all the Flare pickle files and combine into list ... need to revise this now that they are not objects!
            os.chdir(folder)
            pickles=sorted(glob.glob('*.p')) #sorted by name... may want to change that, or change my ID naming convention
            for p in pickles:
                self.list.append(pickle.load(open(str(i)+'OCFlare.p','rb')))
            os.chdir('/Users/wheatley/Documents/Solar/occulted_flares')

        self.Notes=['']*len(self.list) #list of empty strings with the same length as flare list
        #convert into Pandas series
        self.list=pd.Series(self.list, name='Flares')
        self.Notes=pd.Series(self.Notes, name='Notes')        

    def write(self,picklename=False):
        '''Pickle flare list object'''
        import pickle
        if not picklename:
            picklename=raw_input('What do you want to name your file? ')
        pickle.dump(self, open(picklename, 'wb'))
        print picklename + ' was created in '+ os.getcwd()

    def writeIDs(self,picklename=False):
        '''Pickle flare list IDs, to make initiation faster'''
        import pickle
        if not picklename:
            picklename=raw_input('What do you want to name your file? ')
        pickle.dump(self.isolate_att('ID'), open(picklename, 'wb'))
        print picklename + ' was created in '+ os.getcwd()

    def __iter__(self):
        '''Returns a generator that iterates over the object'''
        for attr, value in self.__dict__.iteritems():
            yield attr, value

    def select_outliers(self, angle=90.,threshold=10.): #greater than given angle... need to fix bug/feature with ratio
        '''Select certain flares from the list to be examined in visually'''
        indices=[]
        #newlist=copy.deepcopy(self)
        for i,flare in enumerate(self.ID):
            if np.abs(self.list[i].Properties.ratio) > threshold and self.list[i].Observation.Angle < angle:
                indices.append(i)
                newlist.append(self.list[i])    
        return newlist #this won't have all the plot methods attached to it... maybe I should let an object be initialized from list of objects and just tack on the methods

    def isolate_att(self,att,key=False):
        '''Make a list of all the values of the given attribute. Either input name of attribute, or name of dictionary and keyword ie loops['length1']'''
        alist=[]
        if hasattr(self.list[0].Datetimes,att):
            for i in np.arange(0,len(self.list)): alist.append(getattr(self.list[i].Datetimes,att))
        elif hasattr(self.list[0].Properties,att):
            for i in np.arange(0,len(self.list)):
                if not key:
                    alist.append(getattr(self.list[i].Properties,att))
                else:
                    alist.append(getattr(self.list[i].Properties,att)[key])
        elif hasattr(self.list[0].Files,att):
            for i in np.arange(0,len(self.list)):
                if not key:
                    alist.append(getattr(self.list[i].Files,att))
                else:
                    alist.append(getattr(self.list[i].Files,att)[key])                    
        elif hasattr(self.list[0].Observation,att):
            for i in np.arange(0,len(self.list)):
                if not key:
                    alist.append(getattr(self.list[i].Observation,att))
                else:
                    alist.append(getattr(self.list[i].Observation,att)[key])
        elif att=='ID':
            for i in np.arange(0,len(self.list)):
                alist.append(getattr(self.list[i],'ID'))
        return alist

    def att2df(self, flattened=False):
        '''Store the flare attributes in one large, flattened, data frame instead of several small ones'''
        from pandas.io.json import json_normalize
        dfs=[]
        for f in self.list: #all the DataFrames
            dfs.append(f.att2df())
                   #current problem is that each dictionary is stored in the same index in its dataframe. So concatenating the dataframes gives the wrong indexing. We basically need to build a new, re-formatted data frame out of the old ones.
        nid=[df['ID'][0] for df in dfs]
        ndt=[df['Datetimes'][1] for df in dfs]
        nf=[df['Files'][2] for df in dfs]
        npp=[df['Properties'][3] for df in dfs]
        no=[df['Observation'][4] for df in dfs]
        nn=[df['Notes'][5] for df in dfs]
        
        dfid=pd.DataFrame(nid,index=np.arange(len(self.list)),columns=['ID']).to_dict(orient='series')
        dfdt=pd.DataFrame(ndt,index=np.arange(len(self.list)),columns=ndt[0].keys()).to_dict(orient='series')
        dff=pd.DataFrame(nf,index=np.arange(len(self.list)),columns=nf[0].keys()).to_dict(orient='series')
        dfp=pd.DataFrame(npp,index=np.arange(len(self.list)),columns=npp[0].keys()).to_dict(orient='series')
        dfo=pd.DataFrame(no,index=np.arange(len(self.list)),columns=no[0].keys()).to_dict(orient='series')
        dfn=pd.DataFrame(nn,index=np.arange(len(self.list)),columns=['Notes']).to_dict(orient='series')

        onedict={'ID':dfid,'Datetimes':dfdt,'Files':dff,'Properties':dfp,'Observation':dfo,'Notes':dfn}
        oneframe=pd.DataFrame.from_dict(onedict)#dfid.join(dfdt.join(dff.join(dfp.join(dfo.join(dfn)))))

        if flattened: #flatten it because json won't cooperate
            flatdict={}
            for tk in onedict.keys(): #top level of keys
                if len(onedict[tk]) ==1: #notes or ID
                    flatdict[tk]=onedict[tk][tk].tolist()
                else:
                    for sk in onedict[tk].keys():
                        if type(onedict[tk][sk][0]) != dict:
                            flatdict[tk+'.'+sk]=onedict[tk][sk].tolist()
                        else: #flatten even more
                            for ssk in onedict[tk][sk][0].keys():
                                flatdict[tk+'.'+sk+'.'+ssk]=[onedict[tk][sk][i][ssk] for i in range(0,len(onedict[tk][sk]))]                                
            return flatdict
        else:
            return oneframe
    
    def export2csv(self, filename):
        '''Exports flare list attributes to csv file'''
        from pandas.io.json import json_normalize
        onenorm= self.att2df(flattened=True)
            
        os.chdir(self.list[0].Files.dir+self.list[0].Files.folders['flare_lists'])
        
        import csv
        with open(filename,'wb') as f:
            writer=csv.writer(f)
            writer.writerow(onenorm.keys())
            for i in range(0,len(onenorm['ID'])):
                writer.writerow([onenorm[k][i] for k in onenorm.keys()])        
        os.chdir(self.list[0].Files.dir)
        #return onenorm

    def update(self):
        '''re-write attribute pickles for each flare in list'''
        for f in self.list:
            f.export2pickle()         

####################################################  PLOT METHODS  ##################################################################

    def plot_angle_distribution(self,ymax=10):
        '''make a histogram of the angle distributions'''
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        angle=isolate_att('Angle')
        n, bins, patches = plt.hist(angle, 18, facecolor='green', alpha=0.75)
    
        plt.xlabel('Angle between Earth and Mercury (degrees)')
        plt.ylabel('Number of Events')
        ax1.set_ylim([0,ymax])
        ax1.set_xlim([0,180])

        ax1.plot()

        fig.show()
    
    def plot_goes_ratio(self,title= "",ymin=0, ymax=3, labels=1,ylog=True, scatter = True,cc='GOES',save=False,show=True):
        '''make a plot of the GOEs ratio vs. angle'''
        import matplotlib.patches as mpatches
        mc=self.isolate_att("Messenger_GOES")
        gc=self.isolate_att("GOES_GOES")
        angle=self.isolate_att("Angle")
        chisq=self.isolate_att("chisq")
        ids=self.isolate_att("ID")
        dts=self.isolate_att("chisq")
        ylabel='Messenger_GOES/Observed_GOES'
        colors,labelang,labelratio,coordlabel=[],[],[],[]
        for mval,gval,cs,ID,dt,ang in zip(mc,gc,chisq,ids,dts,angle):
            if cc=='GOES':
                #print np.rint(-np.log10(rval))
                if np.rint(-np.log10(gval)) <= 4.0:colors.append('r')
                elif np.rint(-np.log10(gval)) == 5.0:colors.append('m')
                elif np.rint(-np.log10(gval)) == 6.0:colors.append('y')
                elif np.rint(-np.log10(gval)) == 7.0:colors.append('g')
                elif np.rint(-np.log10(gval)) == 8.0:colors.append('b')
                elif np.rint(-np.log10(gval)) == 9.0:colors.append('k')
            else: #color code by other one
                if np.rint(-np.log10(mval)) <= 4.0:colors.append('r')
                elif np.rint(-np.log10(mval)) == 5.0:colors.append('m')
                elif np.rint(-np.log10(mval)) == 6.0:colors.append('y')
                elif np.rint(-np.log10(mval)) == 7.0:colors.append('g')
                elif np.rint(-np.log10(mval)) == 8.0:colors.append('b')
                elif np.rint(-np.log10(mval)) == 9.0:colors.append('k')

            labelang.append(ang) #?
            labelratio.append(mval/gval) #?
            if labels==1:
                coordlabel.append(ID)
            else: #0 or 2
                coordlabel.append(datetime.strftime(dt,'%D %H:%M'))
            if scatter:
                delta = 50
            elif cs == '':#notes column is empty
                delta.append(5000*10*2**np.rint(np.log10(np.abs(mval-gval)))) #difference in size between the GOES classes
            else: #notes carries chisq value
                delta.append(50*10*2**np.rint(np.log10(float(cs))))

        ratio = np.array(mc/gc)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.scatter(np.array(angle), ratio, s=delta, c=colors,alpha=.75)
        ax1.axhline(y=1,linestyle='dashed',color='k')
        
        if labels != 0: 
            for x,y,t in zip(np.array(labelang),labelratio,coordlabel):
                #print x,y,t
                ax1.annotate('%s' % t, xy=(x,y), textcoords='data')
            plt.grid()

   
        plt.xlabel('Angle between Earth and Mercury (degrees)')
        plt.ylabel(ylabel)

        if ylog:
            ax1.set_yscale('log')
        ax1.set_ylim([ymin,ymax])
        ax1.set_xlim([0,180])
        plt.title(title)
        X = mpatches.Patch(color='red', label='X')
        M = mpatches.Patch(color='magenta', label='M')
        C = mpatches.Patch(color='yellow', label='C')
        B = mpatches.Patch(color='green', label='B')
        A= mpatches.Patch(color='blue', label='A')
        K= mpatches.Patch(color='black', label='<A')
        ax1.legend(handles=[X,M,C,B,A,K],loc='upper left',fontsize='medium')

        ax1.plot()
        if save:
            plt.savefig(save)
        if show:
            fig.show()
            
    def plot_goes_messenger(self,title= "",ymin=10**-9, ymax=10**-3, minangle=False,maxangle=False,labels=0,loglog=True, loc=True, scatter = True,cc='GOES',save=False,show=True, channel='short', line=False, linexlim=False):
        '''make a plot of the Messenger_GOES vs GOES_GOES'''
        import matplotlib.patches as mpatches
        from datetime import datetime

        mc=self.isolate_att("Messenger_GOES")
        if channel == 'short':
            gc=self.isolate_att("GOES_GOES")
        else:
            gc=self.isolate_att("GOES_GOES_long")
            
        angles=self.isolate_att("Angle")
        if loc:
            location=self.isolate_att("source_pos")
        else:
            location=np.arange(len(mc)) #just an array of meaningless numbers
        if labels !=2:
            ids=self.isolate_att("ID")
        else:
            ids=self.isolate_att("Messenger_peak")

        gcgt0,mcgt0=[],[]
        colors,coordlabel=[],[]
        classes=['m','y','g','b']

        def color_code(ang,ulimit=180,llimit=0,nsegments=3):
            arange= ulimit - llimit
            #print arange,llimit+arange/nsegments,llimit+2*arange/nsegments
            if ang >= llimit and ang <= llimit+arange/nsegments:return 'b'
            elif ang <= llimit+2*arange/nsegments: return 'r'
            else : return 'g'
        
        for mval,gval,ID,ang,l in zip(mc,gc,ids,angles,location): 
            app=False
            if not np.isnan(gval) and gval >1e-09 and not np.isnan(mval) and mval > 0:
                if not loc and not minangle and not maxangle:app= True
                elif loc and minangle and l != '[0,0]' and l != '' and type(l) == str and ang > minangle: app = True
                elif loc and maxangle and l != '[0,0]' and l != '' and type(l) == str and ang < maxangle: app = True
                elif loc and not minangle and not maxangle and l != '[0,0]' and l != '' and type(l) == str: app = True
                    
                try:
                    if app: #then append, else go on to the next iteration
                        gcgt0.append(gval)
                        mcgt0.append(mval)
                        #print app, ang,gval,mval,l

                        #now determine details of color coding etc
                        if cc=='GOES':
                            idx=-int(np.floor(np.log10(gval))) # int from <9 to >4
                            if idx <=4: colors.append('r')
                            elif idx >=9: colors.append('k')
                            else: colors.append(classes[idx-5])
                        elif cc=='Messenger':
                            idx=-int(np.floor(np.log10(mval))) # int from <9 to >4
                            if idx <=4: colors.append('r')
                            elif idx >=9: colors.append('k')
                            else: colors.append(classes[idx-5])
                        elif cc == 'angle':
                            if minangle: #color code by angle
                                colors.append(color_code(ang,llimit=minangle))
                            elif maxangle: #color code by angle
                                colors.append(color_code(ang,ulimit=maxangle))
                            else: #color code by angle divided into 3 segments. Should do something about the legend labels too...
                                colors.append(color_code(ang))
                        
                        if labels == 2:
                            coordlabel.append(datetime.strftime(ID,'%D %H:%M'))
                        else: #0 or 2
                            coordlabel.append(ID)
                            
                except UnboundLocalError: continue

        def line_of_best_fit(datax,datay,loglog= True,minx=False,maxx=False,minratio=False,maxratio=False):
            #if list convert to np.array
            fdatax,fdatay=[],[]
            if type(datax) == list: datax=np.array(datax)
            if type(datay) == list: datay=np.array(datay)
            if not minx and not maxx: fdatax=datax
            #if not miny and not maxy: fdatay=datay

            for x,y in zip(datax,datay):
                if minx and x > minx and maxx and x < maxx:
                    fdatax.append(x)
                    fdatay.append(y)
                #if maxx and x < maxx:
                #    fdatax.append(x)
                #    fdatay.append(y)
                    #for y in datay:
            #    if miny and y > miny:
            #        fdatay.append(y)
            #    if maxy and y < maxy:
            #        fdatay.append(y)
            if type(fdatax) == list: fdatax=np.array(fdatax)
            if type(fdatay) == list: fdatay=np.array(fdatay)
            
            if loglog:
                fdatax=np.log10(fdatax)
                fdatay=np.log10(fdatay)
        
            #ratio=datay/datax deal with this later
            print np.shape(fdatax),np.shape(fdatay)
            line=np.polyfit(fdatax,fdatay,1)
            print 'line of best fit: 10**('+str(line[0])+'*log10(x) + ' +str(line[1]) +')'
            x=10**(fdatax) #assume log
            yfit = lambda x: 10**(line[0]*np.log10(x)+line[1])
            return x,yfit(x),fdatax,fdatay

        if line:
            lx,ly,fdatax,fdatay = line_of_best_fit(gcgt0,mcgt0,minx=linexlim[0],maxx=linexlim[1])
            #print l
                    
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.scatter(gcgt0,mcgt0, s=50, c=colors,alpha=.75) #for some reason the masking is not being effective!
        ax1.plot(np.arange(10),np.arange(10),linestyle='dashed',color='k')
        if line:
            ax1.plot(lx,ly, linestyle='dotted',color='darkorange',linewidth=2)

        if labels != 0: 
            for x,y,t in zip(gcgt0,mcgt0,coordlabel):
                #print x,y,t
                ax1.annotate('%s' % t, xy=(x,y), textcoords='data')
            plt.grid()

   
        plt.xlabel('GOES Flux, W m$^{-2}$')
        plt.ylabel('Messenger Flux, W m$^{-2}$')

        if loglog:
            ax1.set_yscale('log')
            ax1.set_xscale('log', nonposx='clip')
            
        ax1.set_aspect('equal')
        ax1.set_ylim([ymin,ymax])
        ax1.set_xlim([ymin,ymax]) #let's keep everything on the same scale...
        plt.title(title)
        if cc != 'angle':
            X = mpatches.Patch(color='red', label='X')
            M = mpatches.Patch(color='magenta', label='M')
            C = mpatches.Patch(color='yellow', label='C')
            B = mpatches.Patch(color='green', label='B')
            A= mpatches.Patch(color='blue', label='A')
            K= mpatches.Patch(color='black', label='<A')
            handles=[X,M,C,B,A,K]
        else:
            if minangle:
                t1=str(minangle+(int(minangle))/3)
                t2=str(minangle+(2*int(minangle))/3)
                mxa='180'
                mna=str(minangle)
            elif maxangle:
                t1=str((int(maxangle))/3)
                t2=str((2*int(maxangle))/3)
                mxa=str(maxangle)
                mna='0'                
            else:
                t1='60'
                t2='120'
                mna='0'
                mxa='180'

            B= mpatches.Patch(color='blue', label='$' +mna+' \leq \\theta \leq '+t1+'\degree$')
            R = mpatches.Patch(color='red', label='$'+t1+'\degree \leq \\theta \leq '+t2+'\degree$')
            G = mpatches.Patch(color='green', label='$'+t2+'\degree \leq \\theta \leq'+mxa+'\degree$')
            handles=[B,R,G]
           
        ax1.legend(handles=handles,loc='lower right',fontsize='medium')

        ax1.plot()
        if save:
            plt.savefig(save)
        if show:
            fig.show()
        return fdatax,fdatay

    def hist_ratio(self,title='All flares',gc='all', save=False,show=True):
        '''Make histogram to see distribution of ratios by goes class'''
        ratio=get_ratio(flare_list)
        ratio = ratio[np.where(ratio != -1)]
        #print ratio
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        goes=np.array(convert_goes2flux(flare_list.Flare_properties['GOES_GOES']))
        #print goes
    
        if gc == 'all':
            n, bins, patches = plt.hist(ratio, np.logspace(-3,3,12), facecolor='orange', alpha=0.75)
        if gc == 'A':
            gt,lt =np.where(goes > 10**-8),np.where(goes < 10**-7)
            all=set(gt[0]) & set(lt[0])
            ratio= ratio[list(all)]
            n, bins, patches = plt.hist(ratio, np.logspace(-3,3,12), facecolor='blue', alpha=0.75)
        if gc == 'B':
            gt,lt =np.where(goes > 10**-7),np.where(goes < 10**-6)
            all=set(gt[0]) & set(lt[0])
            ratio= ratio[list(all)]
            n, bins, patches = plt.hist(ratio,np.logspace(-3,3,12) , facecolor='green', alpha=0.75)
        if gc == 'C':
            gt,lt =np.where(goes > 10**-6),np.where(goes < 10**-5)
            all=set(gt[0]) & set(lt[0])
            ratio= ratio[list(all)]
            n, bins, patches = plt.hist(ratio, np.logspace(-3,3,12), facecolor='yellow', alpha=0.75)
        if gc == 'M':
            gt,lt =np.where(goes > 10**-5),np.where(goes < 10**-4)
            all=set(gt[0]) & set(lt[0])
            ratio= ratio[list(all)]
            n, bins, patches = plt.hist(ratio, np.logspace(-3,3,12), facecolor='magenta', alpha=0.75)
        if gc == 'overlay':
            gt,lt =np.where(goes > 10**-8),np.where(goes < 10**-7)
            all=set(gt[0]) & set(lt[0])
            ratioA= ratio[list(all)]
            plt.hist(ratioA, np.logspace(-3,3,12), facecolor='blue', alpha=0.6)
            gt,lt =np.where(goes > 10**-7),np.where(goes < 10**-6)
            all=set(gt[0]) & set(lt[0])
            ratioB= ratio[list(all)]
            plt.hist(ratioB,np.logspace(-3,3,12) , facecolor='green', alpha=0.6)
            gt,lt =np.where(goes > 10**-6),np.where(goes < 10**-5)
            all=set(gt[0]) & set(lt[0])
            ratioC= ratio[list(all)]
            plt.hist(ratioC, np.logspace(-3,3,12), facecolor='yellow', alpha=0.6)
            gt,lt =np.where(goes > 10**-5),np.where(goes < 10**-4)
            all=set(gt[0]) & set(lt[0])
            ratioM= ratio[list(all)]
            plt.hist(ratioM, np.logspace(-3,3,12), facecolor='magenta', alpha=0.6)
        
            M = mpatches.Patch(color='magenta', label='M')
            C = mpatches.Patch(color='yellow', label='C')
            B = mpatches.Patch(color='green', label='B')
            A= mpatches.Patch(color='blue', label='A')
            ax1.legend(handles=[M,C,B,A],loc='upper left',fontsize='medium')
        
    
        plt.xlabel('Ratio between Messenger GOES and actual GOES')
        plt.ylabel('Number of Events')
        plt.gca().set_xscale("log")
        plt.title(title)
        ax1.set_ylim([0,100])
        #ax1.set_xlim([0,150])

        ax1.plot()
        if save:
            plt.savefig(save)
        if show:
            fig.show()

    def plot_all_stereo(self,AIAmap=False,maps=False,oplotlimb=True,oplotpoint=False,save=True,subplots=True,zoom=False,diff=False,return_limb=False,n=4):
        '''plot all stereo plots for the flares in the list'''
        limbcoords=[]
        for flare in self.list:
            if not return_limb:
                flare.plot_stereo(AIAmap=AIAmap,maps=maps,oplotlimb=oplotlimb,oplotpoint=oplotpoint,save=save,subplots=subplots,zoom=zoom,diff=diff,return_limb=return_limb,n=n)
            else:
                limbcoords.append(flare.plot_stereo(AIAmap=AIAmap,maps=maps,oplotlimb=oplotlimb,oplotpoint=oplotpoint,save=save,subplots=subplots,zoom=zoom,diff=diff,return_limb=return_limb,n=n))
        if return_limb:
            return limbcoords
        
    def plot_all_aia(self,AIAmap=False,zoom=False,save=True, wave='304',all_wave=False):
        '''plot all aia plots for the flares in the list'''
        for flare in self.list:
            flare.plot_aia(AIAmap=AIAmap,zoom=zoom,save=save, wave=wave,all_wave=all_wave)

    def plot_all_rhessi_aia(self,AIAmap=False,erange=False, wave='193',tag=False):
        '''plot all aia plots for the flares in the list'''
        for flare in self.list:
            flare.plot_RHESSI_AIA_IDL(erange=[[4,9],[12,18]],tag='lowE')  
            flare.plot_RHESSI_AIA_IDL(erange=[[4,9],[12,18],[18,30]],tag='3E')      
            flare.plot_RHESSI_AIA_IDL(erange=[[4,9],[12,18],[18,30],[30,80]],tag='4E')

    def plot_all_lightcurves(self,AIAmap=True,maps=True,oplotlimb=False,oplotpoint=False,save=False,subplots=False,zoom=False,diff=False,return_limb=True):
        '''plot all lightcurves for the flares in the list'''
        for flare in self.list:
            flare.get_closest_limbcoord(limbcoords=limbcoords,redef=redef,plot=plot,smap=smap)

    def plot_all_spectrograms(self,filename,ivec):
        '''plot all spectrograms for the flares in the list'''
        for flare,i in zip(self.list,ivec):
            flare.sunpy_spectrogram(filename,i)

    def get_all_limbcoord(self,limbcoords=False,redef=False,plot=True,smap=False):
        '''get closest coord on the limb for the flares in the list'''
        for flare in self.list:
            flare.get_closest_limbcoord(limbcoords=limbcoords,redef=redef,plot=plot,smap=smap)
       
