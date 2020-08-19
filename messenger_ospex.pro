;pro messenger_ospex

;first find the maximum of a given spectrum
function find_messenger_peak,obj=obj, data=data, time_interval= time_interval,erange = erange
    eaxis=obj->getaxis(/ct_energy)
    taxis=obj->getaxis(/mean) ;are these formatted from a different date?
    elist=where( (eaxis ge min(erange)) AND (eaxis le max(erange)) )
    phigh=total(data.data(elist,*),1)
    ;print, anytim(min(taxis),/vms),anytim(max(taxis),/vms),anytim(time_interval[0],/vms),anytim(time_interval[1],/vms)
    tlist=where( (taxis ge time_interval[0]) AND (taxis le time_interval[1]))
    phigh=phigh(tlist)
    this_max=max(phigh,this_max_ind)
    this_tr1=taxis(tlist(this_max_ind))+[-8.,8.]
    this_tr=[max(this_tr1(0)),min(this_tr1(1))]
    return, this_tr
end

pro check_messenger_peak, flare_list,i
    ;plot the energy range and time interval just to check if the max is
                                ;where it should be how do I plot
                                ;things in ospex again???
    search_network,/enable
    erange=[3.2,8.5]
    t=flare_list.datetimes.Messenger_datetimes[i]
    time_interval=[anytim(t),anytim(t)+2400.]
    obj=ospex(/no_gui)
    obj-> set, spex_specfile= 'data/dat_files/' + strtrim(flare_list.data_properties.xrs_files[i],1)+'.dat'
    data=obj->getdata(class='spex_data')
    eaxis=obj->getaxis(/ct_energy)
    taxis=obj->getaxis(/mean) ;are these formatted from a different date?
    elist=where( (eaxis ge min(erange)) AND (eaxis le max(erange)) )
    phigh = total(data.data(elist,*),1)
    ;tlist=where( (taxis ge time_interval[0]) AND (taxis le time_interval[1]))
    ;mti=size(tlist,/dim)
    ;newtaxis = [taxis[tlist[0]-1],taxis[tlist]]
    utplot,taxis, phigh,timerange=[time_interval[0],time_interval[1]],yrange=[min(phigh),max(phigh)]
    print, anytim(time_interval,/vms),max(phigh),'     ',flare_list.datetimes.obs_start_time[i],'     ',flare_list.datetimes.obs_end_time[i]
    ;outfilename='data/'+t+'.png'
                                ;write_png,outfilename,tvrd(/true)
    obj_destroy,obj
end
  
;do the fit for a given spex file ....
function messenger_ospex, obj=o,time_interval= time_interval,erange = erange,mask=mask,comp_params=comp_params,outfilename=filename,quiet=quiet
    o-> set, spex_fit_time_interval= time_interval
    o-> set, spex_eband= [[1.06958, 1.65593], [1.65593, 2.56371], [2.56371, 3.96913], $    
                            [3.96913, 6.14501], [6.14501, 9.51370]] ;is this ok?
    ;o->set,spex_bk_time_interval=this_bkg ;do I need to do this for each one? how to automate...
    o->set, fit_function='vth+vth' ;2 thermal components 'vth+bpow' -> thermal + broken power law
    o-> set, fit_comp_params= comp_params;[1.0, 2.0, 1.0, 0.0, 0.3, 1.0] ;em, T, abund for each component
    o-> set, fit_comp_free_mask= mask ;fit high-energy component first. don't fit abundances
    o->set,spex_fit_manual=0
    set_logenv, 'OSPEX_NOINTERACTIVE', '1' ;I really want to be able to use my mouse!!!
    o->set, spex_erange=erange  ;[2.3,8]
    if quiet eq 'True' then o->set,spex_autoplot_enable=0,spex_fit_progbar=0,spex_fitcomp_plot_resid=0 ;surpress the GUI
    ;endif
    o->dofit,/all
    o->savefit,outfile=filename+'high_a.fits'
    o->set,spex_erange = [1.5,8.5] ;expand this for low energy fit
    o-> set, fit_comp_free_mask= [0,0,0,1,1,0] ;now fit low-energy compmonent
    o-> set, fit_comp_free_mask= [0,0,0,1,1,1] ;fit abundances too
    o->dofit,/all
    o->savefit,outfile=filename+'low_a.fits'
    ;o-> set, fit_comp_free_mask= [1,1,0,1,1,0] ;now fit all
    o-> set, fit_comp_free_mask= [1,1,1,1,1,1] ;let's fit the abundances too now
    o->dofit,/all
    o->savefit,outfile=filename+'all_a.fits'
 end

function messenger_ospex_1T, obj=o,time_interval= time_interval,erange = erange,mask=mask,comp_params=comp_params,outfilename=outfilename,quiet=quiet
    o-> set, spex_fit_time_interval= time_interval
    o-> set, spex_eband= [[1.06958, 1.65593], [1.65593, 2.56371], [2.56371, 3.96913], $    
                            [3.96913, 6.14501], [6.14501, 9.51370]] ;is this ok?
    ;o->set,spex_bk_time_interval=this_bkg ;do I need to do this for each one? how to automate...
    o->set, fit_function='vth' ;2 thermal components 'vth+bpow' -> thermal + broken power law
    o-> set, fit_comp_params= comp_params;[1.0, 2.0, 1.0, 0.0, 0.3, 1.0] ;em, T, abund for each component
    o-> set, fit_comp_free_mask= mask ;fit high-energy component first. don't fit abundances
    o->set,spex_fit_manual=0
    set_logenv, 'OSPEX_NOINTERACTIVE', '1' ;I really want to be able to use my mouse!!!
    o->set, spex_erange=erange  ;[2.3,8]
    if quiet eq 'True' then o->set,spex_autoplot_enable=0,spex_fit_progbar=0,spex_fitcomp_plot_resid=0 ;surpress the GUI
    ;endif
    o->dofit,/all
    ;o-> set, fit_comp_free_mask= [1,1,1] ;now fit all
    ;o->savefit,outfile=filename+'_first.fits'
    ;o-> set, fit_comp_free_mask= [1,1,1,1,1,1] ;let's fit the abundances too now
    ;o->dofit,/all
    ;o->savefit,outfile=filename+'all_a.fits'

    o->savefit,outfile=outfilename
end

function get_GOES_sat_flux, time_interval, quiet=quiet, short=short,flux=flux ;time intervals are strings
    a=ogoes()
    goes_class=-1
    ;loadct,13
    ;pick the appropriate satellite for the date
    if anytim(time_interval[0]) lt anytim('27-Jun-2009 22:51:00.000') then sat='13' else if time_interval[0] lt anytim('04-Mar-2010 23:57:00.000') then sat='14' else sat='15'
    print, 'getting goes sat flux at ',time_interval[0]
    a->set, tstart=time_interval[0],tend=time_interval[1],sat=sat ;this has to be a string now not anytim format...
    catch, Error_status
    d=a->getdata(/struct)          ;don't bother with background subtraction...
    print, Error_status,!ERROR_STATE.MSG
    if Error_status eq 0 then begin
                               
    ;if d ne -1 then begin
    if isa(d.ydata,/array) eq 1 then begin
        max_index_short=where(d.ydata[*,1] eq max(d.ydata[*,1]))
        max_index_long=where(d.ydata[*,0] eq max(d.ydata[*,0]))
        goes_class=d.ydata[max_index_long[0],1]
        if flux ne 'flux' then begin
           if short eq 'short' then goes_class=goes_value2class(d.ydata[max_index_short[0],1])
           goes_class=goes_value2class(d.ydata[max_index_long[0],1])
        endif else if short eq 'short' then goes_class=d.ydata[max_index_short[0],1]
        if quiet eq 'False' then begin
            print, goes_class
            long_time_interval=anytim(time_interval)+[-1200.,1200.]
            a->set, tstart=anytim(long_time_interval[0],/vms),tend=anytim(long_time_interval[1],/vms),sat=sat ;longer time interval for plotting
            dd=a->getdata(/struct) ;don't bother with background subtraction...        
            utplot, dd.tarray,dd.ydata[*,0],dd.utbase,/ylog
            oplot, dd.tarray,dd.ydata[*,1], linestyle=2
         endif
     endif
    ;endif
    catch,/cancel
    endif
    print, goes_class
    return, goes_class
end

function get_GOES_class,T,EM,sat=sat,channel=channel, flux=flux
    TMKelvin = T*11.6045
    goes_fluxes,TMKelvin,EM,flong,fshort,sat=sat
    goesnum=flong
    if flux eq 'flux' then begin
        goesclass=flong
        g2=fshort
    endif else begin
        goesclass=goes_value2class(flong)
        g2=goes_value2class(fshort)
    endelse
    print, 'long: ',goesclass,'short:', g2
    if channel eq 'short' then return, g2 else return, goesclass

end

function inspect_GOES, filename,index,flare_list,flare_list_high,flare_list_low,flare_list_all,channel=channel,abund=abund,flux=flux
  if abund eq 'True' then enders=['high_a.fits','low_a.fits','all_a.fits'] else enders=['high.fits','low.fits','all.fits']
  for j=0,2 do begin
        aa=spex_read_fit_results(filename+enders[j])
        infoarr = aa.spex_summ_params
        chisq=aa.spex_summ_chisq
        print, chisq
        time_interval=flare_list.datetimes.obs_start_time[index]
        if anytim(time_interval) lt anytim('27-Jun-2009 22:51:00.000') then sat='13' else if time_interval lt anytim('04-Mar-2010 23:57:00.000') then sat='14' else sat='15'
        print, 'original: ',flare_list.flare_properties.messenger_goes[index]
        case j of
           0: begin
              flare_list_high.flare_properties.messenger_T[index] = infoarr[1] 
              flare_list_high.flare_properties.messenger_EM1[index] = infoarr[0]
              gc=get_GOES_class(infoarr[1],infoarr[0],sat=sat,channel=channel,flux=flux)
              flare_list_high.flare_properties.messenger_goes[index] = gc
              flare_list_high.notes[index] = string(chisq)
              end
           1: begin
              flare_list_low.flare_properties.messenger_T[index] = infoarr[4] 
              flare_list_low.flare_properties.messenger_EM1[index] = infoarr[3]
              gc=get_GOES_class(infoarr[4],infoarr[3],sat=sat,channel=channel,flux=flux)
              flare_list_low.flare_properties.messenger_goes[index] = gc
              flare_list_low.notes[index] = string(chisq)
              end
           2: begin
              flare_list_all.flare_properties.messenger_T[index] = infoarr[1] 
              flare_list_all.flare_properties.messenger_EM1[index] = infoarr[0] 
              gc=get_GOES_class(infoarr[1],infoarr[0],sat=sat,channel=channel,flux=flux)
              flare_list_all.flare_properties.messenger_goes[index] = gc
              flare_list_all.notes[index] = string(chisq)
              end
        endcase
     endfor
    return, [flare_list_high, flare_list_low,flare_list_all]
end

pro goes_check,savename, outfilename,channel=channel, flux=flux,abund=abund
    cd,'/Users/wheatley/Documents/Solar/occulted_flares/flare_lists'
    ;restore,'list_short_formatted.sav'
    restore, savename
    flare_list_high=flare_list  ;make copies of the structure for modification
    if abund eq 'True' then filenames=string(flare_list.id,format='(I)')+'a' else filenames=string(flare_list.id,format='(I)')
    len=size(filenames,/dim)
    goes_list=strarr(len[0])
    gc=goes_list
    xrsnames = flare_list.data_properties.xrs_files;[122:124];[124:214]
    flare_times = flare_list.Datetimes.messenger_datetimes;[122:124];[124:214] ;will need to format these
    ;cd,'data/dat_files'
    for i=0,len[0]-1 do begin
       if i ne 122 then begin
            ;o=ospex(/no_gui)
           ; print, xrsnames[i]+'.dat'
           ; ;o-> set, spex_specfile= 'data/dat_files/'+xrsnames[i]+'.dat'
                                ;o->set, spex_drmfile=this_dfile ;what's this?
           ; ;data=o->getdata(class='spex_data')
            ;;start_time=anytim(flare_times[i]) ;this is already the 'start time', so just go forward
            ;;end_time=anytim(flare_times[i]) +1200. ;is that too much? Maybe stick with 15 or 20 minutes
            ;;peak_time = find_messenger_peak(obj=o,data=data,time_interval= [start_time, end_time],erange = [3.2,8.5]) ;search the high energy band
                                ;print,'data/dat_files/xrs'+filenames[i]+'high.fits'
            ;;flare_list.datetimes.obs_start_time[i]=anytim(peak_time[0],/vms)
            ;;flare_list.datetimes.obs_end_time[i]=anytim(peak_time[1],/vms)
           if anytim(flare_times[i]) lt anytim('27-Jun-2009 22:51:00.000') then sat='13' else if flare_times[i] lt anytim('04-Mar-2010 23:57:00.000') then sat='14' else sat='15'
            aa=spex_read_fit_results('../data/dat_files/'+strtrim(string(filenames[i]),1)+'high_a.fits')
            infoarr = aa.spex_summ_params
            ;print, flare_list.flare_properties.messenger_T[i],flare_list.flare_properties.messenger_EM1[i]
            ;gc[i]=get_GOES_class(flare_list.flare_properties.messenger_T[i],flare_list.flare_properties.messenger_EM1[i],channel=channel,sat=sat)
            gc[i]=get_GOES_class(infoarr[1],infoarr[0],channel=channel,sat=sat, flux=flux)
            if strmid(flare_list.datetimes.obs_start_time[i],0,1) eq "'" then begin 
               time_interval=[strmid(flare_list.datetimes.obs_start_time[i],1,24),strmid(flare_list.datetimes.obs_end_time[i],1,24)]
            endif else begin
               ;time_interval=[peak_time[0],peak_time[1]]
               time_interval=[flare_list.datetimes.obs_start_time[i],flare_list.datetimes.obs_end_time[i]]               
            endelse
            ;print, time_interval;,'here',anytim(time_interval[0])
            goes_list[i]=get_goes_sat_flux(time_interval,quiet='True',short=channel, flux=flux)
         endif else goes_list[i] = -1
    endfor
    ;print, goes_list[121:123]
                                ;print,
                                ;flare_list.flare_properties.RHESSI_goes[0:15]
    print, 'm',gc;,'g',goes_list
    flare_list.flare_properties.goes_goes=goes_list
    flare_list.flare_properties.messenger_goes=gc
    save, flare_list, filename=outfilename
end
    
pro append_chisq,savename,outfilename
    cd,'/Users/wheatley/Documents/Solar/occulted_flares'
    ;restore,'list_short_formatted.sav'
    restore, savename
    flare_list_high=flare_list ;make copies of the structure for modification
    flare_list_low=flare_list_high
    flare_list_all=flare_list_high
    cd,'data/dat_files'
    fnames=file_search('*high.fits')
  
    for i=0,114 do begin
        lists=inspect_GOES(strmid(fnames[i],0,10),i,flare_list,flare_list_high,flare_list_low,flare_list_all)   
    endfor

    flare_list=lists[0]         ;so far it has to be named flare_list...
    save,flare_list,filename=outfilename+'fit_high.sav'
    flare_list=lists[1]
    save,flare_list,filename=outfilename+'fit_low.sav'
    flare_list=lists[2]
    save,flare_list,filename=outfilename+'fit_all.sav'
    cd,'../../'
end


pro use_short,savename,outfilename
    cd,'/Users/wheatley/Documents/Solar/occulted_flares'
    restore, savename+'high_fixed.sav'
    cd,'data/dat_files'
    restore, savename
    fnames=file_search('*high.fits')
  
    for i=0,(fix(size(flare_times,/dim))-1)[0] do begin
        lists=inspect_GOES(strmid(fnames[i],0,10),i,flare_list,flare_list_high,flare_list_low,flare_list_all,channel='short')   
    endfor

    flare_list=lists[0]         ;so far it has to be named flare_list...
    save,flare_list,filename=outfilename+'fit_high.sav'
    flare_list=lists[1]
    save,flare_list,filename=outfilename+'fit_low.sav'
    flare_list=lists[2]
    save,flare_list,filename=outfilename+'fit_all.sav'
    cd,'../../'
end

;main commands
pro do_ospex_fit, savename,outfilename,quiet=quiet,channel=channel,abund=abund,flux=flux
    search_network,/enable
    cd,'/Users/wheatley/Documents/Solar/occulted_flares/flare_lists'
    ;restore,'list_short_formatted.sav'
    restore, savename
    flare_list_high=flare_list ;make copies of the structure for modification
    flare_list_low=flare_list_high
    flare_list_all=flare_list_high
    ids=flare_list.id;[122:124]
    filenames = flare_list.data_properties.xrs_files;[122:124];[124:214]
    flare_times = flare_list.Datetimes.messenger_datetimes;[122:124];[124:214] ;will need to format these
    goes_list=strarr(fix(size(flare_times,/dim)))
    cd,'../data/dat_files'
    erange=[13,15]   ;GOES long - need to check on the peak selections for this case...
    ;erange=[2,25]                  ;GOES short

    for i=0,(fix(size(flare_times,/dim))-1)[0] do begin
        if file_search(filenames[i]+'.dat') eq '' then continue ;if there's no data files, skip this iteration
        ;if i eq 122 then continue;begin                                                    ;time interval error on this one kills ospex
            ;flare_list_high.datetimes.obs_start_time[i]='07-Sep-2012 17:32:11'
            ;flare_list_high.datetimes.obs_end_time[i]='07-Sep-2012 17:32:11'
            ;flare_list_low.datetimes.obs_start_time[i]='07-Sep-2012 17:32:11'
            ;flare_list_low.datetimes.obs_end_time[i]='07-Sep-2012 17:32:11'
            ;flare_list_all.datetimes.obs_start_time[i]='07-Sep-2012 17:32:11'
                                ;flare_list_all.datetimes.obs_end_time[i]='07-Sep-2012
                                ;17:32:11'
           ;continue
                                ;endif else begin
        print,flare_list_high.datetimes.obs_start_time[i]
        if flare_list_high.datetimes.obs_start_time[i] eq '' then begin
            o=ospex(/no_gui)
            o-> set, spex_specfile= strmid(filenames[i],1,10)+'.dat'
                                ;o->set, spex_drmfile=this_dfile ;what's this?
            data=o->getdata(class='spex_data')
            catch, Error_status
            print, Error_status,!ERROR_STATE.MSG
            if Error_status eq 0 then begin
               start_time=anytim(flare_times[i])           ;this is already the 'start time', so just go forwar
               end_time=anytim(flare_times[i]) +1200.      ;is that too much? Maybe stick with 15 or 20 minuidltes
               if strmid(start_time,0,1) eq "'" then begin
                    start_time=anytim(strmid(start_time,1,24))
                    end_time= anytim(start_time) +1200.
                endif 
                print, anytim(start_time,/vms),anytim(end_time,/vms)
                peak_time = find_messenger_peak(obj=o,data=data,time_interval= [start_time, end_time],erange = [3.2,8.5]) ;search the high energy bands
             endif else begin
                ;continue ;if you can't get the data, skip this iteration
                catch,/cancel
             endelse
           ; endif
        endif else begin ;if the fit has already been done
           peak_time=[anytim(flare_list_high.datetimes.obs_start_time[i]),anytim(flare_list_high.datetimes.obs_end_time[i])]
        endelse
        print, 'PEAK TIME:',anytim(peak_time[0],/vms,/truncate),anytim(peak_time[1],/vms,/truncate)
        flare_list_high.datetimes.obs_start_time[i]=anytim(peak_time[0],/vms,/truncate)
        flare_list_high.datetimes.obs_end_time[i]=anytim(peak_time[1],/vms,/truncate)
        flare_list_low.datetimes.obs_start_time[i]=anytim(peak_time[0],/vms,/truncate)
        flare_list_low.datetimes.obs_end_time[i]=anytim(peak_time[1],/vms,/truncate)
        flare_list_all.datetimes.obs_start_time[i]=anytim(peak_time[0],/vms,/truncate)
        flare_list_all.datetimes.obs_end_time[i]=anytim(peak_time[1],/vms,/truncate)
        
        if abund eq 'True' then begin
            if file_search(strtrim(string((ids[i]),format='(I)'),1)+'ahigh_a.fits') eq '' then begin 
               foo=messenger_ospex(obj=o,time_interval= [peak_time[0], peak_time[1]],erange = [2.3,8.5],mask=[1,1,1,0,0,0],comp_params=[1.0, 2.0, 1.0, 0.3, 0.3, 1.0],outfilename=string(ids[i],format='(I)')+'a', quiet=quiet)
               lists=[flare_list_high,flare_list_low,flare_list_all]
            endif else lists=inspect_GOES(strtrim(string(ids[i],format='(I)'),1)+'a',i,flare_list,flare_list_high,flare_list_low,flare_list_all,channel=channel,abund=abund,flux=flux)
        endif else begin
            if file_search(strtrim(string((ids[i]),format='(I)'),1)+'high.fits') eq '' then begin
                foo=messenger_ospex(obj=o,time_interval=[peak_time[0], peak_time[1]],erange = [2.3,8.5],mask=[1,1,0,0,0,0],comp_params=[1.0,2.0, 1.0, 0.3, 0.3,1.0],outfilename=string(ids[i],format='(I)'),quiet=quiet)
                lists=[flare_list_high,flare_list_low,flare_list_all]
            endif else lists=inspect_GOES(strtrim(string(ids[i],format='(I)'),1),i,flare_list,flare_list_high,flare_list_low,flare_list_all,channel=channel,abund=abund,flux=flux)
        endelse
        print, lists[0].flare_properties.goes_goes[0:10]
        goes_list[i]=get_goes_sat_flux([anytim(peak_time[0],/vms),anytim(peak_time[1],/vms)],quiet='True',short=channel)
    endfor

    lists[0].flare_properties.goes_goes=goes_list
    lists[1].flare_properties.goes_goes=goes_list
    lists[2].flare_properties.goes_goes=goes_list
    
    flare_list=lists[0]         ;so far it has to be named flare_list...
    save,flare_list,filename=outfilename+'_fit_high.sav'
    flare_list=lists[1]
    save,flare_list,filename=outfilename+'_fit_low.sav'
    flare_list=lists[2]
    save,flare_list,filename=outfilename+'_fit_all.sav'
    cd,'../../'
end

 
