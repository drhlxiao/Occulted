pro lightcurves, flare_list,energy_bin,sp_time_int=sp_time_int,mask=mask,extract_data=extract_data,quiet=quiet
    llen = size(flare_list.id,/dimensions)
    for i=0,llen[0]-1 do begin
        end_time=anytim(flare_list.datetimes.messenger_datetimes[i])+1200.
        ft = [flare_list.datetimes.messenger_datetimes[i],anytim(end_time,/vms)] 
        ;drate = one_curve_R(ft,energy_bin,sp_time_int=sp_time_int,mask=mask,quiet=quiet)
        outfilename='data/lightcurves/'+ strtrim(string(flare_list.id[i]),1) +'drate.sav'
        if keyword_set(extract_data) eq 1 then save, drate,filename=outfilename
    endfor
end

function one_curve_R, flare_times,energy_bin,sp_time_int=sp_time_int,mask=mask,quiet=quiet
    o=hsi_spectrum()
    o->set,obs_time_interval=[flare_times[0],flare_times[1]]
    o -> set, sp_energy_binning=[energy_bin[0],energy_bin[1]]
    if keyword_set(sp_time_int) eq 0 then o -> set, sp_time_interval=4. else o -> set, sp_time_interval=sp_time_int
    if keyword_set(mask) eq 0 then o -> set, seg_index_mask=[1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0] else o->set,seg_index_mask=mask
    o -> set, /decimation_correct, /rear_decimation_correct
    drate = o-> getdata(sp_data_unit='rate', /sp_data_struct)   
    if keyword_set(quiet) eq 0 then begin
       o -> plotman, sp_data_unit='Flux'
       utplot,average(drate.ut,1),drate.rate,xrange=[drate(0).ut(0),max(drate.ut(0))],yrange=[0,150]           
       eutplot,average(drate.ut,1),drate.rate-drate.erate,drate.rate+drate.erate,width=1d-4
    endif
    obj_destroy,o
    return, drate
 end

function one_curve_M,flare_times,xrs_file,energy_bin,quiet=quiet
    ;let's use GOES energy ranges
    erange=[energy_bin[0],energy_bin[1]]
    time_interval=[anytim(flare_times[0]),anytim(flare_times[1])]
    obj=ospex(/no_gui)
    obj-> set, spex_specfile= 'data/dat_files/' + strtrim(xrs_file,1)+'.dat'
    obj->set_time,time_interval
    data=obj->getdata(class='spex_data')
    eaxis=obj->getaxis(/ct_energy)
    taxis=obj->getaxis(/mean) ;are these formatted from a different date?
    elist=where( (eaxis ge min(erange)) AND (eaxis le max(erange)) )
    phigh = total(data.data(elist,*),1)
    ;tlist=where( (taxis ge time_interval[0]) AND (taxis le time_interval[1]))
    ;mti=size(tlist,/dim)
                                ;newtaxis = [taxis[tlist[0]-1],taxis[tlist]]
    if keyword_set(quiet) eq 0 then utplot,taxis, phigh,timerange=[time_interval[0],time_interval[1]],yrange=[min(phigh),max(phigh)]
    obj_destroy,obj
    Mtim=taxis
    nMtim = anytim(Mtim,/vms)

    data={taxis:nMtim,phigh:phigh}
    return,data
 end

function one_spec_M,flare_times,xrs_file,energy_bin,quiet=quiet
    ;let's use GOES energy ranges
    erange=[energy_bin[0],energy_bin[1]]
    time_interval=[anytim(flare_times[0]),anytim(flare_times[1])]
    tdelta=anytim(flare_times[1])-anytim(flare_times[0])
    obj=ospex(/no_gui)
    obj-> set, spex_specfile= 'data/dat_files/' + strtrim(xrs_file,1)+'.dat'
    obj->set_time,time_interval
    data=obj->getdata(class='spex_data')
    eaxis=obj->getaxis(/ct_energy)
    ;taxis=obj->getaxis(/mean) ;are these formatted from a different date?
    ;elist=where( (eaxis ge min(erange)) AND (eaxis le max(erange)) )
    ctot = total(data.data,2)/tdelta
    ;tlist=where( (taxis ge time_interval[0]) AND (taxis le time_interval[1]))
    ;mti=size(tlist,/dim)
                                ;newtaxis = [taxis[tlist[0]-1],taxis[tlist]]
    if keyword_set(quiet) eq 0 then utplot,taxis, phigh,timerange=[time_interval[0],time_interval[1]],yrange=[min(phigh),max(phigh)]
    obj_destroy,obj
    ;Mtim=taxis
    ;nMtim = anytim(Mtim,/vms)

    data={eaxis:eaxis,ctot:ctot}
    return,data
 end

function one_curve_G,flare_times,quiet=quiet
    ;GOES lightcurves
    ft=[anytim(flare_times[0],/vms),anytim(flare_times[1],/vms)]
    if anytim(ft[0]) lt anytim('27-Jun-2009 22:51:00.000') then sat='13' else if anytim(ft[0]) lt anytim('04-Mar-2010 23:57:00.000') then sat='14' else sat='15'
    a=ogoes()
    a->set, tstart=ft[0],tend=ft[1],sat=sat
    d=a->getdata(/struct)
    if keyword_set(quiet) eq 0 then utplot,d.tarray, d.ydata[0,*],timerange=[time_interval[0],time_interval[1]],yrange=[min(phigh),max(phigh)]
                                ;print, 'before',d.tarray[0:5]
    if isa(d,'struct') then begin
        bb = string(anytim(d.tarray+anytim(d.UTbase),/vms)) ;fix this to vms?
        data={UTbase:d.UTbase,tarray:bb,ydata:d.ydata,yclean:d.yclean}
     endif else data=-1
    ;print, 'after', data.tarray[0:5]
    obj_destroy,a

    return,data
 end

;restore,'flare_lists/more_occulted_final.sav'
;d=get_data(flare_list,[[4,9],[12,18],[18,30],[30,80]],[[1.55,12.4],[3.099,24.8]],export='/data/test_data.sav',/quiet,/goes)

function format_mask,mask ;this is gonna get ugly...
  mask_formatted=[strmid(mask,1,1),strmid(mask,2,1),strmid(mask,3,1),strmid(mask,4,1),strmid(mask,5,1),strmid(mask,6,1),strmid(mask,7,1),strmid(mask,8,1),strmid(mask,9,1),0,0,0,0,0,0,0,0,0]
  return,float(mask_formatted)
  end

pro testloop
  k=0
  for i=0,2 do begin
  for j=0,4 do begin
     print,k
     k=k+1
  endfor
  ;k=k+1
endfor
  end

function get_data,flare_list, Rbins=Rbins, Mbins=Mbins,export=filename, goes=goes,quiet=quiet
    llen = size(flare_list.id,/dimensions)
    start_time=anytim(flare_list.datetimes.obs_start_time[0]);anytim(flare_list.datetimes.messenger_datetimes[0])-600.
    end_time=anytim(flare_list.datetimes.obs_end_time[0]);anytim(flare_list.datetimes.messenger_datetimes[0])+1200.
    mask=flare_list.data_properties.good_det[0]
    mask_formatted=format_mask(mask)
    ft = [anytim(start_time,/vms),anytim(end_time,/vms)]
    if keyword_set(Rbins) eq 0 then Rbins = [[4,9],[12,18],[18,30],[30,80]]
    if keyword_set(Mbins) eq 0 then Mbins = [[1.55,12.4],[3.099,24.8]]    
    nbins=size(Rbins,/dimensions)
    drate = one_curve_R(ft,[Rbins[0,0],Rbins[1,0]],sp_time_int=sp_time_int,mask=mask_formatted,quiet=quiet)
                                ;newdrate={UT:drate.UT,rate:drate.rate,erate:drate.erate,ltime:drate.ltime,erange:strtrim(string(Rbins[0,0]),1)+'-'+strtrim(string(Rbins[1,0]),1)} ;let's convert UT back to string so I can also use in python if I want ;fix this! make one dictionary entry for each energy range
    newdrate={UT:dblarr(2,600),rate:fltarr(600),erate:fltarr(600),ltime:fltarr(600),erange:strtrim(string(Rbins[0,0]),1)+'-'+strtrim(string(Rbins[1,0]),1)}
    Rdata=replicate(newdrate,(nbins[1])*llen[0])
    k=1
    for i=1,llen[0]-1 do begin
        start_time=anytim(flare_list.datetimes.obs_start_time[i]);anytim(flare_list.datetimes.messenger_datetimes[i])-600.
        end_time=anytim(flare_list.datetimes.obs_end_time[i]);anytim(flare_list.datetimes.messenger_datetimes[i])+1200.
        mask=flare_list.data_properties.good_det[i]
        ft = [anytim(start_time,/vms),anytim(end_time,/vms)]
        for j=0,nbins[1]-1 do begin
           drate = one_curve_R(ft,[Rbins[0,j],Rbins[1,j]],sp_time_int=sp_time_int,mask=mask_formatted,quiet='quiet')
           newdrate={UT:drate.UT,rate:drate.rate,erate:drate.erate,ltime:drate.ltime,erange:strtrim(string(Rbins[0,j]),1)+'-'+strtrim(string(Rbins[1,j]),1)}
           Rdata[k] = newdrate
           k=k+1
        endfor
        ;k=k+1
    endfor

    ;; nbins=size(Mbins,/dimensions)    
    ;; ft = [flare_list.datetimes.messenger_datetimes[0],anytim(end_time,/vms)]
    ;; Mdat = one_curve_M(ft,flare_list.data_properties.xrs_files[0],[Mbins[0,0],Mbins[1,0]],quiet=quiet)
    ;; Mdat = {Taxis:dblarr(1000),phigh:fltarr(1000),len:1000L}
    ;; Mdata=replicate(Mdat,(nbins[1]+1)*llen[0]) ;careful... what if these are different sizes?
    ;; k=0
    ;; for j=0,nbins[1]-1 do begin
    ;;    for i=0,llen[0]-1 do begin
    ;;       ;start_time=anytim(flare_list.datetimes.messenger_datetimes[i])-600.
    ;;       ;end_time=anytim(flare_list.datetimes.messenger_datetimes[i])+1200.
    ;;       start_time=anytim(flare_list.datetimes.obs_start_time[i]);anytim(flare_list.datetimes.messenger_datetimes[i])-600.
    ;;       end_time=anytim(flare_list.datetimes.obs_end_time[i]);anytim(flare_list.datetimes.messenger_datetimes[i])+1200.
    ;;       ft = [anytim(start_time,/vms),anytim(end_time,/vms)]
    ;;       print,ft,k
    ;;       data = one_curve_M(ft,flare_list.data_properties.xrs_files[i],[Mbins[0,j],Mbins[1,j]],quiet=quiet)  
    ;;       Mdata[k].taxis = data.taxis
    ;;       Mdata[k].phigh = data.phigh
    ;;       dlen= size(data.phigh,/dimensions)
    ;;       Mdata[k].len = dlen[0]
    ;;       k=k+1
    ;;    endfor
    ;;    k=k+1
    ;; endfor
    ;; ;data={Rdata:Rdata,Mdata:Mdata}
   
    ;; if keyword_set(goes) eq 1 then begin ;something is wrong with how this is getting proccessed!tarray is weird
    ;;     d={UTbase:'foo',tarray:dblarr(1172),Ydata:fltarr(1172,2),len:1172L}
    ;;     Gdata=replicate(d,llen)
    ;;    for i=0,llen[0]-1 do begin
    ;;       a=ogoes()
    ;;       start_time=anytim(flare_list.datetimes.obs_start_time[i]);anytim(flare_list.datetimes.messenger_datetimes[i])-600.
    ;;       end_time=anytim(flare_list.datetimes.obs_end_time[i]);anytim(flare_list.datetimes.messenger_datetimes[i])+1200.
    ;;       ;start_time=anytim(flare_list.datetimes.messenger_datetimes[i])-600.
    ;;       ;end_time=anytim(flare_list.datetimes.messenger_datetimes[i])+1200.
    ;;       ft = [anytim(start_time,/vms),anytim(end_time,/vms)]
    ;;       if anytim(ft[0]) lt anytim('27-Jun-2009 22:51:00.000') then sat='13' else if anytim(ft[0]) lt anytim('04-Mar-2010 23:57:00.000') then sat='14' else sat='15'
    ;;        a->set, tstart=ft[0],tend=ft[1],sat=sat
    ;;        d=a->getdata(/struct)
    ;;        glen= size(d.tarray,/dimensions)
    ;;        Gdata[i].UTbase = d.UTbase
    ;;        Gdata[i].tarray[0:glen-1]=d.tarray
    ;;        print, size(d.ydata,/dim),glen
    ;;        Gdata[i].ydata[0:glen-1,*]=d.ydata
    ;;        Gdata[i].len = glen[0]
    ;;        obj_destroy,a
    ;;    endfor
    ;;    data={Mdata:Mdata,Gdata:Gdata}
    ;;    ;data={Gdata:Gdata}
    ;; endif
       data={Rdata:Rdata}
    ;; if keyword_set(export) eq 1 then write_csv, filename,data ;,header=;tags of structure
    ;save, data, filename='data/lc_data.sav'
    return, data
 end

pro fix_times,data,filename
    ;fix times so they can be converted to Python datettimes
    Rtim=data.rdata.ut
    nRtim = anytim(Rtim,/vms)    
    nRdata={UT:nRtim,rate:data.rdata.rate,erate:data.rdata.erate,ltime:data.rdata.ltime}
    
    Mtim=data.mdata.taxis
    nMtim = anytim(Mtim,/vms)
    nMdata={taxis:nMtim,phigh:data.mdata.phigh,len:data.mdata.len}
    
    Gtim=data.gdata.tarray
    Gdim = size(Gtim,/dim)
    ut=data.gdata.utbase
    for j=0,Gdim[1]-1 do for i=0,Gdim[0]-1 do Gtim[i,j] = Gtim[i,j] + anytim(ut[j])
    nGtim = anytim(Gtim,/vms)
    nGdata={taxis:nGtim,ydata:data.gdata.ydata,len:data.gdata.len}
    
    data={Rdata:nRdata,Mdata:nMdata,Gdata:nGdata}
    save, data, filename=filename
 end

function get_spectrogram, flare_times,ebins,sp_time_int=sp_time_int,mask=mask,quiet=quiet
    o=hsi_spectrum()
    o->set,obs_time_interval=[flare_times[0],flare_times[1]]
    o -> set, sp_energy_binning=ebins
    if keyword_set(sp_time_int) eq 0 then o -> set, sp_time_interval=4. else o -> set, sp_time_interval=sp_time_int
    if keyword_set(mask) eq 0 then o -> set, seg_index_mask=[1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0] else o->set,seg_index_mask=mask
    o -> set, /decimation_correct, /rear_decimation_correct
    drate = o-> getdata(sp_data_unit='rate', /sp_data_struct)   
    if keyword_set(quiet) eq 0 then begin
       ;loadct,13
       ;o -> plotman, sp_data_unit='Flux'
       o -> plot, /pl_spec, sp_data_unit='flux',/log_scale,/ylog,plotman_obj=p, desc='Unsmoothed spectrogram', xrange=[flare_times[0],flare_times[1]], yrange=[ebins[0],ebins[39]]
    endif
    obj_destroy,o
    return, drate
 end

pro spectrograms,flare_list,filename,quiet=quiet
    llen = size(flare_list.id,/dimensions)
    ee=10^(2.*findgen(40)/40+alog10(3))
    epairs=dindgen(40,2)
    for j=0,38 do epairs[j,*]=[ee[j],ee[j+1]]
    ebins=[epairs[0,*],epairs[1,*],epairs[2,*],epairs[3,*],epairs[4,*],epairs[5,*],epairs[6,*],epairs[7,*],epairs[8,*],epairs[9,*],epairs[10,*],epairs[11,*],epairs[12,*],epairs[13,*],epairs[14,*],epairs[15,*],epairs[16,*],epairs[17,*],epairs[18,*],epairs[19,*],epairs[20,*],epairs[21,*],epairs[22,*],epairs[23,*],epairs[24,*],epairs[25,*],epairs[26,*],epairs[27,*],epairs[28,*],epairs[29,*],epairs[30,*],epairs[31,*],epairs[32,*],epairs[33,*],epairs[34,*],epairs[35,*],epairs[36,*],epairs[37,*],epairs[38,*]] ;why is IDL so dumb
    data={UT:strarr(675),rate:fltarr(39,675),erate:fltarr(39,675),ltime:fltarr(39,675),len:450L} ;need to re-format this for a more flexible format...
    data=replicate(data,llen[0])
    ;spect=dindgen()
    
    for i=0,llen[0]-1 do begin
        start_time=anytim(flare_list.datetimes.obs_start_time[i]);anytim(flare_list.datetimes.messenger_datetimes[i])-600.
        end_time=anytim(flare_list.datetimes.obs_end_time[i]);anytim(flare_list.datetimes.messenger_datetimes[i])+1200.
        ;start_time=;anytim(flare_list.datetimes.messenger_datetimes[i])-600.
        ;end_time=;anytim(flare_list.datetimes.messenger_datetimes[i])+1200.
        mask=flare_list.data_properties.good_det[i]
        mask_formatted=format_mask(mask)
        ft = [anytim(start_time,/vms),anytim(end_time,/vms)]
        aa=get_spectrogram(ft,ebins,mask=mask_formatted,quiet=quiet)
        length=size(aa.rate,/dim)
        data[i].UT[0:length[1]-1]=anytim(aa.UT[0],/vms) ;datetime string
        print, data[i].UT[0],data[i].UT[length[1]-1]
        data[i].rate[*,0:length[1]-1]=aa.rate
        data[i].erate[*,0:length[1]-1]=aa.erate
        data[i].ltime=string(anytim(aa.ltime[0],/vms)) ;fix this later...
        data[i].len=length[1]
    endfor
    save, data, ebins, filename=filename
end
    
pro plot_together,data,quiet=quiet;plot RHESSI and Messenger lightcurves together
    ;[4,9]   in this range we have best sensitivity for faint events (i.e. highly occulted flares)
    ;[12,18]  in this range we should see thermal or non-thermal emissions in occulted flares, depending on flare size
    ;[18,30]  this range is for non-thermal emissions in the larger occulted flares, with additional thermal emission from largest flares
    ;[30,80]  this range is normally background in occulted flares except for the very largest flare, but lets monitor this range
    
    ;deal with RHESSI missing data/change of mode

    ;make nice plot with nice colors-actually do I want to move this to Python? It could be nicer and easier
    ;utplot,
    ;utplot,average(drate.ut,1),drate.rate,xrange=[drate(0).ut(0),max(drate.ut(0))],yrange=[0,150]
    ;oplot, stuff
    ;oplot, more stuff

    ;legend_labels = ;make legend labels
    ;al_legend,legend_labels,colors=[],/top,/right
    ;eutplot,average(drate.ut,1),drate.rate-drate.erate,drate.rate+drate.erate,width=1d-4

    ;save all data
end
