pro get_messenger_peak, flare_list,energy_bin,sp_time_int=sp_time_int,mask=mask,extract_data=extract_data,quiet=quiet
    llen = size(flare_list.id,/dimensions)
    for i=0,llen[0]-1 do begin
        end_time=anytim(flare_list.datetimes.messenger_datetimes[i])+1200.
        ft = [flare_list.datetimes.messenger_datetimes[i],anytim(end_time,/vms)] 
        ;drate = one_curve_R(ft,energy_bin,sp_time_int=sp_time_int,mask=mask,quiet=quiet)
        outfilename='data/lightcurves/'+ strtrim(string(flare_list.id[i]),1) +'drate.sav'
        if keyword_set(extract_data) eq 1 then save, drate,filename=outfilename
    endfor
end

function one_curve_M,flare_times,xrs_file,energy_bin,quiet=quiet
    ;let's use GOES energy ranges
    erange=[energy_bin[0],energy_bin[1]]
    time_interval=[anytim(flare_times[0]),anytim(flare_times[1])]
    obj=ospex(/no_gui)
    obj-> set, spex_specfile= 'data/dat_files/' + strtrim(xrs_file,1)+'.dat'
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
    data={taxis:taxis,phigh:phigh}
    return,data
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

