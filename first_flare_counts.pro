function one_curve_M,flare_times,xrs_file,energy_bin,quiet=quiet
    erange=[energy_bin[0],energy_bin[1]]
    time_interval=[anytim(flare_times[0]),anytim(flare_times[1])]
    obj=ospex(/no_gui)
    obj-> set, spex_specfile= 'data/dat_files/' + strtrim(xrs_file,1)+'.dat'
    obj->set_time,time_interval
    data=obj->getdata(class='spex_data')
    eaxis=obj->getaxis(/ct_energy)
    taxis=obj->getaxis(/mean) 
    elist=where( (eaxis ge min(erange)) AND (eaxis le max(erange)) )
    phigh = total(data.data(elist,*),1)

    if keyword_set(quiet) eq 0 then utplot,taxis,	   phigh,timerange=[time_interval[0],time_interval[1]],yrange=[min(phigh),max(phigh)]
    obj_destroy,obj
    Mtim=taxis
    nMtim = anytim(Mtim,/vms)

    data={taxis:nMtim,phigh:phigh}
    return,data

; To run:    
;ft=['28-May-2007 10:25:00','28-May-2007 10:55:00']
;ebin=[6.3,7.0]
;xrs_file='xrs2007148'
;chigh=one_curve_M(ft,xrs_file,ebin,quiet=0)
;print, max(chigh.phigh)

end
