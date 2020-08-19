function load_flare_list, time_range
    ;load the original flare list  
    flare_list_orig_obj=hsi_flare_list( obs_time_interval = time_range )  
    flare_list_orig=flare_list_orig_obj -> getdata()  
    dataq =flare_list_orig_obj -> Get(flag_name = 'DATA_QUALITY' )
    qlook=flare_list_orig_obj -> Get(flag_name = 'IMAGE_STATUS' )
    dim_orig_flare = n_elements(flare_list_orig)  
 
    ;print, 'Orig Flare List: Found ', dim_orig_flare, anytim2utc(flare_list_orig[0].peak_time,/vms)
    struct = {flare_list_orig:flare_list_orig,dataq:dataq,qlook:qlook}
    return, struct
    end

function pick_radius,flare_list_orig
    x_orig = flare_list_orig.position(0)  
    y_orig = flare_list_orig.position(1)  

    ;first just find any flares that might be occulted
    radius_orig = sqrt( x_orig^2 + y_orig^2 )  
    radius_sun= get_rb0p(flare_list_orig.start_time, /radius)

    ;print, mean(radius_orig/radius_sun)

    limb=where(radius_orig/radius_sun gt 1.0)
    print,limb
    ;anytim2utc(list[limb[0]].peak_time,/vms)

    return, limb                ;,limb_flare_list
    end
    ;print, size(limb_flare_list)

function pick_dataq,dataq, threshold
    ; let's choose data quality <2
    dataq_index = where(dataq gt threshold)
    return, dataq_index
    end

function pick_counts,counts, threshold
    ; let's choose data quality <2
    count_index = where(counts gt threshold)
    return, count_index
    end

function is_qlook,qlook
     ;and those with a quicklook image to make our job easier
    qlook_index = where(qlook eq 1)
    return, qlook_index
    end

function do_quicklook, list_filtered
    ; generate quicklook images and take a look
    times = strmid(ptim(list_filtered[0].peak_time,/ecs),0,17) + ['00','59']
    for i=1,9 do begin
        times = times + strmid(ptim(list_filtered[i].peak_time,/ecs),0,17) + ['00','59']
    endfor
    ;timef=ptim(times,/ecs)
    ;print,timef
    o= hsi_image()
    o->set,im_time_interval=times
    data=o->get_data()
    o->panel_display
    obj_destroy,o
    end

function filter_flares, time_range, qthresh,countthresh
    ;time_range = ['2004/08/20 00:00', '2011/01/28 00:00']  
    data = load_flare_list(time_range)
    list=data.flare_list_orig
    dataq=data.dataq
    qlook=data.qlook
    limb = pick_radius(list)
    print, size(limb,/dim)
    dataq_index=pick_dataq(dataq,qthresh)
    print, size(dataq_index,/dim)
    qlook_index=is_qlook(qlook)
    print, size(qlook_index,/dim)
    count_index=pick_counts(list.total_counts,countthresh)
    print, size(count_index,/dim)
    ;now get the flares where indexes match

    length=size(limb,/dim)
    for i=0, length[0]-1 do begin
       if ((where(dataq_index eq limb[i]) ne -1) AND (where(qlook_index eq limb[i]) ne -1)) then begin
          if where(count_index eq limb[i]) ne -1 then begin
             if isa(list_filtered) eq 0 then list_filtered = list[limb[i]] else list_filtered=[list_filtered,list[limb[i]]]
          endif
      endif
   ;endif
   endfor
    print,'list_filtered ', size(list_filtered,/dim)
    return, list_filtered
end


;delvar, list_filtered ;so that it doesn't just append to the old list
time_range = ['2007/05/28 00:00','2013/08/19 00:00'] ;, '2011/01/28 00:00']
qthresh=2
countthresh=100000
list_filtered=filter_flares(time_range, qthresh,countthresh)
save, list_filtered, filename='list_matches.sav'
aa=ptim(list_filtered.peak_time,/ecs)
bb=list_filtered.ID_NUMBER
cc=list_filtered.GOES_CLASS
write_csv, 'list_matches.csv',bb,aa,cc
;save date and time list to .csv, open in python, use webbrowser to
;open all quicklook images in browser. then figure out how to get the
;urls of all the tabs that are left open at the end of inspection,
;extract the date and time from them, and save this as a csv to use as
;a final filter on the flare list. that's a lot of work...
;foo=do_quicklook(list_filtered[0:10])

end
