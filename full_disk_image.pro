;pro full_disk_image

function extend_time_interval,time_int,delta ;delta in seconds
  new_time_int = [anytim(time_int[0]) - delta,anytim(time_int[1]) + delta]
  return, new_time_int ;in anytim format
end

function do_bproj,time_int,nimages,outfilename=outfilename, quiet=quiet
    o=hsi_image(/nogui) 
    o->set,im_time_interval=time_int
    o->set,im_energy_binning=[6,25]
    if nimages eq 1 then begin
    im=o->getdata( image_alg='back',clean_niter=111,$
pixel_size=[16,16]*1,image_dim=[128,128]*1.5,xyoffset=[0,0],$
           det_index_mask=[0,0,0,0,0,0,0,1,1],$
           use_local_average=1)
    ;if quiet ne 'quiet' then o->plotman
                                ;overplot solar disk
                                ;save file    
    if mean(im) ne 0 then begin
        o->plotman,plotman_obj=p
        p->create_plot_file,/png,filename=outfilename+'.png'
        ;save, im, filename=outfilename+'.sav'
        obj_destroy,p
     endif
 endif else begin
    im=o->getdata( image_alg='back',clean_niter=111,$
pixel_size=[16,16]*1,image_dim=[128,128]*1.5,xyoffset=[0,0],$
           det_index_mask=[0,0,0,0,0,0,0,1,0],$
           use_local_average=1)
    ;if quiet ne 'quiet' then o->plotman
                                ;overplot solar disk
                                ;save file    
    if mean(im) ne 0 then begin
        ;o->plotman,plotman_obj=p
        ;p->create_plot_file,/png,filename=outfilename+'_det8.png'
        save, im, filename=outfilename+'_det8.sav'
        ;obj_destroy,p
     endif
     im=o->getdata( image_alg='back',clean_niter=111,$
pixel_size=[16,16]*1,image_dim=[128,128]*1.5,xyoffset=[0,0],$
           det_index_mask=[0,0,0,0,0,0,0,0,1],$
           use_local_average=1)
    ;if quiet ne 'quiet' then o->plotman
                                ;overplot solar disk
                                ;save file    
    if mean(im) ne 0 then begin
        o->plotman,plotman_obj=p
        p->create_plot_file,/png,filename=outfilename+'_det9.png'
        save, im, filename=outfilename+'_det9.sav'
        ;obj_destroy,p
     endif
    endelse
 obj_destroy,o
end

function do_clean,time_int,pix_size,det,image_dim,xyoffset,erange,outfilename=outfilename, quiet=quiet
    o=hsi_image(/nogui) 
    o->set,im_time_interval=time_int
    o->set,im_energy_binning=erange;[6,25]
    im=o->getdata( image_alg='clean',clean_niter=111,$
pixel_size=pix_size*1,image_dim=image_dim*1.5,xyoffset=xyoffset,$
           det_index_mask=det,$
           use_local_average=1)
    ;if quiet ne 'quiet' then o->plotman
                                ;overplot solar disk
                                ;save file    
    if mean(im) ne 0 then begin
                                ;o->plotman,plotman_obj=p
       o->set,im_out_fits_filename=outfilename
       o->fitswrite
        ;p->create_plot_file,/png,filename=outfilename+'.png'
        save, o, filename=outfilename+'_obj.sav'
        ;obj_destroy,p
     endif
    obj_destroy,o
    return,im
end

function get_coords, im_filename ;the sav file
  restore,im_filename
  mx=max(im,loc)                ;should be the brightest source in most images
  coords=array_indices(im,loc)
                                ;convert to solar coords given the
                                ;resolution...center of image (and
                                ;sun) is at [96,96]
  ;radius_sun= int(get_rb0p(anytim('27-Apr-2013 15:10:20.000'), /radius))/16
  ;16"/pix ~60 pix
  xyoffset = [(coords[0]-96)*16,(coords[1]-96)*16]
  return,xyoffset ;store this in source_vis?
end

function recalc_angle,xyoffset,date
   angle = messenger_flare_angle(xyoffset,date)
   return,angle
   end

pro recalc_angles,flare_list,rhessi_list=rhessi_list ;add in option to see if location is given in rhessi flare list
  len=size(flare_list.id,/dim)
  angles=dblarr(len[0])
  loc=strarr(len[0])
  for i=0,len[0]-1 do begin
          date=strmid(flare_list.datetimes.messenger_datetimes[i],0,10)
          angles[i] = recalc_angle([0.,0.],anytim(date)); flare_list.angle[i]
     if file_search('data/bproj_vis/'+strtrim(string(flare_list.id[i]),1)+'_det9.sav') ne '' then begin
          ;print,'here'
          xyoffset=get_coords('data/bproj_vis/'+strtrim(string(flare_list.id[i]),1)+'_det9.sav')
          angles[i]=recalc_angle(xyoffset,anytim(date))
          loc[i]='['+strtrim(string(xyoffset[0]),1)+':'+strtrim(string(xyoffset[1]),1) +']'
       endif else begin
          if keyword_set(rhessi_list) ne 0 then begin
             ;print,'y'          ;stuff
             tr=[anytim(flare_list.datetimes.messenger_datetimes[i])-1800.,anytim(flare_list.datetimes.messenger_datetimes[i])+1800.]
             a=load_flare_list(tr)
             loc[i]='[0,0]'             
             if typename(a.flare_list_orig) eq 'STRUCT' then begin
                b='['+strtrim(string(a.flare_list_orig.position[0,0]),1)+':'+strtrim(string(a.flare_list_orig.position[1,0]),1) +']'
                loc[i] = b[0]
                angles[i]=recalc_angle(a.flare_list_orig[0].position,anytim(date))
              endif
          endif else begin
              loc[i]='[0:0]'
       endelse
       endelse
       print,loc[i],angles[i]
    endfor

  flare_list.Angle=angles
  flare_list.flare_properties.location=loc
  end

function full_disk_image, quiet=quiet
  ;load the flare list
  restore,'flare_lists/more_occulted.sav'
  cd,'data/bproj_vis'
  llen=size(flare_list.id,/dim)
  for i=10,20 do begin;llen[0] do begin
     time_int = extend_time_interval([flare_list.datetimes.obs_start_time[i],flare_list.datetimes.obs_end_time[i]],60.) ;start with 60 s
     print, time_int
      foo=do_bproj(time_int,2,outfilename=strtrim(string(flare_list.id[i]),1),quiet=quiet)
      print, 'flare is ',flare_list.id[i]
  endfor
  end

function image_flare, flare_list,erange, quiet=quiet
  llen=size(flare_list.id,/dim)
  if keyword_set(quiet) then quiet='True'
  pix_size=[1,1]
  det=[0,0,1,1,1,1,1,1,1] ;will eventually have to look at spectrographs to see which ones should and should not be used
  image_dim=[128,128]
  for i=10,20 do begin;llen[0] do begin
                                ;if
                                ;anytim(flare_list.datetimes.obs_end_time[i])
                                ;-
                                ;anytim(flare_list.datetimes.obs_start_time[i])
                                ;le 30. then do time_int =
                                ;extend_time_interval([flare_list.datetimes.obs_start_time[i],flare_list.datetimes.obs_end_time[i]],60.)
                                ;;only do this if original time
                                ;;interval is short...how long of a
     ;;time interval do i actually want
     ;;though?
     print, flare_list.id[i]
     fname='data/bproj_vis/'+strtrim(string(flare_list.id[i]),1)
     if file_search(fname+'_round1.png') ne '' then continue     
     time_int = [flare_list.datetimes.obs_start_time[i],flare_list.datetimes.obs_end_time[i]]
     print, 'here',fname, file_search(fname+'_det8.sav'),file_search(fname+'_det9.sav')
     if file_search(fname+'_det9.sav') ne '' then begin
        xyoffset=get_coords(fname+'_det9.sav')
        if abs(xyoffset[0]) or abs(xyoffset[1]) gt 1000. then begin
           print, 'here',file_search(fname+'_det8.sav')
           if file_search(fname+'_det8.sav') ne '' then xyoffset=get_coords(fname+'_det8.sav') ;else continue
        endif
     endif else continue
     outfilename='data/bproj_vis/'+strtrim(string(flare_list.id[i]),1)+'_round1'
     im=do_clean(time_int,pix_size,det,image_dim,xyoffset,erange,outfilename=outfilename,quiet=quiet)
    endfor
end

 function image_one, flare,i,erange,xyoffset,quiet=quiet
  if keyword_set(quiet) then quiet='True'
  pix_size=[1,1]
  det=[0,0,1,1,1,1,1,1,1]
  image_dim=[128,128]
  fname='data/bproj_vis/'+strtrim(string(flare.id[i]),1)
  ;if file_search(fname+'_round1.png') ne '' then begin     
  time_int = extend_time_interval([flare.datetimes.obs_start_time[i],flare.datetimes.obs_end_time[i]],60.)
  print, 'here',fname,time_int
  ;endif 
     outfilename='data/bproj_vis/'+strtrim(string(flare.id[i]),1)+'_round1'
     im=do_clean(time_int,pix_size,det,image_dim,xyoffset,erange,outfilename=outfilename,quiet=quiet)
  return,im
end
