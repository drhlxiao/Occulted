;within sswidl
;check aspect solution and get spin period

;tr=['2013/05/01 02:32:00','2013/05/01 02:37:50']
;oa=hsi_aspect_solution()
;oa->set,as_roll_solution='PMT'
;oa->set,/as_no_extra
;oa->set,aspect_cntl_level=6
;as=oa->getdata(obs_time_interval=tr,/no_wdelete)
;obj_destroy,oa

;to make imaegs
;make imaging object
pro make_image,starttime, outfilename,outfitsname,nt,det_index_mask,image_alg
  if not keyword_set(nt) then  nt=22 ; 1.5 minute;nt=11
  if not keyword_set(det_index_mask) then det_index_mask=[0,0,0,0,0,1,0,1,1]
  if not keyword_set(image_alg) then image_alg='clean'
o=hsi_image()
;set time range as multiple of spin period
ttt=4.110168
o->set,im_time_interval=anytim(starttime)+[0,nt*ttt];)'2013/05/01 02:32:00')+[0,11*ttt]
ptim,o->get(/im_time_int)
;set energy range
o->set,im_energy_binning=[13,30]
;set output filename
o->set,im_out_fits=outfilename;'clean_20130501_0232_0233_1min_10-30keV_g689_fd.fits'
;run clean for given parameters
;smaller pixel size is needed if finer subcollimators are used
if image_alg eq 'clean' then begin
  im=o->getdata( image_alg=image_alg,$
           pixel_size=[1,1]*3,image_dim=[128,128]*2,xyoffset=[-950,255],$
           det_index_mask=det_index_mask,$
                 clean_frac=0.1,clean_niter=100,clean_mark=0)
endif else begin
  im=o->getdata( image_alg=image_alg,$
           pixel_size=[1,1]*3,image_dim=[128,128]*2,xyoffset=[-950,255],$
               det_index_mask=det_index_mask,$
               vf_circle=1) ;try changing the center to the center where fwdfit says it is
           ;clean_frac=0.1,clean_niter=100,clean_mark=0)
endelse
;to display image
o->fitswrite, this_out_filename=outfitsname
save,o,outfilename
o->plotman
end

pro put_in_struct
  fns=['clean_low1.fits','clean_low2.fits','clean_low3.fits','clean_low4.fits','clean_low5.fits','clean_low6.fits','clean_low7.fits','clean_low8.fits','clean_low9.fits','clean_low10.fits']
  fits2map,fns[0],map1
  fits2map,fns[1],map2
  fits2map,fns[2],map3
  fits2map,fns[3],map4
  fits2map,fns[4],map5
  fits2map,fns[5],map6
  fits2map,fns[6],map7
  fits2map,fns[7],map8
  fits2map,fns[8],map9
  ;fits2map,fns[9],map10
  
  ;lmst={one: map1,$
  ;      two: map2,$
  ;      three: map3, $
  ;      four: map4, $
  ;      five: map5}
  save, map1,map2,map3,map4,map5,map6,map7,map8,map9,filename='test_st.sav'
end

pro put_high
  fns=['clean_test.fits','clean_test2.fits','clean_test3.fits','clean_test4.fits','clean_test5.fits','clean_test6.fits']
  fits2map,fns[0],map1
  fits2map,fns[1],map2
  fits2map,fns[2],map3
  fits2map,fns[3],map4
  fits2map,fns[4],map5
  fits2map,fns[5],map6

  ;fits2map,fns[9],map10
  
  ;lmst={one: map1,$
  ;      two: map2,$
  ;      three: map3, $
  ;      four: map4, $
  ;      five: map5}
  save, map1,map2,map3,map4,map5,map6,filename='rhsi_hi_45s.sav'
  end
  
