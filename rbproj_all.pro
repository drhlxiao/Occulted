;make back projection images for each detector in each time range in
;order to check which detectors to use for the final CLEAN

;to make imaegs
;make imaging object
pro make_image,starttime, outfilename,outfitsname,nt,det_index_mask,ebins
  if not keyword_set(nt) then  nt=22 ; 1.5 minute;nt=11
  if not keyword_set(det_index_mask) then det_index_mask=[0,0,0,0,0,1,0,1,1]
  if not keyword_set(ebins) then ebins=[10,30]
o=hsi_image()
;set time range as multiple of spin period
ttt=4.110168
o->set,im_time_interval=anytim(starttime)+[0,nt*ttt];)'2013/05/01 02:32:00')+[0,11*ttt]
ptim,o->get(/im_time_int)
;set energy range
;o->set,im_energy_binning=[4,8]
o->set,im_energy_binning=ebins
;set output filename
o->set,im_out_fits=outfilename;'clean_20130501_0232_0233_1min_10-30keV_g689_fd.fits'
;run clean for given parameters
;smaller pixel size is needed if finer subcollimators are used
;im=o->getdata( image_alg='clean',$
;           pixel_size=[1,1]*3,image_dim=[128,128]*2,xyoffset=[-950,255],$
;           det_index_mask=det_index_mask,$
;           clean_frac=0.1,clean_niter=100,clean_mark=0)
im=o->getdata( image_alg='vis_fwdfit',$
           pixel_size=[1,1]*3,image_dim=[128,128]*2,xyoffset=[-950,255],$
               det_index_mask=det_index_mask,$
              vf_circle=1)
;im=o->getdata( image_alg='bproj',$
;           pixel_size=[1,1]*3,image_dim=[128,128]*2,xyoffset=[-950,255],$
;               det_index_mask=det_index_mask)
;           ;clean_frac=0.1,clean_niter=100,clean_mark=0)

;to display image
o->fitswrite, this_out_filename=outfitsname
o->plotman
end

pro fwdfit_ebins
  times=['2013/05/01 02:32:00','2013/05/01 02:33:36','2013/05/01 02:35:12']
  ebins=[[10,13],[13,16],[16,19],[19,22],[22,25],[25,30]]
  ;dets=[0,1,2,3,4,5,6,7,8]
  for t=0,1 do begin
     for i=0,5 do begin
        ebin=ebins[*,i]
        ofname='vff_t'+strtrim(string(t),1)+'_ebin_'+strtrim(string(ebin[0]),1)+'-'+strtrim(string(ebin[1]),1)+'_96s'
        ofitsname='vff_t'+strtrim(string(t),1)+'_ebin_'+strtrim(string(ebin[0]),1)+'-'+strtrim(string(ebin[1]),1)+'_96s.fits'
        ;dmask=make_array(1,9,/integer,value=0)
        dmask=[0,0,0,0,0,0,1,1,1]
        ;print, dmask,ofname,ofitsname,times[t]
        make_image,times[t],ofname,ofitsname,24,dmask,ebin
     endfor
  endfor
end


pro bproj_hi
  ;times=['2013/05/01 02:32:00','2013/05/01 02:32:48','2013/05/01 02:33:36','2013/05/01 02:34:24','2013/05/01 02:35:12','2013/05/01 02:36:00']
  times=['2013/05/01 02:32:00','2013/05/01 02:33:36','2013/05/01 02:35:12']
  ;dets=[0,1,2,3,4,5,6,7,8]
  for t=0,2 do begin
     for i=0,8 do begin
        ofname='bproj_hi_t'+strtrim(string(t),1)+'_det_'+strtrim(string(i+1),1)+'_96s'
        ofitsname='bproj_hi_t'+strtrim(string(t),1)+'_det_'+strtrim(string(i+1),1)+'_96s.fits'
        dmask=make_array(1,9,/integer,value=0)
        dmask[i]=1
        ;print, dmask,ofname,ofitsname,times[t]
        make_image,times[t],ofname,ofitsname,24,dmask
     endfor
  endfor
end

;good ones, time cadence 48s
;t[0]: det 8&9
;t[1]: 6?,7,8,9
;t[2]: 8,9
;t[3]: 8,9
;t[4]: 9 (8 shows different location)

;good ones, time cadence 96s
;t[0]: 6,7 (2Kcounts but different location shown?),8,9
;t[1]: 5?,6,7,8,9
;t[2]: 7? 8,9


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
  
