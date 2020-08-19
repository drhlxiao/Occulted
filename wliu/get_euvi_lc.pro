;
;Purpose: 
;  Integrate STEREO EUVI images (or a given FOV) to get the total flux, for flare class prediction
;
;Inputs:
;  files
;  
; optional inputs: 
;  fov=[ [x0,x1], [y0,y1] ]; default for full disk, or specify FOV, then take submap, then integrate full submap.
;  bkg_interval= '2013/10/11' + ['06:00', '06:10']
;
;References:
;  Nitta et al. (2013), http://adsabs.harvard.edu/abs/2013SoPh..288..241N, Solar Phys (2013) 288:241­-254, DOI 10.1007/s11207-013-0307-7
;    Added (2016/09/26) from below: Note: secchi_prep by default returns data array in photons/s for each pixel
;    So Nitta's paper figure had wrong units, should be all photons/sec, not DN/sec; this was clarified with Nitta back in 2014-2015ish.
;  
;To-do:
;
;History:
; 2014/04/14: written, workspace: 
;  /disk2/scr0/weiliu/sdo/events/01_global-wave/13-10-11_reflect-M15-occult/multiw/stereo/run2_total-flux
; 2014/11/22: 1) added option to re-select background interval, use previous run's result as an input, avoid re-do fits file reading;
;   2) used my new value_locate_range(), instead of value_locate()
; 2015/06/27: added keyword drot and ref_time to do drot_map(), useful for specified fov with long duration (solar rotation important)
;   coding followed my get_aia_time_slice.pro
; 2016/08/24: add bkg_time to find nearest neighbor (one data point) for background selection, rather than an average over bkg_interval
;   changed saving bkg_interval to bkg_times_used (new variable)
; 2016/09/26: updated notes above (in References) on units (photons/sec after secchi_prep)
; 2018/08/02, Thu: added option to /normalize_to_1AU (default=1); so anything before today, was not normalized, but flux at the STEREO location
;================================================================================================================

pro get_euvi_lc, files, fov=fov, bkg_interval=bkg_interval, bkg_time=bkg_time, win_size=win_size, new_window=new_window, outfile=outfile, no_sav=no_sav, _extra=_extra   $
    , input_sav=input_sav, new_bkg_interval=new_bkg_interval   $  ;alternative input, override above inputs
    , time_str=time_str, time_dbl=time_dbl, base_time=base_time, lc_euvi=lc_euvi  $  ; optional outputs
    , drot=drot, ref_time=ref_time, keep_limb=keep_limb	 $	; for drot_map(), 2015/06/27
    , normalize_to_1AU=normalize_to_1AU


;=== params ================================================================
checkvar, outfile, 'get_euvi_lc'
checkvar, win_size, [1200, 700]
checkvar, new_window, 1

checkvar, normalize_to_1AU, 1
AU= 149597870700d	; from https://en.wikipedia.org/wiki/Astronomical_unit

if keyword_set(input_sav) and keyword_set(new_bkg_interval) then begin
  ;=== case 0: sav file for lc already done ============================================================================
  restore, input_sav
  nfile= N_elements(time_str)
  time_dbl= anytim(time_str)
  bkg_interval=new_bkg_interval		; override old one

endif else begin
 ;=== case 1: fits file ============================================================================
  nfile= N_elements(files)
  lc0= dblarr(nfile)
  time_dbl= lc0
  exptime= lc0
  time_str= strarr(nfile)

  ;=== main loop ================================================================
  for i=0, nfile-1 do begin
    secchi_prep, files[i], index, data, /silent

    if keyword_set(fov) then begin
      index2map, index, data, map0

      ;--- 2015/06/27 ---------
      if keyword_set(drot) and keyword_set(ref_time) then begin	;--- drot: differential rotate map ---
        map= drot_map(map0, time=ref_time, keep_limb=keep_limb, _extra=_extra)
        map.time= map0.time		; recover original time after drot_map()
      endif else map=map0

      submap=get_sub_map(map, xrange=fov[*,0], yrange=fov[*,1])
      data= submap.data

    endif

    lc0[i]= total(data) 
    exptime[i]= index.exptime
		; Note: secchi_prep by default returns data array in photons/s for each pixel
		;   http://sohowww.nascom.nasa.gov/solarsoft/stereo/secchi/doc/euvi_prep.html
		; before, index.BUNIT= 'DN'
		; after secchi_prep, index.BUNIT= 'PhotonFlux'
		; don't need to /index.exptime, already done with secchi_prep

    ;--- optional (default=1): normalize_to_1AU, added 2018/08/02, Thu --------------------
    if keyword_set(normalize_to_1AU) then begin
      lc0[i] *= (index.dsun_obs / AU)^2
    endif

    if i eq 0 then begin
      if keyword_set(normalize_to_1AU) then where='1 AU' $
         else where='STEREO, @ dsun_obs= ' + trim(index.dsun_obs) + ' meter'
    endif
    
    ;--- get earth time correction -------
		;anytim( addtime(anytim(index.date_obs, /yohkoh), delta=index.ear_time/60.), /ecs)		; from Nitta's stereo_corr_time.pro
    time_earth= anytim(index.date_obs) + index.exptime/2. + index.ear_time		; now at center of exposure, also converted to Earth time to correct for different light travel time
    time_dbl[i]= time_earth
    time_str[i]= anytim(time_earth, /CCSDS)
  
  endfor

  base_time= time_str[0]
  time_dbl1= time_dbl - time_dbl[0]		;anytim(base_time)
endelse 
 ;=== end cases ============================================================================


;--- take out background ------------------------------------

									;if keyword_set(bkg_interval) then bkg_idx= value_locate(time_dbl, anytim(bkg_interval)) >0   $
	;if keyword_set(bkg_interval) then bkg_idx= value_locate_range(time_dbl, anytim(bkg_interval))   $
	; else bkg_idx=[0, 2]

case 1 of
  keyword_set(bkg_interval): bkg_idx= value_locate_range(time_dbl, anytim(bkg_interval))
  keyword_set(bkg_time):     bkg_idx= value_locate_nearest(time_dbl, anytim(bkg_time)) * [1,1]
  else: bkg_idx=[0, 2]
endcase

lc_bkg= mean(lc0[bkg_idx[0]:bkg_idx[1]])
lc_euvi = lc0 - lc_bkg

bkg_times_used = anytim(time_dbl[bkg_idx], /ecs)
box_message, 'bkg_times_used: ' + bkg_times_used[0] + ' - ' + bkg_times_used[1]


;=== output data ==============================================================

if keyword_set(outfile) then begin
  openw, 1, outfile+ '.dat'
  printf, 1, 'time_dbl1 (sec),      time (at Earth),   flux (at '+where+ ') w/o and w/ background (photons/sec), exptime (s);  bkg_times_used, lc_bkg: ', bkg_times_used, lc_bkg, format='(a, 2(a,x), g14.7)'
  for i=0L, nfile-1  do printf,1, time_dbl1[i], time_str[i], lc_euvi[i], lc0[i], exptime[i], format='(g13.6, 2x, a, 2x, 2(g14.7,x), 1x, g13.6)'
  close,1
endif

if not keyword_set(no_sav) then save, filename=outfile + '.sav' $
               , time_str, time_dbl1, base_time, lc_euvi, lc0, exptime, bkg_times_used, lc_bkg  $			;bkg_interval
  , description='time_str, time_dbl1, base_time, lc_euvi, lc0, exptime, bkg_times_used, lc_bkg' $
  , /compress


;--- plot --------------------------------------
if keyword_set(new_window) then window, 11, xs=win_size[0], ys=win_size[1]
utplot, time_dbl1, lc_euvi, base_time, title='STEREO EUVI Integrated Flux', ytitle='Photon Flux (photons/sec)', charsize=1.4, xstyle=1, _extra=_extra
x2png, outfile+ '.png'

end
