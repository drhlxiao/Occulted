

PRO Goes_teem, goes_data, te, em, itim, dti, excess, bflux, btim, dtb, Quiet=quiet, $
               st_times=st_times, int_times=int_times, bck_time=bck_time, bck_dt=bck_dt, $
               dt_resolution=dt_resolution, auto_bck=auto_bck, auto_tim=auto_tim, $
               nbcks=nbcks, nbsigmas=nbsigmas, channel=channel, alt_tsel=alt_tsel, $
               zzb=zzb, auto_peak=auto_peak, schwartz=schwartz, use_goes_bck=use_goes_bck, $
               xxx=xxx, goes8=goes8, goes9=goes9, goes6=goes6, goes7=goes7, garcia=garcia, $
               satellite=satellite

   
   IF(KEYWORD_SET(quiet)) THEN qt = 1 ELSE qt = 0
;Ok, there are two types of structure that you need to deal with,
;Data from RD_GXD has the structure tag name 'GXD_DATA_REC', data
;from READ_GOES_FITS is an anonymous structure.
   strtype = strupcase(tag_names(goes_data, /str))
   IF(strtype EQ 'GXD_DATA_REC') THEN BEGIN
      tim_arr1 = anytim2ints(goes_data)
      flux = [transpose(goes_data.lo), transpose(goes_data.hi)]
   ENDIF ELSE BEGIN
      index = goes_data.time
      typ = size(index) & typ = typ(typ(0)+1)
;the time may be in seconds from 1970
      IF(typ EQ 5) THEN tim_arr1 = t70tot79(index) ELSE tim_arr1 = anytim2ints(index)
      flux = goes_data.flux
   ENDELSE
   
;no. of channels
   nc0 = 2
   nc01 = 1

   tim_arr = int2secarr(tim_arr1) ;time in seconds
   ndset = N_ELEMENTS(tim_arr1)
   nds1 = ndset-1
   dx = tim_arr(1:nds1)-tim_arr(0:nds1-1) ;interval time
   dx = [dx, dx(nds1-1)]
   
   tp1 = fmt_tim(tim_arr1(0), /msec)
   tp2 = fmt_tim(tim_arr1(nds1), /msec)
   IF(qt EQ 0) THEN BEGIN
      print, 'TIME RANGE'
      print, tp1, ' TILL ', tp2
   ENDIF
;open output file
   IF(KEYWORD_SET(outfile)) THEN ofl = outfile ELSE ofl = 0
;and the possiblility of a plot file
   IF(KEYWORD_SET(pfile)) THEN pfl = pfile ELSE pfl = 0

   time0 = tim_arr1(0)

;If i need to choose intervals by plotting
   IF(KEYWORD_SET(channel)) THEN ich0 = fix(channel(0)) ELSE ich0 = 0
   gzz = flux
   y0 = reform(gzz(ich0, *))    ;goes data comes in energy flux...
   plot_ttl = 'GOES FLUX, CHANNEL'+strcompress(string(ich0))
   FOR j = 0, 1 DO gzz(j, *) = gzz(j, *)*dx ;you want total, not the rate

;choose the background time interval
   IF(KEYWORD_SET(zzb)) THEN BEGIN ;no background
      IF(qt EQ 0) THEN print, 'BACKGROUND SET TO ZERO, FOR TESTING...'
      bcnts = fltarr(nc0)
      dtb = 60.0                ; gotta have it
      btim = anytim2ex(tim_arr1(0)) ; gotta have it
   ENDIF ELSE BEGIN
      IF(KEYWORD_SET(auto_bck)) THEN BEGIN
         r = auto_bck_find(tim_arr1, gzz(ich0, *), b_trange, dt_dat = dx, $
                           quiet = qt, /subscripts, dtb_min = 30.0, $
                           dtb_max = 120.0, /use_minimum)
         ib1 = b_trange(0) & ib2 = b_trange(1)
      ENDIF ELSE BEGIN
         IF(qt EQ 0) THEN print, 'BACKGROUND INTERVAL SELECTION::'
         IF(KEYWORD_SET(bck_time)) THEN BEGIN
            s_bckt = size(bck_time)
            IF(s_bckt(0) GT 2) THEN bck_time0 = bck_time(*, *, 0) ELSE bck_time0 = bck_time ;just in case
            bck_time0 = anytim2ints(bck_time0)
            bck_time0 = bck_time0(0) ;just in case...
            IF(KEYWORD_SET(bck_dt)) THEN dtb = bck_dt(0) ELSE dtb = 120.0
            choose_interval, tim_arr1, ib1, ib2, data = y0, $
              plot_title = plot_ttl+ ' CHOOSE BACKGROUND', $
              st_times = bck_time0, int_times = dtb, dt_resolution = -1, alt_tsel = alt_tsel
         ENDIF ELSE BEGIN
            plot_pr0 = ['CHOOSE ANY INTERVAL']
            choose_interval, tim_arr1, ib1, ib2, data = y0, $
              plot_title = plot_ttl+ ' CHOOSE BACKGROUND', $
              dt_resolution = -1, alt_tsel = alt_tsel, plot_prompt = plot_pr0
         ENDELSE
      ENDELSE

;background start times, and dt
      btim = anytim2ex(anytim2ints(time0, off = tim_arr(ib1)))
;interval dt's
      dtb = accum_counts(dx, ib1, ib2)
      IF(qt EQ 0) THEN BEGIN
         tb1 = fmt_tim(anytim2ints(time0, off = tim_arr(ib1)), /msec)
         tb2 = fmt_tim(anytim2ints(time0, off = tim_arr(ib2)+dx(ib2)), /msec)
         print, ' BACKGROUND FROM ', tb1,' TO  ', tb2
      ENDIF
;ok, accumulate the counts
      IF(KEYWORD_SET(use_goes_bck)) THEN BEGIN
         bcnts = dtb*goes_bck0(tim_arr1, flux, [ib1, ib2], quiet = qt)
      ENDIF ELSE bcnts = accum_counts(gzz, ib1, ib2)
   ENDELSE

;Now we choose the time intervals
   IF(KEYWORD_SET(auto_peak)) THEN BEGIN
      y0max = max(y0, ss_st)   ; that's it
      ss_en = ss_st
      print, 'THE PEAK TIME: ', fmt_tim(tim_arr1(ss_st), /msec)
      IF(auto_peak(0) EQ -1) THEN BEGIN
         utplot, tim_arr1, y0, title = plot_ttl, ytitle='Counts/sec'
         outplot, [tim_arr1(ss_st), tim_arr1(ss_st)], [!Cymin, !Cymax]
         yesnox, 'IS THIS AN OK PEAK TIME?', yn, 'y'
         IF(yn EQ 0) THEN choose_interval, tim_arr1, ss_st, ss_en, data = y0, $
           plot_title = plot_ttl, dt_res = -1
      ENDIF
;to automatically choose the time, use the background rate at the maximum
   ENDIF ELSE IF(KEYWORD_SET(auto_tim)) THEN BEGIN
      brate = bcnts(ich0)/dtb(0) ;be sure about scalars
      IF(KEYWORD_SET(nbsigs)) THEN BEGIN
         ex_rate = float(nbsigs)*sqrt(brate/dx)+brate
      ENDIF ELSE IF(KEYWORD_SET(nbcks)) THEN BEGIN
         ex_rate = (1.0+nbcks)*brate
      ENDIF ELSE BEGIN
         ex_rate = 2.0*brate
      ENDELSE
      choose_interval, tim_arr1, ss_st, ss_en, data = y0, b_level = ex_rate, dt_res = dt_resolution, xxx = xxx
;flare times
   ENDIF ELSE BEGIN
      print, 'FIT INTERVAL SELECTION::'
      choose_interval, tim_arr1, ss_st, ss_en, data = y0, plot_title = plot_ttl, $
        st_times = st_times, int_times = int_times, dt_res = dt_resolution, alt_tsel = alt_tsel
   ENDELSE

   IF(ss_st(0) EQ -1) THEN BEGIN
      IF(qt EQ 0) THEN message, /info, 'NO INTERVAL SELECTED, RETURNING'
      te = 0
      em = 0
      itim = 0
      dti = 0
      excess = 0
      bflux = 0
      btim = 0
      dtb = 0
      RETURN
   ENDIF
   
;ok, accumulate the counts
   cnts = accum_counts(gzz, ss_st, ss_en)
;interval dt's
   dti = accum_counts(dx, ss_st, ss_en)
   IF(N_ELEMENTS(dti) GT 1) THEN dti = reform(dti)

;interval start times
   itim = anytim2ex(anytim2ints(time0, off = tim_arr(ss_st)))

;ok, now, cnts
   ntim = N_ELEMENTS(itim(0, *))
   excess = cnts
   bflux = bcnts/dtb(0)
   te = fltarr(ntim)
   em = te
   FOR j = 0, ntim-1 DO excess(*, j) = cnts(*, j)/dti(j)-bflux
   ok = where(excess(1, *) GT 0.0 AND excess(0, *) GT 0.0)
   IF(ok(0) NE -1) THEN BEGIN
      IF(KEYWORD_SET(schwartz)) THEN BEGIN
         goes_tem, reform(excess(0, ok)), reform(excess(1, ok)), teok, emok, $
           goes8 = goes8, goes9 = goes9, goes6 = goes6, goes7 = goes7, $
           satellite=satellite, date = tim_arr1(0), garcia = garcia
         te(ok) = teok
         em(ok) = 49.0+alog10(emok)
      ENDIF ELSE BEGIN
         ratio = reform(excess(1, ok)/excess(0, ok))
         te(ok) = goest_eqn(ratio)
         em(ok) = goesem_eqn(te(ok), reform(excess(0, ok)))
      ENDELSE
   ENDIF ELSE print, 'NO GOOD INTERVALS'

   RETURN
END



