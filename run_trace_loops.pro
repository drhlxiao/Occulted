;program to calculate and plot loops between two helioprojective solar
;coordinates. Uses plots to do plotting, i.e. assumes that a map is
;plotted before calling this program.

;v1, Lucia Kleint, 21.10.2016
;v2, LK, 4.11.16, removed fac1/2 keyword, added ellipse1/2 to include tilting
;                 of ellipses. Added circle keyword, removed tilt keyword.
;                 By calling /circle it will draw a circle with tilt 0.
;                 circle=[1,30] draws circle with tilt 30. Added
;                 length and coordinates as output.
;v3, LK, 21.11.16 fixed issue with loop length calculation. arcsec do
;                 not make sense (esp. at limb). Use km/Mm instead.
;                 Note: cartesian coords are in units of rsun.
;
;INPUTS:
;date = REQUIRED    date (in usual CCSDS format) to obtain distance of sun, arcsec/px
;dx   = OPTIONAL    length of axes of crosses in arcsec (just for plotting)
;ellipse1 = OPTIONAL    =[fac1, tilt1] (units: R_circle, degrees)
;                   if set, ellipse with height fac1*radius_loop will be calculated
;ellipse2 = OPTIONAL    =[fac2, tilt2]
;                   if set, ellipse with height fac2*radius_loop will be calculated
;circle = OPTIONAL  =[1,tilt]
;                   if set, then circle with tilt will be drawn. If
;                   only /circle is set, then tilt is set to
;                   0. Replaces previous tilt keyword.
;                   tilt of loop wrt to surface normal (degrees)
;pos  = OPTIONAL    position of loop footpoints in 2x2 array ([0,0] and [1,0]
;                   are x and y of first point). If not given, user needs to 
;                   select footpoints by clicking on them manually
;legend = OPTIONAL  if set, then display legend for colors of circles/ellipses
;length = OPTIONAL  if set, length of circle and ellipses is returned
;                   in this variable. Order: circle, ell1, ell2 [km].
;cir_coord=OPTIONAL  output cartesian (HEE) coordinates of the circle
;                    (with tilt as specified)
;ell1_coord=OPTIONAL output cartesian (HEE) coordinates of ellipse1
;                    (with tilt as specified)
;ell2_coord=OPTIONAL output cartesian (HEE) coordinates of ellipse2
;                    (with tilt as specified)
;verbose = OPTIONAL  to show debug info
;
;OUTPUTS:
;length ->          loop length in km returned as [circle,ell1,ell2]. Values 
;                   that were not desired (e.g. if ell1 not specified) are set to zero.
;cir_coord ->       cartesian HEE coordinates of the circle in units
;                   of Rsun (i.e 1 means  6.957d5 km). X tow. Earth, Y
;                   tow. West, Z tow. North
;ell1_coord ->      cartesian HEE coordinates of ellipse1 (and 
;                   similarly for ellipse2)
;
;REQUIRES:
;-xy.pro to manually select points
;-cgcolor.pro to display colors
;-get_sun, wcs_au and wcs_conv_???_??? from SSW
;-tilt_loops_cartesian  for tilting loop calculations
;-al_legend from coyote idl routines
;
;RESTRICTIONS:
;- ellipse length is an approximation whose error bar I don't know
;
;EXAMPLE:
;   cir = [1,-30.]
;   ell1 = [1.5,-15.]
;   plot_map,maparr[ind],charsize=1.6
;   date = maparr[ind].time
;   pos = fltarr(2,2)
;   pos[0,0] = 80.3 & pos[1,0] = 187.6 & pos[0,1] = 109.8 & pos[1,1] = 164.1
;   run_trace_loops,dx=dx,date=date,circle=cir,pos=pos,ellipse1=ell1,/legend


PRO run_trace_loops,dx=dx,date=date,ellipse1=ellipse1,ellipse2=ellipse2,circle=circle,$
                    pos=pos,legend=legend,length=length,cir_coord=cir_coord,$
                    ell1_coord=ell1_coord,ell2_coord=ell2_coord,verbose=verbose,$
                    _extra=_extra
                    
if ~keyword_set(dx) then dx=1. ;length of axes of crosses in arcsec (just for plotting)
if ~keyword_set(date) then message,'Keyword date missing!'
if ~keyword_set(verbose) then verbose=0.
;if ~keyword_set(tilt) then tilt=0. ;no longer needed
rsun = 6.957d5 ;km (solar radius)

;assumes that there is a display with a map in arcsec coordinates
;if user does not give footpoint coordinates, select them interactively
    if ~keyword_set(pos) then begin
       print,'Select 2 footpoints of loop'
       xy,count=2,pos=pos
    endif


;draw footpoints
    plots,[pos[0,0]-dx,pos[0,0]+dx],[pos[1,0],pos[1,0]],color=cgcolor('red')
    plots,[pos[0,0],pos[0,0]],[pos[1,0]-dx,pos[1,0]+dx],color=cgcolor('red')
    plots,[pos[0,1]-dx,pos[0,1]+dx],[pos[1,1],pos[1,1]],color=cgcolor('red')
    plots,[pos[0,1],pos[0,1]],[pos[1,1]-dx,pos[1,1]+dx],color=cgcolor('red')

;calculate baseline of selected points
    baseline = sqrt((pos[0,0]-pos[0,1])^2.+(pos[1,0]-pos[1,1])^2.)
    ;print,'Loop baseline length [arcsec]:',baseline

;convert to spherical coordinates (=heliographic stonyhurst with R=1)
;helioprojective coordinates (as in HMI) have no height info and are
;therefore not unique (projection!).
    sun = get_sun(date)
    dsun_obs = sun[0]*wcs_au()  ;solar distance in AU->meters
    wcs_conv_hpc_hg,pos[0,0],pos[1,0],x1,y1,dsun_obs = dsun_obs,ang_units='arcseconds'
    ;default length units are meters
    wcs_conv_hpc_hg,pos[0,1],pos[1,1],x2,y2,dsun_obs = dsun_obs,ang_units='arcseconds'

;-------------------------------------------------------------------------------
;proper way: cartesian coordinates
;convert [x1,y1,1] and [x2,y2,1] [degrees] into cartesian x,y,z with R=rsun
  th_p = [y1,y2]*!dtor  ;radians
  ph_p = [x1,x2]*!dtor
  xcart = 1.D*cos(th_p)*cos(ph_p)
  ycart = 1.D*cos(th_p)*sin(ph_p)
  zcart = 1.D*sin(th_p)

  if verbose then print,'Cartesian X of endpoints:',xcart
  if verbose then print,'Cartesian Y of endpoints:',ycart
  if verbose then print,'Cartesian Z of endpoints:',zcart

;center of loop is slightly below surface (Aschwanden's STEREO loops
;all had centers above surface, not sure how much it matters).
  center = [avg(xcart),avg(ycart),avg(zcart)] ;middle of 2 points
  if verbose then print,'center point: ',center
  p1 = [xcart[0],ycart[0],zcart[0]] ;footpoint 1
  p2 = [xcart[1],ycart[1],zcart[1]] ;footpoint 2
  uvec = center/norm(center)    ;unit vector to loop apex
  nvec = (p1-p2)/norm(p1-p2)    ;unit vector along p1-p2
  zvec = crossp(uvec,nvec)      ;unit vector perpendicular to baseline and apex

  rad = norm(p2-center)         ;loop radius
  t = dindgen(101)/100*!pi
  cir = dblarr(101,3)
  for i=0,100 do cir[i,*] = rad * cos(t[i])*nvec + rad*sin(t[i]) * uvec + center ;create circle
  xpts = cir[*,0]
  ypts = cir[*,1]
  zpts = cir[*,2]
  if verbose then print,'Circle footpoints (tilt 0): ',reform(cir[0,*]),reform(cir[-1,*])
  ;iplot,xpts,ypts,zpts
  ;convert back to hg then hpc
   heights_aa = sqrt(xpts^2.+ypts^2.+zpts^2.)
   thetas_aa  = atan(zpts/sqrt(xpts^2.+ypts^2.))*!radeg
   phis_aa =    atan(ypts/xpts)*!radeg
   tmp = where(abs(phis_aa[1:*]-phis_aa[0:*]) ge !pi/2.) 
   if tmp[0] ne -1 then phis_aa[tmp+1:*] = phis_aa[tmp+1:*]+!pi
   wcs_conv_hg_hpc,phis_aa,thetas_aa,hpln_aa,hplt_aa,hecr=heights_aa,length_units='solRad',$
                   /arcsec,dsun_obs=dsun_obs
   if keyword_set(circle) then plots,hpln_aa,hplt_aa,color=cgcolor('dodger blue'),thick=3

  ;ellipses
  if keyword_set(ellipse1) then begin
     ;todo: test if keyword set properly (not string etc.)
     fac1 = ellipse1[0]
     ell1 = dblarr(101,3)
     for i=0,100 do ell1[i,*] = rad * cos(t[i])*nvec + fac1*rad*sin(t[i]) * uvec + center
     xell1pts = ell1[*,0]
     yell1pts = ell1[*,1]
     zell1pts = ell1[*,2]
     hts_ell1 = sqrt(ell1[*,0]^2. + ell1[*,1]^2. + ell1[*,2]^2.)
     th_ell1 = atan(ell1[*,2]/sqrt( ell1[*,0]^2. + ell1[*,1]^2.))*!radeg
     ph_ell1 = atan(ell1[*,1]/ell1[*,0])*!radeg
     wcs_conv_hg_hpc,ph_ell1,th_ell1,hpln_cell1,hplt_cell1,hecr=hts_ell1,$
                     length_units='solRad',/arcsec,dsun_obs=dsun_obs
     plots,hpln_cell1,hplt_cell1,color=cgcolor('orange'),thick=3
     ;calc eccentricity, find semimajor and -minor axes
     if fac1 ge 1 then major = fac1 else major = 1.
     if fac1 ge 1 then minor = 1. else minor = fac1   
     ecc1 = sqrt(1.-(minor/major)^2.)
     print,'Eccentricity1: ',ecc1
  endif
  if keyword_set(ellipse2) then begin
     fac2 = ellipse2[0]
     ell2 = dblarr(101,3)
     for i=0,100 do ell2[i,*] = rad * cos(t[i])*nvec + fac2*rad*sin(t[i]) * uvec + center
     xell2pts = ell2[*,0]
     yell2pts = ell2[*,1]
     zell2pts = ell2[*,2]
     hts_ell2 = sqrt(ell2[*,0]^2. + ell2[*,1]^2. + ell2[*,2]^2.)
     th_ell2 = atan(ell2[*,2]/sqrt( ell2[*,0]^2. + ell2[*,1]^2.))*!radeg
     ph_ell2 = atan(ell2[*,1]/ell2[*,0])*!radeg
     wcs_conv_hg_hpc,ph_ell2,th_ell2,hpln_cell2,hplt_cell2,hecr=hts_ell2,$
                     length_units='solRad',/arcsec,dsun_obs=dsun_obs
     plots,hpln_cell2,hplt_cell2,color=cgcolor('orange'),thick=3
     ;calc eccentricity
     if fac2 ge 1 then major = fac2 else major = 1.
     if fac2 ge 1 then minor = 1. else minor = fac2  
     ecc2 = sqrt(1.-(minor/major)^2.)
     print,'Eccentricity2: ',ecc2
  endif

;-------------- projections onto loop plane (top/side/front views) ---

;xcart,ycart,zcart
;uvec,nvec,zvec
;save,center,cir,xcart,ycart,zcart,uvec,nvec,zvec,filename='testdata.sav',/ve
;done in tilt_for_projviews

;--------------- tilted loops -----------------------------------
;tilt circular loop
  if keyword_set(circle) then begin
     if n_elements(circle) eq 2 then tilt = circle[1] else tilt=0.
     tilt_loops_cartesian,xpts,ypts,zpts,phout,thout,hout,tilt,outcartesian=cir
     phis_inc = phout
  ;this should no longer occur, but leave it just in case
     tmp = where(abs(phout[1:*]-phout[0:*]) ge !pi/2.) 
     if tmp[0] ne -1 then phis_inc[tmp+1:*] = phis_inc[tmp+1:*]-!pi
     thetas_inc = thout
     heights_inc = hout
     wcs_conv_hg_hpc,phis_inc,thetas_inc,hpln_inc,hplt_inc,hecr=heights_inc,$
                     length_units='solRad',/arcsec,dsun_obs=dsun_obs
     if tilt ne 0. then plots,hpln_inc,hplt_inc,color=cgcolor('blue'),thick=3,lines=2
     if verbose then print,'Circle footpoints (after tilt): ',reform(cir[0,*]),reform(cir[-1,*])

  endif

;tilt ellipse 1 
  if keyword_set(ellipse1) then begin
     if n_elements(ellipse1) eq 2 then tiltell1 = ellipse1[1] else tiltell1=0.
     tilt_loops_cartesian,xell1pts,yell1pts,zell1pts,phout1,thout1,hout1,tiltell1,$
                          outcartesian = ell1
     phis_inc1 = phout1
  ;this should no longer occur, but leave it just in case
     tmp = where(abs(phout1[1:*]-phout1[0:*]) ge !pi/2.) 
     if tmp[0] ne -1 then phis_inc1[tmp+1:*] = phis_inc1[tmp+1:*]-!pi
     thetas_inc1 = thout1
     heights_inc1 = hout1
     wcs_conv_hg_hpc,phis_inc1,thetas_inc1,hpln_inc1,hplt_inc1,hecr=heights_inc1,$
                     length_units='solRad',/arcsec,dsun_obs=dsun_obs
     ;if tilt is 0 then ellipse is already drawn
     if tiltell1 ne 0. then plots,hpln_inc1,hplt_inc1,color=cgcolor('red'),thick=3,lines=2
  endif

;tilt ellipse 2
  if keyword_set(ellipse2) then begin
     if n_elements(ellipse2) eq 2 then tiltell2 = ellipse2[1] else tiltell2=0.
     tilt_loops_cartesian,xell2pts,yell2pts,zell2pts,phout2,thout2,hout2,tiltell2,$
                          outcartesian = ell2
     phis_inc2 = phout2
  ;this should no longer occur, but leave it just in case
     tmp = where(abs(phout2[1:*]-phout2[0:*]) ge !pi/2.) 
     if tmp[0] ne -1 then phis_inc2[tmp+1:*] = phis_inc2[tmp+1:*]-!pi
     thetas_inc2 = thout2
     heights_inc2 = hout2
     wcs_conv_hg_hpc,phis_inc2,thetas_inc2,hpln_inc2,hplt_inc2,hecr=heights_inc2,$
                     length_units='solRad',/arcsec,dsun_obs=dsun_obs
     ;if tilt is 0 then ellipse is already drawn
     if tiltell2 ne 0. then plots,hpln_inc2,hplt_inc2,color=cgcolor('red'),thick=3,lines=2
  endif

;------ add legend ----------
if keyword_set(legend) then al_legend,['circle, tilt=0','circle, tilted',$
        'ellipse, tilt=0','ellipse, tilted'],thick=[3,3,3,3],lines=[0,2,0,2],$
        color=[cgcolor('dodger blue'),cgcolor('blue'),cgcolor('orange'),cgcolor('red')],$
        _extra=_extra
       ; textcol=cgcolor('white'),_extra=_extra


;---------------- calculate loop lengths ----------------------------

;half-circle
      ;;rloop = baseline/2.*(rsun/sun[1]) ;radius (half baseline) in km (projected!!!)
      ;;rloop_arcs = !pi*baseline/2.      ;loop length in arcsec (projected!!!)
      rloop = sqrt( (xcart[0]-xcart[1])^2. + (ycart[0]-ycart[1])^2. + $
                                               (zcart[0]-zcart[1])^2. )/2.*rsun
 
   if keyword_set(circle) then begin
      lloop = sqrt( (xcart[0]-xcart[1])^2. + (ycart[0]-ycart[1])^2. + $
                                               (zcart[0]-zcart[1])^2. )/2.*!pi*rsun ;km
      ;arcsec do not make sense (especially at limb)
      ;;;print,'Circular loop length [arcsec]:',rloop_arcs
      ;;;print,'Circular loop length [Mm]:',rloop*!pi/1000.
      ;the following lengths are all identical (tested)
      ;print,'Circular loop length [rsun]:',sqrt( (xcart[0]-xcart[1])^2. + (ycart[0]-ycart[1])^2. + $
      ;                                           (zcart[0]-zcart[1])^2. )/2.*!pi
      print,'Circular loop length [Mm]:',lloop/1000.
      ;;try also "integration" (identical to xcart,ycart,zcart method)
      ;ll=0
      ;for i=0,n_elements(xpts)-2 do ll = ll+$
      ;sqrt((xpts[i+1]-xpts[i])^2.+(ypts[i+1]-ypts[i])^2.+(zpts[i+1]-zpts[i])^2.)
      ;print,'Circular loop length [Mm]:',ll*rsun/1000.
   endif
;xpts in cartesian


;circumference of ellipses is an integral. Using Ramanujan's
;approximation (and /2 for half-ellipse):
  aa = rloop
  if keyword_set(ellipse1) then begin
     bb = rloop*fac1
   ; rell1 = !pi*(3*(aa+bb)- sqrt((3*aa+bb)*(aa+3*bb)))/2./(rsun/sun[1]) ;arcsec (projected!!!)
     rell1 = !pi*(3*(aa+bb)- sqrt((3*aa+bb)*(aa+3*bb)))/2. ;km
     print,'Ellipse1 length [Mm]:',rell1/1000.
      ;;try also "integration" (tested and identical to ramanujan)
      ;ll=0
      ;for i=0,n_elements(xell1pts)-2 do ll = ll+$
      ;sqrt((xell1pts[i+1]-xell1pts[i])^2.+(yell1pts[i+1]-yell1pts[i])^2.+(zell1pts[i+1]-zell1pts[i])^2.)
      ;print,'Ellipse1 length [Mm]:',ll*rsun/1000.

  endif
  if keyword_set(ellipse2) then begin
     bb2 = rloop*fac2
     rell2 = !pi*(3*(aa+bb2)- sqrt((3*aa+bb2)*(aa+3*bb2)))/2.
     print,'Ellipse2 length [Mm]:',rell2/1000.
  endif

;---------------- return parameters --------------------------------

;------- loop lengths -----------------
;   PRINT,'saving length in variable provided'
   ;set length to zero if circle/ellipse is not defined
   if n_elements(lloop) eq 0 then lloop=0.
   if n_elements(rell1) eq 0 then rell1=0.
   if n_elements(rell2) eq 0 then rell2=0.
   length = reform([lloop,rell1,rell2])

;---- cartesian coordinates for plotting of side/top/front views ---
   cir_coord = cir
   if keyword_set(ellipse1) then ell1_coord = ell1
   if keyword_set(ellipse2) then ell2_coord = ell2


END
