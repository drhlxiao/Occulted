;+
; PURPOSE:
;	To plot the solar limb on a map as seen from other observatories.
;
; EXAMPLES:
;	fits2map, 'work/20111216_Lovejoy/data/STEREO-B_EUVI/20111215_232030_n4euB.fts', map_stereoB
;	plot_map, map_stereoB, /LIMB
;	map_oplot_limb_from_earth, map_stereoB, color=1, linestyle=2
;
;	map_oplot_limb_from_earth, map_stereoB, color=1, linestyle=2, roll_angle_range=[180,360]
;
;	LOADCT,0 & plot_map, map_stereoB,/LIMB,/LOG
;	linecolors & map_oplot_limb_from_earth, map_stereoB, color=1, linestyle=2, /NEAR, thick=3
;	linecolors & map_oplot_limb_from_earth, map_stereoB, color=1, linestyle=2, /NEAR, thick=3, blr=[0,-90,45.*3600],/HARD
;
;
; KEYWORDS:
;	/PSEUDO: use this for plates (and large angular offsets). Otherwise returns straight angles!
;	METHOD 2 only implements the HARD limb for now.
;
;
; RESTRICTIONS:
;	Need ssw STEREO software to work.
;
;
; HISTORY:
;	PSH, 2013/05/22. Hastily written (<10 mins)!.
;	PSH, 2013/05/22. Added keyword roll_angle_range, NEAR_SIDE_ONLY 
;		(the latter could be improved by drawing to the exact poles as perceived by the observatory...).
;	PSH, 2013/05/28: improvements.
;	PSH, 2013/06/02: added METHOD=2. There is a minute difference (1-2") with METHOD 1 that I should track down...
;
;-
function map_oplot_limb_from_earth_f, map, blr=blr, factor=factor, _extra=_extra, $
	roll_angle_range=roll_angle_range, NEAR_SIDE_ONLY=NEAR_SIDE_ONLY, HARD_SPHERE=HARD_SPHERE, METHOD=METHOD, PSEUDO=PSEUDO

	default, roll_angle_range, [0,360]
	default, npts, 360L
	default, factor, 1d	;;1.0 for photospheric limb, about 1.01 for EUV limb...?
	default, METHOD, 1

	Rs=phys_const('R_sun')/phys_const('AU')		;;Rs in [AU]
	a=psh_binning(roll_angle_range[0],roll_angle_range[1],Npts)/radeg()

	IF NOT KEYWORD_SET(blr) THEN BEGIN
		pbr=pb0r(map.TIME,L0=L0,/ARCSEC)
		blr=[pbr[1],L0,pbr[2]]
	ENDIF

	d_sun = asin(Rs)/(map.RSUN/rad2asec())	;;in [AU]

	xx='bla'
	FOR i=0L, npts-1 DO BEGIN
		takeit=1B		
		IF METHOD EQ 1 THEN BEGIN
			;;Assuming we are aligned along HEEQ X-axis:
				IF KEYWORD_SET(HARD_SPHERE) THEN heeq=factor*Rs*[sin(blr[2]/rad2asec()),cos(blr[2]/rad2asec())*cos(a[i]),cos(blr[2]/rad2asec())*sin(a[i])] $	;;;Hard sphere
				ELSE heeq=[0,cos(a[i]),sin(a[i])]*Rs*factor		;;soft/transparent/translucid limb
	
			;;Now, rotate for L-angle around Z-axis, then B-angle around Y-axis:
				
				th=blr[0]/radeg()
				M = [[cos(th),0,sin(th)],[0,1,0],[-sin(th),0,cos(th)]]
				heeq = M # heeq
	
				th=blr[1]/radeg()
				M = [[cos(th),sin(th),0],[-sin(th),cos(th),0],[0,0,1]]
				heeq = M # heeq
	
				xxyyd=heeq2hpc(heeq,[map.ROLL_ANGLE,map.b0,map.L0,map.RSUN])
	
		ENDIF ELSE BEGIN
			IF METHOD EQ 2 THEN BEGIN
				hcr=[factor*Rs*cos(blr[2]/rad2asec()),a[i],factor*Rs*sin(blr[2]/rad2asec())]
				hcc=hcr2hcc(hcr)

				stony = hcc2stonyhurst(hcc,blr[0:1]/radeg())
				hcc=stonyhurst2hcc(stony,[map.B0,map.L0]/radeg())
				hpc = hcc2hpc(hcc, d_sun)
				xxyyd=[hpc[1]*rad2asec(),hpc[2]*rad2asec(),hpc[0]]
				
				;;The following paragraph seems to produce the exact same as the preceding one.
				
				;stony = hcc2stonyhurst(hcc,blr[0:1]/radeg())
				;heeq=stonyhurst2heeq(stony)
				;xxyyd=heeq2hpc(heeq,[map.ROLL_ANGLE,map.b0,map.L0,map.RSUN])

			ENDIF
		ENDELSE

		IF KEYWORD_SET(PSEUDO) THEN BEGIN
			hpc=[xxyyd[2],xxyyd[0]/rad2asec(),xxyyd[1]/rad2asec()]
			hpc_pseudo = hpc2pseudo(hpc)
			xxyyd=[hpc_pseudo[1]*rad2asec(),hpc_pseudo[2]*rad2asec(),hpc[0]]
		ENDIF

		IF KEYWORD_SET(NEAR_SIDE_ONLY) THEN BEGIN
			IF xxyyd[2] GT d_sun  THEN takeit=0B 
		ENDIF

		IF takeit THEN BEGIN
			IF datatype(xx) EQ 'STR' THEN BEGIN
				xx=xxyyd[0]
				yy=xxyyd[1]
			ENDIF ELSE BEGIN
				xx=[xx,xxyyd[0]]
				yy=[yy,xxyyd[1]]
			ENDELSE
		ENDIF
	ENDFOR;i

	IF KEYWORD_SET(NEAR_SIDE_ONLY) THEN BEGIN
		;;Need to shift around the gap, or it'll plot something going through SUn center too...
		;;Not ideal, but will do for now.
		gaps=get_edges(yy,/W)
		tmp=MAX(abs(gaps),ss)
		IF tmp GT 300. THEN BEGIN	;;Arbitrary 300. would not be ok if obs is VERY close to the Sun... Deal with it later...
			yy=SHIFT(yy,-ss-1)
			xx=SHIFT(xx,-ss-1)
		ENDIF
	ENDIF

        OPLOT, xx, yy, _extra=_extra
        coords={xx:xx,yy:yy}
        return, coords
END
