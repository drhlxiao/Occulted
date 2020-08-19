;EL: Relevant code is at the very end, sorry. All the routines up
;front are Pascal's.

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
FUNCTION map_oplot_limb_from_earth, map, blr=blr, factor=factor, _extra=_extra, $
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
        return, [[xx],[yy]]
     END

;+
; PURPOSE:
;	Given solar B0-angle (degrees), L0 (degrees), and solar radius size (in arcsecs), gives position of observer in HEEQ coordinates (units of AU).
;	
;	blr2heeq and heeq2blr appear to be extremely accurate angle-wise. I have some error (2.5%?) in distances, due to my lack of knowledge of exactly how big the Sun is at 1 AU.
;
; EXAMPLES:
;	PRINT, blr2heeq([0,0,959.6])
;	PRINT, blr2heeq([10,30,959.6])
;
;	pbr=pb0r('2008/01/01',SOHO=0,L0=L0,/ARCSEC)
;	PRINT, blr2heeq([pbr[1],L0,pbr[2]])
;
;	pbr=pb0r('2008/01/01',SOHO=0,L0=L0,/ARCSEC,STEREO='B')
;	PRINT, blr2heeq([pbr[1],L0,pbr[2]])
;	PRINT, get_stereo_coord( '2008-01-01T00:00:00', 'B' ,/NOVEL, SYSTEM='HEEQ')/696d3
;		;;~0.1% difference in absolute distances...?? 
;		;;But proportions are ok within the tenth digit!!!! Yeah!!!
;
; HISTORY:
;	PSH, 2010/12/31 written.
;	PSH, 2013/06/01: changed from "ATAN" (actually, not even!) to "ASIN."
;
;-
FUNCTION blr2heeq, blr
	
	;;0/ distance, in Rs:
		d=phys_const('R_sun')/phys_const('AU')/sin(blr[2]/rad2asec())
	;;1/ coordinates of observer in xyz system centered on Sun, where x points towards observer.
		XYZ=[1d,0,0]
	;;2/ rotate by -B-angle around y, to have new xy in same plane as final xy-plane.
		t=-blr[0]/RADEG()
		M=[[cos(t),0,sin(t)],[0,1,0],[-sin(t),0,cos(t)]]
		XYZ #= M
	;;3/ rotate by -L-angle around z, to have X-axis coincide:
		t=-blr[1]/RADEG()
		M=[[cos(t),sin(t),0],[-sin(t),cos(t),0],[0,0,1]]
		XYZ #= M

	;;Final result should be in HEEQ coordinates!
	RETURN, d*XYZ	
END

;+
; PURPOSE:
;	Converts from HEEQ coordinates to helioprojective-cartesian coordiantes.
;	HEEQ: XYZ, 
;	HCP: arcsecs from Sun center
;
; INPUTS:
;	xyz: HEEQ coordinates to convert (in units of AU)	
;	obs: pblr of observer (degrees/degrees/degrees/arcsecs). ACTUALLY, "P" is the ROLL ANGLE (as defined in SSW)!!!!
;
; OUTPUTS:
;	xxyydd: xx,yy: usual map's arcsecs from Sun-center. dd is distance from obs to obj/feature (same units as heeq input).
;	
;
; EXAMPLES:
;	PRINT, heeq2hpc([1,0,0]*phys_const('R_sun')/phys_const('AU'),[0,0,0,959.63])
;	PRINT, heeq2hpc([0,1,0]*phys_const('R_sun')/phys_const('AU'),[0,0,0,959.63],GR=0)
;	PRINT, heeq2hpc([0,1,0]*phys_const('R_sun')/phys_const('AU'),[0,0,0,959.63],GR=1)
;
;
;	PRINT, heeq2hpc([0,1,0]*phys_const('R_sun')/phys_const('AU'),[0,0,0,959.63],M=1)
;	PRINT, heeq2hpc([0,1,0]*phys_const('R_sun')/phys_const('AU'),[0,0,0,959.63],M=2)
;		;;Exactly the same :-( !! All this work implementing the Thompson stuff for nothing :-( :-( :-(
;		;;One the plus side, I had it right the very first time!
;		;;At least, I have added the /PSEUDO_ANGLE keyword..., useful for large angles... (on photographic plates...)
;
;
; HISTORY:
;	PSH, 2010/12/31 written.
;	PSH, 2011/07/22: added keyword GRcorrect. VERY APPROXIMATIVE!!!
;	PSH, 2011/12/04: de-activated GRcorrect: was very wrong!!!
;	PSH, 2013/06/01: added method 2, and /PSEUDO keyword.
;			For now, METHOD 2 does not do the roll angle correction.!
;
;-
FUNCTION heeq2hpc, heeq, pblr, METHOD=METHOD, PSEUDO=PSEUDO

	default, METHOD, 1	;;My original one
		;;1: My original one
		;;2: Using all the coord. transformations.

	IF METHOD EQ 1 THEN BEGIN
		;;1/Get observer's HEEQ coords.
			obs_heeq=blr2heeq(pblr[1:*])
	
		;;2/Get HEEQ coords of obs->obj vector:
			XYZ=-obs_heeq+heeq
		
		;;The idea is to rotate this vector into S/C (observer) coordinates, where X points away from Sun (towards S/C (observer)), and Z is map's "UP".
	
		;;3/Rotate around Z by L-angle:
			t=pblr[2]/RADEG()
			M=[[cos(t),sin(t),0],[-sin(t),cos(t),0],[0,0,1]]
			XYZ #= M
		;;4/Rotate around Y by B-angle:
			t=pblr[1]/RADEG()
			M=[[cos(t),0,sin(t)],[0,1,0],[-sin(t),0,cos(t)]]
			XYZ #= M
		;;5/Rotate around X by -P-angle:
			t=-pblr[0]/RADEG()
			M=[[1,0,0],[0,cos(t),-sin(t)],[0,sin(t),cos(t)]]
			XYZ #= M
		;;6/Transform into hpc, now that Z is map "UP":
			d=sqrt(TOTAL(XYZ^2.))
			xx=!DPI-ATAN(XYZ[1],XYZ[0])
			IF xx GT !DPI THEN xx-=2*!DPI
			yy=!DPI/2 -ACOS(XYZ[2]/d)
			;IF yy GT 2*!DPI THEN yy-=2*!DPI
			
		;;The above METHOD 1 never used /PSEUDO.
	ENDIF ELSE BEGIN			
		IF METHOD EQ 2 THEN BEGIN
			;;heeq->hcc->hcp->pseudo
			DD=rsun2dist(pblr[3])
			hpc=hcc2hpc(heeq2hcc(heeq,pblr[1]/radeg()),DD)
			xx=hpc[1]
			yy=hpc[2]
			d=hpc[0]				
		ENDIF
	ENDELSE
		
	IF KEYWORD_SET(PSEUDO) THEN BEGIN
		tmp1=[d,xx,yy]	
		tmp2=hpc2pseudo(tmp1)
		xx=tmp2[1]
		yy=tmp2[2]
	ENDIF	

	RETURN, [xx*rad2asec(),yy*rad2asec(),d]
END

;+
; PURPOSE:
;	Inverse of heeq2hcp.pro
;
; INPUTS:
;	xyz: 
;
;	pblr: of observer (degrees/degrees/degrees/arcsecs). "P" is actually SSW ROLL ANGLE!
;
; OUTPUTS:
;
; EXAMPLES:
;	pblr=[-20,-7.3,30,1500] & PRINT, hpc2heeq(heeq2hpc([-4,5,1],pblr),pblr)
;
;	PRINT, hpc2heeq(heeq2hpc([-4,5,1],pblr1),pblr2)
;
;
; HISTORY:
;	PSH, 2010/12/31 created.
;
;-
FUNCTION hpc2heeq, xxyydd, pblr


	;;1/From hpc to coordinate system where Z is map's "UP", X is along Sun-observer axis (away from Sun):
		d=xxyydd[2]
		xx=xxyydd[0]/3600.
		yy=xxyydd[1]/3600.
		phi=(180-xx)/RADEG()
		theta=(90-yy)/RADEG()
		XYZ=[0d,0,0]
		XYZ[0]=d*sin(theta)*cos(phi)
		XYZ[1]=d*sin(theta)*sin(phi)
		XYZ[2]=d*cos(theta)
	;;2/Rotate around X by +P-angle:
		t=pblr[0]/RADEG()
		M=[[1,0,0],[0,cos(t),-sin(t)],[0,sin(t),cos(t)]]
		XYZ #= M
	;;3/Rotate around Y by -B-angle:
		t=-pblr[1]/RADEG()
		M=[[cos(t),0,sin(t)],[0,1,0],[-sin(t),0,cos(t)]]
		XYZ #= M
	;;4/Rotate around Z by -L-angle:
		t=-pblr[2]/RADEG()
		M=[[cos(t),sin(t),0],[-sin(t),cos(t),0],[0,0,1]]
		XYZ #= M
	;;5/Add obs_heeq:
		obs_heeq=blr2heeq(pblr[1:*])
		XYZ+=obs_heeq

	RETURN, XYZ
END

;+
; HISTORY:
;       2011/07/15: finally written!
;
; PURPOSE:
;
;
;       Default values are in SI units.
;
;
; EXAMPLE:
;       PRINT, phys_const('AU',/CGS)
;       PRINT, phys_const('eps0')
;       PRINT, phys_const('GM_sun') / phys_const('R_sun')^3. *3600^2.   ;;in Rs^3/hr^2
;
;       PRINT, phys_const('R_sun',/CGS,/ASCII)+', '+phys_const('R_sun2',/CGS,/ASCII)
;
;-
FUNCTION phys_const, what, CGS=CGS, ASCII=ASCII

        IF KEYWORD_SET(CGS) THEN BEGIN
                CASE STRUPCASE(what) OF
                        'C': res = 29979245800d                         ;;[cm/s]
                        'E': res = 4.803206d-10                         ;;[esu]

                        'M_E': res=9.10939d-28                          ;;[g]
                        'M_P': res=1.67262d-24                          ;;[g]

                        'K_B': res=1.3806503d-16                        ;;[erg/K]
                        'G'  : res=6.6726d-8                            ;;[cm^3 g^-1 s^-2]
                        'H'  : res=6.626068d-27                         ;;[erg s]

                        'AU': res= 14959787069100d                      ;;[cm]
                        'AU1': res= 14959787069100d                     ;;[cm]
                        'AU2': res= 1495979d7                           ;;[cm]

                        'R_SUN':  res= phys_const(what)*1d2		;;[cm]
			'R_SUN1': res= phys_const(what)*1d2		;;[cm]
			'R_SUN2': res= phys_const(what)*1d2		;;[cm]
  
			'PC': res=3.08567758d18				;;[cm]	. Google.
              ENDCASE
        ENDIF ELSE BEGIN
                CASE STRUPCASE(what) OF
                        'C': res = 299792458d                           ;;[m/s]
                        'E': res = 1.60217646d-19                       ;;[C]
                        'EPS0': res=1d / 35950207149.4727056/!DPI       ;;[F/m]
                        'MU0': res=4*!DPI*1d-7                          ;;[H/m] or [(V.s)/(A.m)]

                        'M_E': res=9.10939d-31                          ;;[kg]
                        'M_P': res=1.67262d-27                          ;;[kg]

                        'K_B': res=1.3806503d-23                        ;;[J/K] or [m^2 kg s^-2 K^-1]
                        'G'  : res=6.6726d-11                           ;;[N kg^-2 m^2] or [m^3 kg^-1 s^-2]
                        'H'  : res=6.626068d-34                         ;;[J s] or [m^2 kg s^-1]
                        'H_KEV': res=phys_const('h')/1.6d-9		;;[keV.s]

                        'AU' : res= 149597870691d                       ;;[m]
                        'AU1': res= 149597870691d                       ;;[m]   ;;This number I found in SDO/AIA fits header... Also, used by in wcs_au.pro...
                        'AU2': res= 1495979d5                           ;;[m]   ;;According to Allen 4th ed. p.340

                        'R_SUN':  res= phys_const('R_SUN2')
			'R_SUN1': res= phys_const('R_SUN_ARCSEC1')/3600/RADEG()*phys_const('AU')	;;[m]   ;;~696 Mm.
			'R_SUN2': res= 695.508d6                        ;;[m]   ;;Wikipedia and other scientific sources, including Allen 4th ed (p.340)... Also used by wcs_rsun.pro. This is also what is observed by RHESSI SAS and HMI, looks like...
			'R_SUN3': res= 696d6                        	;;[m]   ;;696 Mm is what is used in SDO fits file.
					;;There is an inconsistency in Allen, 4th (p. 340): it gives Rsun, D, and asec, but they do not work out!
			'R_SUN_OBS_ARCSEC' : res= phys_const('R_SUN_OBS_ARCSEC2')
                        'R_SUN_OBS_ARCSEC1': res= 959.63d	;;[arcsec]	;;Allen 4th, p.340, "as viewed from Earth". Martin says 959.64".
                        'R_SUN_OBS_ARCSEC2': res= asin(phys_const('R_SUN2')/phys_const('AU'))*rad2asec()	;;[arcsec]	;;"asin" for observed, "atan" for true. Difference is about 10 mas.
                        'R_SUN_OBS_ARCSEC3': res= asin(phys_const('R_SUN3')/phys_const('AU'))*rad2asec()	;;[arcsec]
					;;10 mas difference between "observed" (hard sphere) and "true".
			'R_SUN_TRUE_ARCSEC' : res= phys_const('R_SUN_TRUE_ARCSEC2')
                        'R_SUN_TRUE_ARCSEC2': res= atan(phys_const('R_SUN2')/phys_const('AU'))*rad2asec()	;;[arcsec]	;;"asin" for observed, "atan" for true. Difference is about 10 mas.
                        'R_SUN_TRUE_ARCSEC3': res= atan(phys_const('R_SUN3')/phys_const('AU'))*rad2asec()	;;[arcsec]
						
			'R_MOON': res=1737.10d3			;;[m] ;;Ref: Wikipedia...

			'GM_EARTH': res=398600.4418d9		;;[m^3/s^2] .See Wikipedia for all planets and Sun...
			'GM_SUN': res=132712440018d9		;;[m^3/s^2] .See Wikipedia for all planets and Sun...

			'PC': res=3.08567758d16			;;[m]	. Google.

                ENDCASE
        ENDELSE

        IF KEYWORD_SET(ASCII) THEN res=strn(res,f='(e20.12)')
        RETURN, res
END

;+
; HISTORY:
;	PSH, 2005/06/07
;	PSH, 2009/10/01: added keyword powerlawindex
;
;
; EXAMPLES:
;	PRINT,psh_binning(1,100,101)		;;spacing is dE=cste=(E2-E1)/(n-1)
;	PRINT,psh_binning(1,100,101,/LOG)	;;spacing is dE= E1*(E2/E1)^(n/(n-1)) - E1*(E2/E1)^((n-1)/(n-2))
;	PRINT,psh_binning(1,100,101,POWER=3.)	;;spacing such that a power law with this index has equal flux in all bins...
;
;	;;TEST:
;		ebins=psh_binning(1,100,101,POWER=2.)
;		i=0  & x=psh_binning(Ebins[i],Ebins[i+1],10) & PRINT, INT_TABULATED(x,x^(-2.))
;		i=10 & x=psh_binning(Ebins[i],Ebins[i+1],10) & PRINT, INT_TABULATED(x,x^(-2.))
;		i=90 & x=psh_binning(Ebins[i],Ebins[i+1],10) & PRINT, INT_TABULATED(x,x^(-2.))
;
;-
FUNCTION psh_binning, startval, endval, nbrpts, LOG=LOG, powerlawindex=powerlawindex
	s=DOUBLE(startval[0])
	e=DOUBLE(endval[0])
	n=DOUBLE(nbrpts[0])

	IF KEYWORD_SET(LOG) THEN RETURN, s*(e/s)^(DINDGEN(n)/(n-1))
	IF KEYWORD_SET(powerlawindex) THEN BEGIN
		g=powerlawindex		;;assumed negative. Not tested for positive yet...
		a=1d - g
		Ftot = ( e^a - s^a )/a
		Fbin=Ftot/(n-1)
		Ebins=s^a
		FOR i=1L, n-1 DO Ebins=[Ebins,Ebins[N_ELEMENTS(Ebins)-1]+a*Fbin]
		RETURN, Ebins^(1./a)
	ENDIF
	
	RETURN, s+(e-s)*DINDGEN(n)/(n-1)
END


FUNCTION rad2asec, ang

	IF exist(ang) THEN RETURN, ang / !DPI * 180d * 3600 ELSE RETURN, 3600d * 180d /!DPI
	
END

;+
; PURPOSE:
;	To replace IDL's !RADEG, which has only a float precision, with a double precision.
;
;
;
;
;-
FUNCTION radeg

	RETURN, 180d / !DPI

END

;+
; PURPOSE:
;	To convert from Rsun (as returned by pb0r.pro, in arcsecs), into distance to Sun center, in [AU].
;
; EXAMPLE:
;	PRINT, rsun2dist((pb0r('2013/11/28 18:37',/ARCSEC))[2])*phys_const('AU')/phys_const('c')
;
;
;-
FUNCTION rsun2dist, rho

	RETURN, phys_const('R_sun')/SIN(rho/rad2asec())/phys_const('AU')

END


fits2map,'../../data/stereo-aia/20130427_105530_n4eua.fts',map
loadct,9
plot_map,map,center=[-700,0],fov=[2,3],/log
c1=map_oplot_limb_from_earth(map)   
c2=map_oplot_limb_from_earth(map,factor=1.01)

END
