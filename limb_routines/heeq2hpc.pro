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
