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
