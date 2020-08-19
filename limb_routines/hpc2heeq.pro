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
