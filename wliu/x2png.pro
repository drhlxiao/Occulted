;=========================================================================
; Purpose: read an x-window plot and write to .png image
;   cf. ~ssw/.../x2jpeg.pro => output file too large, say, 0.5 MB
; Calling sequence: 
;    IDL> xpng, mypng    ; mypng='mypng.png' output filename
; Author: Wei Liu
; 
; Notes on 2015/03/31:
; 1) tvrd(true=true) returns image of array 
;     [nx,ny] if true=0
;     [3,nx,ny] if true=1
;  either format is compabitable with write_png,
;
; 2) but plot_image needs [nx,ny] or [nx,ny,3];
;  so for the latter (3D true color), need 
;    img1= transpose(img, [2,0,1]) to convert img=array[nx,ny,3] used by plot_image to img1=array[3,nx,ny] used by write_png;
;    img= transpose(img1, [1,2,0]) to reverse it.
;  Ref: /home/weiliu/idl_lib/sdo_wei/Marc-DeRosa/tricolor_script_alan13.pro, line 903
;
; History: 2004/11/08, first written
;  2005/01/05: try to use r,g,b to get color png written; but not working,
;    unless the original x-window displace is in color.
;  2007/04/19: made color work with "tvlct", per Kim's suggestion.
;  2007/04/20: add /true to fix the color problem, per Kim's
;  2007/10/30: added "true" keyword, default=1, but allow to set 0
;=========================================================================

pro x2png, mypng, true=true

checkvar, true, 1  ; 2007/10/30

	;rr=indgen(256)  &  gg=rr  &  bb=rr 	; 2005/01/05

tmp=tvrd(true=true)		; 2007/04/20, add /true
write_png, mypng, tmp

	;--- commented 2007/04/20 ------
	;tvlct, rr, gg, bb, /get
	;tmp=tvrd()
	;write_png, mypng, tmp, rr, gg, bb, /verbose		
	   ; 2007/04/19, made it work with "tvlct"; 2005/01/05, add r,g,b, to get color png written, not working
	   ; write_png, mypng, tmp

end
