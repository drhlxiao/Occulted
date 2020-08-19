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


