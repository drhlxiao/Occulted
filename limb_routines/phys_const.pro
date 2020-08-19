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
