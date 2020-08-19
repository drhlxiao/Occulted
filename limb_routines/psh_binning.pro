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
