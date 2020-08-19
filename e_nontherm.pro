FUNCTION e_nontherm,G,A,Ee0,spex=spex,thin=thin,SI=SI

; Function calculates the energy of nonthermal electrons in the
; thick-target assumption given the spectral index G, the low-energy cutoff
; Ee0, and the value A of the flux at 1keV. If /spex is set, A is at 50 keV.
; If /thin is set, the thin-target assumption is used.
;   RETURNS energy in ergs/sec; if /SI set, in J/sec

En = 9.5d24 * G^2 * (G-1) * beta(G-0.5,1.5,/double) * A * Ee0^(-G+1)

if keyword_set(thin) then En = En / G
if keyword_set(spex) then En = En * 5d1^G
if keyword_set(SI) then En = En * 1d-7

return,En

END
