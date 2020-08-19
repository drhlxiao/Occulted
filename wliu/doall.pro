; @doall

wave='195'

 @call_get_euvi_lc.pro

fmax= max(lc_euvi)
goes_flux= 1.39e-11 * fmax
goes_class= goes_flux/1e-5  ; in M-class

;--- take from commandline ------
print, fmax/1e6, goes_class   ; apply Nitta et al. 2013 Sol Phy, Eq. (1) to Watt m^-2, then /1e-5 for M flare, http://adsabs.harvard.edu/abs/2013SoPh..288..241N
;     15.326733      
;     21.304159   => X2.1 (for flare FOV only, not quite applicable for Nitta's formula - take X2.2 at 11:10 UT from my full-disk version ../195/notes.txt, before CCD snow).

; so this translates to an XXX flare if viewed on the disk.

;--- plot again -------
time_dbl1 = time_dbl - time_dbl[0]
utplot, time_dbl1, lc_euvi, base_time, title='STEREO EUVI-' +wave+ ' Angstrom Flare-region Integrated Flux', ytitle='Photon Flux (photons/sec)', charsize=1.4, ylog=0, ystyle=1, yrange=get_range(lc_euvi)
legend, ['max flux F_EUVI (phs/s)= ' + trim(fmax), 'F_GOES (W/m^2) = 1.39 x 10^{-11} F_EUVI = ' + trim(goes_flux), 'GOES class (M)= '+trim(goes_class) $
        , 'See Nitta et al. (2013), Eq. (1), Solar Phys., 2013SoPh..288..241N'], /right, charsize=1.7, box=0
oplot, !x.crange, [1,1]*0, line=2

x2png, 'get_euvi_lc.png'
