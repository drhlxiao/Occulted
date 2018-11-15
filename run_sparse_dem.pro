; Author: Mark Cheung
; Purpose: Sample script for running sparse DEM inversion (see Cheung
; et al. 2015) for details.
; Revision history: 2015-10-20 First version
;                   2015-10-23 Added note about lgT axis
;                   2017-03-25 For initialization, use basis function
;                   sets [Dirac,sigma=0.1,sigma=0.2] instead of the
;                   default, which would also include sigma=0.6.

function run_sparse_dem,files,submap

    fname=strmid(files[0],8,15)
    headarr=list(length=6)
    ;submap=[-1250,-800,0,400]   ;coords of submap
    imgarr=fltarr(700,667,6) ;find the actual size of the map
    ;get headers and data, make the submap
    for i=0,5 do begin
       ;print, i
       mreadfits,files[i],header,data
       fits2map,files[i],map
                                ;for some reason still not correct header format...
       headf=headfits(files[i])
       hstruct=fitshead2struct(headf)
       sub_map,map,smap,/plot,xrange=[submap[0],submap[1]],yrange=[submap[2],submap[3]]
       ;plot_map,smap,/log
       sdata=smap.data
       ;sdata=data[submap[0]:submap[1]-1,submap[2]:submap[3]-1]
       headarr[i]=hstruct
       ;print,headarr[i].WAVELNTH
       imgarr[*,*,i]=sdata
    endfor

    ;make the structure
    ;** Structure MS_259678195004, 4 tags, length=5771872, data length=5771518:
    ;   IMG             FLOAT     Array[600, 400, 6]
    ;   OINDEX          STRUCT    -> MS_259678195003 Array[6] ;headers
    ;   made using fitshead2struct
    ;   BINNING         INT              1
    ;   FOV             INT       Array[4]  ;coordinates of box if submap?

    ; Note!!! The solver assumes the third dimension (of IMG) is arranged according to [94,131,171,193,211,335]
    s={IMG:imgarr,OINDEX:headarr,BINNING: 1 ,FOV:submap}

    ; Initialize solver.  
    ; This step builds up the response functions (which are time-dependent) and the basis functions.
    ; If running inversions over a set of data spanning over a few hours or even days, it is not necessary to re-initialize

    ;IF (0) THEN BEGIN 
    ; As discussed in the Appendix of Cheung et al. (2015), the inversion
    ; can push some EM into temperature bins if the lgT axis goes below
    ; logT ~ 5.7 (see Fig. 18). It is suggested the user check the
    ; dependence of the solution by varying lgTmin.
    lgTmin = 5.7   ; minimum for lgT axis for inversion 
    dlgT   = 0.1   ; width of lgT bin
    nlgT   = 21    ; number of lgT bins
    CATCH, Error_status
    ;print, Error_status
    if Error_status eq 0 then begin
       aia_sparse_em_init, timedepend = s.oindex[0].date_obs, /evenorm, use_lgtaxis=findgen(nlgT)*dlgT+lgTmin, bases_sigmas=[0.0,0.1,0.2]
       lgtaxis = aia_sparse_em_lgtaxis()
    endif else return,'init'
    ;print,'here'
    ; ENDIF

    ; We will use the data in s.img to invert. s.img is a stack of level
    ; 1.5, exposure time normalized AIA pixel values. 
    exptimestr = '['+strjoin(string(s.oindex[0].exptime,format='(F8.5)'),',')+']'
    ;print,exptimestr
    ; This defines the tolerance function for the inversion. 
    ; y denotes the values in DN/s/pixel (i.e. level 1.5, exposure time
    ; normalized)

    ; aia_bp_estimate_error gives the estimated uncertainty for a given
    ; DN / pixel for a given AIA channel. Since y contains DN/pixel/s,
    ; we have to multiply by exposure time first to pass
    ; aia_bp_estimate_error. Then we will divide the output of
    ; aia_bp_estimate_error by the exposure time again.
    ; If one wants to include uncertainties in atomic data, suggest to add
    ; the /temp keyword to aia_bp_estimate_error()
    tolfunc = 'aia_bp_estimate_error(y*'+exptimestr+', [94,131,171,193,211,335], num_images='+strtrim(string(s.binning^2,format='(I4)'),2)+')/'+exptimestr

    ; Do DEM solve. 
    ; Note!!! The solver assumes the third dimension is arranged according to [94,131,171,193,211,335]
    CATCH, Error_status
    ;print, Error_status
    if Error_status eq 0 then begin
       aia_sparse_em_solve, s.img, tolfunc=tolfunc, tolfac=1.4, oem=emcube, status=status, coeff=coeff
       print, 'done solving'
    endif else return,'solve'
    ; emcube contains the emission measure contained in each lgt bin
    ; status contains a mask indicating whether a solution was
    ; found. status = 0 means a solution was found within the given constriants.
    ; If you want to actual coefficients of the basis functions
    ; (i.e. x_i's of Kx = y), then add coeff=coeff to the call to aia_sparse_em_solve
    CATCH, Error_status
    ;print, Error_status
    if Error_status eq 0 then begin
       ; Now show results
       dispimage= aia_sparse_em_display(em=emcube)
    ;print, 'dispimage made'
       ;window, xs=1440, ys=900
       ;tv,dispimage, /true
    ;print, 'before make image'
       ;make image
    ;CATCH, Error_status
    ;print, Error_status
       aia_sparse_em_em2image,emcube,image=image,status=status    
    ;print, 'after make image'
    ;CATCH, Error_status
    ;print, Error_status
    ;print, fname
       ;save it all
       save, emcube,status,coeff,dispimage,image,lgtaxis,filename=fname+'.sav'
       return,1
    endif else return,'display'
      
end
