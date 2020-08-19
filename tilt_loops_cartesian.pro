;program to tilt loops (coordinates) on solar surface

;************************************************************************
;input in cartesian coordinates [xx,yy,zz], X towards earth, Y towards
;solar west, Z towards solar north (HEE system)
;output is tilted loop in spherical coordinates [phi,theta,height] = [phout,thout,hout]
;outcartesian = cartesian coordinates can be returned in outcartesian variable
;verbose      = prints debug info if /verbose is set
;tilt must be given in degrees, it is converted to radians automatically
;uses pass by reference for output
;************************************************************************

;v1, Lucia Kleint, 21.10.2016

PRO tilt_loops_cartesian,xx,yy,zz,phout,thout,hout,tilt,verbose=verbose,$
                         outcartesian=outcartesian
    if keyword_set(verbose) then verbose=1 else verbose=0

;--- put one endpoint of loop to [x,y,z]=[0,0,0] ---
    xx2 = double(xx - xx[0])
    yy2 = double(yy - yy[0])
    zz2 = double(zz - zz[0])


;X is towards Earth, Z towards north, making Y towards west
;point[0] is at [0,0,0]
;idea: find angle to rotate to so that loop lies along Y axis
;need rotation about Z and about X for that

;--- define unit vectors
    yvec  = [0,1.,0]
    zvec  = [0,0,1.]
    endpt = [xx2[-1],yy2[-1],zz2[-1]]


;--- project endpt onto xy plane to determine angle wrt Y
    proj = endpt - (transpose(endpt)#zvec)[0]*zvec
    alpha = (acos( (transpose(proj)#yvec) / (norm(proj)) ))[0] ;radians

;--- rotate about Z by alpha so that loop is completely in Y direction
    xhee2_rot = xx2*cos(alpha) + yy2*sin(alpha)
    if abs(xhee2_rot[-1]) gt .001 then begin
       alpha = -alpha
       xhee2_rot = xx2*cos(alpha) + yy2*sin(alpha)
    endif
    yhee2_rot = -xx2*sin(alpha) + yy2*cos(alpha)
    zhee2_rot = zz2
    endpt2 = [xhee2_rot[-1],yhee2_rot[-1],zhee2_rot[-1]]
;point is now [0,y,z] -> rot about X axis to get [0,y,0]
;something is wrong with these rotations, alpha sometimes should be
;-alpha, probably some acos ambiguity -> just flipped it if point is [>0,y,z]

;--- get angle wrt Y (no proj. necessary since x=0)
    beta = (acos( (transpose(endpt2)#yvec) / (norm(endpt2)) ))[0] 

;--- rotate about X to make loop horizontal (seen from Earth)
    zhee3_rot = -yhee2_rot*sin(beta) + zhee2_rot*cos(beta)
    if abs(zhee3_rot[-1]) gt .001 then begin
       beta = -beta
       zhee3_rot = -yhee2_rot*sin(beta) + zhee2_rot*cos(beta)
    endif
    xhee3_rot = xhee2_rot
    yhee3_rot =  yhee2_rot*cos(beta) + zhee2_rot*sin(beta)
    endpt3 = [xhee3_rot[-1],yhee3_rot[-1],zhee3_rot[-1]]
;point is now [0,Y,0]
    if verbose then print,'endpt [0,n,0]: ',endpt3
    if verbose then print,'angles:',alpha*!radeg,beta*!radeg

;--------------------- INCLINE LOOP --------------------------------
;incline loop (=rotation about Y axis)
;positive angles seem to be towards south-east on the Sun
    incl = double(tilt*!dtor)       

    xhee_inc = xhee3_rot*cos(incl) + zhee3_rot*sin(incl)
    yhee_inc = yhee3_rot
    zhee_inc = -xhee3_rot*sin(incl) + zhee3_rot*cos(incl)
    endpt4   = [xhee_inc[-1],yhee_inc[-1],zhee_inc[-1]]


;---------------- transform back to original coords ----------------
;apply all rotations from above in reverse order

;rot about X by -beta
    xhee4_rot =  xhee_inc
    yhee4_rot =  yhee_inc*cos(-beta) + zhee_inc*sin(-beta)
    zhee4_rot = -yhee_inc*sin(-beta) + zhee_inc*cos(-beta)
    if verbose then print,'after rot about X'
    if verbose then print,xhee4_rot[-1],yhee4_rot[-1],zhee4_rot[-1]   
    if verbose then print,endpt2


;rot about Z by -alpha
    xhee5_rot =  xhee4_rot*cos(-alpha) + yhee4_rot*sin(-alpha)
    yhee5_rot = -xhee4_rot*sin(-alpha) + yhee4_rot*cos(-alpha)
    zhee5_rot =  zhee4_rot
    if verbose then print,'after rot about z'
    if verbose then print,xhee5_rot[-1],yhee5_rot[-1],zhee5_rot[-1]   
    if verbose then print,endpt

;shift back to original coords
    xhee6 = xhee5_rot + xx[0]
    yhee6 = yhee5_rot + yy[0]
    zhee6 = zhee5_rot + zz[0]

    if verbose then print,'input endpoint of loop'
    if verbose then print,[xx[-1],yy[-1],zz[-1]]
    if verbose then print,'output endpoint of loop (should be identical)'
    if verbose then print,[xhee6[-1],yhee6[-1],zhee6[-1]]


;--- transform loop from HEE to spherical coordinates:
    heights_inc = sqrt(xhee6^2.+yhee6^2.+zhee6^2.)
    thetas_inc  = atan(zhee6/sqrt(xhee6^2.+yhee6^2.))*!radeg ;degrees
    phis_inc    = atan(yhee6/xhee6)*!radeg                   ;degrees


;--- return values
    phout = phis_inc
    thout = thetas_inc
    hout  = heights_inc

;if desired, return cartesian coordinates (not very useful for solar plotting)
    if keyword_set(outcartesian) then outcartesian = [[xhee6],[yhee6],[zhee6]]

END
