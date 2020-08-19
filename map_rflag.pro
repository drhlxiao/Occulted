PRO map_rflag,tst,r=r,value=value,shift_center=shift_center,outside=outside,zero=zero

if not keyword_set(r) then r=get_rb0p(tst.time,/radius)

xp=get_map_xp(tst)
yp=get_map_yp(tst)
if keyword_set(shift_center) then begin
  print,'shifted'
  xp=xp+shift_center(0)
  yp=yp+shift_center(1)
endif

rp=sqrt(xp^2+yp^2)
rlist=where(rp le r)
if keyword_set(outside) then rlist=where(rp gt r)

if rlist(0) ne -1 then begin
  if not keyword_set(zero) then tst.data(rlist)=sqrt(-1) else tst.data(rlist)=0
endif

if keyword_set(value) then begin
   print,value
   tst.data(rlist)=value
endif

END
