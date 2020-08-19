FUNCTION make_lc,map,ind=ind,this_display=this_display,noneg=noneg,stop=stop,norm=norm,sunck=sunck,per_second=per_second,above=above,corner=corner,ylog=ylog,no_para=no_para,center_of_mass=center_of_mass,amap=amap,ref_index=ref_index,this_degree=this_degree,more_maps=more_maps,cmap=cmap,separation=separation,pixlist=pixlist,use_same_pixels=use_same_pixels

if not keyword_set(ind) then begin
   image=total(map.data,3)*0.
   this_dim=n_elements(map)
   ;for i=0,this_dim-1 do image=image+map[i].data/max(map[i].data)
   for i=0,this_dim-1 do begin
      this_image=map[i].data
      glist=where(finite(this_image) eq 0)
      if glist(0) ne -1 then this_image(glist)=0
      this_max=max(this_image)
      glist=where(this_image le 0.15*this_max)
      if glist(0) ne -1 then this_image(glist)=0
      image=image+this_image/this_max
   endfor
endif else begin
   image=map[ind].data
endelse

xxx=get_map_xp(map(0))
yyy=get_map_yp(map(0))
                       
if keyword_set(norm) then begin
  tdim=n_elements(map)                       
  image=map(0).data/max(map(0).data)
  for i=1,tdim-1 do image=image+map(i).data/max(map(i).data)                        
endif
 
 if keyword_set(this_display) then image=this_display
                        
if not keyword_set(use_same_pixels) then pixlist=hsi_select_box(image,separation=separation)

selected_sources=map.data*0.

ndim=n_elements(separation)

dim=n_elements(map)
tstmap=map(0).data
xdim=n_elements(tstmap(*,0))
ydim=n_elements(tstmap(0,*))
t=dblarr(dim)
f=findgen(dim,ndim)
fpeak=findgen(dim,ndim,2)
peakx=findgen(dim,ndim,2)
peaky=findgen(dim,ndim,2)
peakxx=findgen(dim,ndim,2)
peakyy=findgen(dim,ndim,2)
com_x=findgen(dim,ndim,2)
com_y=findgen(dim,ndim,2)
fout=findgen(dim)
sout=findgen(dim)
mout=findgen(dim)
m=findgen(dim)

start=0
for i=1,ndim do begin 
  list=pixlist(start:start+separation(i-1)-1)
  start=start+separation(i-1)
  for j=0,dim-1 do begin 
   hhh=map[j].data
   if (keyword_set(noneg)) then hhh=hhh>0
   if keyword_set(corner) then hhh=hhh-average(hhh(0:50,0:50))
   if keyword_set(above) then hhh(where(hhh le above))=0
   f(j,i-1)=total(hhh(list),/nan)
   t(j)=anytim(map[j].time)
   hhh=hhh*0
   hhh(list)=1
   selected_sources(*,*,j)=selected_sources(*,*,j)+hhh
   olist=where(hhh eq 0)
   hhh=float(map[j].data)
   hhh(olist)=0
   fpeak(j,i-1,0)=max(hhh,ind)
   py=fix(ind/xdim)
   px=ind-py*xdim
   peakx(j,i-1,0)=px
   peaky(j,i-1,0)=py
   peakxx(j,i-1,0)=xxx(px,py)
   peakyy(j,i-1,0)=yyy(px,py)
   com_x(j,i-1,0)=total(xxx(list)*hhh(list))/total(hhh(list))
   com_y(j,i-1,0)=total(yyy(list)*hhh(list))/total(hhh(list))
   if not keyword_set(no_para) then begin
   ;peakpix = locate_val(hhh,MAX(hhh))
   ;peakregion = hhh[peakpix[0]-1:peakpix[0]+1,peakpix[1]-1:peakpix[1]+1]
   ;mapeak = PARAPEAK(peakregion)               ; does parabolic fit to 3x3 array
   ;peakx(j,i-1,1)=px+mapeak(0)
   ;peaky(j,i-1,1)=py+mapeak(1)
   ;peakxx(j,i-1,1)=xxx(px,py)+mapeak(0)*map(0).dx
   ;peakyy(j,i-1,1)=yyy(px,py)+mapeak(1)*map(0).dy
   ;fpeak(j,i-1,1)=mapeak(2)
   endif
  endfor
endfor

;outside selected sources
for j=0,dim-1 do begin 
  hhh=selected_sources(*,*,j)
  olist=where(hhh eq 0)
  print,n_elements(olist)
  hhh=map[j].data
  m(j)=max(hhh)
  fout(j)=total(hhh(olist),/nan)
  sout(j)=sigma(hhh(olist))
  mout(j)=max(hhh(olist))
endfor
;stop
if keyword_set(per_second) then f=f/map.dur

p={time:t,flux:f,flux_max:m, flux_out:fout, flux_max_out:mout, sigma_out:sout, fpeak:fpeak, peakx: peakx, peaky: peaky, peakxx: peakxx, peakyy: peakyy, com_x: com_x, com_y: com_y}
loadct,5
utplot,p.time,p.flux(*,0),psym=-1,yrange=[min(f),max(f)],ylog=ylog
for i=1,ndim-1 do begin
  outplot,p.time,p.flux(*,i),psym=-1,color=255/ndim*(ndim-i)
endfor
tint=anytim(p.time,/c)
ff='/disks/sunck/disk1/krucker/lc_data/lc_'
ff='/mydisks/sunck/disk1/krucker/lc_data/lc_'
ff=''
if keyword_set(sunck) then ff='/mydisks/disk1/krucker/lc_data/lc_'
ff=ff+strmid(tint(0),0,4)+strmid(tint(0),5,2)+strmid(tint(0),8,2)
ff=ff+'_'+strmid(tint(0),11,2)+strmid(tint(0),14,2)+strmid(tint(0),17,2)
if keyword_set(addname) then ff=ff+'_'+addname
oldf=findfile(ff+'*')
if oldf(0) eq '' then ff=ff+'_0' else ff=ff+'_'+strtrim(n_elements(oldf),2)

save,p,image,pixlist,separation,ndim,filename=ff+'.dat'

print,'time profiles are saved in: '
print,ff+'.dat'

if keyword_set(center_of_mass) then begin
  lc0=p
  clearplot
  if not keyword_set(this_degree) then this_degree=5
  
   glist=where(finite(lc0.com_x(*,0,0)+lc0.com_y(*,0,0)) eq 1)
   
  loadct2,5
  !P.multi=[0,3,2]
  utplot,lc0.time,lc0.flux(glist,0,0)

  utplot,lc0.time(glist),lc0.com_x(glist,0,0),ystyle=1
  cx0 = poly_fit(lc0.time(glist)-lc0.time(glist(0)),lc0.com_x(glist,0,0),this_degree,yfit=yfitx0)
  outplot,lc0.time(glist),yfitx0,color=4,thick=th2

  utplot,lc0.time(glist),lc0.com_y(glist,0,0),ystyle=1
  cy0 = poly_fit(lc0.time(glist)-lc0.time(glist(0)),lc0.com_y(glist,0,0),this_degree,yfit=yfity0)
  outplot,lc0.time(glist),yfity0,color=4,thick=th2

  if not keyword_set(ref_index) then ref_index=0
  this_yr=[-15,15]
  plot,yfitx0-yfitx0(ref_index),yfity0-yfity0(ref_index),ystyle=1,color=4,yrange=this_yr,xrange=this_yr,thick=th2
  
  utplot,lc0.time(glist),yfitx0-yfitx0(ref_index),ystyle=1,color=4,yrange=this_yr,thick=th2

  utplot,lc0.time(glist),yfity0-yfity0(ref_index),ystyle=1,color=4,yrange=this_yr,thick=th2

  this_xy=[[yfitx0-yfitx0(ref_index)],[yfity0-yfity0(ref_index)]]
 
  amap=map(glist)
  amap.xc=amap.xc-this_xy(*,0)
  amap.yc=amap.yc-this_xy(*,1)
  
  if keyword_set(more_maps) then begin
   cmap=more_maps(glist,*)
   ndim=n_elements(cmap(0,*))
   for i=0,ndim-1 do cmap(*,i).xc=cmap(*,i).xc-this_xy(*,0)
   for i=0,ndim-1 do cmap(*,i).yc=cmap(*,i).yc-this_xy(*,1) 
   for i=0,ndim-1 do print,this_xy(*,0) 
   for i=0,ndim-1 do print,this_xy(*,1) 
  endif
  
  
endif

return,p



END
