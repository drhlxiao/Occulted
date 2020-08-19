
;load halpha wing files -.8 A
fm8=file_search('minus08/*.fits')
  
;load halpha wing files -.5 A
fm5=file_search('minus05/*.fits')
;fm5=fm5[1:5]

;load halpha center files (short)
f00=file_search('center/*.fits')

;load halpha wing files +.5 A
fp5=file_search('plus05/*.fits')

;load halpha wing files +.8 A
fp8=file_search('plus08/*.fits')
;fp8=fp8[1:5]

;load AIA 304 files
faia=file_search('../AIA/304/mtk*.fts')

!p.multi=[0,5,1]
set_plot,'ps'
;aia_lct,r,g,b,wavelnth=304,/load
;rgb_304=[[r],[g],[b]]
pc=.75
for i = 0,4 do begin
   ;mm8=process_file(fm8[i])
   ;mm5=process_file(fm5[i])
   ;m00=process_file(f00[i])
   ;mp5=process_file(fp5[i])
   ;mp8=process_file(fp8[i])  
   ;maia=process_file(faia[i])
   fits2map,fm8[i],map
   mm8=rot_map(map,-1*map.roll_angle)
   bg=where(mm8.data le pc*mean(mm8.data))
   mm8.data[bg]=0
   map_rflag,mm8
   mm8g=map2gmap(mm8,1.,/inverse)
   
   fits2map,fm5[i],map
   mm5=rot_map(map,-1*map.roll_angle)
   bg=where(mm5.data le pc*mean(mm5.data))
   mm5.data[bg]=0
   map_rflag,mm5
   mm5g=map2gmap(mm5,1.,/inverse)

   fits2map,f00[i],map
   m00=rot_map(map,-1*map.roll_angle)
   bg=where(m00.data le pc*mean(m00.data))
   m00.data[bg]=0
   map_rflag,m00
   m00g=map2gmap(m00,1.,/inverse)
   
   fits2map,fp5[i],map
   mp5=rot_map(map,-1*map.roll_angle)
   bg=where(mp5.data le pc*mean(mp5.data))
   mp5.data[bg]=0
   map_rflag,mp5
   mp5g=map2gmap(mp5,1.,/inverse)

   fits2map,fp8[i],map
   mp8=rot_map(map,-1*map.roll_angle)
   bg=where(mp8.data le pc*mean(mp8.data))
   mp8.data[bg]=0
   map_rflag,mp8
   mp8g=map2gmap(mp8,1.,/inverse)
   
   fits2map,faia[i],maia
   gmap=maia
   map_rflag,gmap
                                ;map2gmap(maia,1.,/inverse)

   fname='mitaka'+strtrim(string(i),1)+'.eps'
   device,filename=fname,/encapsulate,/color,xsize=12.5,ysize=4.5,/inches,/helvetica,bits=8

   loadct,13
   plot_map,mm8g,/limb,title=mm8.time + ' Halpha -.8 A',center=[-900,200],fov=[7] ;,rgb_table=rgb_304
   loadct,0
   plot_map,gmap,/overlay, levels=[10,25],/per,smooth_width=4

   loadct,13
   plot_map,mm5g,/limb,title=mm5.time + ' Halpha -.5 A',center=[-900,200],fov=[8]
   loadct,0
   plot_map,gmap,/overlay, levels=[10,25],/per,smooth_width=4;, levels=[.01,.1,1,10],/per

   loadct,13
   plot_map,m00g,/limb,title=m00.time + ' Halpha line center',center=[-900,200],fov=[8]
   loadct,0
   plot_map,gmap,/overlay, levels=[10,25],/per,smooth_width=4;, levels=[.0001,.001,.01,.1]*max(m00.data)

   loadct,13
   plot_map,mp5g,/limb,title=mp5.time + ' Halpha +.5 A',center=[-900,200],fov=[8]
   loadct,0
   plot_map,gmap,/overlay, levels=[10,25],/per,smooth_width=4;, levels=[.01,.1,1,10],/per

   loadct,13
   plot_map,mp8g,/limb,title=mp8.time + ' Halpha +.8 A',center=[-900,200],fov=[8]
   loadct,0
   plot_map,gmap,/overlay, levels=[10,25],/per,smooth_width=4; levels=[.01,.1,1,10],/per

   device,/close
endfor

end
