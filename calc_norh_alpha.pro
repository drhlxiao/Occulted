pro calc_norh_alpha,indexal,alpha,mvda,f17_c,f34_c
files17=findfile('17g10s_part/ipa130501*')
norh_rd_img,files17[0:173],index17,data17
norh_index2map,index17,data17,map17 
files34=findfile('130501_34g10s/ipz130501*')
norh_rd_img,files34,index34,data34
norh_index2map,index34,data34,map34  

tf17=map_rflag_cube(map17,shift_center=[0,0],r=967.)
nanInd = where(~(FINITE(tf17.data)))
tf17.data[nanInd]=0.0
 f17=norh_tb2flux(tf17.data,index17,/intensity)
;f17=norh_tb2flux(data17,index17,/intensity);,abox=[-1410,-950,-50,400],/intensity)

;tf34=map_rflag_cube(map34,shift_center=[0,0],r=967.)
;nanInd = where(~(FINITE(tf34.data)))
;tf34.data[nanInd]=0.0
f34=norh_tb2flux(data34,index34,/intensity)
; f34=norh_tb2flux(data34,index34,/intensity);,abox=[-1410,-950,-50,400],/intensity)

norh_convol, index34,index17,f17,index17_c,f17_c
norh_convol, index17,index34,f34,index34_c,f34_c
;norh_alpha,index17_c,f17_c,index34_c,f34_c,indexal,alpha,mvda
save, index17_c,f17_c,index34_c,f34_c,f17,f34, filename='norh_convol.sav'
;print, alpha
end
