FUNCTION map_rflag_cube,map,r=r,shift_center=shift_center,outside=outside,value=value,shift_random=shift_random,zero=zero


dim=n_elements(map)

fmap=map

this_dim=n_elements(fmap)
if keyword_set(shift_random) then this_ran=randomn(systim(1),this_dim*2)*shift_random else this_ran=fltarr(this_dim*2)
this_ran_x=this_ran(0:this_dim-1)
this_ran_y=this_ran(this_dim:*)

for i=0,dim-1 do begin
   print,i,dim-1
   this_map=map(i)
   map_rflag,this_map,r=r,shift_center=shift_center+[this_ran_x(i),this_ran_y(i)],outside=outside,value=value,zero=zero
   fmap(i)=this_map
endfor

return,fmap

END

