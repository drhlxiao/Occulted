function process_file,file
  fits2map,file,map
  rmap=rot_map(map,-1*map.roll_angle)
                                ;gmap=map2gmap(map,1.,/inverse)
  filename=strmid(file,39,35,/reverse)+'_rot.fits'
  print,filename
  map2fits,map,filename
end
