pro which_stereo,fpos,mpeak,mlonglat,Alonglat,Blonglat
  llen=size(mpeak,/dim)
  mlonglat=findgen(3,llen[0])
  Alonglat=findgen(3,llen[0])
  Blonglat=findgen(3,llen[0])

  for i=0,llen[0]-1 do begin
     ;need to convert fpos[i] from string to fltarr
     fcomma=strsplit(fpos[i],',')
     fc1=strmid(fpos[i],1,fcomma[1]-2)
     fc2=strmid(fpos[i],fcomma[1],strlen(fpos[i])-(fcomma[1]+1))
     fcoord=[float(fc1),float(fc2)]
     ma=messenger_flare_angle(fcoord,mpeak[i],mess_hg_str) ;get messenger angle, longlat. This is in degrees.
     ;print,fpos[i]
                                ;print,mpeak[i], ma
     z = value_locate( mess_hg_str.date, anytim(mpeak[i]))
                                ;print, mpeak[i],fcoord,z,mess_hg_str[z].xyz
     ml=mess_hg_str[z].xyz
     mlonglat[*,i]=[ma,ml[0:1]]
     aa=get_sunspice_sep_angle(mpeak[i],'A','Earth',system='HEE',/degrees)
     all=get_sunspice_lonlat(mpeak[i],'A',system='HEE',/degrees)
     ba=get_sunspice_sep_angle(mpeak[i],'B','Earth',system='HEE',/degrees)
     bll=get_sunspice_lonlat(mpeak[i],'B',system='HEE',/degrees)
     Alonglat[*,i]=[aa,all[1:*]]
     Blonglat[*,i]=[ba,bll[1:*]]
  endfor

end

