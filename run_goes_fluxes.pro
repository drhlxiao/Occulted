pro run_goes_fluxes,tlist,emlist,flong,fshort,abund
  llen=size(tlist,/dim)
  flong=findgen(llen[0])
  fshort=findgen(llen[0])
  for i=0,llen[0]-1 do begin
        tMK=tlist[i]*11.59
        goes_fluxes,tMK,emlist[i],out1,out2,sat=15,abund=abund
        flong[i]=out1
        fshort[i]=out2
    endfor
end
