;pro get_flare_angles

;gets the angles between Mercury and Earth for the flares in Dennis's messenger list

; read in the data
;filename = '/Users/wheatley/Documents/Solar/MiSolFA/statistics/messenger_list.csv'
;flare_list = read_csv(filename)

;angle=dblarr(656)

;for i=0,655 do begin
;   splitstring = strsplit(flare_list.field2[i],'/',/extract)
                                ;mm/dd/yy -> yy/mm/dd
                                ;89/12/15, 22:00:15.234
;   tstring = strsplit(flare_list.field3[i],' ',/extract)
;   strings = splitstring[2]+'/'+splitstring[0]+'/'+splitstring[1]+ ' '+tstring[0]
;   angle[i]=messenger_flare_angle([-200.,200.], strings)
;   print, tstring,strings,angle[i]
;endfor

;anglest={Angle:angle}
;newst = {Index:flare_list.field1, Date:flare_list.field2, Star_Time:flare_list.field3, End_Time:flare_list.field4, Duration:flare_list.field5, Angle:angle}
;write_csv, 'messenger_angles.csv',newst

;pro calc_goes_class

;calculate the GOES class for each flare, using the largest of EM1 and
;EM2

;read in the data
filename = '/Users/wheatley/Documents/Solar/occulted_flares/messenger_angles.csv'
flare_list = read_csv(filename)

goesnum=dblarr(656)
goesclass=strarr(656)

for i=0,655 do begin
   ;splitstring = strsplit(flare_list.field2[i],'/',/extract)
                                ;mm/dd/yy -> yy/mm/dd
                                ;89/12/15, 22:00:15.234
   ;tstring = strsplit(flare_list.field3[i],' ',/extract)
   ;strings = splitstring[2]+'/'+splitstring[0]+'/'+splitstring[1]+ ' '+tstring[0]
   tstring=flare_list.field6[i]
   EM1=float(flare_list.field5[i])
   EM2=float(flare_list.field7[i])
   TMKelvin = float(tstring)*11.6045 ;convert from keV to MK
   if EM1 gt EM2 then EM=EM1 else EM=EM2
   
   goes_fluxes,TMKelvin,EM,flong,fshort,sat=15
   goesnum[i]=flong   
   goesclass[i]=goes_value2class(flong)
endfor

;plot the histogram
cghistoplot, alog10(goesnum), nbins=7 ,/fill, yran=[0,250], xran=[-10,0],ytitle='Number of events',xtitle='log flux, W/m^2'

end
