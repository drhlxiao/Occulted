;pro trace_all_loops
;use trace loops procedure on all flares, store the locations of the
;footpoints, tilt, loop height, etc. in structure or csv for transport
;to Python
    
function format_datetime,dt,date=date ;cuz IDL sucks at datetimes
    if ~ keyword_set(date) then date=strmid(dt,0,4)+strmid(dt,5,2)+strmid(dt,8,2)+'_*eu*fts*' else date=strmid(dt,0,4)+strmid(dt,5,2)+strmid(dt,8,2)
    return, date
    end

function format_footpoints, datastr ;probably need to output a float array not a string... 
    aa=strsplit(datastr,'=',/extract)  
    x1=float(strmid(aa[1],0,9))              ;remove the y at the end
    y1=float(strmid(aa[2],0,9))   
    x2=float(strmid(aa[3],0,9))    
    y2=float(aa[4])
    pos=[[x1,y1],[x2,y2]]
    return, pos
    end

function format_coord,coord
    aa=strsplit(coord,';',/extract)
    xlen=strlen(aa[0])
    ylen=strlen(aa[1])
    numcoord=[float(strmid(aa[0],1,xlen-1)), float(strmid(aa[1],0,ylen-1))]
    return, numcoord
    end

pro draw_all_loops, csvname, ivec=ivec, quiet=quiet,fov=fov ;might have to do this manually...sadness
    loadct,9
    ;if ~keyword_set(ivec) then ivec=findgen(llen[0])
    if ~keyword_set(fov) then fov=[10]
    stuff=read_csv(csvname,header=header)
    llen = size(stuff.field01,/dim)
    for i=0,llen[0]-1 do begin ;llen[0]-1 do begin
        ;flare=flare_list(ivec[i]) doesn't work
        loadct,9
        if stuff.field02[i] ne '' then begin 
            center=format_coord(stuff.field10[i]) ;these need to be ints or floats
            ell=format_coord(stuff.field04[i]) ;these need to be ints or floats
            pos=format_footpoints(stuff.field02[i]) ;need to do some coordinate conversion on this... go figure
            ffname = format_datetime(anytim(anytim(stuff.field01[i]),/ecs))
            outfilename='/Users/wheatley/Documents/Solar/data/stereo_pfloops/'+format_datetime(anytim(anytim(stuff.field01[i]),/ecs),/date)+'loops.png'
            fitsfile=file_search('/Users/wheatley/Documents/Solar/data/stereo_pfloops/'+ffname)
            fits2map,fitsfile[0],map
            plot_map, map,/log,center=center,fov=fov 
            run_trace_loops,pos=pos,date=map.time,ellipse1=ell
            write_png, outfilename,tvrd(/true) ;save the image somehow
        endif
    endfor
    end

pro time_series_loop, basename, pos=pos,center=center, ell=ell,fov=fov ;might have to do this manually...sadness
    loadct,9
    if ~keyword_set(pos) then pos=[[-694.93088,137.16588],[-565.52995,100.66818]]
    if ~keyword_set(center) then center=[-630,125]
    if ~keyword_set(ell) then ell=[1,-45]
    if ~keyword_set(fov) then fov=[8]
    
    ffname = basename+'*eua.fts'
    fitsfiles=file_search('/Users/wheatley/Documents/Solar/occulted_flares/data/stereo-aia/'+ffname)
    llen=size(fitsfiles,/dim)

    for i=0,llen[0]-1 do begin  
        loadct,9
        outstr=strmid(fitsfiles[i],64,21)
        outfilename='/Users/wheatley/Documents/Solar/occulted_flares/data/stereo_pfloops/'+outstr+'_loops.png'
        fits2map,fitsfiles[i],map
        plot_map, map,/log,center=center,fov=fov 
        run_trace_loops,pos=pos,date=map.time,ellipse1=ell
        write_png, outfilename,tvrd(/true) ;save the image somehow
    endfor
    end
