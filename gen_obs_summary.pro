;pro gen_obs_summary

;Generate observing summary figure for list of RHESSI (and messenger) flares

;get time intervals
function get_time_intervals(list)
    for i=0,size(list.date,/dim) do begin
       datetime[i]=anytim(list.date[i]+list.RHESSI_time[i],'format')
       
    endfor
    ;convert strings to anytim format ... j/k IDL is the worst
    return, times
end

;download data if need be
function get_data(args)
  search_network,/enable
end

function obs_summary(times)
    obj = hsi_obs_summary(/no_gui)
    obj -> set, obs_time_interval=[times.start[0], times.end[0]]
    uncorr_data = obj -> getdata(class_name='hsi_full_rate') ; I think that's right ...
    ;how to plot without GUI?

end

; reference:
; https://hesperia.gsfc.nasa.gov/ssw/packages/spex/doc/ospex_explanation.htm
; https://hesperia.gsfc.nasa.gov/ssw/packages/spex/doc/ospex_params_standard.htm
function plot_messenger(specfiles, drmfiles,)
  
    o=ospex(/no_gui)
    o->set, spex_specfile='d:\ospex\hsi_spectrum_20020220_105002.fits'
    o->set, spex_drmfile='d:\ospex\hsi_srm_20020220_105002.fits'
    o -> set, spex_bk_time_interval=[['20-Feb-2002 10:56:23.040', '20-Feb-2002 10:57:02.009'], $
       ['20-Feb-2002 11:22:13.179', '20-Feb-2002 11:22:47.820'] ]
    o->set, spex_eband=get_edge_products([3,22,43,100,240],/edges_2)
    o->set, spex_fit_time_interval=[ ['20-Feb-2002 11:06:03.259', '20-Feb-2002 11:06:11.919'], $
       ['20-Feb-2002 11:06:11.919', '20-Feb-2002 11:06:24.909'], $
       ['20-Feb-2002 11:06:24.909', '20-Feb-2002 11:06:33.570'] ]
    o->set, spex_erange=[19,190]
    o->set, fit_function='vth+bpow'
    o->set, fit_comp_param=[1.0e-005,1.,1.,  .5, 3., 45., 4.5]
    o->set, fit_comp_free = [0,1,1, 1,1,1,1]
    o->set, spex_fit_manual=0
    o->set, spex_autoplot_enable=1
    o->dofit, /all
end

; main program

  cd, '/Users/wheatley/Documents/Solar/occulted_flares'
restore, 'list_selected.sav'
times = get_time_intervals(list_selected)
if keyword_set(network) == False then begin ;if you already have the data  
   get_data(args)
   endif
else cd, '/path/to/data'        ;go to data directory


end
