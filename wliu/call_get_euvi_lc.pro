;  @call_get_euvi_lc.pro

drot=1
ref_time= '1-Sep-2014 11:00:56.204'			; map1.time
keep_limb=0
fov= [ [285.8166, 700.8098], [185.5374, 494.6057] ]		; from make_map_get_fov.pro, 2015/06/27

.r ~/idl_lib/stereo_wei/get_euvi_lc

file_list= '/disk2/scr0/weiliu/sdo/events/05a_flare_Fermi/14-09-01_Fermi-occult-NE/multiw/stereo/files/b/STb-195_20140901-00_0902-2400_files.txt'
readcol, file_list, files0, format='a'

	;--- try for background ---
	;istart= 72-1
	;iend= 98-1

	;istart= 72-1
	;iend= 527-1		;431-1

istart= 96-1	;20140901_080030_n4euB.fts		; 72-1	;20140901_060030_n4euB.fts
iend= 192-1		;20140901_160030_n4euB.fts

files=files0[istart:iend]

	;bkg_interval= '2014/09/01 ' + ['06:55:00', '07:22:00']		;['09:00:00', '10:00:00']

	;from /14-09-01_Fermi-occult-NE/multiw/stereo/run0_toal-flux/195/final_redo-bkg
new_bkg_interval= '2014/09/01 ' + ['09:40:00', '10:03:00']		;['06:55:00', '07:22:00']		;['09:00:00', '10:00:00']
bkg_interval= new_bkg_interval

get_euvi_lc, files, fov=fov, bkg_interval=bkg_interval, win_size=win_size   $
    , time_str=time_str, time_dbl=time_dbl, base_time=base_time, lc_euvi=lc_euvi $
    , drot=drot, ref_time=ref_time, keep_limb=keep_limb
