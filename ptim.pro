;+
; Name: ptim
; Purpose: Print times in ascii format for up to 10 input arguments
;
; Calling Sequence:  ptim, t1, t2, ...
;
; Input arguments:
;	t1, t2 - up to 10 input times (each can be scalar or array) in any anytim format
;   any keyword accepted by anytim
;
; Method:  anytim is called for each argument passed in.  If no anytim keywords are passed in
;   then /vms is used.
;
; Examples:  ptim,t1         ptim,t1,t2,/ecs,/date
;
; Written: Kim Tolbert 06-Aug-2002
; Modifications:
;-

function ptim, a,b,c,d,e,f,g,h,i,j, help=help, _extra=_extra

on_error, 2

n = n_params()

if keyword_set(help) or n eq 0 then begin
	doc_library,'ptim'
	return,foo
endif

; see if any of the ascii format options was already passed in by seeing if anytim
; returns a string.  If not, default ascii format will be vms - put it in _extra structure.
test = anytim(a, _extra=_extra)
if datatype(test) ne 'STR' then begin
	if exist(_extra) then _extra = add_tag(_extra, 1, 'vms') else _extra = {vms: 1}
endif

if n ge 1 then time= anytim(a, _extra=_extra)
if n ge 2 then time= anytim(b, _extra=_extra)
if n ge 3 then time= anytim(c, _extra=_extra)
if n ge 4 then time= anytim(d, _extra=_extra)
if n ge 5 then time= anytim(e, _extra=_extra)
if n ge 6 then time= anytim(f, _extra=_extra)
if n ge 7 then time= anytim(g, _extra=_extra)
if n ge 8 then time= anytim(h, _extra=_extra)
if n ge 9 then time= anytim(i, _extra=_extra)
if n ge 10 then time= anytim(j, _extra=_extra)
return, time
end
