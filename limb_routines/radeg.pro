;+
; PURPOSE:
;	To replace IDL's !RADEG, which has only a float precision, with a double precision.
;
;
;
;
;-
FUNCTION radeg

	RETURN, 180d / !DPI

END

