
FUNCTION rad2asec, ang

	IF exist(ang) THEN RETURN, ang / !DPI * 180d * 3600 ELSE RETURN, 3600d * 180d /!DPI
	
END

