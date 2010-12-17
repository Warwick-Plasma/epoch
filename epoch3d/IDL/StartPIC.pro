FUNCTION getdata, snapshot, wkdir=wkdir, _EXTRA=extra
  COMMON background, wkdir_global
  ON_ERROR, 2

  IF NOT KEYWORD_SET(wkdir) THEN wkdir = wkdir_global

  IF N_PARAMS() EQ 0 THEN BEGIN
    PRINT, "Usage: result = getdata(snapnumber[, wkdir=<dir>, " + $
        "/empty | /rho, /temp, /vx ...])"
    RETURN, "Usage: result = getdata(snapnumber[, wkdir=<dir>, " + $
        "/empty | /rho, /temp, /vx ...])"
  ENDIF

  file = wkdir + STRING(snapshot, FORMAT='("/",I04.04,".cfd")')
  info = FILE_INFO(file)

  IF info.exists EQ 1 THEN BEGIN
    RETURN, LoadCFDFile(file, _EXTRA=extra)
  ENDIF ELSE BEGIN
    file = wkdir + STRING(snapshot, FORMAT='("/",I04.04,".sdf")')
    RETURN, LoadSDFFile(file, _EXTRA=extra)
  ENDELSE
END
