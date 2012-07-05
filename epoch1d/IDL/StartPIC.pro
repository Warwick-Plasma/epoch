FUNCTION getdata, snapshot, wkdir, retro=retro, _EXTRA=extra

  COMPILE_OPT idl2
  COMMON background, wkdir_global, retro_global
  ON_ERROR, 2

  IF N_ELEMENTS(wkdir) EQ 0 THEN wkdir = wkdir_global
  IF NOT KEYWORD_SET(retro) THEN retro = retro_global

  IF N_PARAMS() EQ 0 THEN BEGIN
    PRINT, "Usage: result = getdata(snapnumber[,<wkdir>, " + $
        "/empty | /rho, /temp, /vx ...])"
    RETURN, "Usage: result = getdata(snapnumber[,<wkdir>, " + $
        "/empty | /rho, /temp, /vx ...])"
  ENDIF

  file = wkdir + STRING(snapshot, FORMAT='("/",I04.04,".cfd")')
  info = FILE_INFO(file)

  IF info.exists EQ 1 THEN BEGIN
    IF (N_ELEMENTS(explorer) NE 0) THEN BEGIN
      RETURN, sdf_explorer(file)
    ENDIF
    RETURN, LoadCFDFile(file, retro=retro, _EXTRA=extra)
  ENDIF ELSE BEGIN
    file = wkdir + STRING(snapshot, FORMAT='("/",I04.04,".sdf")')
    IF (N_ELEMENTS(explorer) NE 0) THEN BEGIN
      RETURN, sdf_explorer(file)
    ENDIF ELSE BEGIN
      RETURN, LoadSDFFile(file, retro=retro, _EXTRA=extra)
    ENDELSE
  ENDELSE
END

; --------------------------------------------------------------------------

FUNCTION explore_data, wkdir, snapshot=snapshot

  COMPILE_OPT idl2
  COMMON background, wkdir_global, retro_global
  ON_ERROR, 2

  IF N_ELEMENTS(wkdir) EQ 0 THEN wkdir = wkdir_global

  RETURN, sdf_explorer(wkdir, snapshot=snapshot)
END

; --------------------------------------------------------------------------

PRO quick_view, wkdir, snapshot=snapshot

  COMPILE_OPT idl2
  COMMON background, wkdir_global
  ON_ERROR, 2

  IF NOT KEYWORD_SET(wkdir) THEN wkdir = wkdir_global

  a = create_sdf_visualizer(wkdir, snapshot=snapshot)
END
