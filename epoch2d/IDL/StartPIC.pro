FUNCTION getdata, snapshot, wkdir=wkdir, explorer=explorer, retro=retro, $
    _EXTRA=extra

  COMMON background, wkdir_global, retro_global
  ON_ERROR, 2

  IF NOT KEYWORD_SET(wkdir) THEN wkdir = wkdir_global
  IF NOT KEYWORD_SET(retro) THEN retro = retro_global

  IF N_PARAMS() EQ 0 THEN BEGIN
    PRINT, "Usage: result = getdata(snapnumber[, wkdir=<dir>, " + $
        "/empty | /rho, /temp, /vx ...])"
    RETURN, "Usage: result = getdata(snapnumber[, wkdir=<dir>, " + $
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

PRO quick_view, snapshot, wkdir=wkdir, get_viewer=get_viewer, $
    set_viewer=set_viewer

  COMMON background, wkdir_global
  ON_ERROR, 2

  IF NOT KEYWORD_SET(wkdir) THEN wkdir = wkdir_global

  IF N_PARAMS() EQ 0 THEN BEGIN
    PRINT, "Usage: result = showdata(snapnumber[, wkdir=<dir>])"
    RETURN
  ENDIF

  file = wkdir + STRING(snapshot, FORMAT='("/",I04.04,".cfd")')
  info = FILE_INFO(file)

  IF info.exists EQ 1 THEN BEGIN
    get_viewer = create_sdf_visualizer(file, set_viewer=set_viewer)
  ENDIF ELSE BEGIN
    file = wkdir + STRING(snapshot, FORMAT='("/",I04.04,".sdf")')
    get_viewer = create_sdf_visualizer(file, set_viewer=set_viewer)
  ENDELSE
END
