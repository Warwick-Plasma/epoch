FUNCTION getdata, snapshot_in, wkdir_in, _EXTRA=extra

  COMPILE_OPT idl2
  COMMON background, wkdir_global, retro_global
  ON_ERROR, 2

  snapshot = -1
  retro = -1
  wkdir = ''

  IF (N_ELEMENTS(snapshot_in) NE 0) THEN BEGIN
    IF (SIZE(snapshot_in, /TYPE) EQ 7) THEN BEGIN
      wkdir = snapshot_in
    ENDIF ELSE BEGIN
      snapshot = snapshot_in
    ENDELSE
  ENDIF

  IF (N_ELEMENTS(wkdir_in) NE 0) THEN BEGIN
    IF (SIZE(wkdir_in, /TYPE) EQ 7) THEN BEGIN
      wkdir = wkdir_in
    ENDIF ELSE BEGIN
      snapshot = wkdir_in
    ENDELSE
  ENDIF

  IF (KEYWORD_SET(extra)) THEN BEGIN
    extra_tags = TAG_NAMES(extra)
    gotextra = 0
    FOR i = 0, N_ELEMENTS(extra_tags)-1 DO BEGIN
      CASE extra_tags[i] OF
        'SNAPSHOT': snapshot = extra.snapshot
        'WKDIR': wkdir = extra.wkdir
        'RETRO': retro = extra.retro
        ELSE: BEGIN
          IF (gotextra EQ 0) THEN BEGIN
            new_extra = CREATE_STRUCT(extra_tags[i], extra.(i))
            gotextra = 1
          ENDIF ELSE BEGIN
            new_extra = CREATE_STRUCT(new_extra, extra_tags[i], extra.(i))
          ENDELSE
        END
      ENDCASE
    ENDFOR
  ENDIF

  IF snapshot EQ -1 THEN BEGIN
    PRINT, "Usage: result = getdata(snapnumber[,<wkdir>, " + $
        "/empty | /rho, /temp, /vx ...])"
    RETURN, "Usage: result = getdata(snapnumber[,<wkdir>, " + $
        "/empty | /rho, /temp, /vx ...])"
  ENDIF

  IF (wkdir EQ '') THEN wkdir = wkdir_global
  IF (retro EQ -1) THEN retro = retro_global

  file = wkdir + STRING(snapshot, FORMAT='("/",I04.04,".sdf")')
  RETURN, LoadSDFFile(file, _retro=retro, _EXTRA=new_extra)
END

; --------------------------------------------------------------------------

FUNCTION getstruct, snapshot_in, wkdir_in, _EXTRA=extra
  COMPILE_OPT idl2
  COMMON background, wkdir_global, retro_global

  data = getdata(snapshot_in, wkdir_in, _retro=0, _EXTRA=extra)

  RETURN, data
END

; --------------------------------------------------------------------------

FUNCTION explore_data, wkdir, snapshot=snapshot

  COMPILE_OPT idl2
  COMMON background, wkdir_global, retro_global
  ON_ERROR, 2

  IF N_ELEMENTS(wkdir) EQ 0 THEN wkdir = wkdir_global

  RETURN, sdf_explorer(wkdir, snapshot=snapshot, _struct=0)
END

; --------------------------------------------------------------------------

FUNCTION explore_struct, wkdir, snapshot=snapshot

  COMPILE_OPT idl2
  COMMON background, wkdir_global, retro_global
  ON_ERROR, 2

  IF N_ELEMENTS(wkdir) EQ 0 THEN wkdir = wkdir_global

  RETURN, sdf_explorer(wkdir, snapshot=snapshot, _struct=1)
END

; --------------------------------------------------------------------------

PRO list_variables, snapshot, wkdir

  COMPILE_OPT idl2
  ON_ERROR, 2

  IF (N_ELEMENTS(snapshot) EQ 0) THEN BEGIN
    PRINT, "Usage: list_variables, snapnumber[, <wkdir>]"
    RETURN
  ENDIF

  q = getdata(snapshot, wkdir, _retro=1, /_variables)
END

; --------------------------------------------------------------------------

PRO quick_view, wkdir, snapshot=snapshot

  COMPILE_OPT idl2
  COMMON background, wkdir_global
  ON_ERROR, 2

  IF (N_ELEMENTS(wkdir_in) EQ 0) THEN wkdir = wkdir_global

  a = create_sdf_visualizer(wkdir, snapshot=snapshot)
END

; --------------------------------------------------------------------------

FUNCTION get_wkdir
  COMPILE_OPT idl2
  COMMON background, wkdir_global, retro_global

  RETURN, wkdir_global
END

; --------------------------------------------------------------------------

PRO set_wkdir, wkdir
  COMPILE_OPT idl2
  COMMON background, wkdir_global, retro_global

  IF (N_ELEMENTS(wkdir) NE 0) THEN wkdir_global = wkdir
END

; --------------------------------------------------------------------------

PRO init_StartPIC
  COMPILE_OPT idl2, hidden
  COMMON background, wkdir_global, retro_global
  COMMON gdlset, gdl
  DEFSYSV, '!GDL', EXISTS=gdl

  init_widget
  init_SDFHelp

  retro_global = 1
  set_wkdir, "Data"

  ;device, true_color=24
  device, decompose=0
  device, retain=2

  !p.charsize = 2
  loadct, 3
END
