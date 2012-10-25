PRO init_widget
  COMPILE_OPT idl2, hidden
  COMMON SDF_View_Internal_data, Loaded_Data, newline, use_retro

  newline = STRING(10B)
END

; --------------------------------------------------------------------------

FUNCTION generate_filename, wkdir, snapshot
  COMPILE_OPT idl2, hidden

  file = wkdir + STRING(snapshot, FORMAT='("/",I04.04,".cfd")')
  info = FILE_INFO(file)

  IF info.exists EQ 1 THEN BEGIN
    RETURN, file
  ENDIF ELSE BEGIN
    file = wkdir + STRING(snapshot, FORMAT='("/",I04.04,".sdf")')
    RETURN, file
  ENDELSE
END

; --------------------------------------------------------------------------

FUNCTION count_files, wkdir
  COMPILE_OPT idl2, hidden

  sdf = FLOOR(TOTAL(FILE_SEARCH(wkdir+'/*.sdf') NE ''))
  cfd = FLOOR(TOTAL(FILE_SEARCH(wkdir+'/*.cfd') NE ''))

  IF (sdf NE 0 AND cfd NE 0) THEN BEGIN
    PRINT ,'Cannot use explorer on directory containing both CFD and SDF files'
    RETURN, 0
  END

  RETURN, MAX([sdf,cfd])
END

; --------------------------------------------------------------------------

FUNCTION load_raw, filename, idstruct, only_md=only_md
  COMPILE_OPT idl2, hidden

  IF (STRUPCASE(STRMID(filename, STRLEN(filename) - 3)) EQ 'CFD') THEN BEGIN
    RETURN, LoadCFDFile(filename, _extra=idstruct, _only_md=only_md)
  ENDIF ELSE BEGIN
    RETURN, LoadSDFFile(filename, _extra=idstruct, _only_md=only_md)
  ENDELSE

END

; --------------------------------------------------------------------------

FUNCTION get_sdf_metatext, viewer, element
  COMPILE_OPT idl2, hidden
  COMMON SDF_View_Internal_data, Loaded_Data, newline, use_retro

  WIDGET_CONTROL, viewer, get_uvalue=obj_data

  namestruct = {_silent:1}
  name = (*obj_data.valid_names)[element]
  namestruct = CREATE_STRUCT(namestruct, name, 1L)
  data = load_raw(*obj_data.filename, namestruct, /only_md)
  a = MAX(TAG_NAMES(data) EQ name, lookup)
  IF (a EQ 0) THEN RETURN, ''
  a = MAX(TAG_NAMES(data) EQ name, lookup)

  str = "Name : " + data.(lookup).metadata.friendlyname + newline
  str = str + "Structure Name : " + (TAG_NAMES(data))[lookup] + newline
  IF (MAX(TAG_NAMES(data.(lookup).metadata) EQ 'DIMS') EQ 1) THEN BEGIN
    ndims = (SIZE(data.(lookup).metadata.dims))[1]
    str = str + "Dimensions : " + STRING(ndims, FORMAT='(I1)') + newline
    str = str + "Gridpoints : "
    FOR i = 0, ndims-2 DO BEGIN
      str = str + STRING(data.(lookup).metadata.dims[i], FORMAT='(I4)') + ' x '
    END
    str = str + STRING(data.(lookup).metadata.dims[ndims-1], FORMAT='(I4)') $
        + newline
  ENDIF

  IF (MAX(TAG_NAMES(data.(lookup).metadata) EQ 'NPOINTS') EQ 1) THEN BEGIN
    npoints = data.(lookup).metadata.npoints
    str = str + "Number of points : " + STRING(npoints, FORMAT='(I7)') + newline
  ENDIF

  IF (MAX(TAG_NAMES(data.(lookup).metadata) EQ 'MESH_ID') EQ 1) THEN BEGIN
    str = str + "Associated Grid : " $
        + STRTRIM(STRING(data.(lookup).metadata.mesh_id), 2) + newline
  END

  RETURN, str
END

; --------------------------------------------------------------------------

PRO viewer_event_handler, event ; event handler

  COMPILE_OPT idl2, hidden
  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error

  WIDGET_CONTROL, event.id, get_uvalue=uvalue ; get the uvalue

  IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_BASE' THEN BEGIN
    data = WIDGET_INFO(uvalue.surface, /geometry)
    ; resize event
    draw_wid = event.x - data.xoffset
    draw_hi = event.y - data.yoffset

    win_wid = event.x
    win_hi = event.y
    IF (draw_wid LT 200 OR draw_hi LT 200) THEN BEGIN
      IF (draw_wid LT 200) THEN BEGIN
        draw_wid = 200
        win_wid = data.xoffset + draw_wid
      ENDIF
      IF (draw_wid LT 200) THEN BEGIN
        draw_hi = 200
        win_hi = data.yoffset + draw_hi
      ENDIF
    ENDIF
    WIDGET_CONTROL, uvalue.viewer, /REALIZE, UPDATE=0
    WIDGET_CONTROL, uvalue.slide, xsize=win_wid
    WIDGET_CONTROL, uvalue.surface, xsize=draw_wid, ysize=draw_hi
    WIDGET_CONTROL, uvalue.label, xsize=win_wid
    WIDGET_CONTROL, uvalue.viewer, xsize=win_wid, ysize=win_hi
    WIDGET_CONTROL, uvalue.viewer, /MAP, /UPDATE
    draw_image, uvalue.viewer
  ENDIF

  IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_KILL_REQUEST' THEN BEGIN
    PTR_FREE, uvalue.valid_names
    PTR_FREE, uvalue.valid_types
    PTR_FREE, uvalue.valid_dims
    WIDGET_CONTROL, uvalue.viewer, /destroy
  ENDIF

  IF (N_ELEMENTS(uvalue) EQ 0) THEN RETURN

  IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_SLIDER' THEN BEGIN
    WIDGET_CONTROL, uvalue.viewer, get_uvalue=vuv
    filename = generate_filename(vuv.wkdir, event.value)
    viewer_load_new_file, uvalue.viewer, filename
  ENDIF

  IF (uvalue.name EQ 'colour_button') THEN BEGIN
    XLOADCT, updatecbdata=uvalue.viewer, updatecallback='xloadct_callback', $
        /silent
    draw_image, uvalue.viewer
  ENDIF

  IF (uvalue.name EQ 'draw_button') THEN BEGIN
    draw_image, uvalue.viewer, /force
  ENDIF

  IF (uvalue.name EQ 'iplot_button') THEN BEGIN
    draw_image, uvalue.viewer, /iplot
  ENDIF


  IF (uvalue.name EQ 'save_button') THEN BEGIN
    SET_PLOT, 'ps'
    filename = dialog_pickfile(dialog_parent=uvalue.viewer, $
        default_extension='eps', $
        filter=[['Encapsulated Postscript'], ['*.eps']])
    DEVICE, /encaps, /color, filename=filename
    draw_image, uvalue.viewer, /force, /nowset
    DEVICE, /close
    SET_PLOT, 'x'
  ENDIF

  IF (uvalue.name EQ 'list') THEN BEGIN
    draw_image, uvalue.viewer
  ENDIF

  IF (uvalue.name EQ 'aupdate') THEN BEGIN
    WIDGET_CONTROL, uvalue.viewer, get_uvalue=obj_data
    obj_data.auto_update = 1 - obj_data.auto_update
    WIDGET_CONTROL, uvalue.viewer, set_uvalue=obj_data
    draw_image, uvalue.viewer
  ENDIF

  IF (uvalue.name EQ 'logscale') THEN BEGIN
    WIDGET_CONTROL, uvalue.viewer, get_uvalue=obj_data
    obj_data.logscale = 1 - obj_data.logscale
    WIDGET_CONTROL, uvalue.viewer, set_uvalue=obj_data
    draw_image, uvalue.viewer
  ENDIF

  IF (uvalue.name EQ 'forceaspect') THEN BEGIN
    WIDGET_CONTROL, uvalue.viewer, get_uvalue=obj_data
    obj_data.iso = 1 - obj_data.iso
    WIDGET_CONTROL, uvalue.viewer, set_uvalue=obj_data
    draw_image, uvalue.viewer
  ENDIF

END

; --------------------------------------------------------------------------

PRO explorer_event_handler, event ; event handler

  COMPILE_OPT idl2, hidden
  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error

  WIDGET_CONTROL, event.id, get_uvalue=uvalue ; get the uvalue

  IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_BASE' THEN BEGIN
    ; resize event
    WIDGET_CONTROL, event.id, xsize=605, ysize=400
  ENDIF

  IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_KILL_REQUEST' THEN BEGIN
    PTR_FREE, uvalue.valid_names
    PTR_FREE, uvalue.valid_types
    PTR_FREE, uvalue.valid_dims
    WIDGET_CONTROL, uvalue.viewer, /destroy
  ENDIF

  IF (N_ELEMENTS(uvalue) EQ 0) THEN RETURN

  IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_LIST' THEN BEGIN
    text = get_sdf_metatext(uvalue.viewer, event.index)
    WIDGET_CONTROL, uvalue.panel, set_value=text
  ENDIF

  IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_SLIDER' THEN BEGIN
    WIDGET_CONTROL, uvalue.viewer, get_uvalue=vuv
    filename = generate_filename(vuv.wkdir, event.value)
    explorer_load_new_file, uvalue.viewer, filename
  ENDIF

  IF (uvalue.name EQ 'explorer_close') THEN BEGIN
    load_data, uvalue.viewer, errcode
    IF (errcode EQ 0) THEN WIDGET_CONTROL, uvalue.viewer, /destroy
  ENDIF

  IF (uvalue.name EQ 'explorer_cancel') THEN BEGIN
    WIDGET_CONTROL, uvalue.viewer, /destroy
  ENDIF

END

; --------------------------------------------------------------------------

PRO xloadct_callback, data=data

  COMPILE_OPT idl2, hidden
  draw_image, data

END

; --------------------------------------------------------------------------

PRO load_data, viewer, errcode

  COMPILE_OPT idl2, hidden
  COMMON SDF_View_Internal_data, Loaded_Data, newline, use_retro

  errcode = 1

  WIDGET_CONTROL, viewer, get_uvalue=obj_data
  element = WIDGET_INFO(obj_data.list, /list_select)

  ; Handle case where no items were selected
  IF (N_ELEMENTS(element) EQ 1 AND element[0] LT 0) THEN RETURN

  struct = {_silent:1,_retro:use_retro}
  FOR i = 0, N_ELEMENTS(element) - 1 DO BEGIN
    name = (*obj_data.valid_names)[element[i]]
    struct = CREATE_STRUCT(struct, name, 1)
  ENDFOR
  WIDGET_CONTROL, /hourglass
  Loaded_Data = load_raw(*obj_data.filename, struct)
  WIDGET_CONTROL, /hourglass
  errcode = 0

END

; --------------------------------------------------------------------------

PRO draw_image, viewer, force=force, nowset=nowset, iplot=iplot

  COMPILE_OPT idl2, hidden
  WIDGET_CONTROL, viewer, get_uvalue=obj_data
  IF (NOT obj_data.auto_update AND N_ELEMENTS(force) EQ 0) THEN BEGIN
    clear_draw_surface, obj_data.viewer, text='Click "Redraw" to show image'
    RETURN
  ENDIF
  element = WIDGET_INFO(obj_data.list, /list_select)

  IF (NOT PTR_VALID(obj_data.valid_names) OR element LT 0) THEN BEGIN
    clear_draw_surface, obj_data.viewer, text='Invalid variable requested'
    RETURN
  ENDIF

  ; Use helvetica true type, since it's as good as IDL gets
  DEVICE, SET_FONT='Helvetica', /TT_FONT

  geom = WIDGET_INFO(obj_data.surface, /geometry)
  bar_height = MIN(geom.ysize / 2.0)
  bar1 = FINDGEN(bar_height) / FLOAT(bar_height)
  bar = FLTARR(20, bar_height)
  FOR i = 0, 19 DO BEGIN
    bar[i,*] = bar1
  ENDFOR
  bar_rel_height = FLOAT(bar_height) / FLOAT(geom.ysize)
  bar_rel_pos = 0.5 - (bar_rel_height) / 2.0

  WIDGET_CONTROL, /hourglass
  namestruct = {_silent:1}
  name = (*obj_data.valid_names)[element]
  namestruct = CREATE_STRUCT(namestruct, name, 1L)
  data = load_raw(*obj_data.filename, namestruct, /only_md)
  a = MAX(TAG_NAMES(data) EQ name, lookup)
  IF (a EQ 0) THEN RETURN

  ndims = FLOOR(total(data.(lookup).metadata.dims GT 1))
  IF (ndims EQ 3) THEN BEGIN
    IF (NOT KEYWORD_SET(force) AND NOT KEYWORD_SET(iplot)) THEN BEGIN
      clear_draw_surface, obj_data.viewer, $
          text='Please click on iPlot button to launch an iVolume renderer'
      RETURN
    ENDIF
  ENDIF

  meshname = STRTRIM(swapchr(STRUPCASE(STRTRIM(STRING( $
      data.(lookup).metadata.mesh_id))), '/', '_'))
  namestruct = CREATE_STRUCT(namestruct, meshname, 1L)
  data = load_raw(*obj_data.filename, namestruct)

  ; Have to lookup variable as well, since there's no guarantee that the
  ; Mesh won't have been loaded to an earlier index. In fact it's quite likely.
  a = MAX(TAG_NAMES(data) EQ name, lookup)
  a = MAX(TAG_NAMES(data) EQ meshname, lookup_mesh)

  mesh = data.(lookup_mesh)
  WIDGET_CONTROL, obj_data.surface, get_value=surface_id
  IF (NOT KEYWORD_SET(nowset)) THEN wset, surface_id

  ndims = SIZE(data.(lookup).data, /n_dimensions)

  ; Do something clever with the axes
  lookup_labels = -1
  a = MAX(TAG_NAMES(data.(lookup_mesh)) EQ "LABELS")
  axis = ['x', 'y', 'z']
  IF (a GT 0) THEN BEGIN
    FOR i = 0, ndims - 1 DO BEGIN
      axis[i] = data.(lookup_mesh).labels[i] + '(' $
          + data.(lookup_mesh).units[i] + ')'
    ENDFOR
  ENDIF

  IF (KEYWORD_SET(iplot)) THEN BEGIN
    clear_draw_surface, obj_data.viewer, $
        text='Using IDL iPlot to render. Press "Redraw" to show native display.'
  ENDIF

  plotter = data.(lookup).data
  range = [MIN(plotter), MAX(plotter)]
  disprange = range
  IF (obj_data.logscale AND range[0] LE 0.0) THEN BEGIN
    range[0] = MIN(plotter + range[1] * (plotter LE 0.0))
    disprange[0] = range[0]
    plotter[WHERE(plotter LE 0)] = range[0]
  ENDIF
  IF (obj_data.logscale) THEN BEGIN
    plotter = ALOG(plotter)
    range = ALOG(range)
    disprange = ALOG(disprange)
  ENDIF
  IF (ndims EQ 1) THEN BEGIN
    sz = SIZE(plotter)
    nx = sz[1]
    IF (KEYWORD_SET(iplot)) THEN BEGIN
      iPLOT, mesh.x[0:nx-1], plotter, xtitle=axis[0], $
          ytitle=data.(lookup).metadata.friendlyname + '(' + $
          STRTRIM(STRING(data.(lookup).metadata.units)) + ')', $
          yrange=range, /disable_splash_screen, /no_saveprompt
    ENDIF ELSE BEGIN
      PLOT, mesh.x[0:nx-1], plotter, xtitle=axis[0], $
          ytitle=data.(lookup).metadata.friendlyname + '(' + $
          STRTRIM(STRING(data.(lookup).metadata.units)) + ')', $
          yrange=range, position=[0.2,0.2,0.95,0.95]
    ENDELSE
  ENDIF
  IF (ndims EQ 2) THEN BEGIN
    sz = SIZE(plotter)
    nx = sz[1]
    ny = sz[2]
    IF (KEYWORD_SET(iplot)) THEN BEGIN
      iCONTOUR, plotter, n_levels=40, /fill, rgb_table=1, /insert_colorbar, $
          mesh.x[0:nx-1], mesh.y[0:ny-1], xtitle=axis[0], ytitle=axis[1], $
          zrange=range, /disable_splash_screen, /no_saveprompt
    ENDIF ELSE BEGIN
      CONTOUR, plotter, nlevels=40, /fill, /xsty, /ysty, mesh.x[0:nx-1], $
          mesh.y[0:ny-1], xtitle=axis[0], ytitle=axis[1], zrange=range, $
          iso=obj_data.iso, chars=1.2, position=[0.1,0.1,0.85,0.95]
      TVSCL, bar, /normal, 0.86, bar_rel_pos
      txtpos = 0.87
      form = '(G10.0)'
      XYOUTS, /normal, txtpos, bar_rel_pos, $
          STRTRIM(STRING(disprange[0], format=form)), chars=1.2
      XYOUTS, /normal, txtpos, bar_rel_pos + bar_rel_height, $
          STRTRIM(STRING(disprange[1], format=form)), chars=1.2
      XYOUTS, /normal, txtpos, bar_rel_pos + bar_rel_height / 2.0, $
          STRTRIM(STRING((disprange[0] + disprange[1]) / 2.0, format=form)), $
          chars=1.2
    ENDELSE
  ENDIF

  IF (ndims EQ 3) THEN BEGIN
   ivolume, plotter, /disable_splash_screen, /no_saveprompt
  ENDIF
  WIDGET_CONTROL, /hourglass

END

; --------------------------------------------------------------------------

PRO load_meta_and_populate_sdf, viewer, accepted_types

  COMPILE_OPT idl2, hidden
  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error

  WIDGET_CONTROL, /hourglass
  WIDGET_CONTROL, viewer, get_uvalue=obj_data
  IF (obj_data.cfdfile EQ 1) THEN BEGIN
    data = LoadCFDFile(*obj_data.filename, /_silent, /_variables, _var_list=v, $
        _block_types=types, _block_dims=dims)
    names = v
  ENDIF ELSE BEGIN
    data = LoadSDFFile(*obj_data.filename, /_silent, /_variables, _var_list=v, $
        _block_types=types, _block_dims=dims, _name_list=names, _retro=0)
  ENDELSE

  IF (N_ELEMENTS(accepted_types) EQ 0) THEN BEGIN
    n_valid = N_ELEMENTS(v)
  ENDIF ELSE BEGIN
    n_valid = 0
    FOR i = 0, N_ELEMENTS(v) - 1 DO BEGIN
     n_valid = n_valid + MAX(types[i] EQ accepted_types)
    END
  ENDELSE
  valid_names = STRARR(n_valid)
  friendly_names = STRARR(n_valid)
  valid_types = INTARR(n_valid)
  valid_dims = INTARR(n_valid)

  IF (N_ELEMENTS(accepted_types) EQ 0) THEN BEGIN
    valid_names = STRUPCASE(v)
    friendly_names = names
    valid_types = types
    valid_dims = dims
  ENDIF ELSE BEGIN
    cur = 0
    FOR i = 0, N_ELEMENTS(v) - 1 DO BEGIN
      IF (MAX(types[i] EQ accepted_types) EQ 1) THEN BEGIN
        valid_names[cur] = STRUPCASE(v[i])
        friendly_names[cur] = names[i]
        valid_types[cur] = types[i]
        valid_dims[cur] = dims[i]
        cur = cur + 1
      ENDIF
    ENDFOR
  ENDELSE

  PTR_FREE, obj_data.valid_names
  PTR_FREE, obj_data.friendly_names
  PTR_FREE, obj_data.valid_types
  PTR_FREE, obj_data.valid_dims

  obj_data.valid_names = PTR_NEW(valid_names)
  obj_data.friendly_names = PTR_NEW(friendly_names)
  obj_data.valid_types = PTR_NEW(valid_types)
  obj_data.valid_dims = PTR_NEW(valid_dims)

  WIDGET_CONTROL, obj_data.list, set_value=friendly_names
  WIDGET_CONTROL, viewer, set_uvalue=obj_data

  ; set the label to be the filename
  WIDGET_CONTROL, obj_data.label, set_value=*obj_data.filename
  WIDGET_CONTROL, /hourglass

END

; --------------------------------------------------------------------------

PRO clear_draw_surface, viewer, nowset=nowset, text=text

  COMPILE_OPT idl2, hidden
  WIDGET_CONTROL, viewer, get_uvalue = obj_data
  WIDGET_CONTROL, obj_data.surface, get_value=surface_id
  data = WIDGET_INFO(obj_data.surface, /geometry)
  IF (NOT KEYWORD_SET(nowset)) THEN wset, surface_id
  ERASE

  IF (N_ELEMENTS(text) NE 0) THEN BEGIN
    data = WIDGET_INFO(viewer, /geometry)
    XYOUTS, 0.0, data.ysize/2, text, /device
  ENDIF

END

; --------------------------------------------------------------------------

FUNCTION sdf_explorer, wkdir, snapshot=snapshot, _struct=struct

  COMPILE_OPT idl2, hidden
  COMMON SDF_View_Internal_data, Loaded_Data, newline, use_retro
  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error
  ; common info for the older CFD file format
  COMMON BlockTypes, TYPE_ADDITIONAL, TYPE_MESH, TYPE_MESH_VARIABLE, $
      TYPE_SNAPSHOT

  Loaded_Data = 'Load Cancelled'

  use_retro = 1
  IF (KEYWORD_SET(struct)) THEN use_retro = 0
  IF (N_ELEMENTS(snapshot) EQ 0) THEN snapshot = 0

  mxcount = count_files(wkdir)
  IF (mxcount EQ 0) THEN RETURN, 'Invalid directory'
  IF (snapshot GT mxcount) THEN snapshot = mxcount
  filename = generate_filename(wkdir, snapshot)
  info = FILE_INFO(filename)
  IF (info.exists NE 1) THEN BEGIN
    PRINT, 'File ' + STRTRIM(filename) + ' does not exist'
    RETURN, 'File ' + STRTRIM(filename) + ' does not exist'
  ENDIF

  main = WIDGET_BASE(title='SDF Explorer', xsize=605, ysize=430, $
      tlb_frame_attr=1, /TLB_KILL_REQUEST_EVENTS)
  label_id = WIDGET_LABEL(main, value='', ysize=20, xsize=605, yoffset=35)
  list_id = WIDGET_LIST(main , yoffset=60, xsize=29, ysize=20, /multiple)
  panel_id = WIDGET_TEXT(main, value='', xoffset=205, yoffset=60, ysize=27, $
      xsize=400, /sensitive)
  slide_id = WIDGET_SLIDER(main, minimum=0, maximum=mxcount-1, $
      value=snapshot, xsize=605, ysize=35, $
      uvalue=CREATE_STRUCT({viewer:main, name:'timeslider'}))

  ; Base data that all objects need access to
  base_data = {type:0, filename:PTR_NEW(filename), wkdir:wkdir, viewer:main, $
      list:list_id, cfdfile:0, label:label_id, panel:panel_id, slide:slide_id}

  IF (STRUPCASE(STRMID(filename, STRLEN(filename) - 3)) EQ 'CFD') THEN BEGIN
    base_data.cfdfile = 1
  ENDIF ELSE BEGIN
    base_data.cfdfile = 0
  ENDELSE

  ct_button_id = WIDGET_BUTTON(main, value='Load Data', $
      xsize=200, ysize=20, xoffset=0, yoffset=385, $
      uvalue=CREATE_STRUCT(base_data, 'name', 'explorer_close'))

  ct_button_id = WIDGET_BUTTON(main, value='Cancel', $
      xsize=200, ysize=20, xoffset=0, yoffset=410, $
      uvalue=CREATE_STRUCT(base_data, 'name', 'explorer_cancel'))

  WIDGET_CONTROL, main, set_uvalue=CREATE_STRUCT(base_data, 'name', filename, $
      'valid_names', PTR_NEW(), 'friendly_names', PTR_NEW(), 'valid_types', $
      PTR_NEW(), 'valid_dims', PTR_NEW())
  WIDGET_CONTROL, list_id, set_uvalue=CREATE_STRUCT(base_data, 'name', 'list')

  IF (base_data.cfdfile) THEN BEGIN
    load_meta_and_populate_sdf, main, INDGEN(12)+1
  ENDIF ELSE BEGIN
    load_meta_and_populate_sdf, main, INDGEN(12)+1
  ENDELSE
  WIDGET_CONTROL, main, /realize
  XMANAGER, filename, main, event_handler='explorer_event_handler'

  RETURN, Loaded_Data
END

; --------------------------------------------------------------------------

PRO explorer_load_new_file, viewer, filename

  COMPILE_OPT idl2, hidden
  COMMON SDF_View_Internal_data, Loaded_Data, newline, use_retro
  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error
  ; common info for the older CFD file format
  COMMON BlockTypes, TYPE_ADDITIONAL, TYPE_MESH, TYPE_MESH_VARIABLE, $
      TYPE_SNAPSHOT

  WIDGET_CONTROL, viewer, get_uvalue=obj_data
  PTR_FREE, obj_data.filename
  obj_data.filename = PTR_NEW(filename)
  WIDGET_CONTROL, viewer, set_uvalue=obj_data

  IF (obj_data.cfdfile) THEN BEGIN
    load_meta_and_populate_sdf, viewer, TYPE_MESH_VARIABLE
  ENDIF ELSE BEGIN
    load_meta_and_populate_sdf, viewer, $
        [SDF_BlockTypes.PLAIN_VARIABLE, SDF_BlockTypes.POINT_VARIABLE]
  ENDELSE

END

; --------------------------------------------------------------------------

FUNCTION create_sdf_visualizer, wkdir, snapshot=snapshot

  COMPILE_OPT idl2, hidden
  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error
  ; common info for the older CFD file format
  COMMON BlockTypes, TYPE_ADDITIONAL, TYPE_MESH, TYPE_MESH_VARIABLE, $
      TYPE_SNAPSHOT
  COMMON SDF_View_Internal_data, Loaded_Data, newline, use_retro

  Loaded_Data = 'No Data Loaded'

  IF (N_ELEMENTS(snapshot) EQ 0) THEN snapshot = 0
  filename = generate_filename(wkdir, snapshot)
  info = FILE_INFO(filename)
  IF (info.exists NE 1) THEN BEGIN
    PRINT, 'File ' + STRTRIM(filename) + ' does not exist'
    RETURN, 'File ' + STRTRIM(filename) + ' does not exist'
  ENDIF

  mxcount = count_files(wkdir)
  IF (mxcount EQ 0) THEN RETURN, 'Invalid directory'
  IF (snapshot GT mxcount) THEN snapshot = mxcount
  main = WIDGET_BASE(title='SDF Quick Data Visualizer', xsize=1000, $
      ysize=855, /TLB_KILL_REQUEST_EVENTS, /TLB_SIZE_EVENTS)
  label_id = WIDGET_LABEL(main, value='', ysize=20, xsize=1000, yoffset=35)
  draw_id = WIDGET_DRAW(main, xsize=800, ysize=800, yoffset=60, xoffset=210)
  list_id = WIDGET_LIST(main , yoffset=60, xsize=29, ysize=20)
  slide_id = WIDGET_SLIDER(main, minimum=0, maximum=mxcount-1, $
      value=snapshot, xsize=1000, ysize=35)

  ; Base data that all objects need access to
  base_data = {type:1, filename:PTR_NEW(filename), viewer:main, list:list_id, $
      surface:draw_id, label:label_id, slide:slide_id, auto_update:1, $
      logscale:0, minval:0.0d, maxval:0.0d, use_min:0, use_max:0, cfdfile:0, $
      iso:1, wkdir:wkdir}

  IF (STRUPCASE(STRMID(filename, STRLEN(filename) - 3)) EQ 'CFD') THEN BEGIN
    base_data.cfdfile = 1
  ENDIF ELSE BEGIN
    base_data.cfdfile = 0
  ENDELSE

  WIDGET_CONTROL, main, set_uvalue=CREATE_STRUCT(base_data, 'name', $
      'base' + filename, 'valid_names', PTR_NEW(), 'friendly_names', $
      PTR_NEW(), 'valid_types', PTR_NEW(), 'valid_dims', PTR_NEW())
  WIDGET_CONTROL, list_id, set_uvalue=CREATE_STRUCT(base_data, 'name', 'list')
  WIDGET_CONTROL, slide_id, set_uvalue=CREATE_STRUCT(base_data, 'name', 'slide')

  ; cur_y is the current y position for the control that I'm defining
  cur_y = 385
  ; Checkbox bit under list
  container_height = 75
  checkbox_container = WIDGET_BASE(main, xsize=150, ysize=container_height, $
      yoffset=cur_y, /nonexclusive)
  ; autoupdate
  autoupdate_button = WIDGET_BUTTON(checkbox_container, value='Autoupdate', $
      uvalue=CREATE_STRUCT(base_data, 'name', 'aupdate'))
  WIDGET_CONTROL, autoupdate_button, /set_button
  ; log scale
  log_button = WIDGET_BUTTON(checkbox_container, value='Log scale', $
      uvalue=CREATE_STRUCT(base_data, 'name', 'logscale'))

  ; log scale
  far_button = WIDGET_BUTTON(checkbox_container, value='Force Aspect Ratio', $
      uvalue=CREATE_STRUCT(base_data, 'name', 'forceaspect'))
  WIDGET_CONTROL, far_button, /set_button

  cur_y = cur_y + container_height + 7

  ; Button for redrawing things
  ct_button_id = WIDGET_BUTTON(main, value='Redraw', $
      xsize=200, ysize=20, xoffset=0, yoffset=cur_y, $
      uvalue=CREATE_STRUCT(base_data, 'name', 'draw_button'))
  cur_y = cur_y + 20

  ; Button for changing colour table
  ct_button_id = WIDGET_BUTTON(main, value='Load colour table', $
      xsize=200, ysize=20, xoffset=0, yoffset=cur_y, $
      uvalue=CREATE_STRUCT(base_data, 'name', 'colour_button'))
  cur_y = cur_y + 20

  ; Button for saving files
  ct_button_id = WIDGET_BUTTON(main, value='Save figure', $
      xsize=200, ysize=20, xoffset=0, yoffset=cur_y, $
      uvalue=CREATE_STRUCT(base_data, 'name', 'save_button'))
  cur_y = cur_y + 20

  ct_button_id = WIDGET_BUTTON(main, value='iPlot', $
      xsize=200, ysize=20, xoffset=0, yoffset=cur_y, $
      uvalue=CREATE_STRUCT(base_data, 'name', 'iplot_button'))
  cur_y = cur_y + 20

  WIDGET_CONTROL, main, /realize
  WIDGET_CONTROL, checkbox_container, /realize

  XMANAGER, filename, main, /no_block, event_handler='viewer_event_handler'
  XMANAGER, filename, checkbox_container, /no_block, $
      event_handler='viewer_event_handler'

  clear_draw_surface, main
  IF (base_data.cfdfile) THEN BEGIN
    load_meta_and_populate_sdf, main, TYPE_MESH_VARIABLE
  ENDIF ELSE BEGIN
    load_meta_and_populate_sdf, main, SDF_BlockTypes.PLAIN_VARIABLE
  ENDELSE

  RETURN, loaded_data
END

; --------------------------------------------------------------------------

PRO viewer_load_new_file, viewer, filename

  COMPILE_OPT idl2, hidden
  COMMON SDF_View_Internal_data, Loaded_Data, newline, use_retro
  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error
  ; common info for the older CFD file format
  COMMON BlockTypes, TYPE_ADDITIONAL, TYPE_MESH, TYPE_MESH_VARIABLE, $
      TYPE_SNAPSHOT

  WIDGET_CONTROL, viewer, get_uvalue=obj_data
  element = WIDGET_INFO(obj_data.list, /list_select)
  IF (element GE 0) THEN BEGIN
    name = (*(obj_data.valid_names))[element]
  ENDIF ELSE BEGIN
    name = "blankblank"
  ENDELSE
  PTR_FREE, obj_data.filename
  obj_data.filename = PTR_NEW(filename)
  WIDGET_CONTROL, viewer, set_uvalue=obj_data

  IF (obj_data.cfdfile) THEN BEGIN
    load_meta_and_populate_sdf, viewer, TYPE_MESH_VARIABLE
  ENDIF ELSE BEGIN
    load_meta_and_populate_sdf, viewer, SDF_BlockTypes.PLAIN_VARIABLE
  ENDELSE

  WIDGET_CONTROL, viewer, get_uvalue=obj_data

  IF (MAX(*(obj_data.valid_names) EQ name,loc) NE 0) THEN BEGIN
    WIDGET_CONTROL,obj_data.list,set_list_select=loc
    draw_image, viewer, /force
  ENDIF ELSE BEGIN
    clear_draw_surface, viewer, text='Please select a variable to view'
  ENDELSE

END
