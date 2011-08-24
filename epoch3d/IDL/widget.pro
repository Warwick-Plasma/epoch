PRO viewer_event_handler, event ; event handler

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

  IF (uvalue.name EQ 'colour_button') THEN BEGIN
    XLOADCT, updatecbdata=uvalue.viewer, updatecallback='xloadct_callback', $
        /silent
    draw_image, uvalue.viewer
  ENDIF

  IF (uvalue.name EQ 'draw_button') THEN BEGIN
    draw_image, uvalue.viewer, /force
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

  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error

  WIDGET_CONTROL, event.id, get_uvalue=uvalue ; get the uvalue

  IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_BASE' THEN BEGIN
    ; resize event
    WIDGET_CONTROL, event.id, xsize=205, ysize=400
  ENDIF

  IF TAG_NAMES(event, /STRUCTURE_NAME) EQ 'WIDGET_KILL_REQUEST' THEN BEGIN
    PTR_FREE, uvalue.valid_names
    PTR_FREE, uvalue.valid_types
    PTR_FREE, uvalue.valid_dims
    WIDGET_CONTROL, uvalue.viewer, /destroy
  ENDIF

  IF (N_ELEMENTS(uvalue) EQ 0) THEN RETURN

  IF (uvalue.name EQ 'explorer_close') THEN BEGIN
    load_data, uvalue.viewer
    WIDGET_CONTROL, uvalue.viewer, /destroy
  ENDIF

  IF (uvalue.name EQ 'explorer_cancel') THEN BEGIN
    WIDGET_CONTROL, uvalue.viewer, /destroy
  ENDIF

END

; --------------------------------------------------------------------------

PRO xloadct_callback, data=data

    draw_image, data

END

; --------------------------------------------------------------------------

FUNCTION load_raw, filename, idstruct

  IF (STRUPCASE(STRMID(filename, STRLEN(filename) - 3)) EQ 'CFD') THEN BEGIN
    RETURN, loadcfdfile(filename, _extra=idstruct)
  ENDIF ELSE BEGIN
    RETURN, loadsdffile(filename, _extra=idstruct)
  ENDELSE

END

; --------------------------------------------------------------------------

PRO load_data, viewer, idstruct

  COMMON SDF_View_Internal_data, Loaded_Data

  WIDGET_CONTROL, viewer, get_uvalue=obj_data
  element = WIDGET_INFO(obj_data.list, /list_select)

  struct = {silent:1}
  FOR i = 0, N_ELEMENTS(element) - 1 DO BEGIN
    name = (*obj_data.valid_names)(element(i))
    struct = CREATE_STRUCT(struct, name, 1)
  ENDFOR
  WIDGET_CONTROL, /hourglass
  Loaded_Data = load_raw(obj_data.filename, struct)
  WIDGET_CONTROL, /hourglass

END

; --------------------------------------------------------------------------

PRO draw_image, viewer, force=force, nowset=nowset

  WIDGET_CONTROL, viewer, get_uvalue=obj_data
  IF (NOT obj_data.auto_update AND N_ELEMENTS(force) EQ 0) THEN RETURN
  element = WIDGET_INFO(obj_data.list, /list_select)

  IF (NOT PTR_VALID(obj_data.valid_names) OR element LT 0) THEN BEGIN
    clear_draw_surface, obj_data.viewer
    RETURN
  ENDIF

  ; Use helvetica true type, since it's as good as IDL gets
  DEVICE, SET_FONT='Helvetica', /TT_FONT

  geom = WIDGET_INFO(obj_data.surface, /geometry)
  bar_height = MIN(geom.ysize / 2.0)
  bar1 = FINDGEN(bar_height) / FLOAT(bar_height)
  bar = FLTARR(20, bar_height)
  FOR i = 0, 19 DO BEGIN
    bar(i,*) = bar1
  ENDFOR
  bar_rel_height = FLOAT(bar_height) / FLOAT(geom.ysize)
  bar_rel_pos = 0.5 - (bar_rel_height) / 2.0

  namestruct = {silent:1}
  name = (*obj_data.valid_names)(element)
  namestruct = CREATE_STRUCT(namestruct, name, 1L)
  data = load_raw(obj_data.filename, namestruct)
  a = MAX(TAG_NAMES(data) EQ name, lookup)
  IF (a EQ 0) THEN RETURN
  meshname = STRTRIM(swapchr(STRUPCASE(STRTRIM(STRING( $
      data.(lookup).metadata.mesh_id))), '/', '_'))
  namestruct = CREATE_STRUCT(namestruct, meshname, 1L)
  data = load_raw(obj_data.filename, namestruct)

  ; Have to lookup variable as well, since there's no guarantee that the
  ; Mesh won't have been loaded to an earlier index. In fact it's quite likely.
  a = MAX(TAG_NAMES(data) EQ name, lookup)
  a = MAX(TAG_NAMES(data) EQ meshname, lookup_mesh)

  mesh = data.(lookup_mesh)
  WIDGET_CONTROL, /hourglass
  WIDGET_CONTROL, obj_data.surface, get_value=surface_id
  IF (NOT KEYWORD_SET(nowset)) THEN wset, surface_id

  ndims = SIZE(data.(lookup).data, /n_dimensions)
  ; (*obj_data.valid_dims)(element)

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

  plotter = data.(lookup).data
  range = [MIN(plotter), MAX(plotter)]
  IF (obj_data.logscale AND range[0] LE 0.0) THEN BEGIN
    range[0] = MIN(plotter + range[1] * (plotter LE 0.0))
    plotter(WHERE(plotter LE 0)) = range[0]
  ENDIF
  IF (obj_data.logscale) THEN BEGIN
    plotter = ALOG(plotter)
    range = ALOG(range)
  ENDIF
  IF (ndims EQ 1) THEN BEGIN
    PLOT, mesh.x, plotter, xtitle=axis[0], ytitle=(tag_names(data))[lookup], $
        yrange=range, position=[0.1,0.1,0.95,0.95]
  ENDIF
  IF (ndims EQ 2) THEN BEGIN
    sz = SIZE(plotter)
    nx = sz(1)
    ny = sz(2)
    CONTOUR, plotter, nlevels=40, /fill, /xsty, /ysty, mesh.x(0:nx-1), $
        mesh.y(0:ny-1), xtitle=axis[0], ytitle=axis[1], zrange=range, $
        iso=obj_data.iso, chars=1.2, position=[0.1,0.1,0.85,0.95]
    TVSCL, bar, /normal, 0.86, bar_rel_pos
    txtpos = 0.87
    form = '(G10.0)'
    XYOUTS, /normal, txtpos, bar_rel_pos, $
        STRTRIM(STRING(range[0], format=form)), chars=1.2
    XYOUTS, /normal, txtpos, bar_rel_pos + bar_rel_height, $
        STRTRIM(STRING(range[1], format=form)), chars=1.2
    XYOUTS, /normal, txtpos, bar_rel_pos + bar_rel_height / 2.0, $
        STRTRIM(STRING((range[0] + range[1]) / 2.0, format=form)), chars=1.2
  ENDIF
  WIDGET_CONTROL, /hourglass

END

; --------------------------------------------------------------------------

PRO load_meta_and_populate_sdf, viewer, accepted_types

  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error

  WIDGET_CONTROL, viewer, get_uvalue=obj_data
  IF (STRUPCASE(STRMID(obj_data.filename, STRLEN(obj_data.filename) - 3)) $
      EQ 'CFD') THEN BEGIN
    data = loadcfdfile(obj_data.filename, /silent, /variables, var_list=v, $
        block_types=types, block_dims=dims)
  ENDIF ELSE BEGIN
    data = loadsdffile(obj_data.filename, /silent, /variables, var_list=v, $
        block_types=types, block_dims=dims)
  ENDELSE

  IF (N_ELEMENTS(accepted_types) EQ 0) THEN BEGIN
    n_valid = N_ELEMENTS(v)
  ENDIF ELSE BEGIN
    n_valid = TOTAL(types EQ accepted_types)
  ENDELSE
  valid_names = STRARR(n_valid)
  valid_types = INTARR(n_valid)
  valid_dims = INTARR(n_valid)

  IF (N_ELEMENTS(accepted_types) EQ 0) THEN BEGIN
    valid_names = STRUPCASE(v)
    valid_types = types
    valid_dims = dims
  ENDIF ELSE BEGIN
    cur = 0
    FOR i = 0, N_ELEMENTS(v) - 1 DO BEGIN
      IF (MAX(types(i) EQ accepted_types) EQ 1) THEN BEGIN
        valid_names(cur) = STRUPCASE(v(i))
        valid_types(cur) = types(i)
        valid_dims(cur) = dims(i)
        cur = cur + 1
      ENDIF
    ENDFOR
  ENDELSE

  PTR_FREE, obj_data.valid_names
  PTR_FREE, obj_data.valid_types
  PTR_FREE, obj_data.valid_dims

  obj_data.valid_names = PTR_NEW(valid_names)
  obj_data.valid_types = PTR_NEW(valid_types)
  obj_data.valid_dims = PTR_NEW(valid_dims)

  WIDGET_CONTROL, obj_data.list, set_value=valid_names
  WIDGET_CONTROL, viewer, set_uvalue=obj_data

  ; set the label to be the filename
  WIDGET_CONTROL, obj_data.label, set_value=obj_data.filename

END

; --------------------------------------------------------------------------

PRO clear_draw_surface, viewer, nowset=nowset

  WIDGET_CONTROL, viewer, get_uvalue = obj_data
  WIDGET_CONTROL, obj_data.surface, get_value=surface_id
  data = WIDGET_INFO(obj_data.surface, /geometry)
  IF (NOT KEYWORD_SET(nowset)) THEN wset, surface_id
  ERASE

END

; --------------------------------------------------------------------------

FUNCTION sdf_explorer, filename

  COMMON SDF_View_Internal_data, Loaded_Data

  Loaded_Data = 'Load Cancelled'

  info = FILE_INFO(filename)
  IF (info.exists NE 1) THEN BEGIN
    PRINT, 'File ' + STRTRIM(filename) + ' does not exist'
    RETURN, 'File ' + STRTRIM(filename) + ' does not exist'
  ENDIF

  main = WIDGET_BASE(title='SDF Explorer', xsize=205, ysize=400, $
      tlb_frame_attr=1, /TLB_KILL_REQUEST_EVENTS)
  label_id = WIDGET_LABEL(main, value='', ysize=20, xsize=205)
  list_id = WIDGET_LIST(main , yoffset=25, xsize=29, ysize=20, /multiple)

  ; Base data that all objects need access to
  base_data = {type:0, filename:filename, viewer:main, list:list_id, $
      label:label_id}

  ct_button_id = WIDGET_BUTTON(main, value='Load Data', $
      xsize=200, ysize=20, xoffset=0, yoffset=350, $
      uvalue=CREATE_STRUCT(base_data, 'name', 'explorer_close'))

  ct_button_id = WIDGET_BUTTON(main, value='Cancel', $
      xsize=200, ysize=20, xoffset=0, yoffset=375, $
      uvalue=CREATE_STRUCT(base_data, 'name', 'explorer_cancel'))

  WIDGET_CONTROL, main, set_uvalue=CREATE_STRUCT(base_data, 'name', filename, $
      'valid_names', PTR_NEW(), 'valid_types', PTR_NEW(), 'valid_dims', $
      PTR_NEW())
  WIDGET_CONTROL, list_id, set_uvalue=CREATE_STRUCT(base_data, 'name', 'list')

  load_meta_and_populate_sdf, main
  WIDGET_CONTROL, main, /realize
  XMANAGER, filename, main, event_handler='explorer_event_handler'

  RETURN, Loaded_Data
END

; --------------------------------------------------------------------------

FUNCTION create_sdf_visualizer, filename, set_viewer=view_in

  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error
  ; common info for the older CFD file format
  COMMON BlockTypes, TYPE_ADDITIONAL, TYPE_MESH, TYPE_MESH_VARIABLE, $
      TYPE_SNAPSHOT

  info = FILE_INFO(filename)
  IF (info.exists NE 1) THEN BEGIN
    PRINT, 'File ' + STRTRIM(filename) + ' does not exist'
    RETURN, 'File ' + STRTRIM(filename) + ' does not exist'
  ENDIF

  IF (N_ELEMENTS(view_in) EQ 0) THEN BEGIN
    main = WIDGET_BASE(title='SDF Quick Data Visualizer', xsize=1000, $
        ysize=820, /TLB_KILL_REQUEST_EVENTS, /TLB_SIZE_EVENTS)
    label_id = WIDGET_LABEL(main, value='', ysize=20, xsize=1000)
    draw_id = WIDGET_DRAW(main, xsize=800, ysize=800, yoffset=25, xoffset=210)
    list_id = WIDGET_LIST(main , yoffset=25, xsize=29, ysize=20)

    ; Base data that all objects need access to
    base_data = {type:1, filename:filename, viewer:main, list:list_id, $
        surface:draw_id, label:label_id, auto_update:1, logscale:0, $
        minval:0.0d, maxval:0.0d, use_min:0, use_max:0, cfdfile:0, iso:1}

    IF (STRUPCASE(STRMID(filename, STRLEN(filename) - 3)) EQ 'CFD') THEN BEGIN
      base_data.cfdfile = 1
    ENDIF ELSE BEGIN
      base_data.cfdfile = 0
    ENDELSE

    WIDGET_CONTROL, main, set_uvalue=CREATE_STRUCT(base_data, 'name', $
        'base' + filename, 'valid_names', PTR_NEW(), 'valid_types', $
        PTR_NEW(), 'valid_dims', PTR_NEW())
    WIDGET_CONTROL, list_id, set_uvalue=CREATE_STRUCT(base_data, 'name', 'list')

    ; cur_y is the current y position for the control that I'm defining
    cur_y = 350
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

    ; Button for changing colour table
    ct_button_id = WIDGET_BUTTON(main, value='Save figure', $
        xsize=200, ysize=20, xoffset=0, yoffset=cur_y, $
        uvalue=CREATE_STRUCT(base_data, 'name', 'save_button'))
    cur_y = cur_y + 20

    WIDGET_CONTROL, main, /realize
    WIDGET_CONTROL, checkbox_container, /realize

    XMANAGER, filename, main, /no_block, event_handler='viewer_event_handler'
    XMANAGER, filename, checkbox_container, /no_block, $
        event_handler='viewer_event_handler'
  ENDIF ELSE BEGIN
    main = view_in
    WIDGET_CONTROL, main, get_uvalue = base_data
    IF (base_data.type EQ 0) THEN BEGIN
      PRINT , "Cannot use a handle to a lister to control a visualizer"
      RETURN, view_in
    ENDIF
    base_data.filename = filename
    WIDGET_CONTROL, main, set_uvalue = base_data
  ENDELSE
  clear_draw_surface, main
  IF (base_data.cfdfile) THEN BEGIN
    load_meta_and_populate_sdf, main, TYPE_MESH_VARIABLE
  ENDIF ELSE BEGIN
    load_meta_and_populate_sdf, main, SDF_BlockTypes.PLAIN_VARIABLE
  ENDELSE
  RETURN, main
END
