FUNCTION readvar, handle, varstruct, offset
  COMMON gdlset, gdl

  IF (gdl) THEN BEGIN
    outvar = varstruct
    POINT_LUN, handle, offset
    READU, handle, outvar
    RETURN, outvar
  ENDIF ELSE BEGIN
    outvar = ASSOC(handle, varstruct, offset, /PACKED)
    RETURN, outvar[0]
  ENDELSE
END

; --------------------------------------------------------------------------

FUNCTION LoadSDFFile, filename, Variables=requestv, request_classes=requestc, $
    var_list=var_list, name_list=name_list, block_types=block_types, block_dims=block_dims, $
    silent=silent, errval=errval, retro=retro, only_md=only_md, _extra=extra

  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error

  ON_ERROR, 2
  display = 1 - N_ELEMENTS(silent)
  errval = SDF_Error.NONE

  IF (KEYWORD_SET(retro) and display) THEN BEGIN
    PRINT, "WARNING! loading files in retro mode. It is recommended that " $
        + "you use modern mode."
  ENDIF

  base_mesh_in_place = 0

  IF N_PARAMS() EQ 0 THEN BEGIN
    IF (display) THEN PRINT, $
        "Usage: result = LoadSDFFile(Filename[, /variables])"
    errval = SDF_Error.ERR_BAD_USAGE
    RETURN, "Usage: result = LoadSDFFile(Filename[, /variables])"
  ENDIF

  ; array of names of parameters
  IF (N_ELEMENTS(extra) NE 0) THEN BEGIN
    name_arr = TAG_NAMES(extra)
    element_block = INTARR(N_ELEMENTS(name_arr))
    element_block(*) = 0
  ENDIF ELSE BEGIN
    name_arr = ""
    element_block = 1
  ENDELSE

  id_length = SDF_Common.ID_LENGTH

  header_struct = {sdf:BYTARR(4), endianness:0L, version:0L, $
      revision:0L, code_name:BYTARR(id_length), first_block_location:0LL, $
      summary_location:0LL, summary_size:0L, nblocks:0L, $
      block_header_length:0L, step:0L, time:0D, jobid1:0L, jobid2:0L, $
      string_length:0L, code_io_version:0L, restart_flag:BYTARR(1), $
      subdomain_file:BYTARR(1)}

  CLOSE, 1
  OPENR, 1, filename
  file_header = readvar(1, header_struct, 0)

  IF (file_header.endianness NE 16911887) THEN BEGIN
    CLOSE, 1
    OPENR, 1, filename, /SWAP_ENDIAN
    file_header = readvar(1, header_struct, 0)
  ENDIF

  ; Whole load of boring tests

  IF (STRING(file_header.sdf) NE SDF_Common.MAGIC) THEN BEGIN
    IF (display) THEN PRINT, "The file ", filename, " is not a valid SDF file"
    errcode = SDF_Error.BAD_MAGIC
    CLOSE, 1
    RETURN, 0
  ENDIF

  IF (file_header.version GT SDF_Common.VERSION) THEN BEGIN
    IF (display) THEN BEGIN
      PRINT, "The file ", filename, $
          " is of a version too high to be read by this program"
      PRINT, "Please contact the CFSA, University of Warwick " $
          + "to obtain a new version"
    END
    errcode = SDF_Error.BAD_VERSION
    CLOSE, 1
    RETURN, 0
  ENDIF

  IF (file_header.revision GT SDF_Common.REVISION) THEN BEGIN
    PRINT, "WARNING : The file ", filename, $
        " has a higher revision number than this reader"
    PRINT, "Not all data in the file will be available"
    PRINT, "Please contact the CFSA, University of Warwick " $
        + "to obtain a new version"
  ENDIF

  IF (file_header.nblocks LE 0) THEN BEGIN
    IF (display) THEN BEGIN
      PRINT, "The file ", filename, " either contains no blocks or is corrupted"
    ENDIF
    errcode = SDF_Error.BAD_BLOCK_COUNT
    CLOSE, 1
    RETURN, 0
  ENDIF

  ; The file seems valid, spool through blocks

  f = {filename: filename, timestep: file_header.step, time: file_header.time}

  IF (N_ELEMENTS(requestv) NE 0 AND display) THEN BEGIN
    PRINT, "Available elements are "
  ENDIF

  var_list = strarr(file_header.nblocks)
  name_list = strarr(file_header.nblocks)
  block_types = intarr(file_header.nblocks)
  block_dims = intarr(file_header.nblocks)

  string_length = file_header.string_length
  offset = file_header.first_block_location
  vBlock = 0
  FOR iBlock = 0, file_header.nblocks - 1 DO BEGIN
    b = readvar(1, {next_block_location:0LL, data_location:0LL, $
        id:BYTARR(id_length), data_length:0LL, blocktype:0L, datatype:0L, $
        ndims:0L, fname:BYTARR(string_length)}, offset)
    b = CREATE_STRUCT(b, 'start', offset)
    SWITCH b.blocktype OF
    SDF_Blocktypes.PLAIN_MESH:
    SDF_Blocktypes.POINT_MESH:
    SDF_Blocktypes.PLAIN_VARIABLE:
    SDF_Blocktypes.POINT_VARIABLE: BEGIN

      b = CREATE_STRUCT(b, 'idname', STRTRIM(STRING(b.id), 2))
      b = CREATE_STRUCT(b, 'name', STRTRIM(STRING(b.fname), 2))
      b = CREATE_STRUCT(b, 'class', "")
      b.idname = swapchr(b.idname, ' ', '_')
      b.idname = swapchr(b.idname, '/', '_')
      b.name = swapchr(b.name, ' ', '_')
      pos = STRPOS(b.name, '/')
      b.name = swapchr(b.fname, '/', '_')
      IF (pos GT 0) THEN BEGIN
        b.class = STRMID(b.name, 0, pos)
        b.name = STRMID(b.name, pos+1)
      ENDIF

      var_list[iBlock] = STRTRIM(b.idname)
      name_list[iBlock] = STRTRIM(b.fname)
      block_types[iBlock] = b.blocktype
      block_dims[iBlock] = b.ndims


      IF (N_ELEMENTS(requestv) NE 0) THEN BEGIN
        IF (display) THEN BEGIN
          PRINT, STRTRIM(STRING(vBlock + 1), 2) + ") " + b.name + " (" $
              + b.class + ") : " + STRTRIM(STRING(b.ndims), 2) + "D " $
              + SDF_Blocktype_names[b.blocktype]
          vBlock = vBlock + 1
          element_block(*) = 1
        END
      ENDIF ELSE BEGIN
        SDFHandleBlock, file_header, b, f, offset, name_arr, element_block, $
            retro=retro, only_md=only_md
      ENDELSE
      END
    ELSE:
    ENDSWITCH

    offset = b.next_block_location
  ENDFOR

  ; Check for any names given which have not been understood
  Errcount = 0
  FOR iEl = 0, N_ELEMENTS(name_arr)-1 DO BEGIN
    IF (element_block(iEl) EQ 0) THEN BEGIN
      IF (display) THEN BEGIN
        PRINT, "WARNING! Unrecognised variable requested (", name_arr(iEl), ")"
      ENDIF
      Errcount = Errcount + 1
    ENDIF
  ENDFOR

  IF (Errcount NE 0) THEN BEGIN
    IF (display) THEN BEGIN
      PRINT, "You have specified nonexistant variables. To list available " $
          + "variables, use the '/variables' switch."
    ENDIF
  ENDIF

  CLOSE, 1

  RETURN, f
END

; --------------------------------------------------------------------------

PRO SDFHandleBlock, file_header, block_header, outputobject, offset, $
    name_arr, element_block, md=md, retro=retro, only_md = only_md

  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error

  NameMatch = SDFCheckName(block_header, name_arr, element_block)

  IF (NameMatch EQ 1) THEN BEGIN
    CASE block_header.blocktype OF
      SDF_Blocktypes.PLAIN_MESH: BEGIN
        SDFGetPlainMesh, file_header, block_header, outputobject, offset, $
            md=md, retro=retro, only_md=only_md
      END
      SDF_Blocktypes.POINT_MESH: BEGIN
        SDFGetPointMesh, file_header, block_header, outputobject, offset, $
            md=md, retro=retro, only_md=only_md
      END
      SDF_Blocktypes.PLAIN_VARIABLE: BEGIN
        SDFGetPlainVar, file_header, block_header, outputobject, offset, $
            md=md, retro=retro, only_md=only_md
      END
      SDF_Blocktypes.POINT_VARIABLE: BEGIN
        SDFGetPointVar, file_header, block_header, outputobject, offset, $
            md=md, retro=retro, only_md=only_md
      END
    ELSE:
    ENDCASE
  ENDIF
END

; --------------------------------------------------------------------------

PRO SDFGetPlainMesh, file_header, block_header, output_struct, offset, md=md, $
    retro=retro, only_md=only_md

  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error

  mdflag=KEYWORD_SET(only_md)
  id_length = SDF_Common.ID_LENGTH
  offset = block_header.start + file_header.block_header_length
  mesh_header = readvar(1, { mults:DBLARR(block_header.ndims), $
      labels:BYTARR(block_header.ndims * id_length), $
      units:BYTARR(block_header.ndims * id_length), geometry:0L, $
      minval:DBLARR(block_header.ndims), maxval:DBLARR(block_header.ndims), $
      dims:LONARR(block_header.ndims)}, offset)

  ndims = block_header.ndims
  labels = STRARR(ndims)
  units = STRARR(ndims)
  FOR iDim = 0, ndims-1 DO BEGIN
    n0 = iDim * id_length
    n1 = n0 + id_length - 1
    labels[iDim] = STRTRIM(STRING(mesh_header.labels[n0:n1]))
    units[iDim] = STRTRIM(STRING(mesh_header.units[n0:n1]))
  ENDFOR
  labelstr = labels
  labels[0] = 'X'
  IF (ndims GT 1) THEN labels[1] = 'Y'
  IF (ndims GT 2) THEN labels[2] = 'Z'

  IF (~mdflag) THEN BEGIN
    CASE block_header.datatype OF
      SDF_Datatypes.REAL4: BEGIN
        datastruct = CREATE_STRUCT(labels[0], FLTARR(mesh_header.dims[0]))
        FOR iDim = 1, ndims-1 DO BEGIN
          datastruct = CREATE_STRUCT(datastruct, labels[iDim], $
              FLTARR(mesh_header.dims[iDim]))
        ENDFOR
      END
      SDF_Datatypes.REAL8: BEGIN
        datastruct = CREATE_STRUCT(labels[0], DBLARR(mesh_header.dims[0]))
        FOR iDim = 1, ndims-1 DO BEGIN
          datastruct = CREATE_STRUCT(datastruct, labels[iDim], $
              DBLARR(mesh_header.dims[iDim]))
        ENDFOR
      END
    ENDCASE
    offset = block_header.data_location
    d = readvar(1, datastruct, offset)
  ENDIF ELSE BEGIN
    d = CREATE_STRUCT('ONLYMD',1)
  ENDELSE

  d = CREATE_STRUCT(d,'LABELS', labelstr)
  d = CREATE_STRUCT(d, 'UNITS', units)
  d = CREATE_STRUCT(d, 'NPTS', mesh_header.dims)
  mesh_header = CREATE_STRUCT(mesh_header, 'FRIENDLYNAME', STRTRIM(STRING(block_header.fname),2))

  obj = CREATE_STRUCT('metadata', mesh_header, d)
  md = mesh_header
  IF (N_ELEMENTS(d) NE 0) THEN BEGIN
    IF (KEYWORD_SET(retro)) THEN BEGIN
      output_struct = CREATE_STRUCT(output_struct, d)
    ENDIF ELSE BEGIN
      output_struct = CREATE_STRUCT(output_struct, block_header.idname, obj)
    ENDELSE
    ; Hack to add a cell centred grid for plotting node-centred values
    IF (block_header.idname EQ 'grid' AND ~mdflag) THEN BEGIN
      iDim = 0
      nx = mesh_header.dims[iDim] - 1
      x = 0.5 * (d.(iDim)[0:nx-1] + d.(iDim)[1:nx])
      xc = CREATE_STRUCT(labels[iDim], x)
      FOR iDim = 1, ndims-1 DO BEGIN
        nx = mesh_header.dims[iDim] - 1
        x = 0.5 * (d.(iDim)[0:nx-1] + d.(iDim)[1:nx])
        xc = CREATE_STRUCT(xc, labels[iDim], x)
      ENDFOR
      output_struct = CREATE_STRUCT(output_struct, xc)
    ENDIF
  ENDIF
END

; --------------------------------------------------------------------------

PRO SDFGetPointMesh, file_header, block_header, output_struct, offset, md=md, $
    retro=retro, only_md=only_md

  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error

  mdflag=KEYWORD_SET(only_md)
  id_length = SDF_Common.ID_LENGTH
  offset = block_header.start + file_header.block_header_length
  mesh_header = readvar(1, { mults:DBLARR(block_header.ndims), $
      labels:BYTARR(block_header.ndims * id_length), $
      units:BYTARR(block_header.ndims * id_length), geometry:0L, $
      minval:DBLARR(block_header.ndims), maxval:DBLARR(block_header.ndims), $
      npoints:0LL}, offset)

  ndims = block_header.ndims
  labels = STRARR(ndims)
  units = STRARR(ndims)
  FOR iDim = 0, ndims-1 DO BEGIN
    n0 = iDim * id_length
    n1 = n0 + id_length - 1
    labels[iDim] = STRTRIM(STRING(mesh_header.labels[n0:n1]))
    units[iDim] = STRTRIM(STRING(mesh_header.units[n0:n1]))
  ENDFOR
  labelstr = labels
  labels[0] = 'X'
  IF (ndims GT 1) THEN labels[1] = 'Y'
  IF (ndims GT 2) THEN labels[2] = 'Z'

  IF (~mdflag) THEN BEGIN
    CASE block_header.datatype OF
      SDF_Datatypes.REAL4: BEGIN
        datastruct = CREATE_STRUCT(labels[0], FLTARR(mesh_header.npoints))
        FOR iDim = 1, ndims-1 DO BEGIN
          datastruct = CREATE_STRUCT(datastruct, labels[iDim], $
              FLTARR(mesh_header.npoints))
        ENDFOR
      END
      SDF_Datatypes.REAL8: BEGIN
        datastruct = CREATE_STRUCT(labels[0], DBLARR(mesh_header.npoints))
        FOR iDim = 1, ndims-1 DO BEGIN
          datastruct = CREATE_STRUCT(datastruct, labels[iDim], $
              DBLARR(mesh_header.npoints))
        ENDFOR
      END
    ENDCASE
    offset = block_header.data_location
    d = readvar(1, datastruct, offset)
  ENDIF ELSE BEGIN
    d=CREATE_STRUCT('ONLYMD',1)
  ENDELSE

  d = CREATE_STRUCT(d, 'LABELS', labelstr)
  d = CREATE_STRUCT(d, 'UNITS', units)
  d = CREATE_STRUCT(d, 'NPART', mesh_header.npoints)
  mesh_header = CREATE_STRUCT(mesh_header, 'FRIENDLYNAME', STRTRIM(STRING(block_header.fname),2))

  obj = CREATE_STRUCT('metadata', mesh_header, d)
  md = mesh_header
  IF (N_ELEMENTS(d) NE 0) THEN BEGIN
    IF (KEYWORD_SET(retro)) THEN BEGIN
      output_struct = CREATE_STRUCT(output_struct, d)
    ENDIF ELSE BEGIN
      output_struct = CREATE_STRUCT(output_struct, block_header.idname, obj)
    ENDELSE
  ENDIF

END

; --------------------------------------------------------------------------

PRO SDFGetPlainVar, file_header, block_header, output_struct, offset, md=md, $
    retro=retro, only_md=only_md

  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error

  mdflag=KEYWORD_SET(only_md)
  id_length = SDF_Common.ID_LENGTH
  offset = block_header.start + file_header.block_header_length
  var_header = readvar(1, {mult:0D, units:BYTARR(id_length), $
      mesh_id:BYTARR(id_length), dims:LONARR(block_header.ndims), $
      stagger:0L}, offset)

  struct_name = 'data'
  IF (KEYWORD_SET(retro)) THEN struct_name = block_header.idname

  IF (~mdflag) THEN BEGIN
    CASE block_header.datatype OF
      SDF_Datatypes.REAL4: BEGIN
        datastruct = CREATE_STRUCT(struct_name, FLTARR(var_header.dims))
      END
      SDF_Datatypes.REAL8: BEGIN
        datastruct = CREATE_STRUCT(struct_name, DBLARR(var_header.dims))
      END
    ENDCASE
    offset = block_header.data_location
    d = readvar(1, datastruct, offset)
  ENDIF ELSE BEGIN
    d=CREATE_STRUCT('ONLYMD',1)
  ENDELSE

  var_header = CREATE_STRUCT(var_header, 'FRIENDLYNAME', STRTRIM(STRING(block_header.fname),2))
  obj = CREATE_STRUCT('metadata', var_header, d)
  md = var_header
  IF (N_ELEMENTS(d) NE 0) THEN BEGIN
    IF (KEYWORD_SET(retro)) THEN BEGIN
      output_struct = CREATE_STRUCT(output_struct, d)
    ENDIF ELSE BEGIN
      output_struct = CREATE_STRUCT(output_struct, block_header.idname, obj)
    ENDELSE
  ENDIF
END

; --------------------------------------------------------------------------

PRO SDFGetPointVar, file_header, block_header, output_struct, offset, $
    md=md, retro=retro, only_md=only_md

  COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
      SDF_Datatypes, SDF_Error

  mdflag=KEYWORD_SET(only_md)
  id_length = SDF_Common.ID_LENGTH
  offset = block_header.start + file_header.block_header_length
  var_header = readvar(1, {mult:0D, units:BYTARR(id_length), $
      mesh_id:BYTARR(id_length), npoints:0LL}, offset)

  struct_name = 'data'
  IF (KEYWORD_SET(retro)) THEN struct_name = block_header.idname

  IF (~mdflag) THEN BEGIN
    CASE block_header.datatype OF
      SDF_Datatypes.REAL4: BEGIN
        datastruct = CREATE_STRUCT(struct_name, FLTARR(var_header.npoints))
      END
      SDF_Datatypes.REAL8: BEGIN
        datastruct = CREATE_STRUCT(struct_name, DBLARR(var_header.npoints))
      END
    ENDCASE

    offset = block_header.data_location
    d = readvar(1, datastruct, offset)
  ENDIF ELSE BEGIN
    d = CREATE_STRUCT('MDONLY',1)
  ENDELSE

  var_header = CREATE_STRUCT(var_header, 'FRIENDLYNAME', STRTRIM(STRING(block_header.fname),2))
  obj = CREATE_STRUCT('metadata', var_header, d)
  md = var_header
  IF (N_ELEMENTS(d) NE 0) THEN BEGIN
    IF (KEYWORD_SET(retro)) THEN BEGIN
      output_struct = CREATE_STRUCT(output_struct, d)
    ENDIF ELSE BEGIN
      output_struct = CREATE_STRUCT(output_struct, block_header.idname, obj)
    ENDELSE
  ENDIF
END
