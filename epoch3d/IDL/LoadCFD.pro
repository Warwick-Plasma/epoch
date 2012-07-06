FUNCTION readvar, handle, varstruct, offset

  COMPILE_OPT idl2, hidden
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

FUNCTION LoadCFDFile, filename, _variables=requestv, _silent=silent, $
    _retro=retro, _only_md=only_md, _errvar=errval, _var_list=var_list, $
    _block_types=block_types, _block_dims=block_dims, _name_list=name_list, $
    _extra=extra

  COMPILE_OPT idl2, hidden
  COMMON BlockTypes, TYPE_ADDITIONAL, TYPE_MESH, TYPE_MESH_VARIABLE, $
      TYPE_SNAPSHOT
  COMMON MeshTypes, MESH_CARTESIAN, MESH_PARTICLE

  ON_ERROR, 2

  display = 1 - N_ELEMENTS(silent)

  IF (KEYWORD_SET(retro) and display) THEN BEGIN
    PRINT, "WARNING! loading files in retro mode. It is recommended that " $
        + "you use modern mode."
  ENDIF

  Version = 1L
  Revision = 1L

  base_mesh_in_place = 0

  offset = 0LL

  IF N_PARAMS() EQ 0 THEN BEGIN
    IF (display) THEN PRINT, $
        "Usage: result = LoadCFDFile(Filename[, /variables])"
    RETURN, "Usage: result = LoadCFDFile(Filename[, /variables])"
  ENDIF

  ; array of names of parameters
  IF (N_ELEMENTS(extra) NE 0) THEN BEGIN
    name_arr = TAG_NAMES(extra)
    element_block = INTARR(N_ELEMENTS(name_arr))
    element_block[*] = 0
  ENDIF ELSE BEGIN
    name_arr = ""
    element_block = 1
  ENDELSE

  CLOSE, 1
  OPENR, 1, filename
  fileheader = readvar(1, {cfd:BYTARR(3), headeroff:0L, blockheaderoff:0L, $
      Version:0L, Revision:0L, MaxString:0L, nblocks:0L, endianness:0L, $
      start_sec:0L, start_millisec:0L, step:0L, time:0D}, 0)

  MaxStringLen = fileheader.MaxString
  offset = fileheader.headeroff

  ; Whole load of boring tests

  IF (STRING(fileheader.cfd) NE "CFD") THEN BEGIN
    IF (display) THEN PRINT, "The file ", filename, " is not a valid CFD file"
    CLOSE, 1
    RETURN, 0
  ENDIF

  IF (fileheader.Version GT Version) THEN BEGIN
    IF (display) THEN BEGIN
      PRINT, "The file ", filename, $
          " is of a version too high to be read by this program"
      PRINT, "Please contact the CFSA, University of Warwick " $
          + "to obtain a new version"
    ENDIF
    CLOSE, 1
    RETURN, 0
  ENDIF

  IF (fileheader.Revision GT Revision) THEN BEGIN
    IF (display) THEN BEGIN
      PRINT, "WARNING : The file ", filename, $
          " has a higher revision number than this reader"
      PRINT, "Not all data in the file will be available"
      PRINT, "Please contact the CFSA, University of Warwick " $
          + "to obtain a new version"
    ENDIF
  ENDIF

  IF (fileheader.nblocks LE 0) THEN BEGIN
    IF (display) THEN BEGIN
      PRINT, "The file ", filename, " either contains no blocks or is corrupted"
    ENDIF
    CLOSE, 1
    RETURN, 0
  ENDIF

  ; The file seems valid, spool through blocks

  var_list = strarr(fileheader.nblocks)
  block_types = intarr(fileheader.nblocks)
  block_dims = intarr(fileheader.nblocks)

  IF ((fileheader.version EQ 1) AND (fileheader.revision EQ 0)) THEN BEGIN
    f = {filename: filename}
  ENDIF ELSE BEGIN
    f = {filename: filename, snapshot: fileheader.step, time: fileheader.time}
  ENDELSE

  IF (N_ELEMENTS(requestv) NE 0) THEN BEGIN
    IF (display) THEN PRINT, "Available elements are "
  ENDIF

  vBlock = 0
  FOR iBlock = 0, fileheader.nblocks-1 DO BEGIN
    blockheader = readvar(1, {Name:BYTARR(MaxStringLen), $
        Class:BYTARR(MaxStringLen), Type:0L, BlockMDLen:0LL, BlockLen:0LL}, $
        offset)

    ; Read what we know of the block header, so skip the rest
    offset = offset + fileheader.blockheaderoff
    IF (N_ELEMENTS(requestv) NE 0) THEN BEGIN
      HandleBlock, fileheader, blockheader, f, offset, name_arr, $
          /only_md, md=md, retro=retro
      q = ReturnIDLUsable(blockheader, md)
      IF (q EQ 1) THEN BEGIN
        IF (display) THEN PRINT, STRTRIM(STRING(vBlock + 1), 2), ") ", $
            STRTRIM(STRING(blockheader.Name), 2), " (" + $
            swapchr(STRTRIM(STRING(blockheader.Class), 2), ' ', '_') + $
            ") : " + ReturnFriendlyTypeName(blockheader, md)
        vBlock = vBlock + 1
      ENDIF
      element_block[*] = 1
    ENDIF ELSE BEGIN
      HandleBlock, fileheader, blockheader, f, offset, name_arr, $
          element_block, md=md, retro=retro, only_md=only_md
    ENDELSE
    var_list[iBlock] = STRTRIM(STRING(blockheader.name))
    block_types[iBlock] = blockheader.type
    IF (N_ELEMENTS(md) GT 0) THEN block_dims[iBlock] = md.nd
    offset = offset + blockheader.BlockLen
  ENDFOR

  ; Check for any names given which have not been understood
  Errcount = 0
  FOR iEl = 0, N_ELEMENTS(name_arr)-1 DO BEGIN
    IF (element_block[iEl] EQ 0) THEN BEGIN
      PRINT, "WARNING! Unrecognised variable requested (", name_arr[iEl], ")"
      Errcount = Errcount + 1
    ENDIF
  ENDFOR

  IF (Errcount NE 0) THEN BEGIN
    PRINT, "You have specified nonexistant variables. To list available"
    PRINT, "variables, use the 'list_variables' routine."
  ENDIF

  CLOSE, 1

  RETURN, f
END

; --------------------------------------------------------------------------

PRO HandleBlock, fileheader, blockheader, outputobject, offset, name_arr, $
    element_block, md=md, only_md=only_md, retro=retro

  COMPILE_OPT idl2, hidden
  COMMON BlockTypes, TYPE_ADDITIONAL, TYPE_MESH, TYPE_MESH_VARIABLE, $
      TYPE_SNAPSHOT

  NameMatch = CheckName(blockheader, name_arr, element_block)
  IF (NameMatch EQ 1 || blockheader.Type EQ TYPE_SNAPSHOT) THEN BEGIN
    IF (blockheader.Type EQ TYPE_MESH) THEN BEGIN
      GetMesh, fileheader, blockheader, outputobject, offset, $
          only_md=only_md, md=md, byname=NameMatch, retro=retro
      RETURN
    ENDIF
    IF (blockheader.Type EQ TYPE_MESH_VARIABLE) THEN  BEGIN
      GetMeshVar, fileheader, blockheader, outputobject, offset, $
          only_md=only_md, md=md, retro=retro
      RETURN
    ENDIF
    IF (blockheader.Type EQ TYPE_SNAPSHOT) THEN  BEGIN
      GetSnapShot, fileheader, blockheader, outputobject, offset, $
          only_md=only_md, md=md, retro=retro
      RETURN
    ENDIF
  ENDIF
END

; --------------------------------------------------------------------------

PRO GetMesh, file_header, block_header, output_struct, offset, $
    only_md=only_md, md=md, byname=byname, retro=retro

  COMPILE_OPT idl2, hidden
  COMMON MeshTypes, MESH_CARTESIAN, MESH_PARTICLE
  COMMON ParticleCoords, PARTICLE_CARTESIAN, PARTICLE_POLAR, $
      PARTICLE_CYLINDRICAL

  mesh_header = readvar(1, {MeshType:0L, nd:0L, sof:0L}, offset)

  mdonly_f = 0
  IF (n_elements(only_md) NE 0) THEN mdonly_f = 1
  byname_f = 0
  IF (n_elements(byname) NE 0) THEN byname_f = byname

  IF (mesh_header.MeshType EQ MESH_CARTESIAN) THEN BEGIN
    ; Read in the cartesian metadata
    ; Read in the actual mesh

    names = ["x", "y", "z", "a", "b", "c", "d", "e", "f"]

    mesh_header = readvar(1, {MeshType:0L, nd:0L, sof:0L, $
        npts:LONARR(mesh_header.nd)}, offset)

    IF (mdonly_f NE 1) THEN BEGIN
      IF (mesh_header.sof EQ 4) THEN BEGIN
        datastruct = CREATE_STRUCT(names[0], FLTARR(mesh_header.npts[0]))
        FOR iDim = 1, mesh_header.nd-1 DO BEGIN
          datastruct = CREATE_STRUCT(datastruct, names[iDim], $
              FLTARR(mesh_header.npts[iDim]))
        ENDFOR
      ENDIF
      IF (mesh_header.sof EQ 8) THEN BEGIN
        datastruct = CREATE_STRUCT(names[0], DBLARR(mesh_header.npts[0]))
        FOR iDim = 1, mesh_header.nd-1 DO BEGIN
          datastruct = CREATE_STRUCT(datastruct, names[iDim], $
              DBLARR(mesh_header.npts[iDim]))
        ENDFOR
      ENDIF
    ENDIF

    mesh_header = CREATE_STRUCT(mesh_header, 'FRIENDLYNAME', $
        STRTRIM(STRING(block_header.Name), 2))

    IF(mdonly_f NE 1) THEN BEGIN
      d = readvar(1, datastruct, offset + block_header.BlockMDLen)
      d = CREATE_STRUCT(d, mesh_header)
    ENDIF
  ENDIF ELSE IF (mesh_header.MeshType EQ MESH_PARTICLE AND $
      byname_f EQ 1) THEN BEGIN
    mesh_header = readvar(1, {MeshType:0L, nd:0L, sof:0L, CoordType:0L, $
        npart:0LL}, offset)
    IF (mdonly_f NE 1) THEN BEGIN
      IF (mesh_header.sof EQ 4) THEN $
          datastruct = CREATE_STRUCT("ParticlePositions", $
              FLTARR(mesh_header.npart, mesh_header.nd))
      IF (mesh_header.sof EQ 8) THEN $
          datastruct = CREATE_STRUCT("ParticlePositions", $
              DBLARR(mesh_header.npart, mesh_header.nd))
      ; PRINT, "Warning, you have loaded a particle mesh. At present, " + $
      ;     "IDL support for particle meshes is very limited"
      ; PRINT, "The raw data has simply been loaded, and you will have to " + $
      ;     "decode it yourself"
      d = readvar(1, datastruct, offset + block_header.BlockMDLen)
      d = CREATE_STRUCT("MetaData", mesh_header, d)
    ENDIF
  ENDIF

  md = mesh_header
  IF(mdonly_f NE 1 AND N_ELEMENTS(d) NE 0) THEN BEGIN
    output_struct = CREATE_STRUCT(output_struct, $
        STRTRIM(STRING(block_header.Name), 2), d)
  ENDIF
END

; --------------------------------------------------------------------------

PRO GetMeshVar, file_header, block_header, output_struct, offset, $
    md=md, only_md=only_md, retro=retro

  COMPILE_OPT idl2, hidden
  COMMON VarTypes, VAR_CARTESIAN, VAR_PARTICLE

  var_header = readvar(1, {VarType:0L, nd:0L, sof:0L}, offset)
  mesh_header = readvar(1, {mesh_id:BYTARR(file_header.MaxString), $
    mesh_class:BYTARR(file_header.MaxString)}, $
    offset + block_header.BlockMDLen - 2 * file_header.MaxString)
  mdonly_f = 1

  IF (N_ELEMENTS(only_md) EQ 0) THEN mdonly_f = 0

  struct_name = 'data'
  IF (KEYWORD_SET(retro)) THEN $
      struct_name = STRTRIM(STRING(block_header.name), 2)

  IF (var_header.VarType EQ VAR_CARTESIAN) THEN BEGIN
    ; Read in the actual variable
    var_header = readvar(1, {VarType:0L, nd:0L, sof:0L, $
        npts:LONARR(var_header.nd)}, offset)
    IF (mdonly_f NE 1) THEN BEGIN
      IF (var_header.sof EQ 4) THEN $
          datastruct = CREATE_STRUCT(struct_name, FLTARR(var_header.npts))
      IF (var_header.sof EQ 8) THEN $
          datastruct = CREATE_STRUCT(struct_name, DBLARR(var_header.npts))
    ENDIF

    var_header = CREATE_STRUCT(var_header, 'FRIENDLYNAME', $
        STRTRIM(STRING(block_header.Name), 2))
    md = var_header
    IF (mdonly_f NE 1) THEN BEGIN
      d = readvar(1, datastruct, offset + block_header.BlockMDLen)
    ENDIF
  ENDIF ELSE IF (var_header.VarType EQ VAR_PARTICLE) THEN BEGIN
    var_header = readvar(1, {VarType:0L, nd:0L, sof:0L, npart:0LL}, offset)
    IF (mdonly_f NE 1) THEN BEGIN
      IF (var_header.sof EQ 4) THEN $
          datastruct = CREATE_STRUCT(struct_name, FLTARR(var_header.npart))
      IF (var_header.sof EQ 8) THEN $
          datastruct = CREATE_STRUCT(struct_name, DBLARR(var_header.npart))
    ENDIF

    md = var_header
    IF (mdonly_f NE 1) THEN BEGIN
      d = readvar(1, datastruct, offset + block_header.BlockMDLen)
    ENDIF
  ENDIF

  ; Attach the meshname and class to the metadata structure
  var_header = CREATE_STRUCT(var_header, mesh_header)

  IF(mdonly_f NE 1 AND N_ELEMENTS(d) NE 0) THEN BEGIN
    obj = CREATE_STRUCT('metadata', var_header, d)
  ENDIF ELSE BEGIN
    obj = CREATE_STRUCT('metadata', var_header)
  ENDELSE

  IF (KEYWORD_SET(retro)) THEN BEGIN
    output_struct = CREATE_STRUCT(output_struct, d)
  ENDIF ELSE BEGIN
    output_struct = CREATE_STRUCT(output_struct, $
        STRTRIM(STRING(block_header.name), 2), obj)
  ENDELSE
END

; --------------------------------------------------------------------------

PRO GetSnapshot, file_header, block_header, output_struct, offset, $
    md=md, only_md=only_md, retro=retro

  COMPILE_OPT idl2, hidden

  snap_header = readvar(1, {Snapshot:0L, Time:0D}, offset)

  mdonly_f = 1
  IF (N_ELEMENTS(only_md) EQ 0) THEN mdonly_f = 0

  IF (mdonly_f NE 1) THEN BEGIN
    output_struct = CREATE_STRUCT(output_struct, snap_header)
  ENDIF
END
