; Takes the string provided by instring and replaces all occurences of
; the character fromchar with the character tochar
FUNCTION swapchr, instring, fromchar, tochar
  stringbytes = BYTE(instring)
  frombytes = BYTE(fromchar)
  occurences = WHERE(stringbytes EQ frombytes(0), nfound)
  IF nfound EQ 0 THEN RETURN, instring
  tobytes = BYTE(tochar)
  stringbytes(occurences) = tobytes(0)
  RETURN, STRING(stringbytes)
END

; --------------------------------------------------------------------------

; Function to test whether a block name appears in a given namelist
; If the namelist is empty then the block is assumed valid
FUNCTION CheckName, block_header, namelist, element_block
  IF (namelist[0] EQ "") THEN RETURN, 1

  FOR i = 0, N_ELEMENTS(namelist)-1 DO BEGIN
    name = STRUPCASE(STRTRIM(STRING(block_header.Name)))
    class = swapchr(STRUPCASE(STRTRIM(STRING(block_header.Class))), ' ', '_')
    IF (namelist[i] EQ name || $
        namelist[i] EQ class || namelist[i] EQ "EMPTY") THEN BEGIN
      element_block[i] = 1
      IF (namelist[i] NE "EMPTY") THEN RETURN, 1
    ENDIF
  ENDFOR
  RETURN, 0
END

; --------------------------------------------------------------------------

; Return whether or not IDL knows how to deal with the block
; Stops people seing VISIT specific blocks etc.
FUNCTION ReturnIDLUsable, block_header, block_metadata
  COMMON BlockTypes, TYPE_ADDITIONAL, TYPE_MESH, TYPE_MESH_VARIABLE, $
      TYPE_SNAPSHOT
  COMMON MeshTypes, MESH_CARTESIAN, MESH_PARTICLE
  COMMON VarTypes, VAR_CARTESIAN, VAR_PARTICLE

  test = block_header.Type
  ; Always see meshes
  IF (test EQ TYPE_MESH) THEN RETURN, 1
  ; Always show fluid variables
  IF (test EQ TYPE_MESH_VARIABLE) THEN RETURN, 1
  ; Don't want to show time types (always load them)
  IF (test EQ TYPE_SNAPSHOT) THEN RETURN, 0

  RETURN, 0
END

; --------------------------------------------------------------------------

; Returns a user friendly name for when people list variables
FUNCTION ReturnFriendlyTypeName, block_header, block_metadata
  COMMON BlockTypes, TYPE_ADDITIONAL, TYPE_MESH, TYPE_MESH_VARIABLE, $
      TYPE_SNAPSHOT
  COMMON MeshTypes, MESH_CARTESIAN, MESH_PARTICLE
  COMMON VarTypes, VAR_CARTESIAN, VAR_PARTICLE


  test = block_header.Type
  IF (test EQ TYPE_MESH) THEN BEGIN
    IF (block_metadata.MeshType EQ MESH_CARTESIAN) THEN BEGIN
      RETURN, STRTRIM(STRING(block_metadata.nd), 2) + "D cartesian grid"
    ENDIF
    IF (block_metadata.MeshType EQ MESH_PARTICLE) THEN BEGIN
      RETURN, STRTRIM(STRING(block_metadata.nd), 2) + "D particle grid"
    ENDIF
  ENDIF

  IF (test EQ TYPE_MESH_VARIABLE) THEN BEGIN
    IF (block_metadata.VarType EQ VAR_CARTESIAN) THEN BEGIN
      RETURN, STRTRIM(STRING(block_metadata.nd), 2) + "D cartesian variable"
    ENDIF
    IF (block_metadata.VarType EQ VAR_PARTICLE) THEN BEGIN
      RETURN, "particle variable"
    ENDIF
  ENDIF
  RETURN, "Unknown block"
END

; --------------------------------------------------------------------------

COMMON BlockTypes, TYPE_ADDITIONAL, TYPE_MESH, TYPE_MESH_VARIABLE, TYPE_SNAPSHOT
COMMON MeshTypes, MESH_CARTESIAN, MESH_PARTICLE
COMMON VarTypes, VAR_CARTESIAN, VAR_PARTICLE
COMMON ParticleCoords, PARTICLE_CARTESIAN, PARTICLE_POLAR, PARTICLE_CYLINDRICAL

TYPE_ADDITIONAL = 0L
TYPE_MESH = 1L
TYPE_MESH_VARIABLE = 2L
TYPE_SNAPSHOT = 3L

MESH_CARTESIAN = 0L
MESH_PARTICLE = 1L

VAR_CARTESIAN = 0L
VAR_PARTICLE = 1L

PARTICLE_CARTESIAN = 0
PARTICLE_POLAR = 1
PARTICLE_CYLINDRICAL = 2

END
