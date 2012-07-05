; Function to test whether a block name appears in a given namelist
; If the namelist is empty then the block is assumed valid
FUNCTION SDFCheckName, blockheader, namelist, element_block

  COMPILE_OPT idl2, hidden

  IF (namelist[0] EQ "") THEN RETURN, 1

  FOR i = 0, N_ELEMENTS(namelist)-1 DO BEGIN
    name = STRUPCASE(blockheader.name)
    class = STRUPCASE(blockheader.class)
    fullname = class + '_' + name
    namecmp = STRUPCASE(STRTRIM(namelist[i]))
    IF (namecmp EQ fullname || namecmp EQ name $
        || namecmp EQ class || namecmp EQ "EMPTY") THEN BEGIN
      element_block[i] = 1
      IF (namecmp NE "EMPTY") THEN RETURN, 1
    ENDIF
  ENDFOR
  RETURN, 0
END

; --------------------------------------------------------------------------

COMPILE_OPT idl2, hidden
COMMON SDF_Common_data, SDF_Common, SDF_Blocktypes, SDF_Blocktype_names, $
    SDF_Datatypes, SDF_Error

SDF_Common = { $
    ID_LENGTH:32L, $
    ENDIANNESS:16911887L, $
    VERSION:1L, $
    REVISION:1L, $
    MAGIC:"SDF1" }

SDF_Blocktypes = { $
    SCRUBBED:-1L, $
    NULL:0L, $
    PLAIN_MESH:1L, $
    POINT_MESH:2L, $
    PLAIN_VARIABLE:3L, $
    POINT_VARIABLE:4L, $
    CONSTANT:5L, $
    ARRAY:6L, $
    RUN_INFO:7L, $
    SOURCE:8L, $
    STITCHED_TENSOR:9L, $
    STITCHED_MATERIAL:10L, $
    STITCHED_MATVAR:11L, $
    STITCHED_SPECIES:12L, $
    FAMILY:13L }

SDF_Blocktype_names = [ $
    "Invalid block", $
    "Plain mesh", $
    "Point mesh", $
    "Plain variable", $
    "Point variable", $
    "Constant", $
    "Simple array", $
    "Run information", $
    "Source code", $
    "Stitched tensor", $
    "Stitched material", $
    "Stitched material variable", $
    "Stitched species", $
    "Particle family" ]

SDF_Datatypes = { $
    NULL:0L, $
    INTEGER4:1L, $
    INTEGER8:2L, $
    REAL4:3L, $
    REAL8:4L, $
    REAL16:5L, $
    CHARACTER:6L, $
    LOGICAL:7L, $
    OTHER:8L }

SDF_Error = { $
    NONE:0L, $
    BAD_USAGE:1L, $
    BAD_MAGIC:2L, $
    BAD_VERSION:3L, $
    BAD_BLOCK_COUNT:4L }

END
