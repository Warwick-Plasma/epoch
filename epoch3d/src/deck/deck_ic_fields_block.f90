MODULE deck_ic_fields_block

  USE shared_data
  USE strings_advanced
  USE shared_parser_data
  USE strings
  USE shunt
  USE evaluator

  IMPLICIT NONE

  SAVE

CONTAINS

  FUNCTION HandleICFieldsDeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER :: HandleICFieldsDeck

    HandleICFieldsDeck=ERR_NONE
    IF (Element .EQ. blank .OR. Value .EQ. blank) RETURN

    IF (StrCmp(Element,"ex")) THEN
       CALL EvaluateStringInSpace(Value,Ex(-1:nx+2,-1:ny+2,-1:nz+2),(/-1,nx+2/),(/-1,ny+2/),(/-1,nz+2/),HandleICFieldsDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"ey")) THEN
       CALL EvaluateStringInSpace(Value,Ex(-1:nx+2,-1:ny+2,-1:nz+2),(/-1,nx+2/),(/-1,ny+2/),(/-1,nz+2/),HandleICFieldsDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"ez")) THEN
       CALL EvaluateStringInSpace(Value,Ex(-1:nx+2,-1:ny+2,-1:nz+2),(/-1,nx+2/),(/-1,ny+2/),(/-1,nz+2/),HandleICFieldsDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"bx")) THEN
       CALL EvaluateStringInSpace(Value,Ex(-1:nx+2,-1:ny+2,-1:nz+2),(/-1,nx+2/),(/-1,ny+2/),(/-1,nz+2/),HandleICFieldsDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"by")) THEN
       CALL EvaluateStringInSpace(Value,Ex(-1:nx+2,-1:ny+2,-1:nz+2),(/-1,nx+2/),(/-1,ny+2/),(/-1,nz+2/),HandleICFieldsDeck)
    ENDIF

    IF (StrCmp(Element,"bz")) THEN
       CALL EvaluateStringInSpace(Value,Ex(-1:nx+2,-1:ny+2,-1:nz+2),(/-1,nx+2/),(/-1,ny+2/),(/-1,nz+2/),HandleICFieldsDeck)
    ENDIF

  END FUNCTION HandleICFieldsDeck

  FUNCTION CheckICFieldsBlock()

    INTEGER :: CheckICFieldsBlock

    !Should do error checking but can't be bothered at the moment
    CheckICFieldsBlock=ERR_NONE

  END FUNCTION CheckICFieldsBlock

END MODULE deck_ic_fields_block
