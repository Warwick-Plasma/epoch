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

  FUNCTION handle_ic_fields_deck(element,value)
    CHARACTER(*),INTENT(IN) :: element,value
    INTEGER :: handle_ic_fields_deck

    handle_ic_fields_deck=ERR_NONE
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element,"ex")) THEN
       CALL evaluate_string_in_space(value,ex(-1:nx+2),(/-1,nx+2/),handle_ic_fields_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"ey")) THEN
       CALL evaluate_string_in_space(value,ey(-1:nx+2),(/-1,nx+2/),handle_ic_fields_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"ez")) THEN
       CALL evaluate_string_in_space(value,ez(-1:nx+2),(/-1,nx+2/),handle_ic_fields_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"bx")) THEN
       CALL evaluate_string_in_space(value,bx(-1:nx+2),(/-1,nx+2/),handle_ic_fields_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"by")) THEN
       CALL evaluate_string_in_space(value,by(-1:nx+2),(/-1,nx+2/),handle_ic_fields_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"bz")) THEN
       CALL evaluate_string_in_space(value,bz(-1:nx+2),(/-1,nx+2/),handle_ic_fields_deck)
       RETURN
    ENDIF

  END FUNCTION handle_ic_fields_deck

  FUNCTION check_ic_fields_block()

    INTEGER :: check_ic_fields_block

    !Should do error checking but can't be bothered at the moment
    check_ic_fields_block=ERR_NONE

  END FUNCTION check_ic_fields_block

END MODULE deck_ic_fields_block
