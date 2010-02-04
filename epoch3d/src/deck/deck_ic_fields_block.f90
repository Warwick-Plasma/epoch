MODULE deck_ic_fields_block

  USE strings_advanced

  IMPLICIT NONE

  SAVE

CONTAINS

  FUNCTION handle_ic_fields_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_ic_fields_deck

    handle_ic_fields_deck = c_err_none
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element, "ex")) THEN
      CALL evaluate_string_in_space(value, ex(-1:nx+2, -1:ny+2, -1:nz+2), &
          (/-1, nx+2/), (/-1, ny+2/), (/-1, nz+2/), handle_ic_fields_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "ey")) THEN
      CALL evaluate_string_in_space(value, ex(-1:nx+2, -1:ny+2, -1:nz+2), &
          (/-1, nx+2/), (/-1, ny+2/), (/-1, nz+2/), handle_ic_fields_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "ez")) THEN
      CALL evaluate_string_in_space(value, ex(-1:nx+2, -1:ny+2, -1:nz+2), &
          (/-1, nx+2/), (/-1, ny+2/), (/-1, nz+2/), handle_ic_fields_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "bx")) THEN
      CALL evaluate_string_in_space(value, ex(-1:nx+2, -1:ny+2, -1:nz+2), &
          (/-1, nx+2/), (/-1, ny+2/), (/-1, nz+2/), handle_ic_fields_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "by")) THEN
      CALL evaluate_string_in_space(value, ex(-1:nx+2, -1:ny+2, -1:nz+2), &
          (/-1, nx+2/), (/-1, ny+2/), (/-1, nz+2/), handle_ic_fields_deck)
    ENDIF

    IF (str_cmp(element, "bz")) THEN
      CALL evaluate_string_in_space(value, ex(-1:nx+2, -1:ny+2, -1:nz+2), &
          (/-1, nx+2/), (/-1, ny+2/), (/-1, nz+2/), handle_ic_fields_deck)
    ENDIF

  END FUNCTION handle_ic_fields_deck



  FUNCTION check_ic_fields_block()

    INTEGER :: check_ic_fields_block

    ! Should do error checking but can't be bothered at the moment
    check_ic_fields_block = c_err_none

  END FUNCTION check_ic_fields_block

END MODULE deck_ic_fields_block
