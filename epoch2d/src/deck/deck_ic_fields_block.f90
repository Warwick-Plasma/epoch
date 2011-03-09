MODULE deck_fields_block

  USE strings_advanced
  USE simple_io

  IMPLICIT NONE

  SAVE
  INTEGER(KIND=MPI_OFFSET_KIND) :: offset

CONTAINS

  SUBROUTINE fields_start

    offset = 0

  END SUBROUTINE fields_start



  FUNCTION handle_fields_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_fields_deck
    CHARACTER(LEN=string_length) :: filename
    INTEGER :: err
    LOGICAL :: got_file

    handle_fields_deck = c_err_none
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element, "offset")) THEN
      offset = as_long_integer_simple(value, handle_fields_deck)
      RETURN
    ENDIF

    CALL get_filename(value, filename, got_file, err)

    IF (str_cmp(element, "ex")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, ex, offset, &
            handle_fields_deck)
      ELSE
        CALL evaluate_string_in_space(value, ex, &
            -2, nx+3, -2, ny+3, handle_fields_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "ey")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, ey, offset, &
            handle_fields_deck)
      ELSE
        CALL evaluate_string_in_space(value, ey, &
            -2, nx+3, -2, ny+3, handle_fields_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "ez")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, ez, offset, &
            handle_fields_deck)
      ELSE
        CALL evaluate_string_in_space(value, ez, &
            -2, nx+3, -2, ny+3, handle_fields_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "bx")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, bx, offset, &
            handle_fields_deck)
      ELSE
        CALL evaluate_string_in_space(value, bx, &
            -2, nx+3, -2, ny+3, handle_fields_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "by")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, by, offset, &
            handle_fields_deck)
      ELSE
        CALL evaluate_string_in_space(value, by, &
            -2, nx+3, -2, ny+3, handle_fields_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "bz")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, bz, offset, &
            handle_fields_deck)
      ELSE
        CALL evaluate_string_in_space(value, bz, &
            -2, nx+3, -2, ny+3, handle_fields_deck)
      ENDIF
      RETURN
    ENDIF

  END FUNCTION handle_fields_deck



  FUNCTION check_fields_block()

    INTEGER :: check_fields_block

    ! Should do error checking but can't be bothered at the moment
    check_fields_block = c_err_none

  END FUNCTION check_fields_block

END MODULE deck_fields_block
