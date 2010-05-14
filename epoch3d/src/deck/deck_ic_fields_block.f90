MODULE deck_ic_fields_block

  USE strings_advanced
  USE simple_io

  IMPLICIT NONE

  SAVE
  INTEGER(KIND=MPI_OFFSET_KIND) :: offset

CONTAINS

  SUBROUTINE fields_start

    offset = 0

  END SUBROUTINE fields_start



  FUNCTION handle_ic_fields_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_ic_fields_deck
    CHARACTER(LEN=string_length) :: filename
    INTEGER :: err
    LOGICAL :: got_file

    handle_ic_fields_deck = c_err_none
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element, "offset")) THEN
      offset = as_long_integer_simple(value, handle_ic_fields_deck)
      RETURN
    ENDIF

    CALL get_filename(value, filename, got_file, err)

    IF (str_cmp(element, "ex")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, ex, offset, &
            handle_ic_fields_deck)
      ELSE
        CALL evaluate_string_in_space(value, ex(-1:nx+2, -1:ny+2, -1:nz+2), &
            (/-1, nx+2/), (/-1, ny+2/), (/-1, nz+2/), handle_ic_fields_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "ey")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, ey, offset, &
            handle_ic_fields_deck)
      ELSE
        CALL evaluate_string_in_space(value, ey(-1:nx+2, -1:ny+2, -1:nz+2), &
            (/-1, nx+2/), (/-1, ny+2/), (/-1, nz+2/), handle_ic_fields_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "ez")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, ez, offset, &
            handle_ic_fields_deck)
      ELSE
        CALL evaluate_string_in_space(value, ez(-1:nx+2, -1:ny+2, -1:nz+2), &
            (/-1, nx+2/), (/-1, ny+2/), (/-1, nz+2/), handle_ic_fields_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "bx")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, bx, offset, &
            handle_ic_fields_deck)
      ELSE
        CALL evaluate_string_in_space(value, bx(-1:nx+2, -1:ny+2, -1:nz+2), &
            (/-1, nx+2/), (/-1, ny+2/), (/-1, nz+2/), handle_ic_fields_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "by")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, by, offset, &
            handle_ic_fields_deck)
      ELSE
        CALL evaluate_string_in_space(value, by(-1:nx+2, -1:ny+2, -1:nz+2), &
            (/-1, nx+2/), (/-1, ny+2/), (/-1, nz+2/), handle_ic_fields_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "bz")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, bz, offset, &
            handle_ic_fields_deck)
      ELSE
        CALL evaluate_string_in_space(value, bz(-1:nx+2, -1:ny+2, -1:nz+2), &
            (/-1, nx+2/), (/-1, ny+2/), (/-1, nz+2/), handle_ic_fields_deck)
      ENDIF
      RETURN
    ENDIF

  END FUNCTION handle_ic_fields_deck



  FUNCTION check_ic_fields_block()

    INTEGER :: check_ic_fields_block

    ! Should do error checking but can't be bothered at the moment
    check_ic_fields_block = c_err_none

  END FUNCTION check_ic_fields_block

END MODULE deck_ic_fields_block
