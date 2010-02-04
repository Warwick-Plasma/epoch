MODULE deck_control_block

  USE shared_data
  USE strings
  USE strings_advanced

  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: control_block_elements = 13
  LOGICAL, DIMENSION(control_block_elements) :: control_block_done = .FALSE.
  CHARACTER(LEN=string_length), DIMENSION(control_block_elements) :: control_block_name = (/"nx", "npart", "nsteps", "t_end", "x_start", "x_end", "dt_multiplier", "dlb", "dlb_threshold", "initial_conditions", "icfile", "restart_snapshot", "neutral_background"/)

CONTAINS

  FUNCTION handle_control_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_control_deck
    INTEGER :: loop, elementselected
    handle_control_deck = c_err_unknown_element

    elementselected = 0

    DO loop = 1, control_block_elements
      IF(str_cmp(element, TRIM(ADJUSTL(control_block_name(loop))))) THEN
        elementselected = loop
        EXIT
      ENDIF
    ENDDO

    IF (elementselected .EQ. 0) RETURN
    IF (control_block_done(elementselected)) THEN
      handle_control_deck = c_err_preset_element
      RETURN
    ENDIF
    control_block_done(elementselected) = .TRUE.
    handle_control_deck = c_err_none

    SELECT CASE (elementselected)
    CASE(1)
      nx_global = as_integer(value, handle_control_deck)
    CASE(2)
      npart_global = as_long_integer(value, handle_control_deck)
    CASE(3)
      nsteps = as_integer(value, handle_control_deck)
    CASE(4)
      t_end = as_real(value, handle_control_deck)
    CASE(5)
      x_start = as_real(value, handle_control_deck)
    CASE(6)
      x_end = as_real(value, handle_control_deck)
    CASE(7)
      dt_multiplier = as_real(value, handle_control_deck)
    CASE(8)
      dlb = as_logical(value, handle_control_deck)
    CASE(9)
      dlb_threshold = as_real(value, handle_control_deck)
    CASE(10)
      ictype = as_integer(value, handle_control_deck)
    CASE(11)
      icfile%value = value(1:MIN(LEN(value), data_dir_max_length))
    CASE(12)
      restart_snapshot = as_integer(value, handle_control_deck)
    CASE(13)
      neutral_background = as_logical(value, handle_control_deck)
    END SELECT

  END FUNCTION handle_control_deck



  FUNCTION check_control_block()

    INTEGER :: check_control_block, index

    check_control_block = c_err_none

    ! npart is not a required variable
    control_block_done(2) = .TRUE.

    ! If not using external load then don't need a file
    IF (IAND(ictype, c_ic_external) .EQ. 0) control_block_done(11) = .TRUE.
    ! If not using restart then don't need a restart number
    IF (IAND(ictype, c_ic_restart)  .EQ. 0) control_block_done(12) = .TRUE.

    ! The neutral background is still beta, so hide it if people don't want it
    control_block_done(13) = .TRUE.

    restart = IAND(ictype, c_ic_restart) .NE. 0

    DO index = 1, control_block_elements
      IF (.NOT. control_block_done(index)) THEN
        IF (rank .EQ. 0) THEN
          PRINT *, "***ERROR***"
          PRINT *, "Required control block element ", TRIM(ADJUSTL(control_block_name(index))), " absent. Please create this entry in the input deck"
          WRITE(40, *) ""
          WRITE(40, *) "***ERROR***"
          WRITE(40, *) "Required control block element ", TRIM(ADJUSTL(control_block_name(index))), " absent. Please create this entry in the input deck"
        ENDIF
        check_control_block = c_err_missing_elements
      ENDIF
    ENDDO

  END FUNCTION check_control_block

END MODULE deck_control_block
