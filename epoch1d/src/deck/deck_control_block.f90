MODULE deck_control_block

  USE strings_advanced
  USE fields

  IMPLICIT NONE

  SAVE
  INTEGER, PARAMETER :: control_block_elements = 12
  LOGICAL, DIMENSION(control_block_elements) :: control_block_done = .FALSE.
  CHARACTER(LEN=string_length), DIMENSION(control_block_elements) :: &
      control_block_name = (/ &
          "nx                ", &
          "npart             ", &
          "nsteps            ", &
          "t_end             ", &
          "x_start           ", &
          "x_end             ", &
          "dt_multiplier     ", &
          "dlb_threshold     ", &
          "icfile            ", &
          "restart_snapshot  ", &
          "neutral_background", &
          "field_order       " /)

CONTAINS

  FUNCTION handle_control_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_control_deck
    INTEGER :: loop, elementselected, field_order

    handle_control_deck = c_err_unknown_element

    elementselected = 0

    DO loop = 1, control_block_elements
      IF (str_cmp(element, TRIM(ADJUSTL(control_block_name(loop))))) THEN
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
      dlb_threshold = as_real(value, handle_control_deck)
      dlb = .TRUE.
    CASE(9)
      icfile%value = value(1:MIN(LEN(value), data_dir_max_length))
    CASE(10)
      restart_snapshot = as_integer(value, handle_control_deck)
      ic_from_restart = .TRUE.
    CASE(11)
      neutral_background = as_logical(value, handle_control_deck)
    CASE(12)
      field_order = as_integer(value, handle_control_deck)
      IF (field_order .NE. 2 .AND. field_order .NE. 4 &
          .AND. field_order .NE. 6) THEN
        handle_control_deck = c_err_bad_value
      ELSE
        CALL set_field_order(field_order)
      ENDIF
    END SELECT

  END FUNCTION handle_control_deck



  FUNCTION check_control_block()

    INTEGER :: check_control_block, index

    check_control_block = c_err_none

    ! npart is not a required variable
    control_block_done(2) = .TRUE.

    ! dlb threshold is optional
    control_block_done(8) = .TRUE.

    ! external input deck is optional
    control_block_done(9) = .TRUE.

    ! restart snapshot is optional
    control_block_done(10) = .TRUE.

    ! The neutral background is still beta, so hide it if people don't want it
    control_block_done(11) = .TRUE.

    ! field_order is optional
    control_block_done(12) = .TRUE.

    DO index = 1, control_block_elements
      IF (.NOT. control_block_done(index)) THEN
        IF (rank .EQ. 0) THEN
          PRINT *, "***ERROR***"
          PRINT *, "Required control block element ", &
              TRIM(ADJUSTL(control_block_name(index))), &
              " absent. Please create this entry in the input deck"
          WRITE(40, *) ""
          WRITE(40, *) "***ERROR***"
          WRITE(40, *) "Required control block element ", &
              TRIM(ADJUSTL(control_block_name(index))), &
              " absent. Please create this entry in the input deck"
        ENDIF
        check_control_block = c_err_missing_elements
      ENDIF
    ENDDO

  END FUNCTION check_control_block

END MODULE deck_control_block
