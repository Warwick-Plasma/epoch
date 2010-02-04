MODULE deck_ic_laser_block

  USE strings_advanced
  USE laser

  IMPLICIT NONE

  SAVE

  TYPE(laser_block), POINTER :: working_laser
  LOGICAL :: direction_set = .FALSE.
  INTEGER :: direction

CONTAINS

  FUNCTION handle_ic_laser_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_ic_laser_deck
    REAL(num) :: dummy
    TYPE(primitive_stack) :: output
    INTEGER :: ix, iy

    handle_ic_laser_deck = c_err_none
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element, "direction")) THEN
      ! If the direction has already been set, simply ignore further calls to it
      IF (direction_set) RETURN
      direction = as_direction(value, handle_ic_laser_deck)
      direction_set = .TRUE.
      CALL init_laser(direction, working_laser)
      RETURN
    ENDIF

    IF (.NOT. direction_set) THEN
      IF (rank .EQ. 0) THEN
        PRINT *, "***ERROR*** Cannot set laser properties before &
            &direction is set"
        WRITE(40, *) "***ERROR*** Cannot set laser properties before &
            &direction is set"
      ENDIF
      extended_error_string = "direction"
      handle_ic_laser_deck = c_err_required_element_not_set
      RETURN
    ENDIF

    IF (str_cmp(element, "amp")) THEN
      working_laser%amp = as_real(value, handle_ic_laser_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "freq")) THEN
      working_laser%freq = as_real(value, handle_ic_laser_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "profile")) THEN
      working_laser%profile = 0.0_num
      output%stack_point = 0
      CALL tokenize(value, output, handle_ic_laser_deck)
      IF (working_laser%direction .EQ. c_bd_left .OR. &
          working_laser%direction .EQ. c_bd_right) THEN
        DO iy = 1, ny
          working_laser%profile(iy) = &
              evaluate_at_point(output, 0, iy, handle_ic_laser_deck)
        ENDDO
      ELSE IF (working_laser%direction .EQ. c_bd_up .OR. &
          working_laser%direction .EQ. c_bd_down) THEN
        DO ix = 1, nx
          working_laser%profile(ix) = &
              evaluate_at_point(output, ix, 0, handle_ic_laser_deck)
        ENDDO
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "phase")) THEN
      working_laser%phase = 0.0_num
      output%stack_point = 0
      CALL tokenize(value, output, handle_ic_laser_deck)
      IF (working_laser%direction .EQ. c_bd_left .OR. &
          working_laser%direction .EQ. c_bd_right) THEN
        DO iy = 1, ny
          working_laser%phase(iy) = &
              evaluate_at_point(output, 0, iy, handle_ic_laser_deck)
        ENDDO
      ELSE IF (working_laser%direction .EQ. c_bd_up .OR. &
          working_laser%direction .EQ. c_bd_down) THEN
        DO ix = 1, nx
          working_laser%phase(ix) = &
              evaluate_at_point(output, ix, 0, handle_ic_laser_deck)
        ENDDO
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "t_start")) THEN
      working_laser%t_start = as_time(value, handle_ic_laser_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "t_end")) THEN
      working_laser%t_end = as_time(value, handle_ic_laser_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "t_profile")) THEN
      working_laser%use_time_function = .TRUE.
      working_laser%time_function%stack_point = 0
      CALL tokenize(value, working_laser%time_function, handle_ic_laser_deck)
      ! evaluate it once to check that it's a valid block
      dummy = evaluate(working_laser%time_function, handle_ic_laser_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "pol")) THEN
      ! Convert from degrees to radians
      working_laser%pol = pi * as_real(value, handle_ic_laser_deck)/180.0_num
      RETURN
    ENDIF

    IF (str_cmp(element, "id")) THEN
      working_laser%id = as_integer(value, handle_ic_laser_deck)
      RETURN
    ENDIF

    handle_ic_laser_deck = c_err_unknown_element

  END FUNCTION handle_ic_laser_deck



  FUNCTION check_ic_laser_block()

    INTEGER :: check_ic_laser_block

    ! Should do error checking but can't be bothered at the moment
    check_ic_laser_block = c_err_none

  END FUNCTION check_ic_laser_block



  SUBROUTINE laser_start

    ! Every new laser uses the internal time function
    ALLOCATE(working_laser)
    working_laser%use_time_function = .FALSE.

  END SUBROUTINE laser_start



  SUBROUTINE laser_end

    CALL attach_laser(working_laser)
    direction_set = .FALSE.

  END SUBROUTINE laser_end



  FUNCTION as_time(value, err)

    CHARACTER(LEN=*), INTENT(IN) :: value
    INTEGER, INTENT(INOUT) :: err
    REAL(num) :: as_time

    IF (str_cmp(value, "start")) THEN
      as_time = 0.0_num
      RETURN
    ENDIF

    IF (str_cmp(value, "end")) THEN
      as_time = t_end
      RETURN
    ENDIF

    as_time = as_real(value, err)

  END FUNCTION as_time

END MODULE deck_ic_laser_block
