MODULE deck_ic_laser_block

  USE strings_advanced
  USE laser

  IMPLICIT NONE

  SAVE

  TYPE(laser_block), POINTER :: working_laser
  LOGICAL :: boundary_set = .FALSE.
  INTEGER :: boundary

CONTAINS

  FUNCTION handle_ic_laser_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_ic_laser_deck
    REAL(num) :: dummy
    TYPE(primitive_stack) :: output
    INTEGER :: ix, iy, ierr

    handle_ic_laser_deck = c_err_none
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element, "boundary") .OR. str_cmp(element, "direction")) THEN
      IF (str_cmp(element, "direction") .AND. rank .EQ. 0) THEN
        WRITE(*, *) '*** WARNING ***'
        WRITE(*, *) 'Element "direction" in the block "laser" is deprecated.'
        WRITE(*, *) 'Please use the element name "boundary" instead.'
      ENDIF
      ! If the boundary has already been set, simply ignore further calls to it
      IF (boundary_set) RETURN
      boundary = as_boundary(value, handle_ic_laser_deck)
      boundary_set = .TRUE.
      CALL init_laser(boundary, working_laser)
      RETURN
    ENDIF

    IF (.NOT. boundary_set) THEN
      IF (rank .EQ. 0) THEN
        WRITE(*, *) '*** ERROR ***'
        WRITE(*, *) 'Cannot set laser properties before boundary is set'
        WRITE(40,*) '*** ERROR ***'
        WRITE(40,*) 'Cannot set laser properties before boundary is set'
        CALL MPI_ABORT(comm, errcode, ierr)
      ENDIF
      extended_error_string = "boundary"
      handle_ic_laser_deck = c_err_required_element_not_set
      RETURN
    ENDIF

    IF (str_cmp(element, "amp")) THEN
      working_laser%amp = as_real(value, handle_ic_laser_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "irradiance")) THEN
      working_laser%amp = SQRT(as_real(value, handle_ic_laser_deck) &
          / (c*epsilon0/2.0_num))
      RETURN
    ENDIF

    IF (str_cmp(element, "omega") .OR. str_cmp(element, "freq")) THEN
      IF (rank .EQ. 0 .AND. str_cmp(element, "freq")) THEN
        WRITE(*,*) '*** WARNING ***'
        WRITE(*,*) 'Element "freq" in the block "laser" is deprecated.'
        WRITE(*,*) 'Please use the element name "omega" instead.'
      ENDIF
      working_laser%omega = as_real(value, handle_ic_laser_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "frequency")) THEN
      working_laser%omega = 2.0_num * pi * as_real(value, handle_ic_laser_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "lambda")) THEN
      working_laser%omega = &
          2.0_num * pi * c / as_real(value, handle_ic_laser_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "profile")) THEN
      working_laser%profile = 0.0_num
      CALL initialise_stack(output)
      CALL tokenize(value, output, handle_ic_laser_deck)
      IF (working_laser%boundary .EQ. c_bd_x_min &
          .OR. working_laser%boundary .EQ. c_bd_x_max) THEN
        DO iy = 1, ny
          working_laser%profile(iy) = &
              evaluate_at_point(output, 0, iy, handle_ic_laser_deck)
        ENDDO
      ELSE IF (working_laser%boundary .EQ. c_bd_y_max &
          .OR. working_laser%boundary .EQ. c_bd_y_min) THEN
        DO ix = 1, nx
          working_laser%profile(ix) = &
              evaluate_at_point(output, ix, 0, handle_ic_laser_deck)
        ENDDO
      ENDIF
      CALL deallocate_stack(output)
      RETURN
    ENDIF

    IF (str_cmp(element, "phase")) THEN
      working_laser%phase = 0.0_num
      CALL initialise_stack(output)
      CALL tokenize(value, output, handle_ic_laser_deck)
      IF (working_laser%boundary .EQ. c_bd_x_min &
          .OR. working_laser%boundary .EQ. c_bd_x_max) THEN
        DO iy = 1, ny
          working_laser%phase(iy) = &
              evaluate_at_point(output, 0, iy, handle_ic_laser_deck)
        ENDDO
      ELSE IF (working_laser%boundary .EQ. c_bd_y_max &
          .OR. working_laser%boundary .EQ. c_bd_y_min) THEN
        DO ix = 1, nx
          working_laser%phase(ix) = &
              evaluate_at_point(output, ix, 0, handle_ic_laser_deck)
        ENDDO
      ENDIF
      CALL deallocate_stack(output)
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
      CALL initialise_stack(working_laser%time_function)
      CALL tokenize(value, working_laser%time_function, handle_ic_laser_deck)
      ! evaluate it once to check that it's a valid block
      dummy = evaluate(working_laser%time_function, handle_ic_laser_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "pol_angle")) THEN
      working_laser%pol_angle = as_real(value, handle_ic_laser_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "pol")) THEN
      ! Convert from degrees to radians
      working_laser%pol_angle = &
          pi * as_real(value, handle_ic_laser_deck) / 180.0_num
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
    TYPE(laser_block), POINTER :: current
    INTEGER :: error

    check_ic_laser_block = c_err_none

    error = 0
    current=>laser_x_min
    DO WHILE(ASSOCIATED(current))
      IF (current%omega .LT. 0.0_num) error = IOR(error, 1)
      IF (current%amp .LT. 0.0_num) error = IOR(error, 2)
      current=>current%next
    ENDDO

    current=>laser_x_max
    DO WHILE(ASSOCIATED(current))
      IF (current%omega .LT. 0.0_num) error = IOR(error, 1)
      IF (current%amp .LT. 0.0_num) error = IOR(error, 2)
      current=>current%next
    ENDDO

    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      IF (current%omega .LT. 0.0_num) error = IOR(error, 1)
      IF (current%amp .LT. 0.0_num) error = IOR(error, 2)
      current=>current%next
    ENDDO

    current=>laser_y_max
    DO WHILE(ASSOCIATED(current))
      IF (current%omega .LT. 0.0_num) error = IOR(error, 1)
      IF (current%amp .LT. 0.0_num) error = IOR(error, 2)
      current=>current%next
    ENDDO

    IF (IAND(error,1) .NE. 0) THEN
      IF (rank .EQ. 0) THEN
        WRITE(*, *) '*** ERROR ***'
        WRITE(*, *) 'Must define a "lambda" or "omega" for every laser.'
        WRITE(40,*) '*** ERROR ***'
        WRITE(40,*) 'Must define a "lambda" or "omega" for every laser.'
      ENDIF
      check_ic_laser_block = c_err_missing_elements
    ENDIF

    IF (IAND(error,2) .NE. 0) THEN
      IF (rank .EQ. 0) THEN
        WRITE(*, *) '*** ERROR ***'
        WRITE(*, *) 'Must define an "amp" or "irradiance" for every laser.'
        WRITE(40,*) '*** ERROR ***'
        WRITE(40,*) 'Must define an "amp" or "irradiance" for every laser.'
      ENDIF
      check_ic_laser_block = c_err_missing_elements
    ENDIF

  END FUNCTION check_ic_laser_block



  SUBROUTINE laser_start

    ! Every new laser uses the internal time function
    ALLOCATE(working_laser)
    working_laser%use_time_function = .FALSE.

  END SUBROUTINE laser_start



  SUBROUTINE laser_end

    CALL attach_laser(working_laser)
    boundary_set = .FALSE.

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
