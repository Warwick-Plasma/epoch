MODULE deck_laser_block

  USE strings_advanced
  USE laser

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: laser_deck_initialise, laser_deck_finalise
  PUBLIC :: laser_block_start, laser_block_end
  PUBLIC :: laser_block_handle_element, laser_block_check

  TYPE(laser_block), POINTER :: working_laser
  LOGICAL :: boundary_set = .FALSE.
  INTEGER :: boundary

CONTAINS

  SUBROUTINE laser_deck_initialise

  END SUBROUTINE laser_deck_initialise



  SUBROUTINE laser_deck_finalise

  END SUBROUTINE laser_deck_finalise



  SUBROUTINE laser_block_start

    IF (deck_state .EQ. c_ds_first) RETURN

    ! Every new laser uses the internal time function
    ALLOCATE(working_laser)
    working_laser%use_time_function = .FALSE.

  END SUBROUTINE laser_block_start



  SUBROUTINE laser_block_end

    IF (deck_state .EQ. c_ds_first) RETURN

    CALL attach_laser(working_laser)
    boundary_set = .FALSE.

  END SUBROUTINE laser_block_end



  FUNCTION laser_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    REAL(num) :: dummy
    INTEGER :: ierr, io

    errcode = c_err_none
    IF (deck_state .EQ. c_ds_first) RETURN
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element, "boundary") .OR. str_cmp(element, "direction")) THEN
      IF (rank .EQ. 0 .AND. str_cmp(element, "direction")) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Element "direction" in the block "laser" is deprecated.'
          WRITE(io,*) 'Please use the element name "boundary" instead.'
        ENDDO
      ENDIF
      ! If the boundary has already been set, simply ignore further calls to it
      IF (boundary_set) RETURN
      boundary = as_boundary(value, errcode)
      boundary_set = .TRUE.
      CALL init_laser(boundary, working_laser)
      RETURN
    ENDIF

    IF (.NOT. boundary_set) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Cannot set laser properties before boundary is set'
        ENDDO
        CALL MPI_ABORT(comm, errcode, ierr)
      ENDIF
      extended_error_string = "boundary"
      errcode = c_err_required_element_not_set
      RETURN
    ENDIF

    IF (str_cmp(element, "amp")) THEN
      working_laser%amp = as_real(value, errcode)
      RETURN
    ENDIF

    ! SI (W/m^2)
    IF (str_cmp(element, "irradiance") .OR. str_cmp(element, "intensity")) THEN
      working_laser%amp = SQRT(as_real(value, errcode) &
          / (c*epsilon0/2.0_num))
      RETURN
    ENDIF

    IF (str_cmp(element, "irradiance_w_cm2") &
        .OR. str_cmp(element, "intensity_w_cm2")) THEN
      working_laser%amp = SQRT(as_real(value, errcode) &
          / (c*epsilon0/2.0_num)) * 100_num 
      RETURN
    ENDIF

    IF (str_cmp(element, "omega") .OR. str_cmp(element, "freq")) THEN
      IF (rank .EQ. 0 .AND. str_cmp(element, "freq")) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Element "freq" in the block "laser" is deprecated.'
          WRITE(io,*) 'Please use the element name "omega" instead.'
        ENDDO
      ENDIF
      working_laser%omega = as_real(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "frequency")) THEN
      working_laser%omega = 2.0_num * pi * as_real(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "lambda")) THEN
      working_laser%omega = 2.0_num * pi * c / as_real(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "profile")) THEN
      working_laser%profile = as_real(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "phase")) THEN
      working_laser%phase = as_real(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "t_start")) THEN
      working_laser%t_start = as_time(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "t_end")) THEN
      working_laser%t_end = as_time(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "t_profile")) THEN
      working_laser%use_time_function = .TRUE.
      CALL initialise_stack(working_laser%time_function)
      CALL tokenize(value, working_laser%time_function, errcode)
      ! evaluate it once to check that it's a valid block
      dummy = evaluate(working_laser%time_function, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "pol_angle")) THEN
      working_laser%pol_angle = as_real(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "pol")) THEN
      ! Convert from degrees to radians
      working_laser%pol_angle = pi * as_real(value, errcode) / 180.0_num
      RETURN
    ENDIF

    IF (str_cmp(element, "id")) THEN
      working_laser%id = as_integer(value, errcode)
      RETURN
    ENDIF

    errcode = c_err_unknown_element

  END FUNCTION laser_block_handle_element



  FUNCTION laser_block_check() RESULT(errcode)

    INTEGER :: errcode
    TYPE(laser_block), POINTER :: current
    INTEGER :: error, io

    errcode = c_err_none

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

    IF (IAND(error,1) .NE. 0) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "lambda" or "omega" for every laser.'
        ENDDO
      ENDIF
      errcode = c_err_missing_elements
    ENDIF

    IF (IAND(error,2) .NE. 0) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define an "amp" or "irradiance" for every laser.'
        ENDDO
      ENDIF
      errcode = c_err_missing_elements
    ENDIF

  END FUNCTION laser_block_check



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

END MODULE deck_laser_block
