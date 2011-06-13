MODULE deck_control_block

  USE strings_advanced
  USE fields

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: control_deck_initialise, control_deck_finalise
  PUBLIC :: control_block_start, control_block_end
  PUBLIC :: control_block_handle_element, control_block_check

  INTEGER, PARAMETER :: control_block_elements = 12 + 4 * c_ndims
  LOGICAL, DIMENSION(control_block_elements) :: control_block_done
  CHARACTER(LEN=string_length), DIMENSION(control_block_elements) :: &
      control_block_name = (/ &
          "nx                ", &
          "ny                ", &
          "x_min             ", &
          "x_max             ", &
          "y_min             ", &
          "y_max             ", &
          "nprocx            ", &
          "nprocy            ", &
          "npart             ", &
          "nsteps            ", &
          "t_end             ", &
          "dt_multiplier     ", &
          "dlb_threshold     ", &
          "icfile            ", &
          "restart_snapshot  ", &
          "neutral_background", &
          "field_order       ", &
          "stdout_frequency  ", &
          "use_random_seed   ", &
          "smooth_currents   " /)
  CHARACTER(LEN=string_length), DIMENSION(control_block_elements) :: &
      alternate_name = (/ &
          "nx                ", &
          "ny                ", &
          "x_start           ", &
          "x_end             ", &
          "y_start           ", &
          "y_end             ", &
          "nprocx            ", &
          "nprocy            ", &
          "npart             ", &
          "nsteps            ", &
          "t_end             ", &
          "dt_multiplier     ", &
          "dlb_threshold     ", &
          "icfile            ", &
          "restart_snapshot  ", &
          "neutral_background", &
          "field_order       ", &
          "stdout_frequency  ", &
          "use_random_seed   ", &
          "smooth_currents   " /)

CONTAINS

  SUBROUTINE control_deck_initialise

    IF (deck_state .NE. c_ds_first) RETURN
    control_block_done = .FALSE.

  END SUBROUTINE control_deck_initialise



  SUBROUTINE control_deck_finalise

  END SUBROUTINE control_deck_finalise



  SUBROUTINE control_block_start

  END SUBROUTINE control_block_start



  SUBROUTINE control_block_end

  END SUBROUTINE control_block_end



  FUNCTION control_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    INTEGER :: loop, elementselected, field_order, ierr, io

    errcode = c_err_none
    IF (deck_state .NE. c_ds_first) RETURN

    errcode = c_err_unknown_element

    elementselected = 0

    DO loop = 1, control_block_elements
      IF (str_cmp(element, TRIM(ADJUSTL(control_block_name(loop)))) &
          .OR. str_cmp(element, TRIM(ADJUSTL(alternate_name(loop))))) THEN
        elementselected = loop
        EXIT
      ENDIF
    ENDDO

    IF (elementselected .EQ. 0) RETURN
    IF (control_block_done(elementselected)) THEN
      errcode = c_err_preset_element
      RETURN
    ENDIF
    control_block_done(elementselected) = .TRUE.
    errcode = c_err_none

    SELECT CASE (elementselected)
    CASE(1)
      nx_global = as_integer(value, errcode)
    CASE(2)
      ny_global = as_integer(value, errcode)
    CASE(c_ndims+1)
      x_min = as_real(value, errcode)
    CASE(c_ndims+2)
      x_max = as_real(value, errcode)
    CASE(c_ndims+3)
      y_min = as_real(value, errcode)
    CASE(c_ndims+4)
      y_max = as_real(value, errcode)
    CASE(3*c_ndims+1)
      nprocx = as_integer(value, errcode)
    CASE(3*c_ndims+2)
      nprocy = as_integer(value, errcode)
    CASE(4*c_ndims+1)
      npart_global = as_long_integer(value, errcode)
    CASE(4*c_ndims+2)
      nsteps = as_integer(value, errcode)
    CASE(4*c_ndims+3)
      t_end = as_real(value, errcode)
    CASE(4*c_ndims+4)
      dt_multiplier = as_real(value, errcode)
    CASE(4*c_ndims+5)
      dlb_threshold = as_real(value, errcode)
      dlb = .TRUE.
    CASE(4*c_ndims+6)
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The "icfile" option is no longer supported.'
          WRITE(io,*) 'Please use the "import" directive instead'
        ENDDO
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    CASE(4*c_ndims+7)
      restart_snapshot = as_integer(value, errcode)
      ic_from_restart = .TRUE.
    CASE(4*c_ndims+8)
      neutral_background = as_logical(value, errcode)
    CASE(4*c_ndims+9)
      field_order = as_integer(value, errcode)
      IF (field_order .NE. 2 .AND. field_order .NE. 4 &
          .AND. field_order .NE. 6) THEN
        errcode = c_err_bad_value
      ELSE
        CALL set_field_order(field_order)
      ENDIF
    CASE(4*c_ndims+10)
      stdout_frequency = as_integer(value, errcode)
    CASE(4*c_ndims+11)
      use_random_seed = as_logical(value, errcode)
    CASE(4*c_ndims+12)
      smooth_currents = as_logical(value, errcode)
    END SELECT

  END FUNCTION control_block_handle_element



  FUNCTION control_block_check() RESULT(errcode)

    INTEGER :: errcode, index, io

    errcode = c_err_none

    ! nprocx/y and npart are optional
    control_block_done(3*c_ndims+1:4*c_ndims+1) = .TRUE.

    ! Only one of nsteps or t_end need be specified
    IF (control_block_done(4*c_ndims+2)) &
        control_block_done(4*c_ndims+3) = .TRUE.
    IF (control_block_done(4*c_ndims+3)) &
        control_block_done(4*c_ndims+2) = .TRUE.

    ! All entries after t_end are optional
    control_block_done(4*c_ndims+4:) = .TRUE.

    DO index = 1, control_block_elements
      IF (.NOT. control_block_done(index)) THEN
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Required control block element ', &
                TRIM(ADJUSTL(control_block_name(index))), &
                ' absent. Please create this entry in the input deck'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF
    ENDDO

    IF (.NOT. neutral_background) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The option "neutral_background=F" is not supported', &
              ' in this version of EPOCH.'
        ENDDO
      ENDIF
      errcode = c_err_terminate
    ENDIF

  END FUNCTION control_block_check

END MODULE deck_control_block
