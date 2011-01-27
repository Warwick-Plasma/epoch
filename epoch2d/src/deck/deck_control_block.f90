MODULE deck_control_block

  USE strings_advanced
  USE fields

  IMPLICIT NONE

  SAVE
  INTEGER, PARAMETER :: control_block_elements = 12 + 4 * c_ndims
  LOGICAL, DIMENSION(control_block_elements) :: control_block_done = .FALSE.
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

  FUNCTION handle_control_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_control_deck
    INTEGER :: loop, elementselected, field_order, ierr

    handle_control_deck = c_err_unknown_element

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
      handle_control_deck = c_err_preset_element
      RETURN
    ENDIF
    control_block_done(elementselected) = .TRUE.
    handle_control_deck = c_err_none

    SELECT CASE (elementselected)
    CASE(1)
      nx_global = as_integer(value, handle_control_deck)
    CASE(2)
      ny_global = as_integer(value, handle_control_deck)
    CASE(c_ndims+1)
      x_min = as_real(value, handle_control_deck)
    CASE(c_ndims+2)
      x_max = as_real(value, handle_control_deck)
    CASE(c_ndims+3)
      y_min = as_real(value, handle_control_deck)
    CASE(c_ndims+4)
      y_max = as_real(value, handle_control_deck)
    CASE(3*c_ndims+1)
      nprocx = as_integer(value, handle_control_deck)
    CASE(3*c_ndims+2)
      nprocy = as_integer(value, handle_control_deck)
    CASE(4*c_ndims+1)
      npart_global = as_long_integer(value, handle_control_deck)
    CASE(4*c_ndims+2)
      nsteps = as_integer(value, handle_control_deck)
    CASE(4*c_ndims+3)
      t_end = as_real(value, handle_control_deck)
    CASE(4*c_ndims+4)
      dt_multiplier = as_real(value, handle_control_deck)
    CASE(4*c_ndims+5)
      dlb_threshold = as_real(value, handle_control_deck)
      dlb = .TRUE.
    CASE(4*c_ndims+6)
      IF (rank .EQ. 0) THEN
        WRITE(*, *) '***ERROR***'
        WRITE(*, *) 'The "icfile" option is no longer supported.'
        WRITE(*, *) 'Please use the "import" directive instead'
        WRITE(40,*) '***ERROR***'
        WRITE(40,*) 'The "icfile" option is no longer supported.'
        WRITE(40,*) 'Please use the "import" directive instead'
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    CASE(4*c_ndims+7)
      restart_snapshot = as_integer(value, handle_control_deck)
      ic_from_restart = .TRUE.
    CASE(4*c_ndims+8)
      neutral_background = as_logical(value, handle_control_deck)
    CASE(4*c_ndims+9)
      field_order = as_integer(value, handle_control_deck)
      IF (field_order .NE. 2 .AND. field_order .NE. 4 &
          .AND. field_order .NE. 6) THEN
        handle_control_deck = c_err_bad_value
      ELSE
        CALL set_field_order(field_order)
      ENDIF
    CASE(4*c_ndims+10)
      stdout_frequency = as_integer(value, handle_control_deck)
    CASE(4*c_ndims+11)
      use_random_seed = as_logical(value, handle_control_deck)
    CASE(4*c_ndims+12)
      smooth_currents = as_logical(value, handle_control_deck)
    END SELECT

  END FUNCTION handle_control_deck



  FUNCTION check_control_block()

    INTEGER :: check_control_block, index

    check_control_block = c_err_none

    ! nprocx/y and npart are optional
    control_block_done(3*c_ndims+1:4*c_ndims+1) = .TRUE.

    ! All entries after t_end are optional
    control_block_done(4*c_ndims+4:) = .TRUE.

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

    IF (.NOT. neutral_background) THEN
      IF (rank .EQ. 0) THEN
        WRITE(*, *) '***ERROR***'
        WRITE(*, *) 'The option "neutral_background=F" is not supported', &
            ' in this version of EPOCH.'
        WRITE(40,*)
        WRITE(40,*) '***ERROR***'
        WRITE(40,*) 'The option "neutral_background=F" is not supported', &
            ' in this version of EPOCH.'
      ENDIF
      check_control_block = c_err_terminate
    ENDIF

  END FUNCTION check_control_block

END MODULE deck_control_block
