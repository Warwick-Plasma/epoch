MODULE deck_control_block

  USE strings_advanced
  USE fields
  USE timer

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: control_deck_initialise, control_deck_finalise
  PUBLIC :: control_block_start, control_block_end
  PUBLIC :: control_block_handle_element, control_block_check

  INTEGER, PARAMETER :: control_block_elements = 25 + 4 * c_ndims
  LOGICAL, DIMENSION(control_block_elements) :: control_block_done
  CHARACTER(LEN=string_length), DIMENSION(control_block_elements) :: &
      control_block_name = (/ &
          'nx                       ', &
          'ny                       ', &
          'x_min                    ', &
          'x_max                    ', &
          'y_min                    ', &
          'y_max                    ', &
          'nprocx                   ', &
          'nprocy                   ', &
          'npart                    ', &
          'nsteps                   ', &
          't_end                    ', &
          'dt_multiplier            ', &
          'dlb_threshold            ', &
          'icfile                   ', &
          'restart_snapshot         ', &
          'neutral_background       ', &
          'field_order              ', &
          'stdout_frequency         ', &
          'use_random_seed          ', &
          'smooth_currents          ', &
          'field_ionisation         ', &
          'use_multiphoton          ', &
          'use_bsi                  ', &
          'particle_tstart          ', &
          'use_migration            ', &
          'migration_interval       ', &
          'use_exact_restart        ', &
          'allow_cpu_reduce         ', &
          'check_stop_file_frequency', &
          'stop_at_walltime         ', &
          'stop_at_walltime_file    ', &
          'simplify_deck            ', &
          'print_constants          ' /)
  CHARACTER(LEN=string_length), DIMENSION(control_block_elements) :: &
      alternate_name = (/ &
          'nx                       ', &
          'ny                       ', &
          'x_start                  ', &
          'x_end                    ', &
          'y_start                  ', &
          'y_end                    ', &
          'nprocx                   ', &
          'nprocy                   ', &
          'npart                    ', &
          'nsteps                   ', &
          't_end                    ', &
          'dt_multiplier            ', &
          'dlb_threshold            ', &
          'icfile                   ', &
          'restart_snapshot         ', &
          'neutral_background       ', &
          'field_order              ', &
          'stdout_frequency         ', &
          'use_random_seed          ', &
          'smooth_currents          ', &
          'use_field_ionise         ', &
          'multiphoton              ', &
          'bsi                      ', &
          'particle_tstart          ', &
          'migrate_particles        ', &
          'migration_interval       ', &
          'use_exact_restart        ', &
          'allow_cpu_reduce         ', &
          'check_stop_frequency     ', &
          'stop_at_walltime         ', &
          'stop_at_walltime_file    ', &
          'simplify_deck            ', &
          'print_constants          ' /)

CONTAINS

  SUBROUTINE control_deck_initialise

    IF (deck_state .EQ. c_ds_first) THEN
      control_block_done = .FALSE.
      use_exact_restart = .FALSE.
      allow_cpu_reduce = .TRUE.
      check_walltime = .FALSE.
      simplify_deck = .TRUE.
      print_deck_constants = .FALSE.
      restart_number = 0
      check_stop_frequency = 10
      stop_at_walltime = -1.0_num
      restart_filename = ''
    ENDIF

  END SUBROUTINE control_deck_initialise



  SUBROUTINE control_deck_finalise

    CHARACTER(LEN=22) :: filename_fmt

    IF (.NOT.ic_from_restart) use_exact_restart = .FALSE.

    IF (ic_from_restart) THEN
      IF (TRIM(restart_filename) .EQ. '') THEN
        WRITE(filename_fmt, '(''(i'', i3.3, ''.'', i3.3, '', ".sdf")'')') &
            n_zeros, n_zeros
        WRITE(restart_filename, filename_fmt) restart_number
      ENDIF
      full_restart_filename = TRIM(filesystem) &
          // TRIM(data_dir) // '/' // TRIM(restart_filename)
    ENDIF

    IF (deck_state .EQ. c_ds_first) RETURN

    IF (stop_at_walltime .GE. 0.0_num) THEN
      check_walltime = .TRUE.
      timer_collect = .TRUE.
    ENDIF

  END SUBROUTINE control_deck_finalise



  SUBROUTINE control_block_start

  END SUBROUTINE control_block_start



  SUBROUTINE control_block_end

  END SUBROUTINE control_block_end



  FUNCTION control_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    CHARACTER(LEN=string_length) :: str_tmp
    CHARACTER(LEN=1) :: c
    INTEGER :: errcode
    INTEGER :: loop, elementselected, field_order, ierr, io, iu, i
    LOGICAL :: isnum

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
      nx_global = as_integer_print(value, element, errcode)
    CASE(2)
      ny_global = as_integer_print(value, element, errcode)
    CASE(c_ndims+1)
      x_min = as_real_print(value, element, errcode)
    CASE(c_ndims+2)
      x_max = as_real_print(value, element, errcode)
    CASE(c_ndims+3)
      y_min = as_real_print(value, element, errcode)
    CASE(c_ndims+4)
      y_max = as_real_print(value, element, errcode)
    CASE(3*c_ndims+1)
      nprocx = as_integer_print(value, element, errcode)
    CASE(3*c_ndims+2)
      nprocy = as_integer_print(value, element, errcode)
    CASE(4*c_ndims+1)
      npart_global = as_long_integer_print(value, element, errcode)
    CASE(4*c_ndims+2)
      nsteps = as_integer_print(value, element, errcode)
    CASE(4*c_ndims+3)
      t_end = as_real_print(value, element, errcode)
    CASE(4*c_ndims+4)
      dt_multiplier = as_real_print(value, element, errcode)
    CASE(4*c_ndims+5)
      dlb_threshold = as_real_print(value, element, errcode)
      use_balance = .TRUE.
    CASE(4*c_ndims+6)
      IF (rank .EQ. 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The "icfile" option is no longer supported.'
          WRITE(io,*) 'Please use the "import" directive instead'
        ENDDO
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    CASE(4*c_ndims+7)
      isnum = .TRUE.
      str_tmp = TRIM(ADJUSTL(value))
      DO i = 1,LEN_TRIM(str_tmp)
        c = str_tmp(i:i)
        IF (c .LT. '0' .OR. c .GT. '9') THEN
          isnum = .FALSE.
          EXIT
        ENDIF
      ENDDO
      IF (isnum) THEN
        restart_number = as_integer_print(value, element, errcode)
      ELSE
        restart_filename = TRIM(str_tmp)
      ENDIF
      ic_from_restart = .TRUE.
    CASE(4*c_ndims+8)
      neutral_background = as_logical_print(value, element, errcode)
    CASE(4*c_ndims+9)
      field_order = as_integer_print(value, element, errcode)
      IF (field_order .NE. 2 .AND. field_order .NE. 4 &
          .AND. field_order .NE. 6) THEN
        errcode = c_err_bad_value
      ELSE
        CALL set_field_order(field_order)
      ENDIF
    CASE(4*c_ndims+10)
      stdout_frequency = as_integer_print(value, element, errcode)
    CASE(4*c_ndims+11)
      use_random_seed = as_logical_print(value, element, errcode)
    CASE(4*c_ndims+12)
      smooth_currents = as_logical_print(value, element, errcode)
    CASE(4*c_ndims+13)
      use_field_ionisation = as_logical_print(value, element, errcode)
    CASE(4*c_ndims+14)
      use_multiphoton = as_logical_print(value, element, errcode)
    CASE(4*c_ndims+15)
      use_bsi = as_logical_print(value, element, errcode)
    CASE(4*c_ndims+16)
      particle_push_start_time = as_real_print(value, element, errcode)
    CASE(4*c_ndims+17)
      use_particle_migration = as_logical_print(value, element, errcode)
    CASE(4*c_ndims+18)
      particle_migration_interval = as_integer_print(value, element, errcode)
    CASE(4*c_ndims+19)
      use_exact_restart = as_logical_print(value, element, errcode)
    CASE(4*c_ndims+20)
      allow_cpu_reduce = as_logical_print(value, element, errcode)
    CASE(4*c_ndims+21)
      check_stop_frequency = as_integer_print(value, element, errcode)
    CASE(4*c_ndims+22)
      stop_at_walltime = as_real_print(value, element, errcode)
    CASE(4*c_ndims+23)
      IF (rank .EQ. 0) THEN
        OPEN(unit=lu, status='OLD', iostat=ierr, &
            file=TRIM(data_dir) // '/' // TRIM(value))
        IF (ierr .EQ. 0) THEN
          READ(lu,*,iostat=ierr) stop_at_walltime
          CLOSE(lu)
        ENDIF
      ENDIF
      CALL MPI_BCAST(stop_at_walltime, 1, mpireal, 0, comm, errcode)
    CASE(4*c_ndims+24)
      simplify_deck = as_logical_print(value, element, errcode)
    CASE(4*c_ndims+25)
      print_deck_constants = as_logical_print(value, element, errcode)
    END SELECT

  END FUNCTION control_block_handle_element



  FUNCTION control_block_check() RESULT(errcode)

    INTEGER :: errcode, index, io, iu

    errcode = c_err_none

    ! nprocx/y/z and npart are optional
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
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
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
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
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
