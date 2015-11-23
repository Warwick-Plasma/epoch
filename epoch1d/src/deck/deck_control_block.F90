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

  INTEGER, PARAMETER :: control_block_elements = 28 + 4 * c_ndims
  LOGICAL, DIMENSION(control_block_elements) :: control_block_done
  CHARACTER(LEN=string_length), DIMENSION(control_block_elements) :: &
      control_block_name = (/ &
          'nx                       ', &
          'x_min                    ', &
          'x_max                    ', &
          'nprocx                   ', &
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
          'print_constants          ', &
          'allow_missing_restart    ', &
          'print_eta_string         ', &
          'n_zeros                  ' /)
  CHARACTER(LEN=string_length), DIMENSION(control_block_elements) :: &
      alternate_name = (/ &
          'nx                       ', &
          'x_start                  ', &
          'x_end                    ', &
          'nprocx                   ', &
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
          'print_constants          ', &
          'allow_missing_restart    ', &
          'print_eta_string         ', &
          'n_zeros                  ' /)

CONTAINS

  SUBROUTINE control_deck_initialise

    IF (deck_state == c_ds_first) THEN
      control_block_done = .FALSE.
      use_exact_restart = .FALSE.
      allow_cpu_reduce = .TRUE.
      check_walltime = .FALSE.
      simplify_deck = .TRUE.
      print_deck_constants = .FALSE.
      allow_missing_restart = .FALSE.
      print_eta_string = .FALSE.
      restart_number = 0
      check_stop_frequency = 10
      stop_at_walltime = -1.0_num
      restart_filename = ''
      n_zeros_control = -1
    ENDIF

  END SUBROUTINE control_deck_initialise



  SUBROUTINE control_deck_finalise

    CHARACTER(LEN=22) :: filename_fmt
    INTEGER :: io, iu

    IF (n_zeros_control > 0) THEN
      IF (n_zeros_control < 4) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'n_zeros was less than 4 and has been ignored'
          ENDDO
        ENDIF
        n_zeros_control = -1
      ELSE
        n_zeros = n_zeros_control
      ENDIF
    ENDIF

    IF (ic_from_restart) THEN
      IF (TRIM(restart_filename) == '') THEN
        WRITE(filename_fmt, '(''(i'', i3.3, ''.'', i3.3, '', ".sdf")'')') &
            n_zeros, n_zeros
        WRITE(restart_filename, filename_fmt) restart_number
      ENDIF
      full_restart_filename = TRIM(filesystem) &
          // TRIM(data_dir) // '/' // TRIM(restart_filename)

      CALL check_valid_restart
    ENDIF

    IF (.NOT.ic_from_restart) use_exact_restart = .FALSE.

    IF (deck_state == c_ds_first) RETURN

    IF (stop_at_walltime >= 0.0_num) THEN
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

    IF (deck_state /= c_ds_first) RETURN

    errcode = c_err_unknown_element
    elementselected = 0

    DO loop = 1, control_block_elements
      IF (str_cmp(element, TRIM(ADJUSTL(control_block_name(loop)))) &
          .OR. str_cmp(element, TRIM(ADJUSTL(alternate_name(loop))))) THEN
        elementselected = loop
        EXIT
      ENDIF
    ENDDO

    IF (elementselected == 0) RETURN

    IF (control_block_done(elementselected)) THEN
      errcode = c_err_preset_element
      RETURN
    ENDIF

    control_block_done(elementselected) = .TRUE.
    errcode = c_err_none

    SELECT CASE (elementselected)
    CASE(1)
      nx_global = as_integer_print(value, element, errcode)
    CASE(c_ndims+1)
      x_min = as_real_print(value, element, errcode)
    CASE(c_ndims+2)
      x_max = as_real_print(value, element, errcode)
    CASE(3*c_ndims+1)
      nprocx = as_integer_print(value, element, errcode)
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
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The "icfile" option is no longer supported.'
          WRITE(io,*) 'Please use the "import" directive instead'
        ENDDO
      ENDIF
      CALL abort_code(c_err_bad_value)
    CASE(4*c_ndims+7)
      isnum = .TRUE.
      str_tmp = TRIM(ADJUSTL(value))
      DO i = 1,LEN_TRIM(str_tmp)
        c = str_tmp(i:i)
        IF (c < '0' .OR. c > '9') THEN
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
      IF (field_order /= 2 .AND. field_order /= 4 &
          .AND. field_order /= 6) THEN
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
      IF (rank == 0) THEN
        OPEN(unit=lu, status='OLD', iostat=ierr, &
            file=TRIM(data_dir) // '/' // TRIM(value))
        IF (ierr == 0) THEN
          READ(lu,*,iostat=ierr) stop_at_walltime
          CLOSE(lu)
        ENDIF
      ENDIF
      CALL MPI_BCAST(stop_at_walltime, 1, mpireal, 0, comm, errcode)
    CASE(4*c_ndims+24)
      simplify_deck = as_logical_print(value, element, errcode)
    CASE(4*c_ndims+25)
      print_deck_constants = as_logical_print(value, element, errcode)
    CASE(4*c_ndims+26)
      allow_missing_restart = as_logical_print(value, element, errcode)
    CASE(4*c_ndims+27)
      print_eta_string = as_logical_print(value, element, errcode)
    CASE(4*c_ndims+28)
      n_zeros_control = as_integer_print(value, element, errcode)
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
        IF (rank == 0) THEN
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
      IF (rank == 0) THEN
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



  SUBROUTINE check_valid_restart

    CHARACTER(LEN=c_id_length) :: code_name
    INTEGER :: step, code_io_version, string_len
    REAL(num) :: time
    LOGICAL :: restart_flag
    TYPE(sdf_file_handle) :: sdf_handle
    LOGICAL :: valid = .TRUE.

    CALL sdf_open(sdf_handle, full_restart_filename, comm, c_sdf_read)

    IF (sdf_handle%error_code == 0) THEN
      CALL sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
          string_len, restart_flag)
      CALL sdf_close(sdf_handle)
    ELSE
      restart_flag = .FALSE.
    ENDIF

    IF (.NOT. restart_flag) THEN
      valid = .FALSE.
      IF (.NOT. allow_missing_restart) THEN
        IF (rank == 0) THEN
          PRINT*, '*** ERROR ***'
          PRINT*, 'SDF file ', TRIM(full_restart_filename), &
              ' is not a restart dump. Unable to continue.'
        ENDIF
        CALL abort_code(c_err_io_error)
        STOP
      ENDIF
    ENDIF

    IF (valid .AND. .NOT.str_cmp(code_name, 'Epoch1d')) THEN
      valid = .FALSE.
      IF (.NOT. allow_missing_restart) THEN
        IF (rank == 0) THEN
          PRINT*, '*** ERROR ***'
          PRINT*, 'SDF restart file was not generated by Epoch1d. Unable to ', &
              'continue.'
        ENDIF
        CALL abort_code(c_err_io_error)
        STOP
      ENDIF
    ENDIF

    IF (.NOT. valid) THEN
      ic_from_restart = .FALSE.
      IF (rank == 0) THEN
        PRINT*, '*** WARNING ***'
        PRINT*, 'No valid restart dump found. Using initial conditions instead.'
      ENDIF
    ENDIF

  END SUBROUTINE check_valid_restart

END MODULE deck_control_block
