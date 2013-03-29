MODULE deck_io_block

  USE strings_advanced

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: io_deck_initialise, io_deck_finalise
  PUBLIC :: io_block_start, io_block_end
  PUBLIC :: io_block_handle_element, io_block_check

  INTEGER, PARAMETER :: io_block_elements = num_vars_to_dump + 23
  INTEGER :: block_number, full_io_block, restart_io_block
  LOGICAL, DIMENSION(io_block_elements) :: io_block_done
  LOGICAL, PRIVATE :: got_name, got_dump_source_code, got_dump_input_decks
  CHARACTER(LEN=string_length), DIMENSION(io_block_elements) :: io_block_name
  CHARACTER(LEN=string_length), DIMENSION(io_block_elements) :: alternate_name
  CHARACTER(LEN=string_length) :: name
  TYPE(io_block_type), POINTER :: io_block

CONTAINS

  SUBROUTINE io_deck_initialise

    INTEGER :: i

    block_number = 0
    IF (deck_state .NE. c_ds_first) RETURN

    alternate_name = ''
    io_block_name(c_dump_part_grid        ) = 'particles'
    io_block_name(c_dump_grid             ) = 'grid'
    io_block_name(c_dump_part_species     ) = 'species_id'
    io_block_name(c_dump_part_weight      ) = 'particle_weight'
    io_block_name(c_dump_part_px          ) = 'px'
    io_block_name(c_dump_part_py          ) = 'py'
    io_block_name(c_dump_part_pz          ) = 'pz'
    io_block_name(c_dump_part_vx          ) = 'vx'
    io_block_name(c_dump_part_vy          ) = 'vy'
    io_block_name(c_dump_part_vz          ) = 'vz'
    io_block_name(c_dump_part_charge      ) = 'charge'
    io_block_name(c_dump_part_mass        ) = 'mass'
    io_block_name(c_dump_part_id          ) = 'id'
    io_block_name(c_dump_part_ek          ) = 'ek'
    alternate_name(c_dump_part_ek         ) = 'particle_energy'

    io_block_name(c_dump_part_opdepth     ) = ''
    io_block_name(c_dump_part_qed_energy  ) = ''
    io_block_name(c_dump_part_opdepth_tri ) = ''
#ifdef PHOTONS
    io_block_name(c_dump_part_opdepth     ) = 'optical_depth'
    io_block_name(c_dump_part_qed_energy  ) = 'qed_energy'
#ifdef TRIDENT_PHOTONS
    io_block_name(c_dump_part_opdepth_tri ) = 'trident_optical_depth'
#endif
#endif

    io_block_name(c_dump_ex               ) = 'ex'
    io_block_name(c_dump_ey               ) = 'ey'
    io_block_name(c_dump_ez               ) = 'ez'
    io_block_name(c_dump_bx               ) = 'bx'
    io_block_name(c_dump_by               ) = 'by'
    io_block_name(c_dump_bz               ) = 'bz'
    io_block_name(c_dump_jx               ) = 'jx'
    io_block_name(c_dump_jy               ) = 'jy'
    io_block_name(c_dump_jz               ) = 'jz'
    io_block_name(c_dump_ekbar            ) = 'ekbar'
    io_block_name(c_dump_mass_density     ) = 'mass_density'
    io_block_name(c_dump_charge_density   ) = 'charge_density'
    io_block_name(c_dump_number_density   ) = 'number_density'
    io_block_name(c_dump_temperature      ) = 'temperature'
    io_block_name(c_dump_dist_fns         ) = 'distribution_functions'
    io_block_name(c_dump_probes           ) = 'particle_probes'
    io_block_name(c_dump_ejected_particles) = 'ejected_particles'
    io_block_name(c_dump_ekflux           ) = 'ekflux'
    io_block_name(c_dump_poynt_flux       ) = 'poynt_flux'
    io_block_name(c_dump_cpml_psi_eyx     ) = 'cpml_psi_eyx'
    io_block_name(c_dump_cpml_psi_ezx     ) = 'cpml_psi_ezx'
    io_block_name(c_dump_cpml_psi_byx     ) = 'cpml_psi_byx'
    io_block_name(c_dump_cpml_psi_bzx     ) = 'cpml_psi_bzx'
    io_block_name(c_dump_cpml_psi_exy     ) = 'cpml_psi_exy'
    io_block_name(c_dump_cpml_psi_ezy     ) = 'cpml_psi_ezy'
    io_block_name(c_dump_cpml_psi_bxy     ) = 'cpml_psi_bxy'
    io_block_name(c_dump_cpml_psi_bzy     ) = 'cpml_psi_bzy'
    io_block_name(c_dump_absorption       ) = 'absorption'

    i = num_vars_to_dump
    io_block_name (i+1 ) = 'dt_snapshot'
    io_block_name (i+2 ) = 'full_dump_every'
    io_block_name (i+3 ) = 'restart_dump_every'
    io_block_name (i+4 ) = 'force_first_to_be_restartable'
    io_block_name (i+5 ) = 'force_final_to_be_restartable'
    alternate_name(i+5 ) = 'force_last_to_be_restartable'
    io_block_name (i+6 ) = 'use_offset_grid'
    io_block_name (i+7 ) = 'extended_io_file'
    io_block_name (i+8 ) = 'dt_average'
    alternate_name(i+8 ) = 'averaging_period'
    io_block_name (i+9 ) = 'nstep_average'
    alternate_name(i+9 ) = 'min_cycles_per_average'
    io_block_name (i+10) = 'nstep_snapshot'
    io_block_name (i+11) = 'dump_source_code'
    io_block_name (i+12) = 'dump_input_decks'
    io_block_name (i+13) = 'dump_first'
    io_block_name (i+14) = 'dump_last'
    alternate_name(i+14) = 'dump_final'
    io_block_name (i+15) = 'restartable'
    io_block_name (i+16) = 'name'
    io_block_name (i+17) = 'time_start'
    io_block_name (i+18) = 'time_stop'
    io_block_name (i+19) = 'nstep_start'
    io_block_name (i+20) = 'nstep_stop'
    io_block_name (i+21) = 'dump_at_nsteps'
    alternate_name(i+21) = 'nsteps_dump'
    io_block_name (i+22) = 'dump_at_times'
    alternate_name(i+22) = 'times_dump'
    io_block_name (i+23) = 'dump_cycle'

    track_ejected_particles = .FALSE.
    averaged_var_block = 0
    new_style_io_block = .FALSE.
    n_io_blocks = 0

  END SUBROUTINE io_deck_initialise



  SUBROUTINE io_deck_finalise

    INTEGER :: i, io, ierr

    n_io_blocks = block_number

    IF (n_io_blocks .GT. 0) THEN
      IF (deck_state .EQ. c_ds_first) THEN
        IF (.NOT.new_style_io_block .AND. n_io_blocks .NE. 1) THEN
          IF (rank .EQ. 0) THEN
            DO io = stdout, du, du - stdout ! Print to stdout and to file
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Cannot have multiple unnamed "output" blocks.'
            ENDDO
          ENDIF
          CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        ENDIF

        ALLOCATE(io_block_list(n_io_blocks))
        DO i = 1, n_io_blocks
          CALL init_io_block(io_block_list(i))
        ENDDO
      ELSE
        DO i = 1, n_io_blocks
          IF (io_block_list(i)%dt_average .GT. t_end) THEN
            IF (rank .EQ. 0) THEN
              DO io = stdout, du, du - stdout ! Print to stdout and to file
                WRITE(io,*) '*** WARNING ***'
                WRITE(io,*) 'Averaging time is longer than t_end, will set', &
                    ' averaging time equal'
                WRITE(io,*) 'to t_end.'
              ENDDO
            ENDIF
            io_block_list(i)%dt_average = t_end
          ENDIF
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE io_deck_finalise



  SUBROUTINE io_block_start

    io_block_done = .FALSE.
    got_name = .FALSE.
    got_dump_source_code = .FALSE.
    got_dump_input_decks = .FALSE.
    block_number = block_number + 1
    IF (deck_state .NE. c_ds_first .AND. block_number .GT. 0) THEN
      io_block => io_block_list(block_number)
      io_block%dump_first = .FALSE.
    ENDIF

  END SUBROUTINE io_block_start



  SUBROUTINE io_block_end

    INTEGER :: io, ierr
    CHARACTER(LEN=c_max_string_length) :: list_filename

    IF (deck_state .EQ. c_ds_first) RETURN

    IF (io_block%dumpmask(c_dump_ejected_particles) .NE. c_io_never) THEN
      track_ejected_particles = .TRUE.
    ENDIF

    IF (.NOT. got_name) THEN
      IF (new_style_io_block) THEN
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Cannot mix old and new style output blocks.'
            WRITE(io,*) 'You can either have multiple, named output blocks ', &
                'or a single unnamed one.'
          ENDDO
        ENDIF
        CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      ENDIF
      io_block%name = 'normal'
    ENDIF

    ! Delete any existing visit file lists
    IF (rank .EQ. 0) THEN
      list_filename = TRIM(ADJUSTL(data_dir)) // '/' &
          // TRIM(io_block%name) // '.visit'
      OPEN(unit=lu, status='UNKNOWN', file=list_filename)
      CLOSE(unit=lu, status='DELETE')
    ENDIF
    io_block%dumpmask(c_dump_jx) = IOR(io_block%dumpmask(c_dump_jx), c_io_field)
    io_block%dumpmask(c_dump_jy) = IOR(io_block%dumpmask(c_dump_jy), c_io_field)
    io_block%dumpmask(c_dump_jz) = IOR(io_block%dumpmask(c_dump_jz), c_io_field)

    IF (.NOT.got_dump_source_code) THEN
      IF (io_block%restart .OR. .NOT.new_style_io_block) &
          io_block%dump_source_code = .TRUE.
    ENDIF

    IF (.NOT.got_dump_input_decks) THEN
      IF (io_block%restart .OR. .NOT.new_style_io_block) &
          io_block%dump_input_decks = .TRUE.
    ENDIF

    CALL set_restart_dumpmasks

  END SUBROUTINE io_block_end



  FUNCTION io_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode, style_error
    INTEGER :: loop, elementselected, mask, mask_element, ierr, io
    INTEGER :: is, subset, n_list
    INTEGER, ALLOCATABLE :: subsets(:)
    LOGICAL :: bad
    INTEGER, PARAMETER :: c_err_new_style_ignore = 1
    INTEGER, PARAMETER :: c_err_new_style_global = 2
    INTEGER, PARAMETER :: c_err_old_style_ignore = 3

    errcode = c_err_none
    IF (deck_state .EQ. c_ds_first) THEN
      IF (str_cmp(element, 'name')) new_style_io_block = .TRUE.
      RETURN
    ENDIF

    errcode = c_err_unknown_element

    elementselected = 0

    DO loop = 1, io_block_elements
      IF (str_cmp(element, TRIM(ADJUSTL(io_block_name(loop)))) &
          .OR. str_cmp(element, TRIM(ADJUSTL(alternate_name(loop))))) THEN
        elementselected = loop
        EXIT
      ENDIF
    ENDDO

    IF (elementselected .EQ. 0) RETURN

    IF (io_block_done(elementselected)) THEN
      errcode = c_err_preset_element
      RETURN
    ENDIF
    io_block_done(elementselected) = .TRUE.
    errcode = c_err_none
    style_error = c_err_none

    SELECT CASE (elementselected-num_vars_to_dump)
    CASE(1)
      io_block%dt_snapshot = as_real(value, errcode)
      IF (io_block%dt_snapshot .LT. 0.0_num) io_block%dt_snapshot = 0.0_num
    CASE(2)
      IF (new_style_io_block) THEN
        style_error = c_err_new_style_ignore
      ELSE
        full_dump_every = as_integer(value, errcode)
        IF (full_dump_every .EQ. 0) full_dump_every = 1
      ENDIF
    CASE(3)
      IF (new_style_io_block) THEN
        style_error = c_err_new_style_ignore
      ELSE
        restart_dump_every = as_integer(value, errcode)
        IF (restart_dump_every .EQ. 0) restart_dump_every = 1
      ENDIF
    CASE(4)
      IF (new_style_io_block) style_error = c_err_new_style_global
      force_first_to_be_restartable = as_logical(value, errcode)
    CASE(5)
      IF (new_style_io_block) style_error = c_err_new_style_global
      force_final_to_be_restartable = as_logical(value, errcode)
    CASE(6)
      IF (new_style_io_block) style_error = c_err_new_style_global
      use_offset_grid = as_logical(value, errcode)
    CASE(7)
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The "extended_io_file" option is no longer supported.'
          WRITE(io,*) 'Please use the "import" directive instead'
        ENDDO
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    CASE(8)
      io_block%dt_average = as_real(value, errcode)
    CASE(9)
      io_block%nstep_average = as_integer(value, errcode)
    CASE(10)
      io_block%nstep_snapshot = as_integer(value, errcode)
      IF (io_block%nstep_snapshot .LT. 0) io_block%nstep_snapshot = 0
    CASE(11)
      io_block%dump_source_code = as_logical(value, errcode)
      got_dump_source_code = .TRUE.
    CASE(12)
      io_block%dump_input_decks = as_logical(value, errcode)
      got_dump_input_decks = .TRUE.
    CASE(13)
      io_block%dump_first = as_logical(value, errcode)
    CASE(14)
      io_block%dump_last = as_logical(value, errcode)
    CASE(15)
      IF (.NOT.new_style_io_block) style_error = c_err_old_style_ignore
      io_block%restart = as_logical(value, errcode)
    CASE(16)
      io_block%name = value
      got_name = .TRUE.
    CASE(17)
      io_block%time_start = as_real(value, errcode)
    CASE(18)
      io_block%time_stop = as_real(value, errcode)
    CASE(19)
      io_block%nstep_start = as_integer(value, errcode)
    CASE(20)
      io_block%nstep_stop = as_integer(value, errcode)
    CASE(21)
      IF (.NOT.new_style_io_block) style_error = c_err_old_style_ignore
      CALL get_allocated_array(value, io_block%dump_at_nsteps, errcode)
    CASE(22)
      IF (.NOT.new_style_io_block) style_error = c_err_old_style_ignore
      CALL get_allocated_array(value, io_block%dump_at_times, errcode)
    CASE(23)
      io_block%dump_cycle = as_integer(value, errcode)
    END SELECT

    IF (style_error .EQ. c_err_old_style_ignore) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Element "' &
              // TRIM(ADJUSTL(io_block_name(elementselected))) &
              // '" not ', 'allowed in an unnamed output block.'
          WRITE(io,*) 'It has been ignored.'
          WRITE(io,*)
        ENDDO
      ENDIF
    ELSE IF (style_error .EQ. c_err_new_style_ignore) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Element "' &
              // TRIM(ADJUSTL(io_block_name(elementselected))) &
              // '" not ', 'allowed in a named output block.'
          WRITE(io,*) 'It has been ignored.'
          WRITE(io,*)
        ENDDO
      ENDIF
    ELSE IF (style_error .EQ. c_err_new_style_global) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Element "' &
              // TRIM(ADJUSTL(io_block_name(elementselected))) &
              // '" should be moved to ', 'an "output_global" block.'
          WRITE(io,*) 'Its value will be applied to all output blocks.'
          WRITE(io,*)
        ENDDO
      ENDIF
    ENDIF

    IF (elementselected .GT. num_vars_to_dump) RETURN

    mask_element = elementselected
    ALLOCATE(subsets(n_subsets+1))
    CALL as_list(value, subsets, n_list, errcode)

    DO is = 1, n_list
      subset = subsets(is)
      IF (is .EQ. 1) THEN
        mask = subset
      ELSE
        mask = subset_list(subset)%mask
      ENDIF

      ! If setting dumpmask for features which haven't been compiled
      ! in then issue a warning
#ifndef PARTICLE_PROBES
      IF (mask_element .EQ. c_dump_probes .AND. mask .NE. c_io_never) THEN
        errcode = c_err_pp_options_wrong
        extended_error_string = '-DPARTICLE_PROBES'
        mask = c_io_never
      ENDIF
#endif

#ifdef PARTICLE_ID4
#define PARTICLE_ID
#endif
#ifndef PARTICLE_ID
      IF (mask_element .EQ. c_dump_part_id .AND. mask .NE. c_io_never) THEN
        errcode = c_err_pp_options_wrong
        extended_error_string = '-DPARTICLE_ID'
        mask = c_io_never
      ENDIF
#endif

      ! Setting some flags like species and average
      ! wastes memory if the parameters make no sense. Do sanity checking.

      IF (IAND(mask, c_io_species) .NE. 0 &
          .OR. IAND(mask, c_io_no_sum) .EQ. 0) THEN
        bad = .TRUE.
        ! Check for sensible per species variables
        IF (mask_element .EQ. c_dump_ekbar) bad = .FALSE.
        IF (mask_element .EQ. c_dump_ekflux) bad = .FALSE.
        IF (mask_element .EQ. c_dump_mass_density) bad = .FALSE.
        IF (mask_element .EQ. c_dump_charge_density) bad = .FALSE.
        IF (mask_element .EQ. c_dump_number_density) bad = .FALSE.
        IF (mask_element .EQ. c_dump_temperature) bad = .FALSE.
        IF (mask_element .EQ. c_dump_jx) bad = .FALSE.
        IF (mask_element .EQ. c_dump_jy) bad = .FALSE.
        IF (mask_element .EQ. c_dump_jz) bad = .FALSE.
        IF (bad) THEN
          IF (rank .EQ. 0 .AND. IAND(mask, c_io_species) .NE. 0) THEN
            DO io = stdout, du, du - stdout ! Print to stdout and to file
              WRITE(io,*) '*** WARNING ***'
              WRITE(io,*) 'Attempting to set per species property for "' &
                  // TRIM(element) // '" which'
              WRITE(io,*) 'does not support this property. Ignoring.'
            ENDDO
          ENDIF
          mask = IAND(mask, NOT(c_io_species))
          mask = IOR(mask, c_io_no_sum)
        ENDIF
      ENDIF

      IF (IAND(mask, c_io_averaged) .NE. 0) THEN
        bad = .TRUE.
        ! Check for sensible averaged variables
        IF (mask_element .EQ. c_dump_ex) bad = .FALSE.
        IF (mask_element .EQ. c_dump_ey) bad = .FALSE.
        IF (mask_element .EQ. c_dump_ez) bad = .FALSE.
        IF (mask_element .EQ. c_dump_bx) bad = .FALSE.
        IF (mask_element .EQ. c_dump_by) bad = .FALSE.
        IF (mask_element .EQ. c_dump_bz) bad = .FALSE.
        IF (mask_element .EQ. c_dump_jx) bad = .FALSE.
        IF (mask_element .EQ. c_dump_jy) bad = .FALSE.
        IF (mask_element .EQ. c_dump_jz) bad = .FALSE.

        ! Unset 'no_sum' dumpmask for grid variables
        IF (.NOT.bad) mask = IAND(mask, NOT(c_io_no_sum))

        IF (mask_element .EQ. c_dump_ekbar) bad = .FALSE.
        IF (mask_element .EQ. c_dump_mass_density) bad = .FALSE.
        IF (mask_element .EQ. c_dump_charge_density) bad = .FALSE.
        IF (mask_element .EQ. c_dump_number_density) bad = .FALSE.
        IF (mask_element .EQ. c_dump_temperature) bad = .FALSE.
        IF (bad) THEN
          IF (rank .EQ. 0) THEN
            DO io = stdout, du, du - stdout ! Print to stdout and to file
              WRITE(io,*) '*** WARNING ***'
              WRITE(io,*) 'Attempting to set average property for "' &
                  // TRIM(element) // '" which'
              WRITE(io,*) 'does not support this property. Ignoring.'
            ENDDO
          ENDIF
          mask = IAND(mask, NOT(c_io_averaged))
        ELSE
          any_average = .TRUE.
          io_block%any_average = .TRUE.
          IF (IAND(mask, c_io_average_single) .NE. 0 .AND. num .NE. r4) THEN
            io_block%averaged_data(mask_element)%dump_single = .TRUE.
          ENDIF
          IF (averaged_var_block(mask_element) .NE. 0) THEN
            PRINT*,'error ',mask_element,block_number
          ENDIF
          averaged_var_block(mask_element) = block_number
        ENDIF
      ENDIF

      IF (is .EQ. 1) THEN
        io_block%dumpmask(mask_element) = mask
      ELSE
        subset_list(subset)%dumpmask(block_number,mask_element) = mask
      ENDIF
    ENDDO

    DEALLOCATE(subsets)

  END FUNCTION io_block_handle_element



  FUNCTION io_block_check() RESULT(errcode)

    INTEGER :: errcode, io, i

    ! Just assume that anything not included except for the compulsory
    ! elements is not wanted
    errcode = c_err_none

    ! Other control parameters are optional
    i = num_vars_to_dump
    io_block_done(i+1:i+7) = .TRUE.
    io_block_done(i+10:io_block_elements) = .TRUE.
    ! Averaging info not compulsory unless averaged variable selected
    IF (.NOT. any_average) io_block_done(i+8:i+9) = .TRUE.

    IF (.NOT. io_block_done(i+8) .AND. .NOT. io_block_done(i+9)) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Required output block element ', &
              TRIM(ADJUSTL(io_block_name(i+8))), &
              ' absent. Please create this entry in the input deck'
        ENDDO
      ENDIF
      errcode = c_err_missing_elements
    ENDIF

    ! Can't check the io_block if it hasn't been allocated.
    IF (n_io_blocks .EQ. 0) RETURN

    IF (io_block%dt_average .GT. io_block%dt_snapshot) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Averaging time is longer than dt_snapshot, will set', &
              ' averaging time equal'
          WRITE(io,*) 'to dt_snapshot.'
        ENDDO
      ENDIF
      io_block%dt_average = io_block%dt_snapshot
    ENDIF

  END FUNCTION io_block_check



  SUBROUTINE init_io_block(io_block)

    TYPE(io_block_type) :: io_block
    INTEGER :: i

    io_block%name = ''
    io_block%dt_snapshot = -1.0_num
    io_block%time_prev = 0.0_num
    io_block%time_first = 0.0_num
    io_block%dt_average = -1.0_num
    io_block%dt_min_average = -1.0_num
    io_block%average_time = -1.0_num
    io_block%nstep_snapshot = -1
    io_block%nstep_prev = 0
    io_block%nstep_first = 0
    io_block%nstep_average = -1
    io_block%restart = .FALSE.
    io_block%dump = .FALSE.
    io_block%any_average = .FALSE.
    io_block%dump_first = .FALSE.
    io_block%dump_last = .TRUE.
    io_block%dump_source_code = .FALSE.
    io_block%dump_input_decks = .FALSE.
    io_block%dumpmask = 0
    io_block%time_start = -1.0_num
    io_block%time_stop  = HUGE(1.0_num)
    io_block%nstep_start = -1
    io_block%nstep_stop  = HUGE(1)
    io_block%dump_cycle  = HUGE(1)
    NULLIFY(io_block%dump_at_nsteps)
    NULLIFY(io_block%dump_at_times)
    DO i = 1, num_vars_to_dump
      io_block%averaged_data(i)%dump_single = .FALSE.
    ENDDO

  END SUBROUTINE init_io_block



  SUBROUTINE set_restart_dumpmasks

    ! Set the dumpmask for variables required to restart

    ! Particles
    io_block%dumpmask(c_dump_part_grid) = &
        IOR(io_block%dumpmask(c_dump_part_grid), c_io_restartable)
    io_block%dumpmask(c_dump_part_species) = &
        IOR(io_block%dumpmask(c_dump_part_species), c_io_restartable)
    io_block%dumpmask(c_dump_part_weight) = &
        IOR(io_block%dumpmask(c_dump_part_weight), c_io_restartable)
    io_block%dumpmask(c_dump_part_px) = &
        IOR(io_block%dumpmask(c_dump_part_px), c_io_restartable)
    io_block%dumpmask(c_dump_part_py) = &
        IOR(io_block%dumpmask(c_dump_part_py), c_io_restartable)
    io_block%dumpmask(c_dump_part_pz) = &
        IOR(io_block%dumpmask(c_dump_part_pz), c_io_restartable)
#ifdef PHOTONS
    io_block%dumpmask(c_dump_part_opdepth) = &
        IOR(io_block%dumpmask(c_dump_part_opdepth), c_io_restartable)
    io_block%dumpmask(c_dump_part_qed_energy) = &
        IOR(io_block%dumpmask(c_dump_part_qed_energy), c_io_restartable)
#ifdef TRIDENT_PHOTONS
    io_block%dumpmask(c_dump_part_opdepth_tri) = &
        IOR(io_block%dumpmask(c_dump_part_opdepth_tri), c_io_restartable)
#endif
#endif
    ! Fields
    io_block%dumpmask(c_dump_grid) = &
        IOR(io_block%dumpmask(c_dump_grid), c_io_restartable)
    io_block%dumpmask(c_dump_ex) = &
        IOR(io_block%dumpmask(c_dump_ex), c_io_restartable)
    io_block%dumpmask(c_dump_ey) = &
        IOR(io_block%dumpmask(c_dump_ey), c_io_restartable)
    io_block%dumpmask(c_dump_ez) = &
        IOR(io_block%dumpmask(c_dump_ez), c_io_restartable)
    io_block%dumpmask(c_dump_bx) = &
        IOR(io_block%dumpmask(c_dump_bx), c_io_restartable)
    io_block%dumpmask(c_dump_by) = &
        IOR(io_block%dumpmask(c_dump_by), c_io_restartable)
    io_block%dumpmask(c_dump_bz) = &
        IOR(io_block%dumpmask(c_dump_bz), c_io_restartable)
    io_block%dumpmask(c_dump_jx) = &
        IOR(io_block%dumpmask(c_dump_jx), c_io_restartable)
    io_block%dumpmask(c_dump_jy) = &
        IOR(io_block%dumpmask(c_dump_jy), c_io_restartable)
    io_block%dumpmask(c_dump_jz) = &
        IOR(io_block%dumpmask(c_dump_jz), c_io_restartable)
    ! CPML boundaries
    io_block%dumpmask(c_dump_cpml_psi_eyx) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_eyx), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_ezx) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_ezx), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_byx) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_byx), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_bzx) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_bzx), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_exy) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_exy), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_ezy) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_ezy), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_bxy) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_bxy), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_bzy) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_bzy), c_io_restartable)

  END SUBROUTINE set_restart_dumpmasks

END MODULE deck_io_block
