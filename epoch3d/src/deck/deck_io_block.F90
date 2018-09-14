! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009      Chris Brady <C.S.Brady@warwick.ac.uk>
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE deck_io_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: io_deck_initialise, io_deck_finalise
  PUBLIC :: io_block_start, io_block_end
  PUBLIC :: io_block_handle_element, io_block_check, copy_io_block

  INTEGER, PARAMETER :: ov = 33
  INTEGER, PARAMETER :: io_block_elements = num_vars_to_dump + ov
  INTEGER :: block_number, nfile_prefixes
  INTEGER :: rolling_restart_io_block
  INTEGER :: o1, o2, o3, o4, o5, o6, o7, o8
  LOGICAL, DIMENSION(io_block_elements) :: io_block_done
  LOGICAL, PRIVATE :: got_name, got_dump_source_code, got_dump_input_decks
  LOGICAL, PRIVATE :: warning_printed
  CHARACTER(LEN=string_length), DIMENSION(io_block_elements) :: io_block_name
  CHARACTER(LEN=string_length), DIMENSION(io_block_elements) :: alternate_name
  CHARACTER(LEN=c_id_length), ALLOCATABLE :: io_prefixes(:)
  TYPE(io_block_type), POINTER :: io_block

CONTAINS

  SUBROUTINE io_deck_initialise

    INTEGER :: i

    block_number = 0
    IF (deck_state /= c_ds_first) RETURN

    any_average = .FALSE.
    warning_printed = .FALSE.
    alternate_name = ''
    io_block_name (c_dump_part_grid        ) = 'particles'
    alternate_name(c_dump_part_grid        ) = 'particle_grid'
    io_block_name (c_dump_grid             ) = 'grid'
    alternate_name(c_dump_grid             ) = 'field_grid'
    io_block_name (c_dump_part_species     ) = 'species_id'
    io_block_name (c_dump_part_weight      ) = 'particle_weight'
    alternate_name(c_dump_part_weight      ) = 'weight'
    io_block_name (c_dump_part_px          ) = 'px'
    io_block_name (c_dump_part_py          ) = 'py'
    io_block_name (c_dump_part_pz          ) = 'pz'
    io_block_name (c_dump_part_vx          ) = 'vx'
    io_block_name (c_dump_part_vy          ) = 'vy'
    io_block_name (c_dump_part_vz          ) = 'vz'
    io_block_name (c_dump_part_charge      ) = 'charge'
    io_block_name (c_dump_part_mass        ) = 'mass'
    alternate_name(c_dump_part_mass        ) = 'rest_mass'
    io_block_name (c_dump_part_id          ) = 'id'
    io_block_name (c_dump_part_ek          ) = 'ek'
    alternate_name(c_dump_part_ek          ) = 'particle_energy'
    io_block_name (c_dump_part_rel_mass    ) = 'relativistic_mass'
    io_block_name (c_dump_part_gamma       ) = 'gamma'

    io_block_name (c_dump_part_opdepth     ) = ''
    io_block_name (c_dump_part_qed_energy  ) = ''
    io_block_name (c_dump_part_opdepth_tri ) = ''
#ifdef PHOTONS
    io_block_name (c_dump_part_opdepth     ) = 'optical_depth'
    io_block_name (c_dump_part_qed_energy  ) = 'qed_energy'
#ifdef TRIDENT_PHOTONS
    io_block_name (c_dump_part_opdepth_tri ) = 'trident_optical_depth'
#endif
#endif
#ifdef WORK_DONE_INTEGRATED
    io_block_name (c_dump_part_work_x      ) = 'work_x'
    io_block_name (c_dump_part_work_y      ) = 'work_y'
    io_block_name (c_dump_part_work_z      ) = 'work_z'
    io_block_name (c_dump_part_work_x_total) = 'work_x_total'
    io_block_name (c_dump_part_work_y_total) = 'work_y_total'
    io_block_name (c_dump_part_work_z_total) = 'work_z_total'
#endif

    io_block_name (c_dump_ex               ) = 'ex'
    io_block_name (c_dump_ey               ) = 'ey'
    io_block_name (c_dump_ez               ) = 'ez'
    io_block_name (c_dump_bx               ) = 'bx'
    io_block_name (c_dump_by               ) = 'by'
    io_block_name (c_dump_bz               ) = 'bz'
    io_block_name (c_dump_jx               ) = 'jx'
    io_block_name (c_dump_jy               ) = 'jy'
    io_block_name (c_dump_jz               ) = 'jz'
    io_block_name (c_dump_ekbar            ) = 'ekbar'
    io_block_name (c_dump_mass_density     ) = 'mass_density'
    io_block_name (c_dump_charge_density   ) = 'charge_density'
    io_block_name (c_dump_number_density   ) = 'number_density'
    io_block_name (c_dump_ppc              ) = 'ppc'
    alternate_name(c_dump_ppc              ) = 'particles_per_cell'
    io_block_name (c_dump_average_weight   ) = 'average_weight'
    io_block_name (c_dump_temperature      ) = 'temperature'
    io_block_name (c_dump_dist_fns         ) = 'distribution_functions'
    io_block_name (c_dump_probes           ) = 'particle_probes'
    io_block_name (c_dump_ejected_particles) = 'ejected_particles'
    io_block_name (c_dump_ekflux           ) = 'ekflux'
    io_block_name (c_dump_poynt_flux       ) = 'poynt_flux'
    io_block_name (c_dump_cpml_psi_eyx     ) = 'cpml_psi_eyx'
    io_block_name (c_dump_cpml_psi_ezx     ) = 'cpml_psi_ezx'
    io_block_name (c_dump_cpml_psi_byx     ) = 'cpml_psi_byx'
    io_block_name (c_dump_cpml_psi_bzx     ) = 'cpml_psi_bzx'
    io_block_name (c_dump_cpml_psi_exy     ) = 'cpml_psi_exy'
    io_block_name (c_dump_cpml_psi_ezy     ) = 'cpml_psi_ezy'
    io_block_name (c_dump_cpml_psi_bxy     ) = 'cpml_psi_bxy'
    io_block_name (c_dump_cpml_psi_bzy     ) = 'cpml_psi_bzy'
    io_block_name (c_dump_cpml_psi_exz     ) = 'cpml_psi_exz'
    io_block_name (c_dump_cpml_psi_eyz     ) = 'cpml_psi_eyz'
    io_block_name (c_dump_cpml_psi_bxz     ) = 'cpml_psi_bxz'
    io_block_name (c_dump_cpml_psi_byz     ) = 'cpml_psi_byz'
    io_block_name (c_dump_absorption       ) = 'absorption'
    io_block_name (c_dump_total_energy_sum ) = 'total_energy_sum'

    i = num_vars_to_dump
    o1 = 1
    io_block_name (i+1 ) = 'dt_snapshot'
    io_block_name (i+2 ) = 'full_dump_every'
    io_block_name (i+3 ) = 'restart_dump_every'
    io_block_name (i+4 ) = 'force_first_to_be_restartable'
    io_block_name (i+5 ) = 'force_final_to_be_restartable'
    alternate_name(i+5 ) = 'force_last_to_be_restartable'
    io_block_name (i+6 ) = 'use_offset_grid'
    o2 = 7
    io_block_name (i+7 ) = 'extended_io_file'
    o3 = 8
    io_block_name (i+8 ) = 'dt_average'
    alternate_name(i+8 ) = 'averaging_period'
    o4 = 9
    io_block_name (i+9 ) = 'nstep_average'
    alternate_name(i+9 ) = 'min_cycles_per_average'
    o5 = 10
    io_block_name (i+10) = 'nstep_snapshot'
    io_block_name (i+11) = 'dump_source_code'
    io_block_name (i+12) = 'dump_input_decks'
    io_block_name (i+13) = 'dump_first'
    io_block_name (i+14) = 'dump_last'
    alternate_name(i+14) = 'dump_final'
    o6 = 15
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
    o7 = 23
    io_block_name (i+23) = 'dump_cycle'
    o8 = 24
    io_block_name (i+24) = 'file_prefix'
    io_block_name (i+25) = 'rolling_restart'
    io_block_name (i+26) = 'dump_cycle_first_index'
    io_block_name (i+27) = 'filesystem'
    io_block_name (i+28) = 'dump_first_after_restart'
    io_block_name (i+29) = 'dump_at_walltimes'
    alternate_name(i+29) = 'walltimes_dump'
    io_block_name (i+30) = 'walltime_interval'
    alternate_name(i+30) = 'walltime_snapshot'
    io_block_name (i+31) = 'walltime_start'
    io_block_name (i+32) = 'walltime_stop'
    io_block_name (i+ov) = 'disabled'

    track_ejected_particles = .FALSE.
    dump_absorption = .FALSE.
    averaged_var_block = 0
    new_style_io_block = .FALSE.
    n_io_blocks = 0
    rolling_restart_io_block = 0

  END SUBROUTINE io_deck_initialise



  SUBROUTINE io_deck_finalise

    INTEGER :: i, io, iu, n_zeros_estimate, n_dumps
    REAL(num) :: dumps
#ifndef NO_IO
    CHARACTER(LEN=c_max_path_length) :: list_filename
#endif

    n_io_blocks = block_number
    block_number = 0

    IF (n_io_blocks > 0) THEN
      ! First pass
      IF (deck_state == c_ds_first) THEN
        IF (.NOT.new_style_io_block .AND. n_io_blocks /= 1) THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Cannot have multiple unnamed "output" blocks.'
            END DO
          END IF
          CALL abort_code(c_err_preset_element)
        END IF

        ALLOCATE(io_prefixes(n_io_blocks+1))
        nfile_prefixes = 1
        io_prefixes(1) = ''

        ALLOCATE(io_block_list(n_io_blocks))
        DO i = 1, n_io_blocks
          CALL init_io_block(io_block_list(i))
          IF (i == rolling_restart_io_block) THEN
            io_block_list(i)%rolling_restart = .TRUE.
            io_block_list(i)%dump_cycle = 1
            io_block_list(i)%restart = .TRUE.
            io_block_list(i)%prefix_index = 2
            nfile_prefixes = 2
            io_prefixes(2) = 'roll'
          END IF
        END DO
      ! Second pass
      ELSE
        dumps = 1.0_num

        DO i = 1, n_io_blocks
          IF (io_block_list(i)%disabled) CYCLE

          IF (io_block_list(i)%dt_average > t_end) THEN
            IF (rank == 0) THEN
              DO iu = 1, nio_units ! Print to stdout and to file
                io = io_units(iu)
                WRITE(io,*) '*** WARNING ***'
                WRITE(io,*) 'Averaging time is longer than t_end, will set', &
                    ' averaging time equal'
                WRITE(io,*) 'to t_end.'
              END DO
            END IF
            io_block_list(i)%dt_average = t_end
          END IF

          IF (io_block_list(i)%dump_cycle_first_index &
                > io_block_list(i)%dump_cycle) THEN
            IF (rank == 0) THEN
              DO iu = 1, nio_units ! Print to stdout and to file
                io = io_units(iu)
                WRITE(io,*) '*** WARNING ***'
                WRITE(io,*) '"dump_cycle_first_index" cannot be greater ', &
                    'than "dump_cycle"'
                WRITE(io,*) 'Resetting to zero.'
              END DO
            END IF
            io_block_list(i)%dump_cycle_first_index = 0
          END IF

          IF (io_block_list(i)%dt_snapshot > 0.0_num) THEN
            dumps = MAX(dumps, t_end / io_block_list(i)%dt_snapshot)
          END IF
        END DO

        n_dumps = FLOOR(dumps)

        IF (io_block%dump_last .AND. dumps - n_dumps > c_tiny) &
            n_dumps = n_dumps + 1

        IF (.NOT. io_block%dump_first) &
            n_dumps = n_dumps - 1

        n_zeros_estimate = MAX(n_zeros, FLOOR(LOG10(REAL(n_dumps))) + 1)

        IF (n_zeros_control > 0) use_accurate_n_zeros = .FALSE.

        IF (.NOT.use_accurate_n_zeros) THEN
          IF (n_zeros_control > 0 &
              .AND. n_zeros_estimate /= n_zeros_control) THEN
            n_zeros = n_zeros_estimate
            IF (n_zeros > n_zeros_control .AND. rank == 0) THEN
              DO iu = 1, nio_units ! Print to stdout and to file
                io = io_units(iu)
                WRITE(io,*) '*** WARNING ***'
                WRITE(io,'(A,I1,A)') ' Estimated value of n_zeros (', n_zeros, &
                    ') has been overidden by input deck'
              END DO
            END IF
            n_zeros = n_zeros_control
          ELSE IF (n_zeros_estimate > n_zeros) THEN
            n_zeros = n_zeros_estimate
            IF (rank == 0) THEN
              DO iu = 1, nio_units ! Print to stdout and to file
                io = io_units(iu)
                WRITE(io,*) '*** WARNING ***'
                WRITE(io,'(A,I1,A)') ' n_zeros changed to ', n_zeros, &
                    ' to accomodate requested number of snapshots'
              END DO
            END IF
          END IF
        END IF

        ALLOCATE(file_prefixes(nfile_prefixes))
        ALLOCATE(file_numbers(nfile_prefixes))
        DO i = 1,nfile_prefixes
          file_prefixes(i) = TRIM(io_prefixes(i))
          file_numbers(i) = 0
        END DO
        DEALLOCATE(io_prefixes)

#ifndef NO_IO
        ! Remove any left-over VisIt file lists
        IF (.NOT.ic_from_restart .AND. rank == 0) THEN
          list_filename = TRIM(ADJUSTL(data_dir)) // '/full.visit'
          OPEN(unit=lu, status='UNKNOWN', file=list_filename)
          CLOSE(unit=lu, status='DELETE')
          list_filename = TRIM(ADJUSTL(data_dir)) // '/normal.visit'
          OPEN(unit=lu, status='UNKNOWN', file=list_filename)
          CLOSE(unit=lu, status='DELETE')
          list_filename = TRIM(ADJUSTL(data_dir)) // '/restart.visit'
          OPEN(unit=lu, status='UNKNOWN', file=list_filename)
          CLOSE(unit=lu, status='DELETE')
        END IF
#endif
      END IF
    END IF

  END SUBROUTINE io_deck_finalise



  SUBROUTINE io_block_start

    io_block_done = .FALSE.
    got_name = .FALSE.
    got_dump_source_code = .FALSE.
    got_dump_input_decks = .FALSE.
    block_number = block_number + 1
    IF (deck_state /= c_ds_first .AND. block_number > 0) THEN
      io_block => io_block_list(block_number)
      IF (io_block%rolling_restart) THEN
        io_block_done(num_vars_to_dump+o6) = .TRUE.
        io_block_done(num_vars_to_dump+o7) = .TRUE.
        io_block_done(num_vars_to_dump+o8) = .TRUE.
      END IF
    END IF

  END SUBROUTINE io_block_start



  SUBROUTINE io_block_end

    INTEGER :: io, iu, mask
#ifndef NO_IO
    CHARACTER(LEN=c_max_path_length) :: list_filename
#endif

    IF (deck_state == c_ds_first) RETURN

    IF (io_block%disabled) RETURN

    mask = io_block%dumpmask(c_dump_ejected_particles)
    IF (mask /= c_io_none .AND. IAND(mask,c_io_never) == 0) &
        track_ejected_particles = .TRUE.

    IF (io_block%dumpmask(c_dump_absorption) /= c_io_never) THEN
      dump_absorption = .TRUE.
    END IF

    IF (.NOT. got_name) THEN
      IF (new_style_io_block) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Cannot mix old and new style output blocks.'
            WRITE(io,*) 'You can either have multiple, named output blocks ', &
                'or a single unnamed one.'
          END DO
        END IF
        CALL abort_code(c_err_bad_value)
      END IF
      io_block%name = 'normal'
    END IF

#ifndef NO_IO
    ! Delete any existing visit file lists
    IF (.NOT.ic_from_restart .AND. rank == 0) THEN
      list_filename = TRIM(ADJUSTL(data_dir)) // '/' &
          // TRIM(io_block%name) // '.visit'
      OPEN(unit=lu, status='UNKNOWN', file=list_filename)
      CLOSE(unit=lu, status='DELETE')
    END IF
#endif
    io_block%dumpmask(c_dump_jx) = IOR(io_block%dumpmask(c_dump_jx), c_io_field)
    io_block%dumpmask(c_dump_jy) = IOR(io_block%dumpmask(c_dump_jy), c_io_field)
    io_block%dumpmask(c_dump_jz) = IOR(io_block%dumpmask(c_dump_jz), c_io_field)

    IF (.NOT.got_dump_source_code) THEN
      IF (io_block%restart .OR. .NOT.new_style_io_block) &
          io_block%dump_source_code = .TRUE.
    END IF

    IF (.NOT.got_dump_input_decks) THEN
      IF (io_block%restart .OR. .NOT.new_style_io_block) &
          io_block%dump_input_decks = .TRUE.
    END IF

    CALL set_restart_dumpmasks

  END SUBROUTINE io_block_end



  FUNCTION io_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode, style_error
    INTEGER :: loop, elementselected, mask, fullmask = 0, mask_element
    INTEGER :: i, is, subset, n_list, io, iu
    INTEGER, ALLOCATABLE :: subsets(:)
    LOGICAL :: bad, found
    INTEGER, PARAMETER :: c_err_new_style_ignore = 1
    INTEGER, PARAMETER :: c_err_new_style_global = 2
    INTEGER, PARAMETER :: c_err_old_style_ignore = 3

    errcode = c_err_none
    IF (value == blank) RETURN

    IF (deck_state == c_ds_first) THEN
      IF (str_cmp(element, 'name')) new_style_io_block = .TRUE.
      IF (str_cmp(element, 'rolling_restart')) THEN
        IF (rolling_restart_io_block > 0) THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Cannot have multiple "rolling_restart" blocks.'
            END DO
          END IF
          CALL abort_code(c_err_preset_element)
        END IF
        rolling_restart_io_block = block_number
      END IF
      RETURN
    END IF

    IF (io_block%disabled) THEN
      errcode = c_err_none
      RETURN
    END IF

    errcode = c_err_unknown_element

    elementselected = 0

    DO loop = 1, io_block_elements
      IF (str_cmp(element, TRIM(ADJUSTL(io_block_name(loop)))) &
          .OR. str_cmp(element, TRIM(ADJUSTL(alternate_name(loop))))) THEN
        elementselected = loop
        EXIT
      END IF
    END DO

    IF (elementselected == 0) RETURN

    IF (io_block_done(elementselected)) THEN
      errcode = c_err_preset_element
      RETURN
    END IF
    io_block_done(elementselected) = .TRUE.
    errcode = c_err_none
    style_error = c_err_none

    SELECT CASE (elementselected-num_vars_to_dump)
    CASE(1)
      io_block%dt_snapshot = as_real_print(value, element, errcode)
      IF (io_block%dt_snapshot < 0.0_num) io_block%dt_snapshot = 0.0_num
    CASE(2)
      IF (new_style_io_block) THEN
        style_error = c_err_new_style_ignore
      ELSE
        full_dump_every = as_integer_print(value, element, errcode)
        IF (full_dump_every == 0) full_dump_every = 1
      END IF
    CASE(3)
      IF (new_style_io_block) THEN
        style_error = c_err_new_style_ignore
      ELSE
        restart_dump_every = as_integer_print(value, element, errcode)
        IF (restart_dump_every == 0) restart_dump_every = 1
      END IF
    CASE(4)
      IF (new_style_io_block) style_error = c_err_new_style_global
      force_first_to_be_restartable = as_logical_print(value, element, errcode)
    CASE(5)
      IF (new_style_io_block) style_error = c_err_new_style_global
      force_final_to_be_restartable = as_logical_print(value, element, errcode)
    CASE(6)
      IF (new_style_io_block) style_error = c_err_new_style_global
      use_offset_grid = as_logical_print(value, element, errcode)
    CASE(7)
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The "extended_io_file" option is no longer supported.'
          WRITE(io,*) 'Please use the "import" directive instead'
        END DO
      END IF
      CALL abort_code(c_err_unknown_element)
    CASE(8)
      io_block%dt_average = as_real_print(value, element, errcode)
    CASE(9)
      io_block%nstep_average = as_integer_print(value, element, errcode)
    CASE(10)
      io_block%nstep_snapshot = as_integer_print(value, element, errcode)
      IF (io_block%nstep_snapshot < 0) io_block%nstep_snapshot = 0
    CASE(11)
      io_block%dump_source_code = as_logical_print(value, element, errcode)
      got_dump_source_code = .TRUE.
    CASE(12)
      io_block%dump_input_decks = as_logical_print(value, element, errcode)
      got_dump_input_decks = .TRUE.
    CASE(13)
      io_block%dump_first = as_logical_print(value, element, errcode)
    CASE(14)
      io_block%dump_last = as_logical_print(value, element, errcode)
    CASE(15)
      IF (.NOT.new_style_io_block) style_error = c_err_old_style_ignore
      io_block%restart = as_logical_print(value, element, errcode)
    CASE(16)
      DO i = 1,block_number
        IF (TRIM(io_block_list(i)%name) == TRIM(value)) THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Output block "' // TRIM(value) &
                  // '" already defined.'
            END DO
          END IF
          CALL abort_code(c_err_preset_element)
        END IF
      END DO
      io_block%name = value
      got_name = .TRUE.
    CASE(17)
      io_block%time_start = as_real_print(value, element, errcode)
    CASE(18)
      io_block%time_stop = as_real_print(value, element, errcode)
    CASE(19)
      io_block%nstep_start = as_integer_print(value, element, errcode)
    CASE(20)
      io_block%nstep_stop = as_integer_print(value, element, errcode)
    CASE(21)
      IF (.NOT.new_style_io_block) style_error = c_err_old_style_ignore
      CALL get_allocated_array(value, io_block%dump_at_nsteps, errcode)
    CASE(22)
      IF (.NOT.new_style_io_block) style_error = c_err_old_style_ignore
      CALL get_allocated_array(value, io_block%dump_at_times, errcode)
    CASE(23)
      io_block%dump_cycle = as_integer_print(value, element, errcode)
    CASE(24)
      found = .FALSE.
      DO i = 1,nfile_prefixes
        IF (TRIM(io_prefixes(i)) == TRIM(value)) THEN
          found = .TRUE.
          io_block%prefix_index = i
          EXIT
        END IF
      END DO
      IF (.NOT.found) THEN
        nfile_prefixes = nfile_prefixes + 1
        io_prefixes(nfile_prefixes) = TRIM(value)
        io_block%prefix_index = nfile_prefixes
      END IF
    CASE(26)
      io_block%dump_cycle_first_index = &
          as_integer_print(value, element, errcode)
    CASE(27)
      filesystem = TRIM(value) // ':'
    CASE(28)
      io_block%dump_first_after_restart = &
          as_logical_print(value, element, errcode)
    CASE(29)
      IF (.NOT.new_style_io_block) style_error = c_err_old_style_ignore
      CALL get_allocated_array(value, io_block%dump_at_walltimes, errcode)
    CASE(30)
      io_block%walltime_interval = as_real_print(value, element, errcode)
    CASE(31)
      io_block%walltime_start = as_real_print(value, element, errcode)
    CASE(32)
      io_block%walltime_stop = as_real_print(value, element, errcode)
    CASE(ov)
      io_block%disabled = as_logical_print(value, element, errcode)
    END SELECT

    IF (style_error == c_err_old_style_ignore) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Element "' &
              // TRIM(ADJUSTL(io_block_name(elementselected))) &
              // '" not ', 'allowed in an unnamed output block.'
          WRITE(io,*) 'It has been ignored.'
          WRITE(io,*)
        END DO
      END IF
    ELSE IF (style_error == c_err_new_style_ignore) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Element "' &
              // TRIM(ADJUSTL(io_block_name(elementselected))) &
              // '" not ', 'allowed in a named output block.'
          WRITE(io,*) 'It has been ignored.'
          WRITE(io,*)
        END DO
      END IF
    ELSE IF (style_error == c_err_new_style_global) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Element "' &
              // TRIM(ADJUSTL(io_block_name(elementselected))) &
              // '" should be moved to ', 'an "output_global" block.'
          WRITE(io,*) 'Its value will be applied to all output blocks.'
          WRITE(io,*)
        END DO
      END IF
    END IF

    IF (elementselected > num_vars_to_dump) RETURN

    mask_element = elementselected
    ALLOCATE(subsets(n_subsets+1))
    CALL as_list(value, subsets, n_list, errcode)

    DO is = 1, n_list
      subset = subsets(is)
      IF (is == 1) THEN
        mask = subset
        fullmask = mask
      ELSE
        mask = IOR(subset_list(subset)%mask,fullmask)
      END IF

      ! If setting dumpmask for features which haven't been compiled
      ! in then issue a warning
#ifdef NO_PARTICLE_PROBES
      IF (mask_element == c_dump_probes &
          .AND. mask /= c_io_none .AND. IAND(mask,c_io_never) == 0) THEN
        errcode = c_err_pp_options_wrong
        extended_error_string = '-NO_DPARTICLE_PROBES'
        mask = c_io_never
      END IF
#endif

#ifdef PARTICLE_ID4
#define PARTICLE_ID
#endif
#ifndef PARTICLE_ID
      IF (mask_element == c_dump_part_id &
          .AND. mask /= c_io_none .AND. IAND(mask,c_io_never) == 0) THEN
        errcode = c_err_pp_options_missing
        extended_error_string = '-DPARTICLE_ID'
        mask = c_io_never
      END IF
#endif

      ! Setting some flags like species and average
      ! wastes memory if the parameters make no sense. Do sanity checking.

      IF (IAND(mask, c_io_species) /= 0 &
          .OR. IAND(mask, c_io_no_sum) == 0) THEN
        bad = .TRUE.
        ! Check for sensible per species variables
        IF (mask_element == c_dump_ekbar) bad = .FALSE.
        IF (mask_element == c_dump_ekflux) bad = .FALSE.
        IF (mask_element == c_dump_mass_density) bad = .FALSE.
        IF (mask_element == c_dump_charge_density) bad = .FALSE.
        IF (mask_element == c_dump_number_density) bad = .FALSE.
        IF (mask_element == c_dump_ppc) bad = .FALSE.
        IF (mask_element == c_dump_average_weight) bad = .FALSE.
        IF (mask_element == c_dump_temperature) bad = .FALSE.
        IF (mask_element == c_dump_jx) bad = .FALSE.
        IF (mask_element == c_dump_jy) bad = .FALSE.
        IF (mask_element == c_dump_jz) bad = .FALSE.
        IF (bad) THEN
          IF (rank == 0 .AND. IAND(mask, c_io_species) /= 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** WARNING ***'
              WRITE(io,*) 'Attempting to set per species property for "' &
                  // TRIM(element) // '" which'
              WRITE(io,*) 'does not support this property. Ignoring.'
            END DO
          END IF
          mask = IAND(mask, NOT(c_io_species))
          mask = IOR(mask, c_io_no_sum)
        END IF
      END IF

      IF (IAND(mask, c_io_averaged) /= 0) THEN
        bad = .TRUE.
        ! Check for sensible averaged variables
        IF (mask_element == c_dump_ex) bad = .FALSE.
        IF (mask_element == c_dump_ey) bad = .FALSE.
        IF (mask_element == c_dump_ez) bad = .FALSE.
        IF (mask_element == c_dump_bx) bad = .FALSE.
        IF (mask_element == c_dump_by) bad = .FALSE.
        IF (mask_element == c_dump_bz) bad = .FALSE.
        IF (mask_element == c_dump_jx) bad = .FALSE.
        IF (mask_element == c_dump_jy) bad = .FALSE.
        IF (mask_element == c_dump_jz) bad = .FALSE.
        IF (mask_element == c_dump_poynt_flux) bad = .FALSE.

        ! Unset 'no_sum' dumpmask for grid variables
        IF (.NOT.bad) mask = IAND(mask, NOT(c_io_no_sum))

        IF (mask_element == c_dump_ekbar) bad = .FALSE.
        IF (mask_element == c_dump_mass_density) bad = .FALSE.
        IF (mask_element == c_dump_charge_density) bad = .FALSE.
        IF (mask_element == c_dump_number_density) bad = .FALSE.
        IF (mask_element == c_dump_ppc) bad = .FALSE.
        IF (mask_element == c_dump_average_weight) bad = .FALSE.
        IF (mask_element == c_dump_temperature) bad = .FALSE.
        IF (mask_element == c_dump_ekflux) bad = .FALSE.
        IF (bad) THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** WARNING ***'
              WRITE(io,*) 'Attempting to set average property for "' &
                  // TRIM(element) // '" which'
              WRITE(io,*) 'does not support this property. Ignoring.'
            END DO
          END IF
          mask = IAND(mask, NOT(c_io_averaged))
        ELSE
          any_average = .TRUE.
          io_block%any_average = .TRUE.
          IF (IAND(mask, c_io_average_single) /= 0 .AND. num /= r4) THEN
            io_block%averaged_data(mask_element)%dump_single = .TRUE.
          END IF
          i = averaged_var_block(mask_element)
          IF (i /= 0 .AND. i /= block_number) THEN
            IF (rank == 0 .AND. .NOT.warning_printed) THEN
              DO iu = 1, nio_units ! Print to stdout and to file
                io = io_units(iu)
                WRITE(io,*) '*** WARNING ***'
                WRITE(io,*) 'Error occurred whilst assigning the averaging ', &
                    'dumpmask of the variable'
                WRITE(io,*) '"' // TRIM(io_block_name(mask_element)) &
                    // '" in ', 'output block number ', block_number
                WRITE(io,*) 'Only one average per variable can be computed'
                WRITE(io,*) 'If multiple were specified the first averaging ', &
                    'time will be used'
                WRITE(io,*)
              END DO
              warning_printed = .TRUE.
            END IF
          ELSE
            averaged_var_block(mask_element) = block_number
          END IF
        END IF
      END IF

      IF (is == 1) THEN
        io_block%dumpmask(mask_element) = mask
      ELSE
        subset_list(subset)%dumpmask(block_number,mask_element) = mask
      END IF
    END DO

    DEALLOCATE(subsets)

  END FUNCTION io_block_handle_element



  FUNCTION io_block_check() RESULT(errcode)

    INTEGER :: errcode, io, iu, i

    ! Just assume that anything not included except for the compulsory
    ! elements is not wanted
    errcode = c_err_none

    ! Other control parameters are optional
    i = num_vars_to_dump
    io_block_done(i+o1:i+o2) = .TRUE.
    io_block_done(i+o5:io_block_elements) = .TRUE.
    ! Averaging info not compulsory unless averaged variable selected
    IF (.NOT. any_average) io_block_done(i+o3:i+o4) = .TRUE.

    IF (.NOT. io_block_done(i+o3) .AND. .NOT. io_block_done(i+o4)) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Required output block element ', &
              TRIM(ADJUSTL(io_block_name(i+8))), &
              ' absent. Please create this entry in the input deck'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    ! Can't check the io_block if it hasn't been allocated.
    IF (n_io_blocks == 0) RETURN

    IF (io_block%dt_average > io_block%dt_snapshot) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Averaging time is longer than dt_snapshot, will set', &
              ' averaging time equal'
          WRITE(io,*) 'to dt_snapshot.'
        END DO
      END IF
      io_block%dt_average = io_block%dt_snapshot
    END IF

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
    io_block%average_time_start = -1.0_num
    io_block%nstep_snapshot = -1
    io_block%nstep_prev = 0
    io_block%nstep_first = 0
    io_block%nstep_average = -1
    io_block%restart = .FALSE.
    io_block%dump = .FALSE.
    io_block%any_average = .FALSE.
    io_block%dump_first = .TRUE.
    io_block%dump_last = .TRUE.
    io_block%dump_source_code = .FALSE.
    io_block%dump_input_decks = .FALSE.
    io_block%dump_first_after_restart = .FALSE.
    io_block%dumpmask = c_io_none
    io_block%time_start = -1.0_num
    io_block%time_stop  = HUGE(1.0_num)
    io_block%nstep_start = -1
    io_block%nstep_stop  = HUGE(1)
    io_block%dump_cycle  = HUGE(1)
    io_block%dump_cycle_first_index = 0
    io_block%prefix_index = 1
    io_block%rolling_restart = .FALSE.
    io_block%disabled = .FALSE.
    io_block%walltime_interval = -1.0_num
    io_block%walltime_prev = 0.0_num
    io_block%walltime_start = -1.0_num
    io_block%walltime_stop  = HUGE(1.0_num)
    NULLIFY(io_block%dump_at_nsteps)
    NULLIFY(io_block%dump_at_times)
    NULLIFY(io_block%dump_at_walltimes)
    DO i = 1, num_vars_to_dump
      io_block%averaged_data(i)%dump_single = .FALSE.
    END DO

  END SUBROUTINE init_io_block



  SUBROUTINE copy_io_block(io_block, io_block_copy)

    TYPE(io_block_type) :: io_block, io_block_copy
    INTEGER :: i

    io_block_copy = io_block
    NULLIFY(io_block%dump_at_nsteps)
    NULLIFY(io_block%dump_at_times)
    NULLIFY(io_block%dump_at_walltimes)
    DO i = 1, num_vars_to_dump
      io_block_copy%averaged_data(i) = io_block%averaged_data(i)
    END DO

  END SUBROUTINE copy_io_block



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
    io_block%dumpmask(c_dump_cpml_psi_exz) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_exz), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_eyz) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_eyz), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_bxz) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_bxz), c_io_restartable)
    io_block%dumpmask(c_dump_cpml_psi_byz) = &
        IOR(io_block%dumpmask(c_dump_cpml_psi_byz), c_io_restartable)

  END SUBROUTINE set_restart_dumpmasks

END MODULE deck_io_block
