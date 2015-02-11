MODULE diagnostics

  USE calc_df
  USE sdf
  USE deck
  USE dist_fn
  USE epoch_source_info
  USE iterators
  USE probes
  USE version_data
  USE setup
  USE random_generator
  USE strings
  USE window
  USE timer

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: output_routines, create_full_timestring
  PUBLIC :: cleanup_stop_files, check_for_stop_condition
  PUBLIC :: deallocate_file_list

  CHARACTER(LEN=*), PARAMETER :: stop_file = 'STOP'
  CHARACTER(LEN=*), PARAMETER :: stop_file_nodump = 'STOP_NODUMP'

  TYPE(sdf_file_handle) :: sdf_handle
  INTEGER(i8), ALLOCATABLE :: species_offset(:)
  INTEGER(i8), ALLOCATABLE :: ejected_offset(:)
  LOGICAL :: reset_ejected, done_species_offset_init, done_subset_init
  LOGICAL :: restart_flag, dump_source_code, dump_input_decks
  LOGICAL :: dump_field_grid
  LOGICAL, ALLOCATABLE :: dump_point_grid(:)
  LOGICAL, ALLOCATABLE, SAVE :: prefix_first_call(:)
  INTEGER :: isubset
  INTEGER, DIMENSION(num_vars_to_dump) :: iomask
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: iodumpmask, dumped_skip_dir
  INTEGER, SAVE :: last_step = -1
  INTEGER, SAVE :: sdf_max_string_length, max_string_length

  ! Data structures for tracking the list of strings written to '.visit' files
  TYPE string_entry
    CHARACTER(LEN=string_length) :: text
    TYPE(string_entry), POINTER :: next
  END TYPE string_entry

  TYPE string_list
    TYPE(string_entry), POINTER :: head, tail
    INTEGER :: count
  END TYPE string_list

  TYPE(string_list), POINTER :: file_list(:)

  INTERFACE write_particle_variable
    MODULE PROCEDURE &
#if defined(PARTICLE_ID4) || defined(PARTICLE_DEBUG)
        write_particle_variable_i4, &
#endif
#if defined(PARTICLE_ID)
        write_particle_variable_i8, &
#endif
        write_particle_variable_num
  END INTERFACE write_particle_variable

CONTAINS

  SUBROUTINE output_routines(step, force_write)   ! step = step index

    INTEGER, INTENT(INOUT) :: step
    LOGICAL, INTENT(IN), OPTIONAL :: force_write
    CHARACTER(LEN=22) :: filename_fmt
    CHARACTER(LEN=5+n_zeros+c_id_length) :: filename
    CHARACTER(LEN=6+data_dir_max_length+n_zeros+c_id_length) :: full_filename
    CHARACTER(LEN=c_max_string_length) :: dump_type, temp_name
    CHARACTER(LEN=c_id_length) :: temp_block_id
    REAL(num) :: elapsed_time
    REAL(num), DIMENSION(:), ALLOCATABLE :: x_reduced
    REAL(num), DIMENSION(:), ALLOCATABLE :: array
    INTEGER :: code, i, ii, io, ispecies, iprefix, mask, rn, dir, dumped
    INTEGER :: random_state(4)
    INTEGER, ALLOCATABLE :: random_states_per_proc(:)
    INTEGER, DIMENSION(c_ndims) :: dims
    LOGICAL :: convert, force, any_written, restart_id, print_arrays
    LOGICAL, SAVE :: first_call = .TRUE.
    TYPE(particle_species), POINTER :: species
    TYPE(subset), POINTER :: sub
    CHARACTER(LEN=16) :: timestring, eta_timestring
    CHARACTER(LEN=1), DIMENSION(3) :: dim_tags = (/'x', 'y', 'z'/)
    CHARACTER(LEN=5), DIMENSION(6) :: dir_tags = &
        (/'x_max', 'y_max', 'z_max', 'x_min', 'y_min', 'z_min'/)
    INTEGER, DIMENSION(6) :: fluxdir = &
        (/c_dir_x, c_dir_y, c_dir_z, -c_dir_x, -c_dir_y, -c_dir_z/)

#ifdef NO_IO
    RETURN
#endif

    timer_walltime = -1.0_num
    IF (step /= last_step) THEN
      last_step = step
      IF (rank == 0 .AND. stdout_frequency > 0 &
          .AND. MOD(step, stdout_frequency) == 0) THEN
        timer_walltime = MPI_WTIME()
        elapsed_time = timer_walltime - walltime_start
        CALL create_timestring(elapsed_time, timestring)
        IF (print_eta_string) THEN
          eta_timestring = ''
          IF (time .GT. 0.0_num) THEN
            elapsed_time = (t_end - time) * elapsed_time / time
            CALL create_timestring(elapsed_time, eta_timestring)
          ENDIF
          WRITE(*, '(''Time'', g14.6, '' and iteration'', i9, '' after'', &
              & a, ''ETA'',a)') time, step, timestring, eta_timestring
        ELSE
          WRITE(*, '(''Time'', g20.12, '' and iteration'', i12, '' after'', &
              & a)') time, step, timestring
        ENDIF
      ENDIF
    ENDIF

    IF (n_io_blocks <= 0) RETURN

    IF (first_call) THEN
      ALLOCATE(dumped_skip_dir(c_ndims,n_subsets))
      ALLOCATE(file_list(n_io_blocks+2))
      ALLOCATE(prefix_first_call(SIZE(file_prefixes)))
      prefix_first_call = first_call
      DO i = 1,n_io_blocks+2
        file_list(i)%count = 0
      ENDDO
      first_call = .FALSE.
      ! Setting a large output buffer for point data can often make
      ! output much faster.
      ! The default value is set in deck_io_global_block.F90
      CALL sdf_set_point_array_size(sdf_buffer_size)
      sdf_max_string_length = sdf_get_max_string_length()
      max_string_length = MIN(sdf_max_string_length, c_max_string_length)
    ENDIF

    force = .FALSE.
    IF (PRESENT(force_write)) force = force_write

    dims = (/nx_global/)

    reset_ejected = .FALSE.
    any_written = .FALSE.

    IF (n_species > 0) THEN
      ALLOCATE(dump_point_grid(n_species))
    ELSE
      ALLOCATE(dump_point_grid(1))
    ENDIF

    DO iprefix = 1,SIZE(file_prefixes)
      CALL io_test(iprefix, step, print_arrays, force, prefix_first_call)

      IF (.NOT.print_arrays) CYCLE

      IF (.NOT.any_written) THEN
        ALLOCATE(array(-2:nx+3))
        CALL create_subtypes
        any_written = .TRUE.
        IF (timer_collect) THEN
          IF (timer_walltime < 0.0_num) THEN
            CALL timer_start(c_timer_io)
          ELSE
            CALL timer_start(c_timer_io, .TRUE.)
          ENDIF
        ENDIF
      ENDIF

      ! Allows a maximum of 10^999 output dumps, should be enough for anyone
      ! (feel free to laugh when this isn't the case)
      WRITE(filename_fmt, '(''(a, i'', i3.3, ''.'', i3.3, '', ".sdf")'')') &
          n_zeros, n_zeros
      WRITE(filename, filename_fmt) TRIM(file_prefixes(iprefix)), &
          file_numbers(iprefix)
      full_filename = TRIM(filesystem) &
          // TRIM(data_dir) // '/' // TRIM(filename)

      ! Always dump the variables with the 'Every' attribute
      code = c_io_always

      ! Only dump variables with the 'FULL' attributre on full dump intervals
      IF (MOD(file_numbers(iprefix), full_dump_every) == 0 &
          .AND. full_dump_every > -1) code = IOR(code, c_io_full)
      IF (restart_flag) code = IOR(code, c_io_restartable)
      dump_field_grid = .FALSE.

      ! open the file
      CALL sdf_open(sdf_handle, full_filename, comm, c_sdf_write)
      CALL sdf_set_string_length(sdf_handle, c_max_string_length)
      CALL sdf_write_header(sdf_handle, 'Epoch1d', 1, step, time, &
          restart_flag, jobid)
      CALL sdf_write_run_info(sdf_handle, c_version, c_revision, c_minor_rev, &
          c_commit_id, epoch_bytes_checksum, c_compile_machine, &
          c_compile_flags, defines, c_compile_date, run_date)
      CALL sdf_write_cpu_split(sdf_handle, 'cpu_rank', 'CPUs/Original rank', &
          cell_x_max)

      file_numbers(iprefix) = file_numbers(iprefix) + 1

      IF (restart_flag) THEN
        CALL sdf_write_srl(sdf_handle, 'dt', 'Time increment', dt)
        CALL sdf_write_srl(sdf_handle, 'dt_plasma_frequency', &
            'Plasma frequency timestep restriction', dt_plasma_frequency)
        IF (move_window .AND. window_started) THEN
          CALL sdf_write_srl(sdf_handle, 'window_shift_fraction', &
              'Window Shift Fraction', window_shift_fraction)
        ENDIF

        DO io = 1, n_io_blocks
          CALL sdf_write_srl(sdf_handle, &
              'time_prev/'//TRIM(io_block_list(io)%name), &
              'time_prev/'//TRIM(io_block_list(io)%name), &
              io_block_list(io)%time_prev)
          CALL sdf_write_srl(sdf_handle, &
              'nstep_prev/'//TRIM(io_block_list(io)%name), &
              'nstep_prev/'//TRIM(io_block_list(io)%name), &
              io_block_list(io)%nstep_prev)
        ENDDO

        DO ispecies = 1, n_species
          species => io_list(ispecies)
          CALL sdf_write_srl(sdf_handle, 'nppc/' // TRIM(species%name), &
              'Particles/Particles Per Cell/' // TRIM(species%name), &
              species%npart_per_cell)
        ENDDO

        CALL sdf_write_srl(sdf_handle, 'file_prefixes', &
            'Output File Stem Names', file_prefixes)

        CALL sdf_write_srl(sdf_handle, 'file_numbers', &
            'Output File Sequence Numbers', file_numbers)

        IF (need_random_state) THEN
          CALL get_random_state(random_state)
          ALLOCATE(random_states_per_proc(4*nproc))
          CALL MPI_GATHER(random_state, 4, MPI_INTEGER, &
              random_states_per_proc, 4, MPI_INTEGER, 0, comm, errcode)
          CALL sdf_write_srl(sdf_handle, 'random_states', 'Random States', &
              random_states_per_proc)
          DEALLOCATE(random_states_per_proc)
        ENDIF
      ENDIF

      iomask = iodumpmask(1,:)

      CALL write_field(c_dump_ex, code, 'ex', 'Electric Field/Ex', 'V/m', &
          c_stagger_ex, ex)
      CALL write_field(c_dump_ey, code, 'ey', 'Electric Field/Ey', 'V/m', &
          c_stagger_ey, ey)
      CALL write_field(c_dump_ez, code, 'ez', 'Electric Field/Ez', 'V/m', &
          c_stagger_ez, ez)

      CALL write_field(c_dump_bx, code, 'bx', 'Magnetic Field/Bx', 'T', &
          c_stagger_bx, bx)
      CALL write_field(c_dump_by, code, 'by', 'Magnetic Field/By', 'T', &
          c_stagger_by, by)
      CALL write_field(c_dump_bz, code, 'bz', 'Magnetic Field/Bz', 'T', &
          c_stagger_bz, bz)

      CALL write_field(c_dump_jx, code, 'jx', 'Current/Jx', 'A/m^2', &
          c_stagger_jx, jx)
      CALL write_field(c_dump_jy, code, 'jy', 'Current/Jy', 'A/m^2', &
          c_stagger_jy, jy)
      CALL write_field(c_dump_jz, code, 'jz', 'Current/Jz', 'A/m^2', &
          c_stagger_jz, jz)

      IF (cpml_boundaries) THEN
        CALL sdf_write_srl(sdf_handle, 'boundary_thickness', &
            'Boundary thickness', cpml_thickness)

        CALL write_field(c_dump_cpml_psi_eyx, code, 'cpml_psi_eyx', &
            'CPML/Ey_x', 'A/m^2', c_stagger_cell_centre, cpml_psi_eyx)
        CALL write_field(c_dump_cpml_psi_ezx, code, 'cpml_psi_ezx', &
            'CPML/Ez_x', 'A/m^2', c_stagger_cell_centre, cpml_psi_ezx)
        CALL write_field(c_dump_cpml_psi_byx, code, 'cpml_psi_byx', &
            'CPML/By_x', 'A/m^2', c_stagger_cell_centre, cpml_psi_byx)
        CALL write_field(c_dump_cpml_psi_bzx, code, 'cpml_psi_bzx', &
            'CPML/Bz_x', 'A/m^2', c_stagger_cell_centre, cpml_psi_bzx)
      ENDIF

      IF (n_subsets > 0) THEN
        DO i = 1, n_species
          CALL create_empty_partlist(io_list_data(i)%attached_list)
        ENDDO
      ENDIF

      DO isubset = 1, n_subsets + 1
        done_species_offset_init = .FALSE.
        done_subset_init = .FALSE.
        dump_point_grid = .FALSE.
        IF (isubset > 1) io_list => io_list_data
        iomask = iodumpmask(isubset,:)

#ifndef PER_SPECIES_WEIGHT
        CALL write_particle_variable(c_dump_part_weight, code, 'Weight', '', &
            it_output_real)
#else
        IF (IAND(iomask(c_dump_part_weight), code) /= 0) THEN
          CALL build_species_subset

          DO ispecies = 1, n_species
            species => io_list(ispecies)
            IF (IAND(species%dumpmask, code) /= 0 &
                .OR. IAND(code, c_io_restartable) /= 0) THEN
              CALL sdf_write_srl(sdf_handle, 'weight/' // TRIM(species%name), &
                  'Particles/Weight/' // TRIM(species%name), species%weight)
            ENDIF
          ENDDO
        ENDIF
#endif
        CALL write_particle_variable(c_dump_part_px, code, &
            'Px', 'kg.m/s', it_output_real)
        CALL write_particle_variable(c_dump_part_py, code, &
            'Py', 'kg.m/s', it_output_real)
        CALL write_particle_variable(c_dump_part_pz, code, &
            'Pz', 'kg.m/s', it_output_real)

        CALL write_particle_variable(c_dump_part_vx, code, &
            'Vx', 'm/s', it_output_real)
        CALL write_particle_variable(c_dump_part_vy, code, &
            'Vy', 'm/s', it_output_real)
        CALL write_particle_variable(c_dump_part_vz, code, &
            'Vz', 'm/s', it_output_real)

        CALL write_particle_variable(c_dump_part_charge, code, &
            'Q', 'C', it_output_real)
        CALL write_particle_variable(c_dump_part_mass, code, &
            'Mass', 'kg', it_output_real)
        CALL write_particle_variable(c_dump_part_ek, code, &
            'Ek', 'J', it_output_real)
        CALL write_particle_variable(c_dump_part_rel_mass, code, &
            'Relativistic Mass', 'kg', it_output_real)
        CALL write_particle_variable(c_dump_part_gamma, code, &
            'Gamma', '', it_output_real)
#ifdef PARTICLE_DEBUG
        CALL write_particle_variable(c_dump_part_proc, code, &
            'Processor', '', it_output_integer4)
        CALL write_particle_variable(c_dump_part_proc0, code, &
            'Processor_at_t0', '', it_output_integer4)
#endif
#if defined(PARTICLE_ID)
        CALL write_particle_variable(c_dump_part_id, code, &
            'ID', '#', it_output_integer8)
#elif defined(PARTICLE_ID4)
        CALL write_particle_variable(c_dump_part_id, code, &
            'ID', '#', it_output_integer4)
#endif
#ifdef PHOTONS
        CALL write_particle_variable(c_dump_part_opdepth, code, &
            'Optical depth', '', it_output_real)
        CALL write_particle_variable(c_dump_part_qed_energy, code, &
            'QED energy', 'J', it_output_real)
#ifdef TRIDENT_PHOTONS
        CALL write_particle_variable(c_dump_part_opdepth_tri, code, &
            'Trident Depth', '', it_output_real)
#endif
#endif
        CALL write_particle_grid(code)

        ! These are derived variables from the particles
        CALL write_nspecies_field(c_dump_ekbar, code, &
            'ekbar', 'EkBar', 'J', &
            c_stagger_cell_centre, calc_ekbar, array)

        CALL write_nspecies_field(c_dump_mass_density, code, &
            'mass_density', 'Mass_Density', 'kg/m^3', &
            c_stagger_cell_centre, calc_mass_density, array)

        CALL write_nspecies_field(c_dump_charge_density, code, &
            'charge_density', 'Charge_Density', 'C/m^3', &
            c_stagger_cell_centre, calc_charge_density, array)

        CALL write_nspecies_field(c_dump_number_density, code, &
            'number_density', 'Number_Density', '1/m^3', &
            c_stagger_cell_centre, calc_number_density, array)

        CALL write_nspecies_field(c_dump_temperature, code, &
            'temperature', 'Temperature', 'K', &
            c_stagger_cell_centre, calc_temperature, array)

        CALL write_nspecies_field(c_dump_jx, code, &
            'jx', 'Jx', 'A/m^2', &
            c_stagger_cell_centre, calc_per_species_jx, array)

        CALL write_nspecies_field(c_dump_jy, code, &
            'jy', 'Jy', 'A/m^2', &
            c_stagger_cell_centre, calc_per_species_jy, array)

        CALL write_nspecies_field(c_dump_jz, code, &
            'jz', 'Jz', 'A/m^2', &
            c_stagger_cell_centre, calc_per_species_jz, array)

        CALL write_nspecies_flux(c_dump_ekflux, code, &
            'ekflux', 'EkFlux', 'W/m^2', &
            c_stagger_cell_centre, calc_ekflux, array, fluxdir, dir_tags)

        CALL write_field_flux(c_dump_poynt_flux, code, &
            'poynt_flux', 'Poynting Flux', 'W/m^2', &
            c_stagger_cell_centre, calc_poynt_flux, array, fluxdir(1:3), &
            dim_tags)

        IF (isubset /= 1) THEN
          DO i = 1, n_species
            CALL append_partlist(species_list(i)%attached_list, &
                io_list(i)%attached_list)
          ENDDO
          DO i = 1, n_species
            CALL create_empty_partlist(io_list(i)%attached_list)
          ENDDO
        ENDIF
      ENDDO

      io_list => species_list
      iomask = iodumpmask(1,:)
      mask = iomask(c_dump_grid)
      restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
      convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

      ! Write the cartesian mesh
      IF (IAND(mask, code) /= 0) dump_field_grid = .TRUE.
      IF (IAND(mask, c_io_never) /= 0) dump_field_grid = .FALSE.

      IF (dump_field_grid) THEN
        IF (.NOT. use_offset_grid) THEN
          CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid', 'Grid/Grid', &
              xb_global, convert)
        ELSE
          CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid', 'Grid/Grid', &
              xb_offset_global, convert)
          CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid_full', &
              'Grid/Grid_Full', xb_global, convert)
        ENDIF
      ENDIF

      dumped_skip_dir = 0
      dumped = 1

      DO io = 1, n_subsets
        sub => subset_list(io)
        IF (.NOT.sub%dump_field_grid) CYCLE

        DO i = 1, io - 1
          dumped = dumped + SUM(dumped_skip_dir(:,i) - sub%skip_dir)
        ENDDO
        IF (dumped == 0) CYCLE
        dumped = 0

        dumped_skip_dir(:,io) = sub%skip_dir

        dir = 1
        rn = sub%n_local(dir) + 1
        ALLOCATE(x_reduced(rn))

        ii = sub%n_start(dir) + 1
        DO i = 1, rn
          x_reduced(i) = xb_global(ii)
          ii = ii + sub%skip_dir(dir)
        ENDDO
        x_reduced(rn) = xb_global(nx+1)

        IF (.NOT. use_offset_grid) THEN
          temp_block_id = 'grid/r_' // TRIM(sub%name)
          temp_name = 'Grid/Reduced_' // TRIM(sub%name)

          CALL check_name_length('subset', &
              'Grid/Reduced_' // TRIM(sub%name))

          CALL sdf_write_srl_plain_mesh(sdf_handle, TRIM(temp_block_id), &
              TRIM(temp_name), x_reduced, convert)
        ELSE
          temp_block_id = 'grid_full/r_' // TRIM(sub%name)
          temp_name = 'Grid_Full/Reduced_' // TRIM(sub%name)

          CALL check_name_length('subset', &
              'Grid_Full/Reduced_' // TRIM(sub%name))

          CALL sdf_write_srl_plain_mesh(sdf_handle, TRIM(temp_block_id), &
              TRIM(temp_name), x_reduced, convert)

          temp_block_id = 'grid/r_' // TRIM(sub%name)
          temp_name = 'Grid/Reduced_' // TRIM(sub%name)

          CALL check_name_length('subset', &
              'Grid/Reduced_' // TRIM(sub%name))

          dir = 1
          rn = sub%n_local(dir)
          ALLOCATE(x_reduced(rn))

          ii = sub%n_start(dir) + 1
          DO i = 1, rn
            x_reduced(i) = xb_offset_global(ii)
            ii = ii + sub%skip_dir(dir)
          ENDDO

          CALL sdf_write_srl_plain_mesh(sdf_handle, TRIM(temp_block_id), &
              TRIM(temp_name), x_reduced, convert)
        ENDIF

        DEALLOCATE(x_reduced)
        sub%dump_field_grid = .FALSE.
      ENDDO

      IF (IAND(iomask(c_dump_dist_fns), code) /= 0) THEN
        CALL write_dist_fns(sdf_handle, code)
      ENDIF

#ifndef NO_PARTICLE_PROBES
      IF (IAND(iomask(c_dump_probes), code) /= 0) THEN
        CALL write_probes(sdf_handle, code)
      ENDIF
#endif

      IF (dump_input_decks) CALL write_input_decks(sdf_handle)
      IF (dump_source_code .AND. epoch_bytes_len > 0) &
          CALL write_source_info(sdf_handle)

      IF (IAND(iomask(c_dump_absorption), code) /= 0) THEN
        CALL MPI_ALLREDUCE(laser_absorb_local, laser_absorbed, 1, mpireal, &
            MPI_SUM, comm, errcode)
        CALL MPI_ALLREDUCE(laser_inject_local, laser_injected, 1, mpireal, &
            MPI_SUM, comm, errcode)
        IF (laser_injected > 0.0_num) THEN
          laser_absorbed = laser_absorbed / laser_injected
        ELSE
          laser_absorbed = 0.0_num
        ENDIF
        CALL sdf_write_srl(sdf_handle, 'laser_enTotal', &
            'Absorption/Total Laser Energy Injected (J)', laser_injected)
        CALL sdf_write_srl(sdf_handle, 'abs_frac', &
            'Absorption/Fraction of Laser Energy Absorbed (%)', laser_absorbed)
      ENDIF

      IF (IAND(iomask(c_dump_total_energy_sum), code) /= 0) THEN
        CALL calc_total_energy_sum

        CALL sdf_write_srl(sdf_handle, 'total_particle_energy', &
            'Total Particle Energy in Simulation (J)', total_particle_energy)
        CALL sdf_write_srl(sdf_handle, 'total_field_energy', &
            'Total Field Energy in Simulation (J)', total_field_energy)
      ENDIF

      ! close the file
      CALL sdf_close(sdf_handle)

      IF (rank == 0) THEN
        DO io = 1, n_io_blocks
          IF (io_block_list(io)%dump) THEN
            dump_type = TRIM(io_block_list(io)%name)
            CALL append_filename(dump_type, filename, io)
          ENDIF
        ENDDO
        IF (IAND(code, c_io_restartable) /= 0) THEN
          dump_type = 'restart'
          CALL append_filename(dump_type, filename, n_io_blocks+1)
        ENDIF
        IF (IAND(code, c_io_full) /= 0) THEN
          dump_type = 'full'
          CALL append_filename(dump_type, filename, n_io_blocks+2)
        ENDIF
        IF (iprefix > 1) dump_type = TRIM(file_prefixes(iprefix))
        WRITE(stat_unit, '(''Wrote '', a7, '' dump number'', i5, '' at time'', &
          & g20.12, '' and iteration'', i7)') dump_type, &
          file_numbers(iprefix)-1, time, step
        CALL flush_stat_file()
      ENDIF
    ENDDO

    DEALLOCATE(dump_point_grid)

    IF (.NOT.any_written) RETURN

    DEALLOCATE(array)
    IF (ALLOCATED(species_offset))  DEALLOCATE(species_offset)
    IF (ALLOCATED(ejected_offset))  DEALLOCATE(ejected_offset)
    CALL free_subtypes()

    IF (reset_ejected) THEN
      DO i = 1, n_species
        CALL destroy_partlist(ejected_list(i)%attached_list)
      ENDDO
    ENDIF

    IF (timer_collect) CALL timer_stop(c_timer_io)

  END SUBROUTINE output_routines



  SUBROUTINE check_name_length(shorten, string)

    CHARACTER(LEN=*), INTENT(IN) :: shorten, string
    CHARACTER(LEN=c_max_string_length) :: len_string
    INTEGER :: length

    length = LEN_TRIM(string)

    IF (length > max_string_length) THEN
      IF (rank == 0) THEN
        CALL integer_as_string(length, len_string)
        PRINT*, '*** WARNING ***'
        PRINT*, 'Output block name ', TRIM(string), ' is truncated.'
        IF (length > c_max_string_length &
            .AND. length > sdf_max_string_length) THEN
          PRINT*, 'Either shorten the ', TRIM(shorten), ' name or increase ', &
              'the size of "c_max_string_length" ', 'in both EPOCH and the ', &
              'SDF library to at least ', TRIM(len_string)
        ELSE IF (length > sdf_max_string_length) THEN
          PRINT*, 'Either shorten the ', TRIM(shorten), ' name or increase ', &
              'the size of "c_max_string_length" ', 'in the ', &
              'SDF library to at least ', TRIM(len_string)
        ELSE
          PRINT*, 'Either shorten the ', TRIM(shorten), ' name or increase ', &
              'the size of "c_max_string_length" ', 'to at least ', &
              TRIM(len_string)
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE check_name_length



  SUBROUTINE append_filename(listname, filename, list_index)
    ! This routine updates a list of each output type (eg. full, dump, normal)
    ! that can be passed to VisIt to filter the output.

    CHARACTER(LEN=*), INTENT(IN) :: listname, filename
    INTEGER, INTENT(IN) :: list_index
    CHARACTER(LEN=data_dir_max_length+c_id_length+8) :: listfile
    TYPE(string_list), POINTER :: list
    TYPE(string_entry), POINTER :: lcur
    INTEGER :: ierr, i
    LOGICAL :: exists

    list => file_list(list_index)

    listfile = TRIM(data_dir) // '/' // TRIM(listname) // '.visit'
    INQUIRE(file=listfile, exist=exists)

    IF (list%count == 0) THEN
      ALLOCATE(list%head)
      list%tail => list%head
      IF (ic_from_restart .AND. exists) CALL setup_file_list(listfile, list)
    ELSE
      ALLOCATE(list%tail%next)
      list%tail => list%tail%next
    ENDIF

    IF (list%count > 0) THEN
      lcur  => list%head
      IF (TRIM(lcur%text) == TRIM(filename)) RETURN
      DO i = 2,list%count
        lcur  => lcur%next
        IF (TRIM(lcur%text) == TRIM(filename)) RETURN
      ENDDO
    ENDIF

    list%count = list%count + 1
    list%tail%text = TRIM(filename)

    IF (exists) THEN
      OPEN(unit=lu, status='OLD', position='APPEND', file=listfile, iostat=ierr)
    ELSE
      OPEN(unit=lu, status='NEW', file=listfile, iostat=errcode)
    ENDIF

    WRITE(lu,'(a)') TRIM(filename)
    CLOSE(lu)

  END SUBROUTINE append_filename



  SUBROUTINE setup_file_list(listfile, list)
    ! This routine initialises the file list contained in the given
    ! *.visit file which is required when restarting.

    CHARACTER(LEN=*), INTENT(IN) :: listfile
    TYPE(string_list), POINTER :: list
    CHARACTER(LEN=string_length) :: str
    INTEGER :: ierr

    OPEN(unit=lu, status='OLD', file=listfile, iostat=ierr)
    DO
      READ(lu, '(A)', iostat=ierr) str
      IF (ierr < 0) EXIT
      list%tail%text = TRIM(str)
      list%count = list%count + 1
      ALLOCATE(list%tail%next)
      list%tail => list%tail%next
    ENDDO
    CLOSE(lu)

  END SUBROUTINE setup_file_list



  SUBROUTINE deallocate_file_list

    INTEGER :: i, n, nlist, stat
    TYPE(string_entry), POINTER :: current, next

    IF (ASSOCIATED(file_list)) THEN
      DO i = 1, n_io_blocks+2
        nlist = file_list(i)%count
        IF (nlist > 0) THEN
          current => file_list(i)%head
          DO n = 1, nlist
            next => current%next
            DEALLOCATE(current, STAT=stat)
            current => next
          ENDDO
        ENDIF
      ENDDO
      DEALLOCATE(file_list, STAT=stat)
    ENDIF
    DEALLOCATE(iodumpmask, STAT=stat)
    DEALLOCATE(dumped_skip_dir, STAT=stat)
    DEALLOCATE(prefix_first_call, STAT=stat)

  END SUBROUTINE deallocate_file_list



  SUBROUTINE io_test(iprefix, step, print_arrays, force, first_call)

    INTEGER, INTENT(IN) :: iprefix, step
    LOGICAL, INTENT(OUT) :: print_arrays
    LOGICAL, INTENT(IN) :: force
    LOGICAL, DIMENSION(:), INTENT(INOUT) :: first_call
    INTEGER :: id, io, is, nstep_next = 0
    REAL(num) :: t0, t1, time_first
    LOGICAL :: last_call, dump

    IF (.NOT.ALLOCATED(iodumpmask)) &
        ALLOCATE(iodumpmask(n_subsets+1,num_vars_to_dump))

    restart_flag = .FALSE.
    dump_source_code = .FALSE.
    dump_input_decks = .FALSE.
    print_arrays = .FALSE.
    iomask = c_io_none
    iodumpmask = c_io_none

    IF (time >= t_end .OR. step == nsteps) THEN
      last_call = .TRUE.
    ELSE
      last_call = .FALSE.
    ENDIF

    DO io = 1, n_io_blocks
      io_block_list(io)%dump = .FALSE.

      IF (io_block_list(io)%disabled) CYCLE
      IF (io_block_list(io)%prefix_index /= iprefix) CYCLE

      IF (last_call .AND. io_block_list(io)%dump_last) &
          io_block_list(io)%dump = .TRUE.
      IF (first_call(iprefix) .AND. io_block_list(io)%dump_first) &
          io_block_list(io)%dump = .TRUE.

      IF (force) THEN
        io_block_list(io)%dump = .TRUE.
        restart_flag = .TRUE.
      ENDIF

      IF (ASSOCIATED(io_block_list(io)%dump_at_nsteps)) THEN
        DO is = 1, SIZE(io_block_list(io)%dump_at_nsteps)
          IF (step >= io_block_list(io)%dump_at_nsteps(is)) THEN
            io_block_list(io)%dump = .TRUE.
            io_block_list(io)%dump_at_nsteps(is) = HUGE(1)
          ENDIF
        ENDDO
      ENDIF

      IF (ASSOCIATED(io_block_list(io)%dump_at_times)) THEN
        DO is = 1, SIZE(io_block_list(io)%dump_at_times)
          IF (time >= io_block_list(io)%dump_at_times(is)) THEN
            io_block_list(io)%dump = .TRUE.
            io_block_list(io)%dump_at_times(is) = HUGE(1.0_num)
          ENDIF
        ENDDO
      ENDIF

      ! Work out the time that the next dump will occur based on the
      ! current timestep
      t0 = HUGE(1.0_num)
      t1 = HUGE(1.0_num)
      IF (io_block_list(io)%dt_snapshot >= 0.0_num) &
          t0 = io_block_list(io)%time_prev + io_block_list(io)%dt_snapshot
      IF (io_block_list(io)%nstep_snapshot >= 0) THEN
        nstep_next = io_block_list(io)%nstep_prev &
            + io_block_list(io)%nstep_snapshot
        t1 = time + dt * (nstep_next - step)
      ENDIF

      IF (t0 < t1) THEN
        ! Next I/O dump based on dt_snapshot
        time_first = t0
        IF (io_block_list(io)%dt_snapshot > 0 .AND. time >= t0) THEN
          ! Store the most recent output time that qualifies
          DO
            t0 = io_block_list(io)%time_prev + io_block_list(io)%dt_snapshot
            IF (t0 > time) EXIT
            io_block_list(io)%time_prev = t0
          ENDDO
          dump = .TRUE.
          IF (dump .AND. time < io_block_list(io)%time_start)  dump = .FALSE.
          IF (dump .AND. time > io_block_list(io)%time_stop)   dump = .FALSE.
          IF (dump .AND. step < io_block_list(io)%nstep_start) dump = .FALSE.
          IF (dump .AND. step > io_block_list(io)%nstep_stop)  dump = .FALSE.
          IF (dump .AND. time < time_start)  dump = .FALSE.
          IF (dump .AND. time > time_stop)   dump = .FALSE.
          IF (dump .AND. step < nstep_start) dump = .FALSE.
          IF (dump .AND. step > nstep_stop)  dump = .FALSE.
          IF (dump) io_block_list(io)%dump = .TRUE.
        ENDIF
      ELSE
        ! Next I/O dump based on nstep_snapshot
        time_first = t1
        IF (io_block_list(io)%nstep_snapshot > 0 &
            .AND. step >= nstep_next) THEN
          ! Store the most recent output step that qualifies
          DO
            nstep_next = io_block_list(io)%nstep_prev &
                + io_block_list(io)%nstep_snapshot
            IF (nstep_next > step) EXIT
            io_block_list(io)%nstep_prev = nstep_next
          ENDDO
          dump = .TRUE.
          IF (dump .AND. time < io_block_list(io)%time_start)  dump = .FALSE.
          IF (dump .AND. time > io_block_list(io)%time_stop)   dump = .FALSE.
          IF (dump .AND. step < io_block_list(io)%nstep_start) dump = .FALSE.
          IF (dump .AND. step > io_block_list(io)%nstep_stop)  dump = .FALSE.
          IF (dump .AND. time < time_start)  dump = .FALSE.
          IF (dump .AND. time > time_stop)   dump = .FALSE.
          IF (dump .AND. step < nstep_start) dump = .FALSE.
          IF (dump .AND. step > nstep_stop)  dump = .FALSE.
          IF (dump) io_block_list(io)%dump = .TRUE.
        ENDIF
      ENDIF

      io_block_list(io)%average_time_start = &
          time_first - io_block_list(io)%average_time

      IF (io_block_list(io)%dump) THEN
        print_arrays = .TRUE.
        IF (io_block_list(io)%restart) restart_flag = .TRUE.
        IF (io_block_list(io)%dump_source_code) dump_source_code = .TRUE.
        IF (io_block_list(io)%dump_input_decks) dump_input_decks = .TRUE.
        IF (file_numbers(iprefix) > io_block_list(io)%dump_cycle) &
            file_numbers(iprefix) = io_block_list(io)%dump_cycle_first_index
        iomask = IOR(iomask, io_block_list(io)%dumpmask)
        IF (n_subsets /= 0) THEN
          DO is = 1, n_subsets
            iodumpmask(1+is,:) = &
                IOR(iodumpmask(1+is,:), subset_list(is)%dumpmask(io,:))
          ENDDO
        ENDIF
      ENDIF
    ENDDO

    DO io = 1, n_io_blocks
      IF (.NOT. io_block_list(io)%any_average) CYCLE

      IF (time >= io_block_list(io)%average_time_start) THEN
        DO id = 1, num_vars_to_dump
          IF (IAND(io_block_list(io)%dumpmask(id), c_io_averaged) /= 0) THEN
            CALL average_field(id, io_block_list(io)%averaged_data(id))
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    IF (MOD(file_numbers(1), restart_dump_every) == 0 &
        .AND. restart_dump_every > -1) restart_flag = .TRUE.
    IF (first_call(iprefix) .AND. force_first_to_be_restartable) &
        restart_flag = .TRUE.
    IF ( last_call .AND. force_final_to_be_restartable) restart_flag = .TRUE.
    IF (force) restart_flag = .TRUE.

    IF (.NOT.restart_flag .AND. .NOT.new_style_io_block) THEN
      dump_source_code = .FALSE.
      dump_input_decks = .FALSE.
    ENDIF

    IF (first_call(iprefix)) first_call(iprefix) = .FALSE.

    iodumpmask(1,:) = iomask

  END SUBROUTINE io_test



  SUBROUTINE average_field(ioutput, avg)

    INTEGER, INTENT(IN) :: ioutput
    TYPE(averaged_data_block) :: avg
    INTEGER :: n_species_local, ispecies
    REAL(num), DIMENSION(:), ALLOCATABLE :: array

    avg%real_time = avg%real_time + dt
    avg%started = .TRUE.

    n_species_local = avg%n_species + avg%species_sum

    IF (n_species_local <= 0) RETURN

    IF (avg%dump_single) THEN
      SELECT CASE(ioutput)
      CASE(c_dump_ex)
        avg%r4array(:,1) = avg%r4array(:,1) + REAL(ex * dt, r4)
      CASE(c_dump_ey)
        avg%r4array(:,1) = avg%r4array(:,1) + REAL(ey * dt, r4)
      CASE(c_dump_ez)
        avg%r4array(:,1) = avg%r4array(:,1) + REAL(ez * dt, r4)
      CASE(c_dump_bx)
        avg%r4array(:,1) = avg%r4array(:,1) + REAL(bx * dt, r4)
      CASE(c_dump_by)
        avg%r4array(:,1) = avg%r4array(:,1) + REAL(by * dt, r4)
      CASE(c_dump_bz)
        avg%r4array(:,1) = avg%r4array(:,1) + REAL(bz * dt, r4)
      CASE(c_dump_jx)
        avg%r4array(:,1) = avg%r4array(:,1) + REAL(jx(1-ng:nx+ng) * dt, r4)
      CASE(c_dump_jy)
        avg%r4array(:,1) = avg%r4array(:,1) + REAL(jy(1-ng:nx+ng) * dt, r4)
      CASE(c_dump_jz)
        avg%r4array(:,1) = avg%r4array(:,1) + REAL(jz(1-ng:nx+ng) * dt, r4)
      CASE(c_dump_ekbar)
        ALLOCATE(array(-2:nx+3))
        DO ispecies = 1, n_species_local
          CALL calc_ekbar(array, ispecies-avg%species_sum)
          avg%r4array(:,ispecies) = avg%r4array(:,ispecies) &
              + REAL(array * dt, r4)
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_mass_density)
        ALLOCATE(array(-2:nx+3))
        DO ispecies = 1, n_species_local
          CALL calc_mass_density(array, ispecies-avg%species_sum)
          avg%r4array(:,ispecies) = avg%r4array(:,ispecies) &
              + REAL(array * dt, r4)
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_charge_density)
        ALLOCATE(array(-2:nx+3))
        DO ispecies = 1, n_species_local
          CALL calc_charge_density(array, ispecies-avg%species_sum)
          avg%r4array(:,ispecies) = avg%r4array(:,ispecies) &
              + REAL(array * dt, r4)
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_number_density)
        ALLOCATE(array(-2:nx+3))
        DO ispecies = 1, n_species_local
          CALL calc_number_density(array, ispecies-avg%species_sum)
          avg%r4array(:,ispecies) = avg%r4array(:,ispecies) &
              + REAL(array * dt, r4)
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_temperature)
        ALLOCATE(array(-2:nx+3))
        DO ispecies = 1, n_species_local
          CALL calc_temperature(array, ispecies-avg%species_sum)
          avg%r4array(:,ispecies) = avg%r4array(:,ispecies) &
              + REAL(array * dt, r4)
        ENDDO
        DEALLOCATE(array)
      END SELECT
    ELSE
      SELECT CASE(ioutput)
      CASE(c_dump_ex)
        avg%array(:,1) = avg%array(:,1) + ex * dt
      CASE(c_dump_ey)
        avg%array(:,1) = avg%array(:,1) + ey * dt
      CASE(c_dump_ez)
        avg%array(:,1) = avg%array(:,1) + ez * dt
      CASE(c_dump_bx)
        avg%array(:,1) = avg%array(:,1) + bx * dt
      CASE(c_dump_by)
        avg%array(:,1) = avg%array(:,1) + by * dt
      CASE(c_dump_bz)
        avg%array(:,1) = avg%array(:,1) + bz * dt
      CASE(c_dump_jx)
        avg%array(:,1) = avg%array(:,1) + jx(1-ng:nx+ng) * dt
      CASE(c_dump_jy)
        avg%array(:,1) = avg%array(:,1) + jy(1-ng:nx+ng) * dt
      CASE(c_dump_jz)
        avg%array(:,1) = avg%array(:,1) + jz(1-ng:nx+ng) * dt
      CASE(c_dump_ekbar)
        ALLOCATE(array(-2:nx+3))
        DO ispecies = 1, n_species_local
          CALL calc_ekbar(array, ispecies-avg%species_sum)
          avg%array(:,ispecies) = avg%array(:,ispecies) + array * dt
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_mass_density)
        ALLOCATE(array(-2:nx+3))
        DO ispecies = 1, n_species_local
          CALL calc_mass_density(array, ispecies-avg%species_sum)
          avg%array(:,ispecies) = avg%array(:,ispecies) + array * dt
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_charge_density)
        ALLOCATE(array(-2:nx+3))
        DO ispecies = 1, n_species_local
          CALL calc_charge_density(array, ispecies-avg%species_sum)
          avg%array(:,ispecies) = avg%array(:,ispecies) + array * dt
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_number_density)
        ALLOCATE(array(-2:nx+3))
        DO ispecies = 1, n_species_local
          CALL calc_number_density(array, ispecies-avg%species_sum)
          avg%array(:,ispecies) = avg%array(:,ispecies) + array * dt
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_temperature)
        ALLOCATE(array(-2:nx+3))
        DO ispecies = 1, n_species_local
          CALL calc_temperature(array, ispecies-avg%species_sum)
          avg%array(:,ispecies) = avg%array(:,ispecies) + array * dt
        ENDDO
        DEALLOCATE(array)
      END SELECT
    ENDIF

  END SUBROUTINE average_field



  SUBROUTINE write_field(id, code, block_id, name, units, stagger, array)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: block_id, name, units
    INTEGER, INTENT(IN) :: stagger
    REAL(num), DIMENSION(-2:), INTENT(IN) :: array
    REAL(num), DIMENSION(:), ALLOCATABLE :: reduced
    INTEGER :: io, mask, dumped
    INTEGER :: i, ii, rnx
    INTEGER :: subtype, subarray, rsubtype, rsubarray
    INTEGER, DIMENSION(c_ndims) :: dims
    LOGICAL :: convert, dump_skipped, restart_id, normal_id, unaveraged_id
    CHARACTER(LEN=c_id_length) :: temp_block_id, temp_grid_id
    CHARACTER(LEN=c_max_string_length) :: temp_name
    TYPE(averaged_data_block), POINTER :: avg
    TYPE(io_block_type), POINTER :: iob
    TYPE(subset), POINTER :: sub

    mask = iomask(id)
    IF (IAND(mask, code) == 0) RETURN
    IF (IAND(mask, c_io_never) /= 0) RETURN

    ! This is a normal dump and normal output variable
    normal_id = IAND(IAND(code, mask), IOR(c_io_always, c_io_full)) /= 0
    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    ! The variable is either averaged or has snapshot specified
    unaveraged_id = IAND(mask, c_io_averaged) == 0 &
        .OR. IAND(mask, c_io_snapshot) /= 0

    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    dims = (/nx_global/)

    IF (convert) THEN
      subtype  = subtype_field_r4
      IF (id == c_dump_jx .OR. id == c_dump_jy .OR. id == c_dump_jz) THEN
        subarray = subarray_field_big_r4
      ELSE
        subarray = subarray_field_r4
      ENDIF
    ELSE
      subtype  = subtype_field
      IF (id == c_dump_jx .OR. id == c_dump_jy .OR. id == c_dump_jz) THEN
        subarray = subarray_field_big
      ELSE
        subarray = subarray_field
      ENDIF
    ENDIF

    ! Output unaveraged data if:
    !  1. This is a restart dump and a restart variable
    !  2. This is a normal dump and the variable
    !      a) is not averaged
    !      b) is averaged and also has snapshot specified

    dumped_skip_dir = 0
    dumped = 1
    dump_skipped = .FALSE.

    DO io = 1, n_subsets
      IF (IAND(iodumpmask(io+1,id), code) == 0) CYCLE

      sub => subset_list(io)
      IF (.NOT.sub%skip) CYCLE

      ! This should prevent a reduced variable from being dumped multiple
      ! times in the same output file
      DO i = 1, io - 1
        dumped = dumped + SUM(dumped_skip_dir(:,i) - sub%skip_dir)
      ENDDO
      IF (dumped == 0) CYCLE
      dumped = 0

      dumped_skip_dir(:,io) = sub%skip_dir

      rnx = sub%n_local(1)

      ALLOCATE(reduced(rnx))

      ii = sub%n_start(1) + 1
      DO i = 1, rnx
        reduced(i) = array(ii)
        ii = ii + sub%skip_dir(1)
      ENDDO

      IF (convert) THEN
        rsubtype  = sub%subtype_r4
        rsubarray = sub%subarray_r4
      ELSE
        rsubtype  = sub%subtype
        rsubarray = sub%subarray
      ENDIF

      temp_grid_id = 'grid/r_' // TRIM(sub%name)
      temp_block_id = TRIM(block_id) // '/r_' // TRIM(sub%name)
      temp_name = TRIM(name) // '/Reduced_' // TRIM(sub%name)

      CALL check_name_length('subset', &
          TRIM(name) // '/Reduced_' // TRIM(sub%name))

      CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
          TRIM(temp_name), TRIM(units), sub%n_global, stagger, &
          TRIM(temp_grid_id), reduced, rsubtype, rsubarray, convert)

      dump_skipped = .TRUE.
      sub%dump_field_grid = .TRUE.
      DEALLOCATE(reduced)
    ENDDO

    IF (restart_id .OR. (.NOT.dump_skipped .AND. unaveraged_id)) THEN
      CALL sdf_write_plain_variable(sdf_handle, TRIM(block_id), &
          TRIM(name), TRIM(units), dims, stagger, 'grid', array, &
          subtype, subarray, convert)
      dump_field_grid = .TRUE.
    ENDIF

    ! Dump averages
    DO io = 1, n_io_blocks
      iob => io_block_list(io)
      IF (.NOT.iob%dump) CYCLE

      IF (IAND(mask, c_io_averaged) == 0) CYCLE

      avg => iob%averaged_data(id)
      IF (.NOT.avg%started) CYCLE

      IF (avg%dump_single) THEN
        avg%r4array = avg%r4array / REAL(avg%real_time, r4)

        CALL sdf_write_plain_variable(sdf_handle, &
            TRIM(block_id) // '_averaged', TRIM(name) // '_averaged', &
            TRIM(units), dims, stagger, 'grid', &
            avg%r4array(:,1), subtype_field_r4, subarray_field_r4)

        avg%r4array = 0.0_num
      ELSE
        avg%array = avg%array / avg%real_time

        CALL sdf_write_plain_variable(sdf_handle, &
            TRIM(block_id) // '_averaged', TRIM(name) // '_averaged', &
            TRIM(units), dims, stagger, 'grid', &
            avg%array(:,1), subtype_field, subarray_field)

        avg%array = 0.0_num
      ENDIF

      dump_field_grid = .TRUE.
      avg%real_time = 0.0_num
      avg%started = .FALSE.
    ENDDO

  END SUBROUTINE write_field



  SUBROUTINE write_nspecies_field(id, code, block_id, name, units, stagger, &
      func, array)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: block_id, name, units
    INTEGER, INTENT(IN) :: stagger
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    REAL(num), DIMENSION(:), ALLOCATABLE :: reduced
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER :: ispecies, io, mask
    INTEGER :: i, ii, rnx
    INTEGER :: subtype, subarray, rsubtype, rsubarray
    CHARACTER(LEN=c_id_length) :: temp_block_id, temp_grid_id
    CHARACTER(LEN=c_max_string_length) :: temp_name
    LOGICAL :: convert, dump_sum, dump_species, dump_skipped
    LOGICAL :: normal_id, restart_id, unaveraged_id
    TYPE(averaged_data_block), POINTER :: avg
    TYPE(io_block_type), POINTER :: iob
    TYPE(subset), POINTER :: sub

    INTERFACE
      SUBROUTINE func(data_array, current_species)
        USE constants
        REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
        INTEGER, INTENT(IN) :: current_species
      END SUBROUTINE func
    END INTERFACE

    mask = iomask(id)
    IF (IAND(mask, code) == 0) RETURN
    IF (IAND(mask, c_io_never) /= 0) RETURN

    ! This is a normal dump and normal output variable
    normal_id = IAND(IAND(code, mask), IOR(c_io_always, c_io_full)) /= 0

    IF (.NOT.normal_id) RETURN

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    ! The variable is either averaged or has snapshot specified
    unaveraged_id = IAND(mask, c_io_averaged) == 0 &
        .OR. IAND(mask, c_io_snapshot) /= 0

    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (convert) THEN
      subtype  = subtype_field_r4
      subarray = subarray_field_r4
    ELSE
      subtype  = subtype_field
      subarray = subarray_field
    ENDIF

    dims = (/nx_global/)

    dump_sum = unaveraged_id &
        .AND. IAND(mask, c_io_no_sum) == 0 .AND. IAND(mask, c_io_field) == 0
    dump_species = unaveraged_id .AND. IAND(mask, c_io_species) /= 0

    IF (isubset == 1) THEN
      dump_skipped = .FALSE.
    ELSE
      sub => subset_list(isubset-1)
      dump_skipped = sub%skip
    ENDIF

    IF (dump_sum .OR. dump_species) CALL build_species_subset

    IF (dump_sum) THEN
      IF (isubset == 1) THEN
        temp_block_id = TRIM(block_id)
        temp_name = 'Derived/' // TRIM(name)
      ELSE
        CALL check_name_length('subset', 'Derived/' // TRIM(name) &
            // '/Subset_' // TRIM(sub%name))

        temp_block_id = TRIM(block_id) &
            // '/s_' // TRIM(sub%name)
        temp_name = 'Derived/' // TRIM(name) &
            // '/Subset_' // TRIM(sub%name)
      ENDIF

      CALL func(array, 0)

      IF (dump_skipped) THEN
        rnx = sub%n_local(1)

        ALLOCATE(reduced(rnx))

        ii = sub%n_start(1) + 1
        DO i = 1, rnx
          reduced(i) = array(ii)
          ii = ii + sub%skip_dir(1)
        ENDDO

        IF (convert) THEN
          rsubtype  = sub%subtype_r4
          rsubarray = sub%subarray_r4
        ELSE
          rsubtype  = sub%subtype
          rsubarray = sub%subarray
        ENDIF

        CALL check_name_length('subset', &
            TRIM(temp_name) // '/Reduced_' // TRIM(sub%name))

        temp_grid_id = 'grid/r_' // TRIM(sub%name)
        temp_block_id = TRIM(temp_block_id) // '/r_' // TRIM(sub%name)
        temp_name = TRIM(temp_name) // '/Reduced_' // TRIM(sub%name)

        CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
            TRIM(temp_name), TRIM(units), sub%n_global, stagger, &
            TRIM(temp_grid_id), reduced, rsubtype, rsubarray, convert)

        sub%dump_field_grid = .TRUE.
      ELSE
        CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
            TRIM(temp_name), TRIM(units), dims, stagger, 'grid', array, &
            subtype, subarray, convert)
        dump_field_grid = .TRUE.
      ENDIF
    ENDIF

    IF (dump_species .AND. dump_skipped) THEN
      rnx = sub%n_local(1)

      IF (.NOT.ALLOCATED(reduced)) ALLOCATE(reduced(rnx))

      IF (convert) THEN
        rsubtype  = sub%subtype_r4
        rsubarray = sub%subarray_r4
      ELSE
        rsubtype  = sub%subtype
        rsubarray = sub%subarray
      ENDIF

      temp_grid_id = 'grid/r_' // TRIM(sub%name)

      DO ispecies = 1, n_species
        IF (IAND(io_list(ispecies)%dumpmask, code) == 0) CYCLE

        CALL check_name_length('species', &
            'Derived/' // TRIM(name) // '/' // TRIM(io_list(ispecies)%name))

        temp_block_id = TRIM(block_id) // '/' // TRIM(io_list(ispecies)%name)
        temp_name = &
            'Derived/' // TRIM(name) // '/' // TRIM(io_list(ispecies)%name)

        CALL check_name_length('subset', &
            TRIM(temp_name) // '/Reduced_' // TRIM(sub%name))

        temp_block_id = TRIM(temp_block_id) // '/r_' // TRIM(sub%name)
        temp_name = TRIM(temp_name) // '/Reduced_' // TRIM(sub%name)

        CALL func(array, ispecies)

        ii = sub%n_start(1) + 1
        DO i = 1, rnx
          reduced(i) = array(ii)
          ii = ii + sub%skip_dir(1)
        ENDDO

        CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
            TRIM(temp_name), TRIM(units), sub%n_global, stagger, &
            TRIM(temp_grid_id), reduced, rsubtype, rsubarray, convert)

        sub%dump_field_grid = .TRUE.
      ENDDO
    ELSEIF (dump_species) THEN
      DO ispecies = 1, n_species
        IF (IAND(io_list(ispecies)%dumpmask, code) == 0) CYCLE

        CALL check_name_length('species', &
            'Derived/' // TRIM(name) // '/' // TRIM(io_list(ispecies)%name))

        temp_block_id = TRIM(block_id) // '/' // TRIM(io_list(ispecies)%name)
        temp_name = &
            'Derived/' // TRIM(name) // '/' // TRIM(io_list(ispecies)%name)

        CALL func(array, ispecies)

        CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
            TRIM(temp_name), TRIM(units), dims, stagger, 'grid', array, &
            subtype, subarray, convert)
        dump_field_grid = .TRUE.
      ENDDO
    ENDIF

    IF (ALLOCATED(reduced)) DEALLOCATE(reduced)

    IF (isubset /= 1) RETURN

    ! Write averaged data
    DO io = 1, n_io_blocks
      iob => io_block_list(io)
      IF (.NOT.iob%dump) CYCLE

      mask = iob%dumpmask(id)
      IF (IAND(mask, c_io_averaged) == 0) CYCLE

      avg => iob%averaged_data(id)
      IF (.NOT.avg%started) CYCLE

      IF (avg%dump_single) THEN
        avg%r4array = avg%r4array / REAL(avg%real_time, r4)

        IF (avg%species_sum > 0 .AND. IAND(mask, c_io_field) == 0) THEN
          CALL sdf_write_plain_variable(sdf_handle, &
              TRIM(block_id) // '_averaged', &
              'Derived/' // TRIM(name) // '_averaged', &
              TRIM(units), dims, stagger, 'grid', &
              avg%r4array(:,1), subtype_field_r4, subarray_field_r4)
          dump_field_grid = .TRUE.
        ENDIF

        IF (avg%n_species > 0) THEN
          DO ispecies = 1, avg%n_species
            IF (IAND(io_list(ispecies)%dumpmask, code) == 0) CYCLE

            CALL check_name_length('species', 'Derived/' // TRIM(name) &
                // '_averaged/' // TRIM(io_list(ispecies)%name))

            temp_block_id = TRIM(block_id) &
                // '_averaged/' // TRIM(io_list(ispecies)%name)
            temp_name = 'Derived/' // TRIM(name) &
                // '_averaged/' // TRIM(io_list(ispecies)%name)

            CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
                TRIM(temp_name), TRIM(units), dims, stagger, 'grid', &
                avg%r4array(:,ispecies+avg%species_sum), &
                subtype_field_r4, subarray_field_r4)
            dump_field_grid = .TRUE.
          ENDDO
        ENDIF

        avg%r4array = 0.0_num
      ELSE
        avg%array = avg%array / avg%real_time

        IF (avg%species_sum > 0 .AND. IAND(mask, c_io_field) == 0) THEN
          CALL sdf_write_plain_variable(sdf_handle, &
              TRIM(block_id) // '_averaged', &
              'Derived/' // TRIM(name) // '_averaged', &
              TRIM(units), dims, stagger, 'grid', &
              avg%array(:,1), subtype_field, subarray_field)
          dump_field_grid = .TRUE.
        ENDIF

        IF (avg%n_species > 0) THEN
          DO ispecies = 1, avg%n_species
            IF (IAND(io_list(ispecies)%dumpmask, code) == 0) CYCLE

            CALL check_name_length('species', 'Derived/' // TRIM(name) &
                // '_averaged/' // TRIM(io_list(ispecies)%name))

            temp_block_id = TRIM(block_id) &
                // '_averaged/' // TRIM(io_list(ispecies)%name)
            temp_name = 'Derived/' // TRIM(name) &
                // '_averaged/' // TRIM(io_list(ispecies)%name)

            CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
                TRIM(temp_name), TRIM(units), dims, stagger, 'grid', &
                avg%array(:,ispecies+avg%species_sum), &
                subtype_field, subarray_field)
            dump_field_grid = .TRUE.
          ENDDO
        ENDIF

        avg%array = 0.0_num
      ENDIF

      avg%real_time = 0.0_num
      avg%started = .FALSE.
    ENDDO

  END SUBROUTINE write_nspecies_field



  SUBROUTINE write_field_flux(id, code, block_id, name, units, stagger, &
      func, array, fluxdir, dir_tags)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: block_id, name, units
    INTEGER, INTENT(IN) :: stagger
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, DIMENSION(:), INTENT(IN) :: fluxdir
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: dir_tags
    REAL(num), DIMENSION(:), ALLOCATABLE :: reduced
    INTEGER :: ndirs, idir, mask
    INTEGER :: i, ii, rnx
    INTEGER :: subtype, subarray, rsubtype, rsubarray
    INTEGER, DIMENSION(c_ndims) :: dims
    CHARACTER(LEN=c_id_length) :: temp_block_id, temp_grid_id
    CHARACTER(LEN=c_max_string_length) :: temp_name
    LOGICAL :: convert, normal_id, dump_skipped
    TYPE(subset), POINTER :: sub

    INTERFACE
      SUBROUTINE func(data_array, direction)
        USE constants
        REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
        INTEGER, INTENT(IN) :: direction
      END SUBROUTINE func
    END INTERFACE

    mask = iomask(id)
    IF (IAND(mask, code) == 0) RETURN
    IF (IAND(mask, c_io_never) /= 0) RETURN

    ! This is a normal dump and normal output variable
    normal_id = IAND(IAND(code, mask), IOR(c_io_always, c_io_full)) /= 0

    IF (.NOT.normal_id) RETURN

    convert = IAND(mask, c_io_dump_single) /= 0

    IF (convert) THEN
      subtype  = subtype_field_r4
      subarray = subarray_field_r4
    ELSE
      subtype  = subtype_field
      subarray = subarray_field
    ENDIF

    ndirs = SIZE(fluxdir)
    dims = (/nx_global/)

    IF (isubset == 1) THEN
      dump_skipped = .FALSE.
    ELSE
      sub => subset_list(isubset-1)
      dump_skipped = sub%skip
    ENDIF

    IF (dump_skipped) THEN
      rnx = sub%n_local(1)

      ALLOCATE(reduced(rnx))

      IF (convert) THEN
        rsubtype  = sub%subtype_r4
        rsubarray = sub%subarray_r4
      ELSE
        rsubtype  = sub%subtype
        rsubarray = sub%subarray
      ENDIF

      temp_grid_id = 'grid/r_' // TRIM(sub%name)

      DO idir = 1, ndirs
        temp_block_id = TRIM(block_id) // '/' // TRIM(dir_tags(idir))
        temp_name = &
            'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir))

        CALL check_name_length('subset', &
            TRIM(temp_name) // '/Reduced_' // TRIM(sub%name))

        temp_block_id = TRIM(temp_block_id) // '/r_' // TRIM(sub%name)
        temp_name = TRIM(temp_name) // '/Reduced_' // TRIM(sub%name)

        CALL func(array, fluxdir(idir))

        ii = sub%n_start(1) + 1
        DO i = 1, rnx
          reduced(i) = array(ii)
          ii = ii + sub%skip_dir(1)
        ENDDO

        CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
            TRIM(temp_name), TRIM(units), sub%n_global, stagger, &
            TRIM(temp_grid_id), reduced, rsubtype, rsubarray, convert)
      ENDDO

      sub%dump_field_grid = .TRUE.
      DEALLOCATE(reduced)
    ELSE
      idir = 1
      CALL check_name_length('dir tag', &
          'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)))

      DO idir = 1, ndirs
        temp_block_id = TRIM(block_id) // '/' // TRIM(dir_tags(idir))
        temp_name = &
            'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir))

        CALL func(array, fluxdir(idir))

        CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
            TRIM(temp_name), TRIM(units), dims, stagger, 'grid', array, &
            subtype, subarray, convert)
      ENDDO
      dump_field_grid = .TRUE.
    ENDIF

    ! Flux variables not currently averaged

  END SUBROUTINE write_field_flux



  SUBROUTINE write_nspecies_flux(id, code, block_id, name, units, stagger, &
      func, array, fluxdir, dir_tags)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: block_id, name, units
    INTEGER, INTENT(IN) :: stagger
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, DIMENSION(:), INTENT(IN) :: fluxdir
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: dir_tags
    REAL(num), DIMENSION(:), ALLOCATABLE :: reduced
    INTEGER :: ispecies, ndirs, idir, mask
    INTEGER :: i, ii, rnx
    INTEGER :: subtype, subarray, rsubtype, rsubarray
    INTEGER, DIMENSION(c_ndims) :: dims
    CHARACTER(LEN=c_id_length) :: temp_block_id, temp_grid_id
    CHARACTER(LEN=c_max_string_length) :: temp_name
    LOGICAL :: convert, dump_sum, dump_species, dump_skipped
    LOGICAL :: normal_id
    TYPE(subset), POINTER :: sub

    INTERFACE
      SUBROUTINE func(data_array, current_species, direction)
        USE constants
        REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
        INTEGER, INTENT(IN) :: current_species, direction
      END SUBROUTINE func
    END INTERFACE

    mask = iomask(id)
    IF (IAND(mask, code) == 0) RETURN
    IF (IAND(mask, c_io_never) /= 0) RETURN

    ! This is a normal dump and normal output variable
    normal_id = IAND(IAND(code, mask), IOR(c_io_always, c_io_full)) /= 0

    IF (.NOT.normal_id) RETURN

    convert = IAND(mask, c_io_dump_single) /= 0

    IF (convert) THEN
      subtype  = subtype_field_r4
      subarray = subarray_field_r4
    ELSE
      subtype  = subtype_field
      subarray = subarray_field
    ENDIF

    ndirs = SIZE(fluxdir)
    dims = (/nx_global/)

    dump_sum = IAND(mask, c_io_no_sum) == 0 .AND. IAND(mask, c_io_field) == 0
    dump_species = IAND(mask, c_io_species) /= 0

    IF (.NOT.dump_sum .AND. .NOT.dump_species) RETURN

    IF (isubset == 1) THEN
      dump_skipped = .FALSE.
    ELSE
      sub => subset_list(isubset-1)
      dump_skipped = sub%skip
    ENDIF

    CALL build_species_subset

    IF (dump_sum .AND. dump_skipped) THEN
      rnx = sub%n_local(1)

      ALLOCATE(reduced(rnx))

      IF (convert) THEN
        rsubtype  = sub%subtype_r4
        rsubarray = sub%subarray_r4
      ELSE
        rsubtype  = sub%subtype
        rsubarray = sub%subarray
      ENDIF

      idir = 1
      CALL check_name_length('dir tag', &
          'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)))

      temp_grid_id = 'grid/r_' // TRIM(sub%name)

      DO idir = 1, ndirs
        temp_block_id = TRIM(block_id) // '/' // TRIM(dir_tags(idir))
        temp_name = &
            'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir))

        CALL check_name_length('subset', &
            TRIM(temp_name) // '/Reduced_' // TRIM(sub%name))

        temp_block_id = TRIM(temp_block_id) // '/r_' // TRIM(sub%name)
        temp_name = TRIM(temp_name) // '/Reduced_' // TRIM(sub%name)

        CALL func(array, 0, fluxdir(idir))

        ii = sub%n_start(1) + 1
        DO i = 1, rnx
          reduced(i) = array(ii)
          ii = ii + sub%skip_dir(1)
        ENDDO

        CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
            TRIM(temp_name), TRIM(units), sub%n_global, stagger, &
            TRIM(temp_grid_id), reduced, rsubtype, rsubarray, convert)
      ENDDO

      sub%dump_field_grid = .TRUE.
    ELSEIF (dump_sum) THEN
      idir = 1
      CALL check_name_length('dir tag', &
          'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir)))

      DO idir = 1, ndirs
        temp_block_id = TRIM(block_id) // '/' // TRIM(dir_tags(idir))
        temp_name = &
            'Derived/' // TRIM(name) // '/' // TRIM(dir_tags(idir))

        CALL func(array, 0, fluxdir(idir))

        CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
            TRIM(temp_name), TRIM(units), dims, stagger, 'grid', array, &
            subtype, subarray, convert)
      ENDDO
      dump_field_grid = .TRUE.
    ENDIF

    IF (dump_species .AND. dump_skipped) THEN
      rnx = sub%n_local(1)

      IF (.NOT.ALLOCATED(reduced)) ALLOCATE(reduced(rnx))

      IF (convert) THEN
        rsubtype  = sub%subtype_r4
        rsubarray = sub%subarray_r4
      ELSE
        rsubtype  = sub%subtype
        rsubarray = sub%subarray
      ENDIF

      temp_grid_id = 'grid/r_' // TRIM(sub%name)

      DO ispecies = 1, n_species
        IF (IAND(io_list(ispecies)%dumpmask, code) == 0) CYCLE

        idir = 1
        CALL check_name_length('species', 'Derived/' // TRIM(name) &
            // '_' // TRIM(dir_tags(idir)) // '/' &
            // TRIM(io_list(ispecies)%name))

        DO idir = 1, ndirs
          temp_block_id = TRIM(block_id) &
              // '_' // TRIM(dir_tags(idir)) // '/' &
              // TRIM(io_list(ispecies)%name)
          temp_name = 'Derived/' // TRIM(name) &
              // '_' // TRIM(dir_tags(idir)) // '/' &
              // TRIM(io_list(ispecies)%name)

          CALL check_name_length('subset', &
              TRIM(temp_name) // '/Reduced_' // TRIM(sub%name))

          temp_block_id = TRIM(temp_block_id) // '/r_' // TRIM(sub%name)
          temp_name = TRIM(temp_name) // '/Reduced_' // TRIM(sub%name)

          CALL func(array, ispecies, fluxdir(idir))

          ii = sub%n_start(1) + 1
          DO i = 1, rnx
            reduced(i) = array(ii)
            ii = ii + sub%skip_dir(1)
          ENDDO

          CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
              TRIM(temp_name), TRIM(units), sub%n_global, stagger, &
              TRIM(temp_grid_id), reduced, rsubtype, rsubarray, convert)
        ENDDO
        sub%dump_field_grid = .TRUE.
      ENDDO
    ELSEIF (dump_species) THEN
      DO ispecies = 1, n_species
        IF (IAND(io_list(ispecies)%dumpmask, code) == 0) CYCLE

        idir = 1
        CALL check_name_length('species', 'Derived/' // TRIM(name) &
            // '_' // TRIM(dir_tags(idir)) // '/' &
            // TRIM(io_list(ispecies)%name))

        DO idir = 1, ndirs
          temp_block_id = TRIM(block_id) &
              // '_' // TRIM(dir_tags(idir)) // '/' &
              // TRIM(io_list(ispecies)%name)
          temp_name = 'Derived/' // TRIM(name) &
              // '_' // TRIM(dir_tags(idir)) // '/' &
              // TRIM(io_list(ispecies)%name)

          CALL func(array, ispecies, fluxdir(idir))

          CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
              TRIM(temp_name), TRIM(units), dims, stagger, 'grid', array, &
              subtype, subarray, convert)
        ENDDO
        dump_field_grid = .TRUE.
      ENDDO
    ENDIF

    IF (ALLOCATED(reduced)) DEALLOCATE(reduced)

    ! Flux variables not currently averaged

  END SUBROUTINE write_nspecies_flux



  SUBROUTINE build_species_subset

    INTEGER :: i, l
    TYPE(particle), POINTER :: current, next
    LOGICAL :: use_particle
    REAL(num) :: gamma, random_num, part_mc

    IF (done_subset_init) RETURN
    done_subset_init = .TRUE.

    IF (isubset == 1) THEN
      io_list => species_list
      RETURN
    ENDIF

    io_list => io_list_data

    l = isubset - 1
    DO i = 1, n_species
      io_list(i) = species_list(i)
      io_list(i)%count = 0
      io_list(i)%name = 'subset_' // TRIM(subset_list(l)%name) // '/' // &
          TRIM(species_list(i)%name)
      CALL create_empty_partlist(io_list(i)%attached_list)

      IF (.NOT. subset_list(l)%use_species(i)) THEN
        io_list(i)%dumpmask = c_io_never
        CYCLE
      ENDIF

      part_mc = c * species_list(i)%mass

      current => species_list(i)%attached_list%head
      DO WHILE (ASSOCIATED(current))
        next => current%next
        use_particle = .TRUE.
        IF (subset_list(l)%use_gamma) THEN
#ifdef PER_PARTICLE_CHARGE_MASS
          part_mc = c * current%mass
#endif
          gamma = SQRT(SUM((current%part_p / part_mc)**2) + 1.0_num)
          IF (subset_list(l)%use_gamma_min &
              .AND. gamma < subset_list(l)%gamma_min) use_particle = .FALSE.
          IF (subset_list(l)%use_gamma_max &
              .AND. gamma > subset_list(l)%gamma_max) use_particle = .FALSE.
        ENDIF

        IF (subset_list(l)%use_x_min &
            .AND. current%part_pos < subset_list(l)%x_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_x_max &
            .AND. current%part_pos > subset_list(l)%x_max) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_px_min &
            .AND. current%part_p(1) < subset_list(l)%px_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_px_max &
            .AND. current%part_p(1) > subset_list(l)%px_max) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_py_min &
            .AND. current%part_p(2) < subset_list(l)%py_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_py_max &
            .AND. current%part_p(2) > subset_list(l)%py_max) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_pz_min &
            .AND. current%part_p(3) < subset_list(l)%pz_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_pz_max &
            .AND. current%part_p(3) > subset_list(l)%pz_max) &
                use_particle = .FALSE.

#ifndef PER_SPECIES_WEIGHT
        IF (subset_list(l)%use_weight_min &
            .AND. current%weight < subset_list(l)%weight_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_weight_max &
            .AND. current%weight > subset_list(l)%weight_max) &
                use_particle = .FALSE.

#endif
#ifdef PER_PARTICLE_CHARGE_MASS
        IF (subset_list(l)%use_charge_min &
            .AND. current%charge < subset_list(l)%charge_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_charge_max &
            .AND. current%charge > subset_list(l)%charge_max) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_mass_min &
            .AND. current%mass < subset_list(l)%mass_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_mass_max &
            .AND. current%mass > subset_list(l)%mass_max) &
                use_particle = .FALSE.

#endif
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
        IF (subset_list(l)%use_id_min &
            .AND. current%id < subset_list(l)%id_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_id_max &
            .AND. current%id > subset_list(l)%id_max) &
                use_particle = .FALSE.
#endif

        IF (subset_list(l)%use_random) THEN
          random_num = random()
          IF (random_num > subset_list(l)%random_fraction) &
              use_particle = .FALSE.
        ENDIF

        IF (use_particle) THEN
          ! Move particle to io_list
          CALL remove_particle_from_partlist(species_list(i)%attached_list, &
              current)
          CALL add_particle_to_partlist(io_list(i)%attached_list, current)
        ENDIF
        current => next
      ENDDO
    ENDDO

  END SUBROUTINE build_species_subset



  SUBROUTINE species_offset_init

    INTEGER(i8) :: species_count
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: npart_species_per_proc
    INTEGER :: i, ispecies
    TYPE(particle_species), POINTER :: spec

    IF (done_species_offset_init) RETURN
    done_species_offset_init = .TRUE.

    IF (.NOT.ALLOCATED(species_offset))  ALLOCATE(species_offset(n_species))

    CALL build_species_subset

    ALLOCATE(npart_species_per_proc(nproc))
    species_offset = 0

    npart_global = 0
    DO ispecies = 1, n_species
      spec => io_list(ispecies)

      CALL MPI_ALLGATHER(spec%attached_list%count, 1, MPI_INTEGER8, &
          npart_species_per_proc, 1, MPI_INTEGER8, comm, errcode)
      species_count = 0
      DO i = 1, nproc
        IF (rank == i-1) species_offset(ispecies) = species_count
        species_count = species_count + npart_species_per_proc(i)
      ENDDO
      spec%count = species_count
      spec%count_update_step = step
      npart_global = npart_global + species_count

      CALL sdf_write_cpu_split(sdf_handle, 'cpu/' // TRIM(spec%name), &
          'CPU split/' // TRIM(spec%name), npart_species_per_proc)
    ENDDO

    IF (track_ejected_particles &
        .AND. .NOT.ALLOCATED(ejected_offset)) THEN
      ALLOCATE(ejected_offset(n_species))
      ejected_offset = 0

      DO ispecies = 1, n_species
        spec => ejected_list(ispecies)

        CALL MPI_ALLGATHER(spec%attached_list%count, 1, MPI_INTEGER8, &
            npart_species_per_proc, 1, MPI_INTEGER8, comm, errcode)
        species_count = 0
        DO i = 1, nproc
          IF (rank == i-1) ejected_offset(ispecies) = species_count
          species_count = species_count + npart_species_per_proc(i)
        ENDDO
        spec%count = species_count
        spec%count_update_step = step

        CALL sdf_write_cpu_split(sdf_handle, 'cpu/' // TRIM(spec%name), &
            'CPU split/' // TRIM(spec%name), npart_species_per_proc)
      ENDDO
    ENDIF

    DEALLOCATE(npart_species_per_proc)

  END SUBROUTINE species_offset_init



  SUBROUTINE write_particle_grid(code)

    INTEGER, INTENT(IN) :: code
    INTEGER :: ispecies, id, mask
    LOGICAL :: convert, dump_grid, restart_id

    id = c_dump_part_grid
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (IAND(mask, c_io_never) == 0 &
        .AND. (IAND(mask, code) /= 0 .OR. ANY(dump_point_grid))) THEN
      CALL build_species_subset

      DO ispecies = 1, n_species
        current_species => io_list(ispecies)

        mask = current_species%dumpmask
        dump_grid = dump_point_grid(ispecies)

        IF (IAND(mask, code) /= 0) dump_grid = .TRUE.
        IF (IAND(mask, c_io_never) /= 0) dump_grid = .FALSE.
        IF (IAND(code, c_io_restartable) /= 0) dump_grid = .TRUE.

        IF (dump_grid) THEN
          CALL species_offset_init()
          IF (npart_global == 0) RETURN

          CALL sdf_write_point_mesh(sdf_handle, &
              'grid/' // TRIM(current_species%name), &
              'Grid/Particles/' // TRIM(current_species%name), &
              TRIM(current_species%name), &
              io_list(ispecies)%count, c_dimension_1d, &
              it_output_position, species_offset(ispecies), convert)
        ENDIF
      ENDDO
    ENDIF

    IF (isubset /= 1) RETURN

    id = c_dump_ejected_particles
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (IAND(mask, c_io_never) == 0 .AND. IAND(mask, code) /= 0) THEN
      reset_ejected = .TRUE.

      DO ispecies = 1, n_species
        CALL species_offset_init()
        IF (npart_global == 0) RETURN

        current_species => ejected_list(ispecies)
        CALL sdf_write_point_mesh(sdf_handle, &
            'grid/' // TRIM(current_species%name), &
            'Grid/Particles/' // TRIM(current_species%name), &
            TRIM(current_species%name), &
            ejected_list(ispecies)%count, c_dimension_1d, &
            it_output_position, ejected_offset(ispecies), convert)
      ENDDO
    ENDIF

  END SUBROUTINE write_particle_grid



  SUBROUTINE write_particle_variable_num(id_in, code, name, units, iterator)

    INTEGER, INTENT(IN) :: id_in, code
    CHARACTER(LEN=*), INTENT(IN) :: name, units
    CHARACTER(LEN=c_id_length) :: temp_block_id
    INTEGER :: ispecies, id, mask
    LOGICAL :: convert, found, restart_id

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, param)
        USE constants
        REAL(num) :: iterator
        REAL(num), DIMENSION(:), INTENT(OUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    id = id_in
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (IAND(mask, c_io_never) == 0 .AND. IAND(mask, code) /= 0) THEN
      CALL build_species_subset

      DO ispecies = 1, n_species
        current_species => io_list(ispecies)

        IF (IAND(current_species%dumpmask, code) /= 0 &
            .OR. IAND(code, c_io_restartable) /= 0) THEN
          CALL species_offset_init()
          IF (npart_global == 0) RETURN

          found = sdf_get_block_id(sdf_handle, &
              'grid/' // TRIM(current_species%name), temp_block_id)
          CALL sdf_write_point_variable(sdf_handle, &
              lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
              'Particles/' // TRIM(name) // '/' // TRIM(current_species%name), &
              TRIM(current_species%name), &
              TRIM(units), io_list(ispecies)%count, temp_block_id, &
              iterator, id_in, species_offset(ispecies), convert)
          dump_point_grid(ispecies) = .TRUE.
        ENDIF
      ENDDO
    ENDIF

    id = c_dump_ejected_particles
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (IAND(mask, c_io_never) == 0 .AND. IAND(mask, code) /= 0) THEN
      reset_ejected = .TRUE.

      DO ispecies = 1, n_species
        CALL species_offset_init()
        IF (npart_global == 0) RETURN

        current_species => ejected_list(ispecies)
        found = sdf_get_block_id(sdf_handle, &
            'grid/' // TRIM(current_species%name), temp_block_id)
        CALL sdf_write_point_variable(sdf_handle, &
            lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
            'Particles/' // TRIM(name) // '/' // TRIM(current_species%name), &
            TRIM(current_species%name), &
            TRIM(units), ejected_list(ispecies)%count, temp_block_id, &
            iterator, id_in, ejected_offset(ispecies), convert)
      ENDDO
    ENDIF

  END SUBROUTINE write_particle_variable_num



#if defined(PARTICLE_ID4) || defined(PARTICLE_DEBUG)
  SUBROUTINE write_particle_variable_i4(id_in, code, name, units, iterator)

    INTEGER, INTENT(IN) :: id_in, code
    CHARACTER(LEN=*), INTENT(IN) :: name, units
    CHARACTER(LEN=c_id_length) :: temp_block_id
    INTEGER :: ispecies, id, mask
    LOGICAL :: convert, found, restart_id

    INTERFACE
      FUNCTION iterator(array, npart_it, start, param)
        USE constants
        INTEGER(i4) :: iterator
        INTEGER(i4), DIMENSION(:), INTENT(OUT) :: array
        INTEGER, INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    id = id_in
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (IAND(mask, c_io_never) == 0 .AND. IAND(mask, code) /= 0) THEN
      CALL build_species_subset

      DO ispecies = 1, n_species
        current_species => io_list(ispecies)

        IF (IAND(current_species%dumpmask, code) /= 0 &
            .OR. IAND(code, c_io_restartable) /= 0) THEN
          CALL species_offset_init()
          IF (npart_global == 0) RETURN

          found = sdf_get_block_id(sdf_handle, &
              'grid/' // TRIM(current_species%name), temp_block_id)
          CALL sdf_write_point_variable(sdf_handle, &
              lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
              'Particles/' // TRIM(name) // '/' // TRIM(current_species%name), &
              TRIM(current_species%name), &
              TRIM(units), io_list(ispecies)%count, temp_block_id, &
              iterator, id_in, species_offset(ispecies), convert)
          dump_point_grid(ispecies) = .TRUE.
        ENDIF
      ENDDO
    ENDIF

    id = c_dump_ejected_particles
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (IAND(mask, c_io_never) == 0 .AND. IAND(mask, code) /= 0) THEN
      reset_ejected = .TRUE.

      DO ispecies = 1, n_species
        CALL species_offset_init()
        IF (npart_global == 0) RETURN

        current_species => ejected_list(ispecies)
        found = sdf_get_block_id(sdf_handle, &
            'grid/' // TRIM(current_species%name), temp_block_id)
        CALL sdf_write_point_variable(sdf_handle, &
            lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
            'Particles/' // TRIM(name) // '/' // TRIM(current_species%name), &
            TRIM(current_species%name), &
            TRIM(units), ejected_list(ispecies)%count, temp_block_id, &
            iterator, id_in, ejected_offset(ispecies), convert)
      ENDDO
    ENDIF

  END SUBROUTINE write_particle_variable_i4
#endif



#if defined(PARTICLE_ID)
  SUBROUTINE write_particle_variable_i8(id_in, code, name, units, iterator)

    INTEGER, INTENT(IN) :: id_in, code
    CHARACTER(LEN=*), INTENT(IN) :: name, units
    CHARACTER(LEN=c_id_length) :: temp_block_id
    INTEGER :: ispecies, id, mask
    LOGICAL :: convert, found, restart_id

    INTERFACE
      FUNCTION iterator(array, npart_it, start, param)
        USE constants
        INTEGER(i8) :: iterator
        INTEGER(i8), DIMENSION(:), INTENT(OUT) :: array
        INTEGER, INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    id = id_in
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (IAND(mask, c_io_never) == 0 .AND. IAND(mask, code) /= 0) THEN
      CALL build_species_subset

      DO ispecies = 1, n_species
        current_species => io_list(ispecies)

        IF (IAND(current_species%dumpmask, code) /= 0 &
            .OR. IAND(code, c_io_restartable) /= 0) THEN
          CALL species_offset_init()
          IF (npart_global == 0) RETURN

          found = sdf_get_block_id(sdf_handle, &
              'grid/' // TRIM(current_species%name), temp_block_id)
          CALL sdf_write_point_variable(sdf_handle, &
              lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
              'Particles/' // TRIM(name) // '/' // TRIM(current_species%name), &
              TRIM(current_species%name), &
              TRIM(units), io_list(ispecies)%count, temp_block_id, &
              iterator, id_in, species_offset(ispecies), convert)
          dump_point_grid(ispecies) = .TRUE.
        ENDIF
      ENDDO
    ENDIF

    id = c_dump_ejected_particles
    mask = iomask(id)

    ! This is a restart dump and a restart variable
    restart_id = IAND(IAND(code, mask), c_io_restartable) /= 0
    convert = IAND(mask, c_io_dump_single) /= 0 .AND. .NOT.restart_id

    IF (IAND(mask, c_io_never) == 0 .AND. IAND(mask, code) /= 0) THEN
      reset_ejected = .TRUE.

      DO ispecies = 1, n_species
        CALL species_offset_init()
        IF (npart_global == 0) RETURN

        current_species => ejected_list(ispecies)
        found = sdf_get_block_id(sdf_handle, &
            'grid/' // TRIM(current_species%name), temp_block_id)
        CALL sdf_write_point_variable(sdf_handle, &
            lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
            'Particles/' // TRIM(name) // '/' // TRIM(current_species%name), &
            TRIM(current_species%name), &
            TRIM(units), ejected_list(ispecies)%count, temp_block_id, &
            iterator, id_in, ejected_offset(ispecies), convert)
      ENDDO
    ENDIF

  END SUBROUTINE write_particle_variable_i8
#endif



  FUNCTION lowercase(string_in) RESULT(string_out)

    CHARACTER(LEN=*), PARAMETER :: lwr = 'abcdefghijklmnopqrstuvwxyz'
    CHARACTER(LEN=*), PARAMETER :: upr = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    CHARACTER(LEN=*), INTENT(IN) :: string_in
    CHARACTER(LEN=LEN(string_in)) :: string_out
    INTEGER :: i, idx

    string_out = string_in

    DO i = 1, LEN(string_out)
      idx = INDEX(upr, string_out(i:i))
      IF (idx /= 0) string_out(i:i) = lwr(idx:idx)
    ENDDO

  END FUNCTION lowercase



  SUBROUTINE create_timestring(time, timestring)

    REAL(num), INTENT(IN) :: time
    CHARACTER(LEN=*), INTENT(INOUT) :: timestring ! length at least 15
    INTEGER :: days, hours, minutes, seconds, frac_seconds

    days = INT(time) / 60 / 60  / 24
    hours = INT(time) / 60 / 60 - days * 24
    minutes = INT(time) / 60 - (days * 24 + hours) * 60
    seconds = INT(time) - ((days * 24 + hours) * 60 + minutes) * 60
    frac_seconds = FLOOR((time - INT(time)) * 100)

    WRITE(timestring, '(i2,'':'',i2.2,'':'',i2.2,'':'',i2.2,''.'',i2.2)') &
        days, hours, minutes, seconds, frac_seconds

  END SUBROUTINE create_timestring



  SUBROUTINE create_full_timestring(time, timestring)

    REAL(num), INTENT(IN) :: time
    CHARACTER(LEN=*), INTENT(INOUT) :: timestring ! length at least 48
    INTEGER :: days, hours, minutes, seconds, frac_seconds, var
    CHARACTER(LEN=8) :: varstring
    CHARACTER(LEN=4) :: intstring, fracstring
    LOGICAL :: string_started

    days = INT(time) / 60 / 60  / 24
    hours = INT(time) / 60 / 60 - days * 24
    minutes = INT(time) / 60 - (days * 24 + hours) * 60
    seconds = INT(time) - ((days * 24 + hours) * 60 + minutes) * 60
    frac_seconds = FLOOR((time - INT(time)) * 100)

    timestring = ''
    string_started = .FALSE.

    var = days
    varstring = ' day'
    IF (var > 0) THEN
      CALL integer_as_string(var, intstring)
      IF (string_started) THEN
        timestring = TRIM(timestring) // ', ' // TRIM(intstring) &
            // TRIM(varstring)
      ELSE
        timestring = TRIM(timestring) // TRIM(intstring) // TRIM(varstring)
      ENDIF
      IF (var > 1) timestring = TRIM(timestring) // 's'
      string_started = .TRUE.
    ENDIF

    var = hours
    varstring = ' hour'
    IF (var > 0) THEN
      CALL integer_as_string(var, intstring)
      IF (string_started) THEN
        timestring = TRIM(timestring) // ', ' // TRIM(intstring) &
            // TRIM(varstring)
      ELSE
        timestring = TRIM(timestring) // TRIM(intstring) // TRIM(varstring)
      ENDIF
      IF (var > 1) timestring = TRIM(timestring) // 's'
      string_started = .TRUE.
    ENDIF

    var = minutes
    varstring = ' minute'
    IF (var > 0) THEN
      CALL integer_as_string(var, intstring)
      IF (string_started) THEN
        timestring = TRIM(timestring) // ', ' // TRIM(intstring) &
            // TRIM(varstring)
      ELSE
        timestring = TRIM(timestring) // TRIM(intstring) // TRIM(varstring)
      ENDIF
      IF (var > 1) timestring = TRIM(timestring) // 's'
      string_started = .TRUE.
    ENDIF

    var = seconds
    varstring = ' seconds'
    IF (var > 0 .OR. frac_seconds > 0 .OR. .NOT.string_started) THEN
      CALL integer_as_string(var, intstring)
      WRITE(fracstring, '(i2.2)') frac_seconds
      IF (string_started) THEN
        timestring = TRIM(timestring) // ', ' // TRIM(intstring) // '.' &
            // TRIM(fracstring) // TRIM(varstring)
      ELSE
        timestring = TRIM(timestring) // TRIM(intstring) // '.' &
            // TRIM(fracstring) // TRIM(varstring)
      ENDIF
    ENDIF

  END SUBROUTINE create_full_timestring



  SUBROUTINE cleanup_stop_files()

    INTEGER :: ierr

    IF (rank /= 0) RETURN

    OPEN(unit=lu, status='OLD', &
        file=TRIM(data_dir) // '/' // TRIM(stop_file), iostat=ierr)
    IF (ierr == 0) CLOSE(lu, status='DELETE')

    OPEN(unit=lu, status='OLD', &
        file=TRIM(data_dir) // '/' // TRIM(stop_file_nodump), iostat=ierr)
    IF (ierr == 0) CLOSE(lu, status='DELETE')

  END SUBROUTINE cleanup_stop_files



  SUBROUTINE check_for_stop_condition(halt, force_dump)

    LOGICAL, INTENT(OUT) :: halt, force_dump
    INTEGER :: ierr
    INTEGER, SAVE :: check_counter = 0
    LOGICAL :: buffer(2), got_stop_condition, got_stop_file
    REAL(num) :: walltime

    IF (check_stop_frequency <= 0 .AND. .NOT.check_walltime) RETURN

    walltime = -1.0_num
    IF (check_walltime) &
        CALL check_walltime_auto(walltime, halt)

    IF (halt) THEN
      force_dump = .TRUE.
      RETURN
    ENDIF

    IF (check_stop_frequency < 0) RETURN

    got_stop_condition = .FALSE.
    force_dump = .FALSE.
    check_counter = check_counter + 1

    IF (check_counter < check_stop_frequency) RETURN
    check_counter = 0

    IF (rank == 0) THEN
      ! Since we're checking for a STOP file, we might as well check the
      ! walltime as well
      IF (check_walltime) THEN
        IF (walltime < 0.0_num) walltime = MPI_WTIME()
        IF (walltime - real_walltime_start >= stop_at_walltime) THEN
          got_stop_condition = .TRUE.
          force_dump = .TRUE.
          PRINT*,'Stopping because "stop_at_walltime" has been exceeded.'
        ENDIF
      ENDIF

      ! Next check if stop file exists
      OPEN(unit=lu, status='OLD', iostat=ierr, &
          file=TRIM(data_dir) // '/' // TRIM(stop_file))
      IF (ierr == 0) THEN
        got_stop_file = .TRUE.
        got_stop_condition = .TRUE.
        force_dump = .TRUE.
        CLOSE(lu, status='DELETE')
      ELSE
        OPEN(unit=lu, status='OLD', iostat=ierr, &
            file=TRIM(data_dir) // '/' // TRIM(stop_file_nodump))
        IF (ierr == 0) THEN
          got_stop_file = .TRUE.
          got_stop_condition = .TRUE.
          force_dump = .FALSE.
          CLOSE(lu, status='DELETE')
        ELSE
          got_stop_file = .FALSE.
        ENDIF
      ENDIF

      IF (got_stop_file) PRINT*,'Stopping because "STOP" file has been found.'

      buffer(1) = got_stop_condition
      buffer(2) = force_dump
    ENDIF

    CALL MPI_BCAST(buffer, 2, MPI_LOGICAL, 0, comm, errcode)
    got_stop_condition = buffer(1)
    force_dump = buffer(2)

    IF (got_stop_condition) halt = .TRUE.

  END SUBROUTINE check_for_stop_condition



  SUBROUTINE check_walltime_auto(walltime, halt)

    REAL(num), INTENT(INOUT) :: walltime
    LOGICAL, INTENT(OUT) :: halt
    INTEGER, PARAMETER :: tag = 2001
    INTEGER :: msg, request, i
    INTEGER :: status_ignore(MPI_STATUS_SIZE)
    LOGICAL :: flag
    LOGICAL, ALLOCATABLE, SAVE :: completed(:)
    LOGICAL, SAVE :: yet_to_sync = .TRUE.
    LOGICAL, SAVE :: first = .TRUE.
    LOGICAL, SAVE :: all_completed = .FALSE.
    REAL(num), SAVE :: wall0
    REAL(num), PARAMETER :: frac = 1.0_num
    REAL(num) :: timeout

    halt = all_completed
    IF (all_completed) RETURN

    IF (walltime < 0) walltime = MPI_WTIME()
    IF ((walltime + timer_average(c_timer_step) + timer_average(c_timer_io) &
        + timer_average(c_timer_balance) - real_walltime_start) &
        < frac * stop_at_walltime) RETURN

    IF (rank == 0) THEN
      IF (first) THEN
        ALLOCATE(completed(nproc-1))
        completed = .FALSE.
        first = .FALSE.
      ENDIF
      wall0 = walltime
      timeout = 2.0_num * timer_average(c_timer_step)
      DO
        all_completed = .TRUE.
        DO i = 1,nproc-1
          IF (.NOT.completed(i)) THEN
            ! The Platform-MPI interface for MPI_IPROBE is broken and will not
            ! allow us to use MPI_STATUS_IGNORE
            CALL MPI_IPROBE(i, tag, comm, flag, status_ignore, errcode)
            completed(i) = flag
            IF (flag) THEN
              CALL MPI_RECV(msg, 0, MPI_INTEGER, i, tag, comm, &
                  MPI_STATUS_IGNORE, errcode)
            ELSE
              all_completed = .FALSE.
            ENDIF
          ENDIF
        ENDDO
        IF (all_completed) THEN
          DEALLOCATE(completed)
          msg = -1
          DO i = 1,nproc-1
            CALL MPI_ISEND(msg, 1, MPI_INTEGER, i, tag, comm, request, errcode)
            CALL MPI_REQUEST_FREE(request, errcode)
          ENDDO
          halt = all_completed
          PRINT*,'Stopping because "stop_at_walltime" has been exceeded.'
          RETURN
        ENDIF
        walltime = MPI_WTIME()
        IF (walltime - wall0 > timeout) THEN
          msg = 1
          DO i = 1,nproc-1
            IF (completed(i)) THEN
              CALL MPI_ISEND(msg, 1, MPI_INTEGER, i, tag, comm, request, &
                  errcode)
              CALL MPI_REQUEST_FREE(request, errcode)
            ENDIF
          ENDDO
          RETURN
        ENDIF
      ENDDO
      RETURN
    ENDIF

    IF (yet_to_sync) THEN
      CALL MPI_ISEND(0, 0, MPI_INTEGER, 0, tag, comm, request, errcode)
      CALL MPI_REQUEST_FREE(request, errcode)
      yet_to_sync = .FALSE.
    ENDIF
    CALL MPI_RECV(msg, 1, MPI_INTEGER, 0, tag, comm, &
        MPI_STATUS_IGNORE, errcode)
    IF (msg < 0) all_completed = .TRUE.
    halt = all_completed

  END SUBROUTINE check_walltime_auto



  SUBROUTINE epoch_write_source_info(h)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=c_id_length) :: stitched_ids(3)
    CHARACTER(LEN=c_id_length) :: time_string
    CHARACTER(LEN=512) :: string_array(6)
    INTEGER :: n

    n = 0

    IF (SIZE(epoch_bytes) > 1 .OR. (TRIM(epoch_bytes_checksum_type) /= '' &
        .AND. ICHAR(epoch_bytes_checksum_type(1:1)) /= 0)) THEN
      n = n + 1
      CALL sdf_safe_copy_id(h, 'epoch_source/source', stitched_ids(n))
      CALL sdf_write_datablock(h, stitched_ids(n), &
          'EPOCH source code', epoch_bytes, &
          epoch_bytes_padding, epoch_bytes_mimetype, &
          epoch_bytes_checksum_type, epoch_bytes_checksum)
    ENDIF

    IF (SIZE(epoch_bytes) == 1 .AND. SIZE(epoch_diff_bytes) > 1) THEN
      n = n + 1
      CALL sdf_safe_copy_id(h, 'epoch_source/diff', stitched_ids(n))
      CALL sdf_write_datablock(h, stitched_ids(n), &
          'EPOCH repository differences', epoch_diff_bytes, &
          epoch_diff_bytes_padding, epoch_diff_bytes_mimetype, &
          epoch_diff_bytes_checksum_type, epoch_diff_bytes_checksum)
    ENDIF

    n = n + 1
    CALL sdf_safe_copy_id(h, 'epoch_source/info', stitched_ids(n))
    WRITE(time_string, '(I20)') epoch_bytes_compile_date

    string_array(1) = TRIM(epoch_bytes_git_version)
    string_array(2) = TRIM(epoch_bytes_compile_date_string)
    string_array(3) = TRIM(ADJUSTL(time_string))
    string_array(4) = TRIM(epoch_bytes_compile_machine_info)
    string_array(5) = TRIM(epoch_bytes_compiler_info)
    string_array(6) = TRIM(epoch_bytes_compiler_flags)

    CALL sdf_write_namevalue(h, stitched_ids(n), &
        'EPOCH repository information', &
        (/'git_version         ', &
          'compile_date_string ', &
          'compile_date_seconds', &
          'compile_machine_info', &
          'compiler_info       ', &
          'compiler_flags      '/), string_array)

    CALL sdf_write_stitched(h, 'epoch_source', 'EPOCH source', &
        stitched_ids(1), c_stagger_cell_centre, stitched_ids, n)

  END SUBROUTINE epoch_write_source_info



  SUBROUTINE write_source_info(h)

    TYPE(sdf_file_handle) :: h

    CALL sdf_write_source_info(h)
    CALL epoch_write_source_info(h)
    !CALL write_input_decks(h)

  END SUBROUTINE write_source_info

END MODULE diagnostics
