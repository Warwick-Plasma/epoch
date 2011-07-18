MODULE diagnostics

  USE calc_df
  USE sdf
  USE deck
  USE dist_fn
  USE encoded_source
  USE iterators
  USE mpi_subtype_control
  USE probes
  USE shared_data
  USE version_data
  USE setup

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_dt, output_routines, iterate_charge

  TYPE(sdf_file_handle) :: sdf_handle
  INTEGER(KIND=8), ALLOCATABLE :: species_offset(:)

CONTAINS

  SUBROUTINE output_routines(step)   ! step = step index

    INTEGER, INTENT(INOUT) :: step
    LOGICAL :: print_arrays, last_call
    CHARACTER(LEN=9+data_dir_max_length+n_zeros) :: filename, filename_desc
    CHARACTER(LEN=8) :: dump_type
    REAL(num), DIMENSION(:), ALLOCATABLE :: array
    INTEGER :: code, ispecies
    INTEGER, DIMENSION(c_ndims) :: dims
    LOGICAL :: restart_flag
    TYPE(particle_species), POINTER :: species

    CHARACTER(LEN=1), DIMENSION(3) :: dim_tags = (/'x', 'y', 'z'/)
    CHARACTER(LEN=2), DIMENSION(6) :: dir_tags = &
        (/'x-', 'x+', 'y-', 'y+', 'z-', 'z+'/)
    INTEGER, DIMENSION(6) :: fluxdir = &
        (/-c_dir_x, c_dir_x, -c_dir_y, c_dir_y, -c_dir_z, c_dir_z/)

#ifdef NO_IO
    RETURN
#endif

    IF (rank .EQ. 0 .AND. stdout_frequency .GT. 0 &
        .AND. MOD(step, stdout_frequency) .EQ. 0) THEN
      WRITE(*, '("Time", g20.12, " and iteration", i7, " after ", &
          & f8.1, " seconds")') time, step, MPI_WTIME() - walltime_start
    ENDIF

    CALL io_test(step, print_arrays, last_call)

    IF (.NOT.print_arrays) RETURN

    dims = (/nx_global/)

    ! Allows a maximum of 10^999 output dumps, should be enough for anyone
    ! (feel free to laugh when this isn't the case)
    WRITE(filename_desc, '("(a, ''/'', i", i3.3, ".", i3.3, ", ''.sdf'')")') &
        n_zeros, n_zeros
    WRITE(filename, filename_desc) TRIM(data_dir), output_file

    ! Always dump the variables with the "Every" attribute
    code = c_io_always

    ! Only dump variables with the "FULL" attributre on full dump intervals
    IF (MOD(output_file, full_dump_every) .EQ. 0 &
        .AND. full_dump_every .GT. -1) code = IOR(code, c_io_full)
    IF (MOD(output_file, restart_dump_every) .EQ. 0 &
        .AND. restart_dump_every .GT. -1) code = IOR(code, c_io_restartable)
    IF (last_call .AND. force_final_to_be_restartable) &
        code = IOR(code, c_io_restartable)

    CALL create_subtypes(code)

    ! Set a restart_flag to pass to the file header
    IF (IAND(code, c_io_restartable) .NE. 0) THEN
      restart_flag = .TRUE.
    ELSE
      restart_flag = .FALSE.
    ENDIF

    ALLOCATE(array(-2:nx+3))

    ! open the file
    CALL sdf_open(sdf_handle, filename, rank, comm, c_sdf_write)
    CALL sdf_write_header(sdf_handle, 'Epoch1d', 1, step, time, restart_flag, &
        jobid)
    CALL sdf_write_run_info(sdf_handle, c_version, c_revision, c_commit_id, &
        sha1sum, c_compile_machine, c_compile_flags, defines, c_compile_date, &
        run_date)

    CALL write_particle_grid(code)

    ! Write the cartesian mesh
    IF (IAND(dumpmask(c_dump_grid), code) .NE. 0) THEN
      IF (.NOT. use_offset_grid) THEN
        CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid', 'Grid/Grid', &
            xb_global)
      ELSE
        CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid', 'Grid/Grid', &
            xb_offset_global)
        CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid_full', &
            'Grid/Grid_Full', xb_global)
      ENDIF
    ENDIF

#ifdef PER_PARTICLE_WEIGHT
    CALL write_particle_variable(c_dump_part_weight, code, 'Weight', &
        iterate_weight)
#else
    IF (IAND(dumpmask(c_dump_part_weight), code) .NE. 0) THEN
      DO ispecies = 1, n_species
        species => species_list(ispecies)
        IF (IAND(species%dumpmask, code) .NE. 0 &
            .OR. IAND(code, c_io_restartable) .NE. 0) THEN
          CALL sdf_write_srl(sdf_handle, 'weight/' // TRIM(species%name), &
              'Particles/Weight/' // TRIM(species%name), species%weight)
        ENDIF
      ENDDO
    ENDIF
#endif
    CALL write_particle_variable(c_dump_part_px, code, 'Px', iterate_px)
    CALL write_particle_variable(c_dump_part_py, code, 'Py', iterate_py)
    CALL write_particle_variable(c_dump_part_pz, code, 'Pz', iterate_pz)

    CALL write_particle_variable(c_dump_part_vx, code, 'Vx', iterate_vx)
    CALL write_particle_variable(c_dump_part_vy, code, 'Vy', iterate_vy)
    CALL write_particle_variable(c_dump_part_vz, code, 'Vz', iterate_vz)

    CALL write_particle_variable(c_dump_part_charge, code, 'Q', iterate_charge)
    CALL write_particle_variable(c_dump_part_mass, code, 'Mass', iterate_mass)
#ifdef PARTICLE_DEBUG
    CALL write_particle_variable(c_dump_part_grid, code, 'Processor', &
        iterate_processor)
    CALL write_particle_variable(c_dump_part_grid, code, 'Processor_at_t0', &
        iterate_processor0)
#endif

    CALL write_field(c_dump_ex, code, 'ex', 'Electric Field/Ex', 'Volt/m', &
        c_stagger_ex, ex)
    CALL write_field(c_dump_ey, code, 'ey', 'Electric Field/Ey', 'Volt/m', &
        c_stagger_ey, ey)
    CALL write_field(c_dump_ez, code, 'ez', 'Electric Field/Ez', 'Volt/m', &
        c_stagger_ez, ez)

    CALL write_field(c_dump_bx, code, 'bx', 'Magnetic Field/Bx', 'Tesla', &
        c_stagger_bx, bx)
    CALL write_field(c_dump_by, code, 'by', 'Magnetic Field/By', 'Tesla', &
        c_stagger_by, by)
    CALL write_field(c_dump_bz, code, 'bz', 'Magnetic Field/Bz', 'Tesla', &
        c_stagger_bz, bz)

    CALL write_field(c_dump_jx, code, 'jx', 'Current/Jx', 'Amp', &
        c_stagger_jx, jx)
    CALL write_field(c_dump_jy, code, 'jy', 'Current/Jy', 'Amp', &
        c_stagger_jy, jy)
    CALL write_field(c_dump_jz, code, 'jz', 'Current/Jz', 'Amp', &
        c_stagger_jz, jz)

    ! These are derived variables from the particles
    CALL write_nspecies_field(c_dump_ekbar, code, 'ekbar', &
        'Derived/EkBar', '?', c_stagger_cell_centre, &
        calc_ekbar, array)

    CALL write_nspecies_field(c_dump_mass_density, code, 'mass_density', &
        'Derived/Mass_Density', '?', c_stagger_cell_centre, &
        calc_mass_density, array)

    CALL write_nspecies_field(c_dump_charge_density, code, 'charge_density', &
        'Derived/Charge_Density', '?', c_stagger_cell_centre, &
        calc_charge_density, array)

    CALL write_nspecies_field(c_dump_number_density, code, 'number_density', &
        'Derived/Number_Density', '?', c_stagger_cell_centre, &
        calc_number_density, array)

    CALL write_nspecies_field(c_dump_temperature, code, 'temperature', &
        'Derived/Temperature', '?', c_stagger_cell_centre, &
        calc_temperature, array)

    CALL write_nspecies_flux(c_dump_ekflux, code, 'ekflux', &
        'Derived/EkFlux', '?', c_stagger_cell_centre, &
        calc_ekflux, array, fluxdir, dir_tags)

    CALL write_nspecies_flux(c_dump_poynt_flux, code, 'poynt_flux', &
        'Derived/Poynting Flux', '?', c_stagger_cell_centre, &
        calc_poynt_flux, array, fluxdir(1:3), dim_tags)

#ifdef FIELD_DEBUG
    array = rank
    CALL sdf_write_plain_variable(sdf_handle, 'rank', 'Processor/Rank', '', &
        dims, c_stagger_cell_centre, 'grid', array, subtype_field, &
        subarray_field)
#endif

    IF (IAND(dumpmask(c_dump_dist_fns), code) .NE. 0) THEN
      CALL write_dist_fns(sdf_handle, code)
    ENDIF

#ifdef PARTICLE_PROBES
    IF (IAND(dumpmask(c_dump_probes), code) .NE. 0) THEN
      CALL write_probes(sdf_handle, code)
    ENDIF
#endif

    IF (restart_flag) THEN
      IF (dump_input_decks) CALL write_input_decks(sdf_handle)
      IF (dump_source_code .AND. SIZE(source_code) .GT. 0) &
          CALL sdf_write_source_code(sdf_handle, "code", &
              "base64_packed_source_code", source_code, last_line, 0)
    ENDIF

    ! close the file
    CALL sdf_close(sdf_handle)

    IF (rank .EQ. 0) THEN
      dump_type = 'normal'
      CALL append_filename(dump_type, output_file)
      IF (IAND(code, c_io_restartable) .NE. 0) THEN
        dump_type = 'restart'
        CALL append_filename(dump_type, output_file)
      ENDIF
      IF (IAND(code, c_io_full) .NE. 0) THEN
        dump_type = 'full'
        CALL append_filename(dump_type, output_file)
      ENDIF
      WRITE(stat_unit, '("Wrote ", a7, " dump number", i5, " at time", g20.12, &
          & " and iteration", i7)') dump_type, output_file, time, step
      CALL flush_stat_file()
    ENDIF

    output_file = output_file + 1

    DEALLOCATE(array)
    IF (ALLOCATED(species_offset)) DEALLOCATE(species_offset)
    CALL free_subtypes()

  END SUBROUTINE output_routines



  SUBROUTINE append_filename(listname, output_file)
    ! This routine updates a list of each output type (eg. full, dump, normal)
    ! that can be passed to VisIt to filter the output.

    CHARACTER(LEN=*), INTENT(IN) :: listname
    INTEGER, INTENT(IN) :: output_file
    CHARACTER(LEN=data_dir_max_length+64) :: listfile, filename_desc
    INTEGER :: ierr
    LOGICAL :: exists

    WRITE(filename_desc, '("(i", i3.3, ".", i3.3, ", ''.sdf'')")') &
        n_zeros, n_zeros
    listfile = TRIM(data_dir) // '/' // TRIM(listname) // '.visit'

    INQUIRE(file=listfile, exist=exists)
    IF (exists) THEN
      OPEN(unit=lu, status='OLD', position='APPEND', file=listfile, iostat=ierr)
    ELSE
      OPEN(unit=lu, status='NEW', file=listfile, iostat=errcode)
    ENDIF

    WRITE(lu,filename_desc) output_file
    CLOSE(lu)

  END SUBROUTINE append_filename



  SUBROUTINE io_test(step, print_arrays, last_call)

    INTEGER, INTENT(IN) :: step
    LOGICAL, INTENT(OUT) :: print_arrays, last_call
    INTEGER :: id
    REAL(num) :: t0, t1, time_first

    print_arrays = .FALSE.
    last_call = .FALSE.

    ! Work out the time that the next dump will occur based on the
    ! current timestep
    t0 = HUGE(1.0_num)
    t1 = HUGE(1.0_num)
    IF (dt_snapshot .GE. 0.0_num) t0 = time_next
    IF (nstep_snapshot .GE. 0) t1 = time + dt * (nstep_next - step)

    IF (t0 .LT. t1) THEN
      ! Next I/O dump based on dt_snapshot
      time_first = t0
      IF (dt_snapshot .GT. 0 .AND. time .GE. time_next) THEN
        time_next  = time_next + dt_snapshot
        print_arrays = .TRUE.
      ENDIF
    ELSE
      ! Next I/O dump based on nstep_snapshot
      time_first = t1
      IF (nstep_snapshot .GT. 0 .AND. step .GE. nstep_next) THEN
        nstep_next = nstep_next + nstep_snapshot
        print_arrays = .TRUE.
      ENDIF
    ENDIF

    DO id = 1, num_vars_to_dump
      IF (IAND(dumpmask(id), c_io_averaged) .NE. 0) THEN
        IF (time .GE. time_first - average_time) THEN
          CALL average_field(id)
        ENDIF
      ENDIF
    ENDDO

    IF ((time .GE. t_end .OR. step .EQ. nsteps) .AND. dump_last) THEN
      last_call = .TRUE.
      print_arrays = .TRUE.
    ENDIF

  END SUBROUTINE io_test



  SUBROUTINE average_field(ioutput)

    INTEGER, INTENT(IN) :: ioutput
    INTEGER :: n_species_local, ispecies, species_sum
    REAL(num), DIMENSION(:), ALLOCATABLE :: array

    averaged_data(ioutput)%real_time = averaged_data(ioutput)%real_time + dt

    species_sum = 0
    n_species_local = 0
    IF (IAND(dumpmask(ioutput), c_io_no_sum) .EQ. 0) &
        species_sum = 1
    IF (IAND(dumpmask(ioutput), c_io_species) .NE. 0) &
        n_species_local = n_species

    n_species_local = n_species_local + species_sum

    SELECT CASE(ioutput)
    CASE(c_dump_ex)
      averaged_data(ioutput)%array(:,1) = &
          averaged_data(ioutput)%array(:,1) + ex * dt
    CASE(c_dump_ey)
      averaged_data(ioutput)%array(:,1) = &
          averaged_data(ioutput)%array(:,1) + ey * dt
    CASE(c_dump_ez)
      averaged_data(ioutput)%array(:,1) = &
          averaged_data(ioutput)%array(:,1) + ez * dt
    CASE(c_dump_bx)
      averaged_data(ioutput)%array(:,1) = &
          averaged_data(ioutput)%array(:,1) + bx * dt
    CASE(c_dump_by)
      averaged_data(ioutput)%array(:,1) = &
          averaged_data(ioutput)%array(:,1) + by * dt
    CASE(c_dump_bz)
      averaged_data(ioutput)%array(:,1) = &
          averaged_data(ioutput)%array(:,1) + bz * dt
    CASE(c_dump_jx)
      averaged_data(ioutput)%array(:,1) = &
          averaged_data(ioutput)%array(:,1) + jx * dt
    CASE(c_dump_jy)
      averaged_data(ioutput)%array(:,1) = &
          averaged_data(ioutput)%array(:,1) + jy * dt
    CASE(c_dump_jz)
      averaged_data(ioutput)%array(:,1) = &
          averaged_data(ioutput)%array(:,1) + jz * dt
    CASE(c_dump_ekbar)
      ALLOCATE(array(-2:nx+3))
      DO ispecies = 1, n_species_local
        CALL calc_ekbar(array, ispecies-species_sum)
        averaged_data(ioutput)%array(:,ispecies) = &
            averaged_data(ioutput)%array(:,ispecies) + array * dt
      ENDDO
      DEALLOCATE(array)
    CASE(c_dump_mass_density)
      ALLOCATE(array(-2:nx+3))
      DO ispecies = 1, n_species_local
        CALL calc_mass_density(array, ispecies-species_sum)
        averaged_data(ioutput)%array(:,ispecies) = &
            averaged_data(ioutput)%array(:,ispecies) + array * dt
      ENDDO
      DEALLOCATE(array)
    CASE(c_dump_charge_density)
      ALLOCATE(array(-2:nx+3))
      DO ispecies = 1, n_species_local
        CALL calc_charge_density(array, ispecies-species_sum)
        averaged_data(ioutput)%array(:,ispecies) = &
            averaged_data(ioutput)%array(:,ispecies) + array * dt
      ENDDO
      DEALLOCATE(array)
    CASE(c_dump_number_density)
      ALLOCATE(array(-2:nx+3))
      DO ispecies = 1, n_species_local
        CALL calc_number_density(array, ispecies-species_sum)
        averaged_data(ioutput)%array(:,ispecies) = &
            averaged_data(ioutput)%array(:,ispecies) + array * dt
      ENDDO
      DEALLOCATE(array)
    CASE(c_dump_temperature)
      ALLOCATE(array(-2:nx+3))
      DO ispecies = 1, n_species_local
        CALL calc_temperature(array, ispecies-species_sum)
        averaged_data(ioutput)%array(:,ispecies) = &
            averaged_data(ioutput)%array(:,ispecies) + array * dt
      ENDDO
      DEALLOCATE(array)
    END SELECT

  END SUBROUTINE average_field



  SUBROUTINE set_dt        ! sets CFL limited step

    dt = dx / c
    IF (dt_plasma_frequency .NE. 0.0_num) dt = MIN(dt, dt_plasma_frequency)
    IF (dt_laser .NE. 0.0_num) dt = MIN(dt, dt_laser)
    dt = dt_multiplier * dt

    IF (.NOT. any_average) RETURN

    average_time = MAX(dt_average, dt * nstep_average)

    IF (dt_min_average .GT. 0) THEN
      IF (dt_min_average .LT. dt) THEN
        IF (rank .EQ. 0) THEN
          PRINT*,'*** WARNING ***'
          PRINT*,'Time step is too small to satisfy "nstep_average"'
          PRINT*,'Averaging will occur over fewer time steps than specified'
          PRINT*,'Set "dt_multiplier" less than ', &
              dt_multiplier * dt_min_average / dt, &
              ' to fix this'
        ENDIF
        dt_min_average = -1
      ENDIF
    ENDIF

  END SUBROUTINE set_dt



  SUBROUTINE energy_account()

  END SUBROUTINE energy_account



  SUBROUTINE write_field(id, code, block_id, name, units, stagger, array)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: block_id, name, units
    INTEGER, INTENT(IN) :: stagger
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER :: should_dump

    IF (IAND(dumpmask(id), code) .EQ. 0) RETURN

    dims = (/nx_global/)

    ! Want the code to output unaveraged data if either restarting or
    ! requested by the user
    should_dump = IOR(c_io_snapshot, IAND(code,c_io_restartable))
    should_dump = IOR(should_dump, NOT(c_io_averaged))

    IF (IAND(dumpmask(id), should_dump) .NE. 0) THEN
      CALL sdf_write_plain_variable(sdf_handle, &
          TRIM(block_id), TRIM(name), &
          TRIM(units), dims, stagger, 'grid', &
          array, subtype_field, subarray_field)
    ENDIF

    IF (IAND(dumpmask(id), c_io_averaged) .NE. 0) THEN
      averaged_data(id)%array = averaged_data(id)%array &
          / averaged_data(id)%real_time

      CALL sdf_write_plain_variable(sdf_handle, &
          TRIM(block_id) // '_averaged', TRIM(name) // '_averaged', &
          TRIM(units), dims, stagger, 'grid', &
          averaged_data(id)%array(:,1), subtype_field, subarray_field)

      averaged_data(id)%real_time = 0.0_num
      averaged_data(id)%array = 0.0_num
    ENDIF

  END SUBROUTINE write_field



  SUBROUTINE write_nspecies_field(id, code, block_id, name, units, stagger, &
      func, array)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: block_id, name, units
    INTEGER, INTENT(IN) :: stagger
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER :: ispecies, should_dump, species_sum
    CHARACTER(LEN=c_max_string_length) :: temp_block_id, temp_name

    INTERFACE
      SUBROUTINE func(data_array, current_species)
        USE shared_data
        REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
        INTEGER, INTENT(IN) :: current_species
      END SUBROUTINE func
    END INTERFACE

    IF (IAND(dumpmask(id), code) .EQ. 0) RETURN

    dims = (/nx_global/)

    ! Want the code to output unaveraged data if either restarting or
    ! requested by the user
    should_dump = IOR(c_io_snapshot, IAND(code,c_io_restartable))
    should_dump = IOR(should_dump, NOT(c_io_averaged))

    IF (IAND(dumpmask(id), should_dump) .NE. 0) THEN
      IF (IAND(dumpmask(id), c_io_no_sum) .EQ. 0) THEN
        CALL func(array, 0)
        CALL sdf_write_plain_variable(sdf_handle, &
            TRIM(ADJUSTL(block_id)), TRIM(ADJUSTL(name)), &
            TRIM(units), dims, stagger, 'grid', &
            array, subtype_field, subarray_field)
      ENDIF

      IF (IAND(dumpmask(id), c_io_species) .NE. 0) THEN
        DO ispecies = 1, n_species
#ifdef TRACER_PARTICLES
          IF (species_list(ispecies)%tracer) CYCLE
#endif
          WRITE(temp_block_id, '(a, "_", a)') TRIM(block_id), &
              TRIM(species_list(ispecies)%name)
          WRITE(temp_name, '(a, "_", a)') TRIM(name), &
              TRIM(species_list(ispecies)%name)
          CALL func(array, ispecies)
          CALL sdf_write_plain_variable(sdf_handle, &
              TRIM(ADJUSTL(temp_block_id)), TRIM(ADJUSTL(temp_name)), &
              TRIM(units), dims, stagger, 'grid', &
              array, subtype_field, subarray_field)
        ENDDO
      ENDIF
    ENDIF

    ! Write averaged data
    IF (IAND(dumpmask(id), c_io_averaged) .NE. 0) THEN
      averaged_data(id)%array = averaged_data(id)%array &
          / averaged_data(id)%real_time

      species_sum = 0
      IF (IAND(dumpmask(id), c_io_no_sum) .EQ. 0) THEN
        species_sum = 1
        CALL sdf_write_plain_variable(sdf_handle, &
            TRIM(block_id) // '_averaged', TRIM(name) // '_averaged', &
            TRIM(units), dims, stagger, 'grid', &
            averaged_data(id)%array(:,1), subtype_field, subarray_field)
      ENDIF

      IF (IAND(dumpmask(id), c_io_species) .NE. 0) THEN
        DO ispecies = 1, n_species
          WRITE(temp_block_id, '(a, "_", a, "_averaged")') TRIM(block_id), &
              TRIM(species_list(ispecies)%name)
          WRITE(temp_name, '(a, "_", a, "_averaged")') TRIM(name), &
              TRIM(species_list(ispecies)%name)
          CALL sdf_write_plain_variable(sdf_handle, &
              TRIM(ADJUSTL(temp_block_id)), TRIM(ADJUSTL(temp_name)), &
              TRIM(units), dims, stagger, 'grid', &
              averaged_data(id)%array(:,ispecies+species_sum), &
              subtype_field, subarray_field)
        ENDDO
      ENDIF

      averaged_data(id)%real_time = 0.0_num
      averaged_data(id)%array = 0.0_num
    ENDIF

  END SUBROUTINE write_nspecies_field



  SUBROUTINE write_nspecies_flux(id, code, block_id, name, units, stagger, &
      func, array, fluxdir, dir_tags)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: block_id, name, units
    INTEGER, INTENT(IN) :: stagger
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, DIMENSION(:), INTENT(IN) :: fluxdir
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: dir_tags
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER :: ispecies, should_dump, ndirs, idir
    CHARACTER(LEN=c_max_string_length) :: temp_block_id, temp_name

    INTERFACE
      SUBROUTINE func(data_array, current_species, direction)
        USE shared_data
        REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
        INTEGER, INTENT(IN) :: current_species, direction
      END SUBROUTINE func
    END INTERFACE

    IF (IAND(dumpmask(id), code) .EQ. 0) RETURN

    ndirs = SIZE(fluxdir)
    dims = (/nx_global/)

    ! Want the code to output unaveraged data if either restarting or
    ! requested by the user
    should_dump = IOR(c_io_snapshot, IAND(code,c_io_restartable))
    should_dump = IOR(should_dump, NOT(c_io_averaged))

    IF (IAND(dumpmask(id), should_dump) .NE. 0) THEN
      IF (IAND(dumpmask(id), c_io_no_sum) .EQ. 0) THEN
        DO idir = 1, ndirs
          WRITE(temp_block_id, '(a, "_", a)') TRIM(block_id), &
              TRIM(dir_tags(idir))
          WRITE(temp_name, '(a, "_", a)') TRIM(name), &
              TRIM(dir_tags(idir))
          CALL func(array, 0, fluxdir(idir))
          CALL sdf_write_plain_variable(sdf_handle, &
              TRIM(ADJUSTL(temp_block_id)), TRIM(ADJUSTL(temp_name)), &
              TRIM(units), dims, stagger, 'grid', &
              array, subtype_field, subarray_field)
        ENDDO
      ENDIF

      IF (IAND(dumpmask(id), c_io_species) .NE. 0) THEN
        DO ispecies = 1, n_species
#ifdef TRACER_PARTICLES
          IF (species_list(ispecies)%tracer) CYCLE
#endif
          DO idir = 1, ndirs
            WRITE(temp_block_id, '(a, "_", a, "_", a)') TRIM(block_id), &
                TRIM(species_list(ispecies)%name), TRIM(dir_tags(idir))
            WRITE(temp_name, '(a, "_", a, "_", a)') TRIM(name), &
                TRIM(species_list(ispecies)%name), TRIM(dir_tags(idir))
            CALL func(array, ispecies, fluxdir(idir))
            CALL sdf_write_plain_variable(sdf_handle, &
                TRIM(ADJUSTL(temp_block_id)), TRIM(ADJUSTL(temp_name)), &
                TRIM(units), dims, stagger, 'grid', &
                array, subtype_field, subarray_field)
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    ! Flux variables not currently averaged

  END SUBROUTINE write_nspecies_flux



  SUBROUTINE species_offset_init

    INTEGER(KIND=8) :: species_count
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: npart_species_per_proc
    INTEGER :: i, ispecies

    IF (ALLOCATED(species_offset)) RETURN

    ALLOCATE(npart_species_per_proc(nproc))
    ALLOCATE(species_offset(n_species))
    species_offset = 0

    DO ispecies = 1, n_species
      CALL MPI_ALLGATHER(species_list(ispecies)%attached_list%count, 1, &
          MPI_INTEGER8, npart_species_per_proc, 1, MPI_INTEGER8, comm, errcode)
      species_count = 0
      DO i = 1, nproc
        IF (rank .EQ. i-1) species_offset(ispecies) = species_count
        species_count = species_count + npart_species_per_proc(i)
      ENDDO
      species_list(ispecies)%count = species_count
    ENDDO

    DEALLOCATE(npart_species_per_proc)

  END SUBROUTINE species_offset_init



  SUBROUTINE write_particle_grid(code)

    INTEGER, INTENT(IN) :: code
    INTEGER :: ispecies

    IF (IAND(dumpmask(c_dump_part_grid), code) .EQ. 0) RETURN

    CALL species_offset_init()

    CALL start_particle_species_only(current_species)

    DO ispecies = 1, n_species
      IF (IAND(current_species%dumpmask, code) .NE. 0 &
          .OR. IAND(code, c_io_restartable) .NE. 0) THEN
        CALL sdf_write_point_mesh(sdf_handle, &
            'grid/' // TRIM(current_species%name), &
            'Grid/Point/' // TRIM(current_species%name), &
            species_list(ispecies)%count, c_dimension_1d, &
            iterate_particles, species_offset(ispecies))
      ENDIF

      CALL advance_particle_species_only(current_species)
    ENDDO

  END SUBROUTINE write_particle_grid



  SUBROUTINE write_particle_variable(id, code, name, iterator)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: name

    INTERFACE
      FUNCTION iterator(array, npart_it, start)
        USE shared_data
        REAL(num) :: iterator
        REAL(num), DIMENSION(:), INTENT(OUT) :: array
        INTEGER, INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
      END FUNCTION iterator
    END INTERFACE

    INTEGER :: ispecies

    IF (IAND(dumpmask(id), code) .EQ. 0) RETURN

    CALL species_offset_init()

    CALL start_particle_species_only(current_species)

    DO ispecies = 1, n_species
      IF (IAND(current_species%dumpmask, code) .NE. 0 &
          .OR. IAND(code, c_io_restartable) .NE. 0) THEN
        CALL sdf_write_point_variable(sdf_handle, &
            lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
            'Particles/' // TRIM(current_species%name) // '/' // &
            TRIM(name), '', &
            species_list(ispecies)%count, &
            'grid/' // TRIM(current_species%name), &
            iterator, species_offset(ispecies))
      ENDIF

      CALL advance_particle_species_only(current_species)
    ENDDO

  END SUBROUTINE write_particle_variable



  FUNCTION lowercase(string_in) RESULT(string_out)

    CHARACTER(LEN=*), PARAMETER :: lwr = 'abcdefghijklmnopqrstuvwxyz'
    CHARACTER(LEN=*), PARAMETER :: upr = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    CHARACTER(LEN=*), INTENT(IN) :: string_in
    CHARACTER(LEN=LEN(string_in)) :: string_out
    INTEGER :: i, idx

    string_out = string_in

    DO i = 1, LEN(string_out)
      idx = INDEX(upr, string_out(i:i))
      IF (idx .NE. 0) string_out(i:i) = lwr(idx:idx)
    ENDDO

  END FUNCTION lowercase

END MODULE diagnostics
