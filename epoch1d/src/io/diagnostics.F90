MODULE diagnostics

  USE calc_df
  USE cfd
  USE deck
  USE dist_fn
  USE encoded_source
  USE iterators
  USE mpi_subtype_control
  USE probes
  USE shared_data
  USE version_data

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_dt, output_routines, iterate_charge

  TYPE(cfd_file_handle) :: cfd_handle
  INTEGER(KIND=8), ALLOCATABLE :: species_offset(:)

CONTAINS

  SUBROUTINE output_routines(i)   ! i = step index

    INTEGER, INTENT(INOUT) :: i
    LOGICAL :: print_arrays, last_call
    CHARACTER(LEN=9+data_dir_max_length+n_zeros) :: filename, filename_desc
    CHARACTER(LEN=c_max_string_length) :: temp_name
    CHARACTER(LEN=8) :: dump_type
    REAL(num), DIMENSION(:), ALLOCATABLE :: array
    REAL(num), DIMENSION(c_ndims) :: stagger = 0.0_num
    INTEGER :: ispecies, code
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER :: restart_flag

    IF (rank .EQ. 0 .AND. stdout_frequency .GT. 0 &
        .AND. MOD(i, stdout_frequency) .EQ. 0) THEN
      WRITE(*, '("Time", g20.12, " and iteration", i7, " after ", &
          & f8.1, " seconds")') time, i, MPI_WTIME() - walltime_start
    ENDIF

    CALL io_test(i, print_arrays, last_call)

    IF (.NOT.print_arrays) RETURN

    dims = (/nx_global/)

    ! Allows a maximum of 10^999 output dumps, should be enough for anyone
    ! (feel free to laugh when this isn't the case)
    WRITE(filename_desc, '("(a, ''/'', i", i3.3, ".", i3.3, ", ''.cfd'')")') &
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

    CALL create_subtypes(IAND(code, c_io_restartable) .NE. 0)

    ! Set a restart_flag to pass to the file header
    IF (IAND(code, c_io_restartable) .NE. 0) THEN
      restart_flag = 1
    ELSE
      restart_flag = 0
    ENDIF

    ALLOCATE(array(-2:nx+3))

    ! open the file
    ! (filename, rank of current process, MPI communicator (can be
    ! MPI_COMM_WORLD), file mode (c_cfd_read or c_cfd_write),
    ! cycle number, simulation time, job id)
    CALL cfd_open(cfd_handle, filename, rank, comm, c_cfd_write, i, time, jobid)
    CALL cfd_write_job_info(cfd_handle, c_code_io_version, c_version, &
        c_revision, defines, c_compile_date, run_date, restart_flag, &
        c_code_name, c_commit_id, sha1sum, c_compile_machine, &
        c_compile_flags, 0)

    CALL write_particle_grid(code)

    ! Write the cartesian mesh
    ! (mesh name, mesh class, x_array, y_array, rank used for writing)
    IF (IAND(dumpmask(c_dump_grid), code) .NE. 0) THEN
      IF (.NOT. use_offset_grid) THEN
        CALL cfd_write_1d_cartesian_grid(cfd_handle, "Grid", "Grid", &
            x_global(1:nx_global), 0)
      ELSE
        CALL cfd_write_1d_cartesian_grid(cfd_handle, "Grid", "Grid", &
            x_offset_global(1:nx_global), 0)
        CALL cfd_write_1d_cartesian_grid(cfd_handle, "Grid_Full", "Grid", &
            x_global(1:nx_global), 0)
      ENDIF
    ENDIF

#ifdef PER_PARTICLE_WEIGHT
    CALL write_particle_variable(c_dump_part_weight, code, 'Weight', &
        iterate_weight)
#else
    IF (IAND(dumpmask(c_dump_part_weight), code) .NE. 0) &
        CALL cfd_write_real_constant("Weight", "Particles", weight, 0)
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

    CALL write_field(c_dump_ex, code, 'Ex', 'Electric Field', ex)
    CALL write_field(c_dump_ey, code, 'Ey', 'Electric Field', ey)
    CALL write_field(c_dump_ez, code, 'Ez', 'Electric Field', ez)

    CALL write_field(c_dump_bx, code, 'Bx', 'Magnetic Field', bx)
    CALL write_field(c_dump_by, code, 'By', 'Magnetic Field', by)
    CALL write_field(c_dump_bz, code, 'Bz', 'Magnetic Field', bz)

    CALL write_field(c_dump_jx, code, 'Jx', 'Current', jx)
    CALL write_field(c_dump_jy, code, 'Jy', 'Current', jy)
    CALL write_field(c_dump_jz, code, 'Jz', 'Current', jz)

    ! These are derived variables from the particles
    CALL write_nspecies_field(c_dump_ekbar, code, &
        'EkBar', 'EkBar', calc_ekbar, array)

    CALL write_nspecies_field(c_dump_mass_density, code, &
        'Mass_density', 'Derived', calc_mass_density, array)

    CALL write_nspecies_field(c_dump_charge_density, code, &
        'Charge_density', 'Derived', calc_charge_density, array)

    CALL write_nspecies_field(c_dump_number_density, code, &
        'Number_density', 'Derived', calc_number_density, array)

    CALL write_nspecies_field(c_dump_temperature, code, &
        'Temperature', 'Derived', calc_temperature, array)

#ifdef FIELD_DEBUG
    array = rank
    CALL cfd_write_1d_cartesian_variable_parallel(cfd_handle, "Rank", &
        "Processor", dims, stagger, "Grid", "Grid", array, subtype_field, &
        subarray_field)
#endif

    IF (IAND(dumpmask(c_dump_dist_fns), code) .NE. 0) THEN
      CALL write_dist_fns(cfd_handle, code)
    ENDIF

#ifdef PARTICLE_PROBES
    IF (IAND(dumpmask(c_dump_probes), code) .NE. 0) THEN
      CALL write_probes(cfd_handle, code)
    ENDIF
#endif

    IF (restart_flag .EQ. 1 .AND. LEN(source_code) .GT. 0) THEN
      CALL write_input_decks(cfd_handle)
      CALL cfd_write_source_code(cfd_handle, "Code", &
          "base64_packed_source_code", source_code, last_line, 0)
    ENDIF

    ! close the file
    CALL cfd_close(cfd_handle)

    output_file = output_file + 1
    IF (rank .EQ. 0) THEN
      IF (IAND(code, c_io_restartable) .NE. 0) THEN
        dump_type = "restart"
      ELSE IF (IAND(code, c_io_full) .NE. 0) THEN
        dump_type = "full"
      ELSE
        dump_type = "normal"
      ENDIF
      WRITE(20, '("Wrote ", a7, " dump number", i5, " at time", g20.12, &
          & " and iteration", i7)') dump_type, output_file-1, time, i
      CALL FLUSH(20)
    ENDIF

    DEALLOCATE(array)
    IF (ALLOCATED(species_offset)) DEALLOCATE(species_offset)
    CALL free_subtypes()

  END SUBROUTINE output_routines



  SUBROUTINE io_test(i, print_arrays, last_call)

    INTEGER, INTENT(IN) :: i
    LOGICAL, INTENT(OUT) :: print_arrays, last_call
    INTEGER :: ioutput

    REAL(num), SAVE :: t1 = 0.0_num
    LOGICAL, SAVE :: first = .TRUE.

    IF (first) THEN
      IF (ic_from_restart) t1 = time
      first = .FALSE.
    ENDIF

    print_arrays = .FALSE.
    last_call = .FALSE.

    DO ioutput = 1, num_vars_to_dump
      IF (IAND(dumpmask(ioutput), c_io_averaged) .NE. 0 &
          .AND. (time &
          .GE. t1 - averaged_data(ioutput)%average_over_real_time)) THEN
        CALL average_field(ioutput)
      ENDIF
    ENDDO

    IF (time .GE. t1) THEN
      print_arrays = .TRUE.
      t1 = t1 + dt_snapshots
    ENDIF

    IF (time .GE. t_end .OR. i .EQ. nsteps) THEN
      last_call = .TRUE.
      print_arrays = .TRUE.
    ENDIF

  END SUBROUTINE io_test



  SUBROUTINE average_field(ioutput)

    INTEGER, INTENT(IN) :: ioutput
    INTEGER :: n_species_local, ispecies
    REAL(num), DIMENSION(:), ALLOCATABLE :: array

    averaged_data(ioutput)%real_time_after_average = &
        averaged_data(ioutput)%real_time_after_average + dt

    n_species_local = 1
    IF (IAND(dumpmask(ioutput), c_io_species) .NE. 0) &
        n_species_local = n_species + 1

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
        CALL calc_ekbar(array, ispecies-1)
        averaged_data(ioutput)%array(:,ispecies) = &
            averaged_data(ioutput)%array(:,ispecies) + array * dt
      ENDDO
      DEALLOCATE(array)
    CASE(c_dump_mass_density)
      ALLOCATE(array(-2:nx+3))
      DO ispecies = 1, n_species_local
        CALL calc_mass_density(array, ispecies-1)
        averaged_data(ioutput)%array(:,ispecies) = &
            averaged_data(ioutput)%array(:,ispecies) + array * dt
      ENDDO
      DEALLOCATE(array)
    CASE(c_dump_charge_density)
      ALLOCATE(array(-2:nx+3))
      DO ispecies = 1, n_species_local
        CALL calc_charge_density(array, ispecies-1)
        averaged_data(ioutput)%array(:,ispecies) = &
            averaged_data(ioutput)%array(:,ispecies) + array * dt
      ENDDO
      DEALLOCATE(array)
    CASE(c_dump_number_density)
      ALLOCATE(array(-2:nx+3))
      DO ispecies = 1, n_species_local
        CALL calc_number_density(array, ispecies-1)
        averaged_data(ioutput)%array(:,ispecies) = &
            averaged_data(ioutput)%array(:,ispecies) + array * dt
      ENDDO
      DEALLOCATE(array)
    CASE(c_dump_temperature)
      ALLOCATE(array(-2:nx+3))
      DO ispecies = 1, n_species_local
        CALL calc_temperature(array, ispecies-1)
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
    IF (dt_min_average .GT. 0.0_num) dt = MIN(dt, dt_min_average)
    dt = dt_multiplier * dt

  END SUBROUTINE set_dt



  SUBROUTINE energy_account()

  END SUBROUTINE energy_account



  SUBROUTINE write_field(id, code, name, class, array)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    REAL(num), DIMENSION(c_ndims) :: stagger = 0.0_num
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER :: should_dump

    IF (IAND(dumpmask(id), code) .EQ. 0) RETURN

    dims = (/nx_global/)

    ! Want the code to output unaveraged data if either restarting or
    ! requested by the user
    should_dump = IOR(c_io_snapshot, IAND(code,c_io_restartable))
    should_dump = IOR(should_dump, NOT(c_io_averaged))

    ! (variable name, variable class, global grid dimensions,
    ! grid stagger, mesh name, mesh class, variable,
    ! mpi type describing data distribution)
    IF (IAND(dumpmask(id), should_dump) .NE. 0) THEN
      CALL cfd_write_1d_cartesian_variable_parallel(cfd_handle, &
          TRIM(name), TRIM(class), dims, stagger, &
          'Grid', 'Grid', array, subtype_field, subarray_field)
    ENDIF

    IF (IAND(dumpmask(id), c_io_averaged) .NE. 0) THEN
      averaged_data(id)%array = averaged_data(id)%array &
          / averaged_data(id)%real_time_after_average

      CALL cfd_write_1d_cartesian_variable_parallel(cfd_handle, &
          TRIM(name) // '_averaged', TRIM(class), dims, stagger, &
          'Grid', 'Grid', averaged_data(id)%array(:,1), &
          subtype_field, subarray_field)

      averaged_data(id)%real_time_after_average = 0.0_num
      averaged_data(id)%array = 0.0_num
    ENDIF

  END SUBROUTINE write_field



  SUBROUTINE write_nspecies_field(id, code, name, class, func, array)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    REAL(num), DIMENSION(c_ndims) :: stagger = 0.0_num
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER :: ispecies, should_dump
    CHARACTER(LEN=c_max_string_length) :: temp_name

    INTERFACE
      SUBROUTINE func(data_array, current_species)
        USE shared_data
        REAL(num), DIMENSION(-2:), INTENT(INOUT) :: data_array
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
      IF (IAND(dumpmask(id), c_io_no_intrinsic) .EQ. 0) THEN
        CALL func(array, 0)
        CALL cfd_write_1d_cartesian_variable_parallel(cfd_handle, &
            TRIM(name), TRIM(class), dims, stagger, &
            'Grid', 'Grid', array, subtype_field, subarray_field)
      ENDIF

      IF (IAND(dumpmask(id), c_io_species) .NE. 0) THEN
        DO ispecies = 1, n_species
          CALL func(array, ispecies)
          WRITE(temp_name, '(a, "_", a)') TRIM(name), &
              TRIM(particle_species(ispecies)%name)
          CALL cfd_write_1d_cartesian_variable_parallel(cfd_handle, &
              TRIM(ADJUSTL(temp_name)), TRIM(class), dims, stagger, &
              'Grid', 'Grid', array, subtype_field, subarray_field)
        ENDDO
      ENDIF
    ENDIF

    ! Write averaged data
    IF (IAND(dumpmask(id), c_io_averaged) .NE. 0) THEN
      averaged_data(id)%array = averaged_data(id)%array &
          / averaged_data(id)%real_time_after_average

      IF (IAND(dumpmask(id), c_io_no_intrinsic) .EQ. 0) THEN
        CALL cfd_write_1d_cartesian_variable_parallel(cfd_handle, &
            TRIM(name) // '_averaged', TRIM(class), dims, stagger, &
            'Grid', 'Grid', averaged_data(id)%array(:,1), &
            subtype_field, subarray_field)
      ENDIF

      IF (IAND(dumpmask(id), c_io_species) .NE. 0) THEN
        DO ispecies = 1, n_species
          WRITE(temp_name, '(a, "_", a, "_averaged")') TRIM(name), &
              TRIM(particle_species(ispecies)%name)
          CALL cfd_write_1d_cartesian_variable_parallel(cfd_handle, &
              TRIM(ADJUSTL(temp_name)), TRIM(class), dims, stagger, &
              'Grid', 'Grid', averaged_data(id)%array(:,ispecies+1), &
              subtype_field, subarray_field)
        ENDDO
      ENDIF

      averaged_data(id)%real_time_after_average = 0.0_num
      averaged_data(id)%array = 0.0_num
    ENDIF

  END SUBROUTINE write_nspecies_field



  SUBROUTINE species_offset_init

    INTEGER(KIND=8) :: species_count
    INTEGER(KIND=8), ALLOCATABLE :: npart_species_per_proc(:)
    INTEGER :: i, ispecies

    IF (ALLOCATED(species_offset)) RETURN

    ALLOCATE(npart_species_per_proc(nproc))
    ALLOCATE(species_offset(n_species))
    species_offset = 0

    DO ispecies = 1, n_species
      CALL MPI_ALLGATHER(particle_species(ispecies)%attached_list%count, 1, &
          MPI_INTEGER8, npart_species_per_proc, 1, MPI_INTEGER8, comm, errcode)
      species_count = 0
      DO i = 1, nproc
        IF (rank .EQ. i-1) species_offset(ispecies) = species_count
        species_count = species_count + npart_species_per_proc(i)
      ENDDO
      particle_species(ispecies)%count = species_count
    ENDDO

    DEALLOCATE(npart_species_per_proc)

  END SUBROUTINE species_offset_init



  SUBROUTINE write_particle_grid(code)

    INTEGER, INTENT(IN) :: code
    INTEGER :: ispecies

    IF (IAND(dumpmask(c_dump_part_grid), code) .EQ. 0) RETURN

    CALL species_offset_init()

    CALL start_particle_family_only(current_family)

    DO ispecies = 1, n_species
      IF (current_family%dump .OR. IAND(code, c_io_restartable) .NE. 0) THEN
        CALL cfd_write_nd_particle_grid_with_iterator_all(cfd_handle, &
            'Particles_' // TRIM(current_family%name), 'Grid', &
            iterate_particles, c_dimension_1d, &
            particle_species(ispecies)%count, npart_per_it, &
            c_particle_cartesian, species_offset(ispecies))
      ENDIF

      CALL advance_particle_family_only(current_family)
    ENDDO

  END SUBROUTINE write_particle_grid



  SUBROUTINE write_particle_variable(id, code, name, iterator)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: name

    INTERFACE
      SUBROUTINE iterator(data, npart_it, start)
        USE shared_data
        REAL(num), DIMENSION(:), INTENT(INOUT) :: data
        INTEGER(8), INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
      END SUBROUTINE iterator
    END INTERFACE

    INTEGER :: ispecies

    IF (IAND(dumpmask(id), code) .EQ. 0) RETURN

    CALL species_offset_init()

    CALL start_particle_family_only(current_family)

    DO ispecies = 1, n_species
      IF (current_family%dump .OR. IAND(code, c_io_restartable) .NE. 0) THEN
        CALL cfd_write_nd_particle_variable_with_iterator_all(cfd_handle, &
            TRIM(name), 'Particles_' // TRIM(current_family%name), iterator, &
            particle_species(ispecies)%count, npart_per_it, &
            'Particles_' // TRIM(current_family%name), 'Grid', &
            species_offset(ispecies))
      ENDIF

      CALL advance_particle_family_only(current_family)
    ENDDO

  END SUBROUTINE write_particle_variable

END MODULE diagnostics
