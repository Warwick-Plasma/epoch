MODULE diagnostics

  USE calc_df
  USE output_cartesian
  USE output_particle
  USE iocontrol
  USE dist_fn
  USE probes
  USE mpi_subtype_control
  USE encoded_source
  USE deck
  USE iterators

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_dt, output_routines, iterate_charge

CONTAINS

  SUBROUTINE output_routines(i)   ! i = step index

    INTEGER, INTENT(INOUT) :: i
    LOGICAL :: print_arrays, last_call
    CHARACTER(LEN=9+data_dir_max_length+n_zeros) :: filename, filename_desc
    CHARACTER(LEN=50) :: temp_name
    CHARACTER(LEN=8) :: dump_type
    REAL(num), DIMENSION(:), ALLOCATABLE :: data
    REAL(num) :: stagger = 0.0_num
    INTEGER(8) :: npart_local, npart_dump_global
    INTEGER :: ispecies, code
    INTEGER :: dims
    INTEGER :: restart_flag

    IF (rank .EQ. 0 .AND. stdout_frequency .GT. 0 &
        .AND. MOD(i, stdout_frequency) .EQ. 0) THEN
      WRITE(*, '("Time", g20.12, " and iteration", i7, " after ", &
          & f8.1, " seconds")') time, i, MPI_WTIME() - walltime_start
    ENDIF

    CALL io_test(i, print_arrays, last_call)

    IF (.NOT.print_arrays) RETURN

    dims = nx_global

    ! Allows a maximum of 10^999 output dumps, should be enough for anyone
    ! (feel free to laugh when this isn't the case)
    WRITE(filename_desc, '("(a, ''/'', i", i3.3, ".", i3.3, ", ''.cfd'')")') &
        n_zeros, n_zeros
    WRITE(filename, filename_desc) TRIM(data_dir), output_file

    ! Always dump the variables with the "Every" attribute
    code = c_io_always

    ! Only dump variables with the "FULL" attributre on full dump intervals
    IF (MOD(output_file, full_dump_every) .EQ. 0)  code = IOR(code, c_io_full)
    IF (MOD(output_file, restart_dump_every) .EQ. 0 &
        .AND. restart_dump_every .GT. -1) code = IOR(code, c_io_restartable)
    IF (last_call .AND. force_final_to_be_restartable) &
        code = IOR(code, c_io_restartable)

    npart_local = &
        get_total_local_dumped_particles(IAND(code, c_io_restartable) .NE. 0)

    CALL MPI_ALLREDUCE(npart_local, npart_dump_global, 1, MPI_INTEGER8, &
        MPI_SUM, comm, errcode)
    CALL create_subtypes(IAND(code, c_io_restartable) .NE. 0)

    ! If the code is doing a restart dump then tell the iterators that this
    ! is a restart dump
    IF (IAND(code, c_io_restartable) .NE. 0) THEN
      iterator_settings%restart = .TRUE.
      restart_flag = 1
    ELSE
      iterator_settings%restart = .FALSE.
      restart_flag = 0
    ENDIF

    ALLOCATE(data(-2:nx+3))

    ! open the file
    ! (filename, rank of current process, MPI communicator (can be
    ! MPI_COMM_WORLD), file mode (c_cfd_read or c_cfd_write),
    ! cycle number, simulation time, job id)
    CALL cfd_open(filename, rank, comm, c_cfd_write, i, time, jobid)
    CALL cfd_write_job_info(restart_flag, sha1sum, 0)

    IF (IAND(dumpmask(c_dump_part_grid), code) .NE. 0) &
        CALL cfd_write_nd_particle_grid_with_iterator_all("Particles", &
            "Part_Grid", iterate_particles, c_dimension_1d, npart_local, &
            npart_dump_global, npart_per_it, c_particle_cartesian, &
            subtype_particle_var)

    ! Write the cartesian mesh
    ! (mesh name, mesh class, x_array, y_array, rank used for writing)
    IF (IAND(dumpmask(c_dump_grid), code) .NE. 0) THEN
      IF (.NOT. use_offset_grid) THEN
        CALL cfd_write_1d_cartesian_grid("Grid", "Grid", &
            x_global(1:nx_global), 0)
      ELSE
        CALL cfd_write_1d_cartesian_grid("Grid", "Grid", &
            x_offset_global(1:nx_global), 0)
        CALL cfd_write_1d_cartesian_grid("Grid_Full", "Grid", &
            x_global(1:nx_global), 0)
      ENDIF
    ENDIF

    ! (variable name, variable class, iterator function,
    ! global number of particles, number of particles to write per iteration,
    ! mesh name, mesh class, mpi type describing data distribution)
    IF (IAND(dumpmask(c_dump_part_species), code) .NE. 0) &
        CALL cfd_write_nd_particle_variable_with_iterator_all("Species", &
            "Particles", iterate_species, npart_dump_global, npart_per_it, &
            "Particles", "Part_Grid", subtype_particle_var)
    IF (IAND(dumpmask(c_dump_part_weight), code) .NE. 0) &
#ifdef PER_PARTICLE_WEIGHT
        CALL cfd_write_nd_particle_variable_with_iterator_all("Weight", &
            "Particles", iterate_weight, npart_dump_global, npart_per_it, &
            "Particles", "Part_Grid", subtype_particle_var)
#else
        CALL cfd_write_real_constant("Weight", "Particles", weight, 0)
#endif
    IF (IAND(dumpmask(c_dump_part_px), code) .NE. 0) &
        CALL cfd_write_nd_particle_variable_with_iterator_all("Px", &
            "Particles", iterate_px, npart_dump_global, npart_per_it, &
            "Particles", "Part_Grid", subtype_particle_var)
    IF (IAND(dumpmask(c_dump_part_py), code) .NE. 0) &
        CALL cfd_write_nd_particle_variable_with_iterator_all("Py", &
            "Particles", iterate_py, npart_dump_global, npart_per_it, &
            "Particles", "Part_Grid", subtype_particle_var)
    IF (IAND(dumpmask(c_dump_part_pz), code) .NE. 0) &
        CALL cfd_write_nd_particle_variable_with_iterator_all("Pz", &
            "Particles", iterate_pz, npart_dump_global, npart_per_it, &
            "Particles", "Part_Grid", subtype_particle_var)
    IF (IAND(dumpmask(c_dump_part_vx), code) .NE. 0) &
        CALL cfd_write_nd_particle_variable_with_iterator_all("Vx", &
            "Particles", iterate_vx, npart_dump_global, npart_per_it, &
            "Particles", "Part_Grid", subtype_particle_var)
    IF (IAND(dumpmask(c_dump_part_vy), code) .NE. 0) &
        CALL cfd_write_nd_particle_variable_with_iterator_all("Vy", &
            "Particles", iterate_vy, npart_dump_global, npart_per_it, &
            "Particles", "Part_Grid", subtype_particle_var)
    IF (IAND(dumpmask(c_dump_part_vz), code) .NE. 0) &
        CALL cfd_write_nd_particle_variable_with_iterator_all("Vz", &
            "Particles", iterate_vz, npart_dump_global, npart_per_it, &
            "Particles", "Part_Grid", subtype_particle_var)
    IF (IAND(dumpmask(c_dump_part_charge), code) .NE. 0) &
        CALL cfd_write_nd_particle_variable_with_iterator_all("Q", &
            "Particles", iterate_charge, npart_dump_global, npart_per_it, &
            "Particles", "Part_Grid", subtype_particle_var)
    IF (IAND(dumpmask(c_dump_part_mass), code) .NE. 0) &
        CALL cfd_write_nd_particle_variable_with_iterator_all("Mass", &
            "Particles", iterate_mass, npart_dump_global, npart_per_it, &
            "Particles", "Part_Grid", subtype_particle_var)
#ifdef PARTICLE_DEBUG
    CALL cfd_write_nd_particle_variable_with_iterator_all("Processor", &
        "Particles", iterate_processor, npart_dump_global, npart_per_it, &
        "Particles", "Part_Grid", subtype_particle_var)
    CALL cfd_write_nd_particle_variable_with_iterator_all("Processor_at_t0", &
        "Particles", iterate_processor0, npart_dump_global, npart_per_it, &
        "Particles", "Part_Grid", subtype_particle_var)
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
    CALL write_nspecies_field(c_dump_ekbar, code, 'EkBar', &
        'EkBar', calc_ekbar, data)

    CALL write_nspecies_field(c_dump_mass_density, code, 'Mass_density', &
        'Derived', calc_mass_density, data)

    CALL write_nspecies_field(c_dump_charge_density, code, 'Charge_density', &
        'Derived', calc_charge_density, data)

    CALL write_nspecies_field(c_dump_number_density, code, 'Number_density', &
        'Derived', calc_number_density, data)

    CALL write_nspecies_field(c_dump_temperature, code, 'Temperature', &
        'Derived', calc_temperature, data)

#ifdef FIELD_DEBUG
    data = rank
    CALL cfd_write_1d_cartesian_variable_parallel("Rank", "Processor", &
        dims, stagger, "Grid", "Grid", data, subtype_field, subarray_field)
#endif

    IF (IAND(dumpmask(c_dump_dist_fns), code) .NE. 0) THEN
      CALL write_dist_fns(code)
    ENDIF

#ifdef PARTICLE_PROBES
    IF (IAND(dumpmask(c_dump_probes), code) .NE. 0) THEN
      CALL write_probes(code)
    ENDIF
#endif

    IF (restart_flag .EQ. 1 .AND. LEN(source_code) .GT. 0) THEN
      CALL write_input_decks
      CALL cfd_write_source_code("Code", "base64_packed_source_code", &
          source_code, last_line, 0)
    ENDIF

    ! close the file
    CALL cfd_close()

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

    DEALLOCATE(data)

  END SUBROUTINE output_routines



  SUBROUTINE io_test(i, print_arrays, last_call)

    INTEGER, INTENT(IN) :: i
    LOGICAL, INTENT(OUT) :: print_arrays, last_call

    REAL(num), SAVE :: t1 = 0.0_num
    LOGICAL, SAVE :: first = .TRUE.

    IF (first) THEN
      IF (ic_from_restart) t1 = time
      first = .FALSE.
    ENDIF

    print_arrays = .FALSE.
    last_call = .FALSE.

    IF (time .GE. t1) THEN
      print_arrays = .TRUE.
      t1 = t1 + dt_snapshots
    ENDIF

    IF (time .GE. t_end .OR. i .EQ. nsteps) THEN
      last_call = .TRUE.
      print_arrays = .TRUE.
    ENDIF

  END SUBROUTINE io_test



  SUBROUTINE set_dt        ! sets CFL limited step

    REAL(num) :: dtx

    dtx = dx/c
    dt = dtx
    ! dt = MIN(dtx**2, dty**2)/SQRT(dtx**2+dty**2)
    IF (dt_laser .NE. 0.0_num) dt = MIN(dt, dt_laser)
    IF (dt_plasma_frequency .NE. 0.0_num) dt = MIN(dt, dt_plasma_frequency)
    dt = dt_multiplier * dt

  END SUBROUTINE set_dt



  SUBROUTINE energy_account()

  END SUBROUTINE energy_account



  SUBROUTINE write_field(id, code, name, class, array)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    REAL(num) :: stagger = 0.0_num
    INTEGER :: dims

    IF (IAND(dumpmask(id), code) .EQ. 0) RETURN

    dims = nx_global

    ! (variable name, variable class, global grid dimensions,
    ! grid stagger, mesh name, mesh class, variable,
    ! mpi type describing data distribution)
    CALL cfd_write_1d_cartesian_variable_parallel( &
        TRIM(name), TRIM(class), dims, stagger, &
        'Grid', 'Grid', array, subtype_field, subarray_field)

  END SUBROUTINE write_field



  SUBROUTINE write_nspecies_field(id, code, name, class, func, data)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    REAL(num) :: stagger = 0.0_num
    INTEGER :: dims
    INTEGER :: ispecies
    CHARACTER(LEN=50) :: temp_name

    INTERFACE
      SUBROUTINE func(data_array, current_species)
        USE shared_data
        REAL(num), DIMENSION(-2:), INTENT(INOUT) :: data_array
        INTEGER, INTENT(IN) :: current_species
      END SUBROUTINE func
    END INTERFACE

    IF (IAND(dumpmask(id), code) .EQ. 0) RETURN

    dims = nx_global

    IF (IAND(dumpmask(id), c_io_no_intrinsic) .EQ. 0) THEN
      CALL func(data, 0)
      CALL cfd_write_1d_cartesian_variable_parallel( &
          TRIM(name), TRIM(class), dims, stagger, &
          'Grid', 'Grid', data, subtype_field, subarray_field)
    ENDIF

    IF (IAND(dumpmask(id), c_io_species) .NE. 0) THEN
      DO ispecies = 1, n_species
        CALL func(data, ispecies)
        WRITE(temp_name, '(a, "_", a)') TRIM(name), &
            TRIM(particle_species(ispecies)%name)
        CALL cfd_write_1d_cartesian_variable_parallel( &
            TRIM(ADJUSTL(temp_name)), TRIM(class), dims, stagger, &
            'Grid', 'Grid', data, subtype_field, subarray_field)
      ENDDO
    ENDIF

  END SUBROUTINE write_nspecies_field

END MODULE diagnostics
