MODULE diagnostics

  USE calc_df
  USE output_cartesian
  USE output_particle
  USE iocontrol
  USE dist_fn
  USE probes
  USE mpi_subtype_control
  USE iterators
  !USE particle_pointer_advance

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_dt, output_routines, iterate_charge

CONTAINS

  SUBROUTINE output_routines(i)   ! i = step index

    INTEGER, INTENT(IN) :: i
    LOGICAL :: print_arrays, last_call
    CHARACTER(LEN=9+data_dir_max_length+n_zeros) :: filename, filename_desc
    CHARACTER(LEN=50) :: temp_name
    REAL(num), DIMENSION(:), ALLOCATABLE :: data
    REAL(num) :: stagger = 0.0_num
    INTEGER(KIND=8) :: n_part_per_it = 100000, npart_local, npart_dump_global
    INTEGER :: ispecies, code
    INTEGER :: dims

    dims = nx_global

    CALL io_test(i, print_arrays, last_call)
    ! Allows a maximum of 10^999 output dumps, should be enough for anyone
    ! (feel free to laugh when this isn't the case)
    WRITE(filename_desc, '("(''nfs:'', a, ''/'', i", i3.3, ".", i3.3, &
        &", ''.cfd'')")'), n_zeros, n_zeros
    WRITE(filename, filename_desc) TRIM(data_dir), output_file

    IF (print_arrays) THEN
      ! Always dump the variables with the "Every" attribute
      code = c_io_always

      ! Only dump variables with the "FULL" attributre on full dump intervals
      IF (MOD(output_file, full_dump_every) .EQ. 0)  code = IOR(code, c_io_full)
      IF (MOD(output_file, restart_dump_every) .EQ. 0 .AND. &
          restart_dump_every .GT. -1) code = IOR(code, c_io_restartable)
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
      ELSE
        iterator_settings%restart = .FALSE.
      ENDIF

      ALLOCATE(data(-2:nx+3))

      ! open the file
      ! (filename, rank_of_current_process, MPI_COMMUNICATOR (can be
      ! MPI_COMM_WORLD), MPI_FILE_MODE (passed straight to MPI_FILE_OPEN))
      CALL cfd_open(filename, rank, comm, MPI_MODE_CREATE + MPI_MODE_WRONLY)
      ! Write the snapshot information
      ! If you prefer the VisIt cycles to display the dump number, change i
      ! for output_file (code_time, n_iterations, rank to write)
      CALL cfd_write_snapshot_data(time, i, 0)

      IF (IAND(dumpmask(1), code) .NE. 0) &
          CALL cfd_write_nd_particle_grid_with_iterator_all("Particles", &
              "Part_Grid", iterate_particles, 1, npart_local, &
              npart_dump_global, npart_per_it, c_particle_cartesian, &
              subtype_particle_var)

      ! Write the cartesian mesh
      ! (Mesh_Name, Mesh_Class, x_array, y_array, rank to write)
      IF (IAND(dumpmask(2), code) .NE. 0) THEN
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

      ! (Variable_Name, Variable_Class, array, global_npart, Mesh_Name,
      ! Mesh_Class, MPI_TYPE describing data distribution)
      IF (IAND(dumpmask(3), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Px", &
              "Particles", iterate_px, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(4), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Py", &
              "Particles", iterate_py, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(5), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Pz", &
              "Particles", iterate_pz, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(6), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Vx", &
              "Particles", iterate_vx, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(7), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Vy", &
              "Particles", iterate_vy, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(8), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Vz", &
              "Particles", iterate_vz, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
#ifdef PART_DEBUG
      CALL cfd_write_nd_particle_variable_with_iterator_all("Processor", &
          "Particles", iterate_processor, npart_dump_global, n_part_per_it, &
          "Particles", "Part_Grid", subtype_particle_var)
      CALL cfd_write_nd_particle_variable_with_iterator_all("Processor_at_t0", &
          "Particles", iterate_processor0, npart_dump_global, n_part_per_it, &
          "Particles", "Part_Grid", subtype_particle_var)
#endif

      ! Serial Cartesian write because shared_grid
      ! parallel equivalent is cfd_Write_1D_Cartesian_Variable_All
      ! (Variable_Name, Variable_Class, Grid_Stagger, Mesh_Name, Mesh_Class,
      ! array, rank of process to write data)
      ! Note that serial writes must be called on each node to keep file
      ! pointer updated
      IF (IAND(dumpmask(9), code) .NE. 0) &
          CALL cfd_write_1d_cartesian_variable_parallel("Ex", &
              "Electric Field", dims, stagger, "Grid", "Grid", &
              ex(1:nx), subtype_field)
      IF (IAND(dumpmask(10), code) .NE. 0) &
          CALL cfd_write_1d_cartesian_variable_parallel("Ey", &
              "Electric Field", dims, stagger, "Grid", "Grid", &
              ey(1:nx), subtype_field)
      IF (IAND(dumpmask(11), code) .NE. 0) &
          CALL cfd_write_1d_cartesian_variable_parallel("Ez", &
              "Electric Field", dims, stagger, "Grid", "Grid", &
              ez(1:nx), subtype_field)

      IF (IAND(dumpmask(12), code) .NE. 0) &
          CALL cfd_write_1d_cartesian_variable_parallel("Bx", &
              "Magnetic Field", dims, stagger, "Grid", "Grid", &
              bx(1:nx), subtype_field)
      IF (IAND(dumpmask(13), code) .NE. 0) &
          CALL cfd_write_1d_cartesian_variable_parallel("By", &
              "Magnetic Field", dims, stagger, "Grid", "Grid", &
              by(1:nx), subtype_field)
      IF (IAND(dumpmask(14), code) .NE. 0) &
          CALL cfd_write_1d_cartesian_variable_parallel("Bz", &
              "Magnetic Field", dims, stagger, "Grid", "Grid", &
              bz(1:nx), subtype_field)

      IF (IAND(dumpmask(15), code) .NE. 0) &
          CALL cfd_write_1d_cartesian_variable_parallel("Jx", &
              "Current", dims, stagger, "Grid", "Grid", &
              jx(1:nx), subtype_field)
      IF (IAND(dumpmask(16), code) .NE. 0) &
          CALL cfd_write_1d_cartesian_variable_parallel("Jy", &
              "Current", dims, stagger, "Grid", "Grid", &
              jy(1:nx), subtype_field)
      IF (IAND(dumpmask(17), code) .NE. 0) &
          CALL cfd_write_1d_cartesian_variable_parallel("Jz", &
              "Current", dims, stagger, "Grid", "Grid", &
              jz(1:nx), subtype_field)

      ! Since these use species lookup tables, have to use the iterator
      ! functions (Variable_Name, Variable_Class, Iterator_Function,
      ! global_npart, npart_per_iteration, Mesh_Name, Mesh_Class,
      ! MPI_TYPE describing data distribution)
      IF (IAND(dumpmask(18), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Q", &
              "Particles", iterate_charge, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(19), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Mass", &
              "Particles", iterate_mass, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)

      IF (IAND(dumpmask(20), code) .NE. 0) THEN
        CALL calc_ekbar(data, 0)
        IF (IAND(dumpmask(20), c_io_no_intrinsic) .EQ. 0) &
            CALL cfd_write_1d_cartesian_variable_parallel("EkBar", "EkBar", &
                dims, stagger, "Grid", "Grid", data(1:nx), subtype_field)
        IF (IAND(dumpmask(20), c_io_species) .NE. 0) THEN
          DO ispecies = 1, n_species
            CALL calc_ekbar(data, ispecies)
            WRITE(temp_name, '("EkBar_", a)') &
                TRIM(particle_species(ispecies)%name)
            CALL cfd_write_1d_cartesian_variable_parallel( &
                TRIM(ADJUSTL(temp_name)), "EkBar", dims, stagger, &
                "Grid", "Grid", data(1:nx), subtype_field)
          ENDDO
        ENDIF
      ENDIF

      ! These are derived variables from the particles
      ! Since you only dump after several particle updates it's actually
      ! quicker to
      IF (IAND(dumpmask(21), code) .NE. 0) THEN
        CALL calc_mass_density(data, 0)
        CALL cfd_write_1d_cartesian_variable_parallel("Mass_Density", &
            "Derived", dims, stagger, "Grid", "Grid", data(1:nx), subtype_field)
        IF (IAND(dumpmask(21), c_io_species) .NE. 0) THEN
          DO ispecies = 1, n_species
            CALL calc_mass_density(data, ispecies)
            WRITE(temp_name, '("Mass_Density_", a)') &
                TRIM(particle_species(ispecies)%name)
            CALL cfd_write_1d_cartesian_variable_parallel( &
                TRIM(ADJUSTL(temp_name)), "Derived", dims, stagger, &
                "Grid", "Grid", data(1:nx), subtype_field)
          ENDDO
        ENDIF
      ENDIF

      IF (IAND(dumpmask(22), code) .NE. 0) THEN
        CALL calc_charge_density(data, 0)
        CALL cfd_write_1d_cartesian_variable_parallel("Charge_Density", &
            "Derived", dims, stagger, "Grid", "Grid", data(1:nx), subtype_field)
        IF (IAND(dumpmask(22), c_io_species) .NE. 0) THEN
          DO ispecies = 1, n_species
            CALL calc_charge_density(data, ispecies)
            WRITE(temp_name, '("Charge_Density_", a)') &
                TRIM(particle_species(ispecies)%name)
            CALL cfd_write_1d_cartesian_variable_parallel( &
                TRIM(ADJUSTL(temp_name)), "Derived", dims, stagger, &
                "Grid", "Grid", data(1:nx), subtype_field)
          ENDDO
        ENDIF
      ENDIF

      IF (IAND(dumpmask(23), code) .NE. 0) THEN
        CALL calc_number_density(data, 0)
        CALL cfd_write_1d_cartesian_variable_parallel("Number_Density", &
            "Derived", dims, stagger, "Grid", "Grid", data(1:nx), subtype_field)
        IF (IAND(dumpmask(23), c_io_species) .NE. 0) THEN
          DO ispecies = 1, n_species
            CALL calc_number_density(data, ispecies)
            WRITE(temp_name, '("Number_Density_", a)') &
                TRIM(particle_species(ispecies)%name)
            CALL cfd_write_1d_cartesian_variable_parallel( &
                TRIM(ADJUSTL(temp_name)), "Derived", dims, stagger, &
                "Grid", "Grid", data(1:nx), subtype_field)
          ENDDO
        ENDIF
      ENDIF

      IF (IAND(dumpmask(24), code) .NE. 0) THEN
#ifdef PER_PARTICLE_WEIGHT
        CALL cfd_write_nd_particle_variable_with_iterator_all("Weight", &
            "Particles", iterate_weight, npart_dump_global, n_part_per_it, &
            "Particles", "Part_Grid", subtype_particle_var)
#else
        CALL cfd_write_real_constant("Weight", "Particles", weight, 0)
#endif
      ENDIF

      IF (IAND(dumpmask(25), code) .NE. 0) THEN
        CALL cfd_write_nd_particle_variable_with_iterator_all("Species", &
            "Particles", iterate_species, npart_dump_global, n_part_per_it, &
            "Particles", "Part_Grid", subtype_particle_var)
      ENDIF

#ifdef FIELD_DEBUG
      data = rank
      CALL cfd_write_1d_cartesian_variable_parallel("Rank", "Processor", &
          dims, stagger, "Grid", "Grid", data(1:nx), subtype_field)
#endif

      IF (IAND(dumpmask(26), code) .NE. 0) THEN
        CALL write_dist_fns(code)
      ENDIF

#ifdef PARTICLE_PROBES
      IF (IAND(dumpmask(27), code) .NE. 0) THEN
        CALL write_probes(code)
      ENDIF
#endif

      IF (IAND(dumpmask(28), code) .NE. 0) THEN
        CALL calc_temperature(data, 0)
        CALL cfd_write_1d_cartesian_variable_parallel("Temperature", &
            "Derived", dims, stagger, "Grid", "Grid", data(1:nx), subtype_field)
        IF (IAND(dumpmask(28), c_io_species) .NE. 0) THEN
          DO ispecies = 1, n_species
            CALL calc_temperature(data, ispecies)
            WRITE(temp_name, '("Temperature_", a)') &
                TRIM(particle_species(ispecies)%name)
            CALL cfd_write_1d_cartesian_variable_parallel( &
                TRIM(ADJUSTL(temp_name)), "Derived", dims, stagger, &
                "Grid", "Grid", data(1:nx), subtype_field)
          ENDDO
        ENDIF
      ENDIF

      ! CLOSE the file
      CALL cfd_close()

      output_file = output_file + 1
      IF (rank .EQ. 0) THEN
        WRITE(20, *) "Dumped data at", time, "at iteration", i, "for dump", &
            output_file-1
        CALL FLUSH(20)
      ENDIF

      DEALLOCATE(data)
    ENDIF

  END SUBROUTINE output_routines



  SUBROUTINE io_test(i, print_arrays, last_call)

    INTEGER, INTENT(IN) :: i
    LOGICAL, INTENT(OUT) :: print_arrays, last_call

    REAL(num), SAVE :: t1 = 0.0_num

    IF (restart) THEN
      t1 = time
      restart = .FALSE.
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

!!$ dtx = dx/(6.0_num * c)
!!$ dty = dy/(6.0_num * c)
    dtx = dx/c
    dt = dtx
    ! dt = MIN(dtx**2, dty**2)/SQRT(dtx**2+dty**2)
    IF (dt_laser .NE. 0.0_num) dt = MIN(dt, dt_laser)
    IF (dt_plasma_frequency .NE. 0.0_num) dt = MIN(dt, dt_plasma_frequency)
#ifdef NEWTONIAN
    dtx = dx/max_part_v
    dt = MIN(dt, dtx)
#endif
    dt = dt_multiplier * dt

  END SUBROUTINE set_dt



  SUBROUTINE energy_account()

  END SUBROUTINE energy_account

END MODULE diagnostics
