MODULE diagnostics

  USE calc_df
  USE output_cartesian
  USE output_particle
  USE iocontrol
  USE dist_fn
  USE probes
  USE mpi_subtype_control
  !USE iterators
  USE particle_pointer_advance

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_dt, output_routines, iterate_charge

CONTAINS

  SUBROUTINE output_routines(i)   ! i = step index

    INTEGER, INTENT(IN) :: i
    LOGICAL :: print_arrays, last_call
    CHARACTER(LEN=9+data_dir_max_length+n_zeros) :: filename, filename_desc
    CHARACTER(LEN=50) :: temp_name
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: data
    REAL(num), DIMENSION(3) :: stagger = 0.0_num
    INTEGER(KIND=8) :: n_part_per_it = 100000, npart_local, npart_dump_global
    INTEGER :: ispecies, code
    INTEGER, DIMENSION(3) :: dims

    dims = (/nx_global, ny_global, nz_global/)

    CALL io_test(i, print_arrays, last_call)
    ! Allows a maximum of 10^999 output dumps, should be enough for anyone
    ! (feel free to laugh when this isn't the case)
    WRITE(filename_desc, '("(''nfs:'', a, ''/'', i", i3.3, ".", i3.3, &
        &"''.cfd'')")'), n_zeros, n_zeros
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

      npart_local = get_total_local_dumped_particles(&
          IAND(code, c_io_restartable) .NE. 0)
      CALL MPI_ALLREDUCE(npart_local, npart_dump_global, 1, MPI_INTEGER8, &
          MPI_SUM, comm, errcode)
      CALL create_subtypes(IAND(code, c_io_restartable) .NE. 0)
      ALLOCATE(data(-2:nx+3, -2:ny+3, -2:nz+3))
      ! open the file
      ! (filename, rank_of_current_process, MPI_COMMUNICATOR
      ! (can be MPI_COMM_WORLD), MPI_FILE_MODE (passed straight to
      ! MPI_FILE_OPEN))
      CALL cfd_open(filename, rank, comm, MPI_MODE_CREATE + MPI_MODE_WRONLY)
      ! Write the snapshot information
      ! If you prefer the VisIT cycles to display the dump number, change i
      ! for output_file (code_time, n_iterations, rank to write)
      CALL cfd_write_snapshot_data(time, i, 0)

      IF (IAND(dumpmask(1), code) .NE. 0) &
          CALL cfd_write_nd_particle_grid_with_iterator_all("Particles", &
              "Part_Grid", iterate_particles, 3, npart_local, &
              npart_dump_global, npart_per_it, c_particle_cartesian, &
              subtype_particle_var)
      ! Write the cartesian mesh
      ! (Mesh_Name, Mesh_Class, x_array, y_array, rank to write)
      IF (IAND(dumpmask(2), code) .NE. 0) THEN
        IF (.NOT. use_offset_grid) THEN
          CALL cfd_write_3d_cartesian_grid("Grid", "Grid", &
              x_global(1:nx_global), y_global(1:ny_global), &
              z_global(1:nz_global), 0)
        ELSE
          CALL cfd_write_3d_cartesian_grid("Grid", "Grid", &
              x_offset_global(1:nx_global), y_offset_global(1:ny_global), &
              z_offset_global(1:nz_global), 0)
          CALL cfd_write_3d_cartesian_grid("Grid_Full", "Grid", &
              x_global(1:nx_global), y_global(1:ny_global), &
              z_global(1:nz_global), 0)
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

      IF (IAND(dumpmask(9), code) .NE. 0) &
          CALL cfd_write_3d_cartesian_variable_parallel("Ex", &
              "Electric Field", dims, stagger, "Grid", "Grid", &
              ex(1:nx, 1:ny, 1:nz), subtype_field)
      IF (IAND(dumpmask(10), code) .NE. 0) &
          CALL cfd_write_3d_cartesian_variable_parallel("Ey", &
              "Electric Field", dims, stagger, "Grid", "Grid", &
              ey(1:nx, 1:ny, 1:nz), subtype_field)
      IF (IAND(dumpmask(11), code) .NE. 0) &
          CALL cfd_write_3d_cartesian_variable_parallel("Ez", &
              "Electric Field", dims, stagger, "Grid", "Grid", &
              ez(1:nx, 1:ny, 1:nz), subtype_field)

      IF (IAND(dumpmask(12), code) .NE. 0) &
          CALL cfd_write_3d_cartesian_variable_parallel("Bx", &
              "Magnetic Field", dims, stagger, "Grid", "Grid", &
              bx(1:nx, 1:ny, 1:nz), subtype_field)
      IF (IAND(dumpmask(13), code) .NE. 0) &
          CALL cfd_write_3d_cartesian_variable_parallel("By", &
              "Magnetic Field", dims, stagger, "Grid", "Grid", &
              by(1:nx, 1:ny, 1:nz), subtype_field)
      IF (IAND(dumpmask(14), code) .NE. 0) &
          CALL cfd_write_3d_cartesian_variable_parallel("Bz", &
              "Magnetic Field", dims, stagger, "Grid", "Grid", &
              bz(1:nx, 1:ny, 1:nz), subtype_field)

      IF (IAND(dumpmask(15), code) .NE. 0) &
          CALL cfd_write_3d_cartesian_variable_parallel("Jx", &
              "Current", dims, stagger, "Grid", "Grid", &
              jx(1:nx, 1:ny, 1:nz), subtype_field)
      IF (IAND(dumpmask(16), code) .NE. 0) &
          CALL cfd_write_3d_cartesian_variable_parallel("Jy", &
              "Current", dims, stagger, "Grid", "Grid", &
              jy(1:nx, 1:ny, 1:nz), subtype_field)
      IF (IAND(dumpmask(17), code) .NE. 0) &
          CALL cfd_write_3d_cartesian_variable_parallel("Jz", &
              "Current", dims, stagger, "Grid", "Grid", &
              jz(1:nx, 1:ny, 1:nz), subtype_field)

      ! Since these use species lookup tables, have to use the iterator
      ! functions (Variable_Name, Variable_Class, Iterator_Function,
      ! global_npart, npart_per_iteration, Mesh_Name, Mesh_Class,
      ! MPI_TYPE describing data distribution)
      IF (IAND(dumpmask(18), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Q", &
              "Particles", iterate_charge, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(19), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("mass", &
              "Particles", iterate_mass, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)

      IF (IAND(dumpmask(20), code) .NE. 0) THEN
        IF (IAND(dumpmask(20), c_io_no_intrinsic) .EQ. 0) THEN
          CALL calc_ekbar(data, 0)
          CALL cfd_write_3d_cartesian_variable_parallel("EkBar", "EkBar", &
              dims, stagger, "Grid", "Grid", data(1:nx,1:ny,1:nz), &
              subtype_field)
        ENDIF
        IF (IAND(dumpmask(20), c_io_species) .NE. 0) THEN
          DO ispecies = 1, n_species
            CALL calc_ekbar(data, ispecies)
            WRITE(temp_name, '("EkBar_", a)') &
                TRIM(particle_species(ispecies)%name)
            CALL cfd_write_3d_cartesian_variable_parallel(&
                TRIM(ADJUSTL(temp_name)), "EkBar", dims, stagger, "Grid", &
                "Grid", data(1:nx, 1:ny, 1:nz), subtype_field)
          ENDDO
        ENDIF
      ENDIF

      ! These are derived variables from the particles
      ! Since you only dump after several particle updates it's actually
      ! quicker to
      IF (IAND(dumpmask(21), code) .NE. 0) THEN
        CALL calc_mass_density(data, 0)
        CALL cfd_write_3d_cartesian_variable_parallel("Mass_Density", &
            "Derived", dims, stagger, "Grid", "Grid", data(1:nx, 1:ny, 1:nz), &
            subtype_field)
        IF (IAND(dumpmask(21), c_io_species) .NE. 0) THEN
          DO ispecies = 1, n_species
            CALL calc_mass_density(data, ispecies)
            WRITE(temp_name, '("Mass_Density_", a)') &
                TRIM(particle_species(ispecies)%name)
            CALL cfd_write_3d_cartesian_variable_parallel(&
                TRIM(ADJUSTL(temp_name)), "Derived", dims, stagger, "Grid", &
                "Grid", data(1:nx, 1:ny, 1:nz), subtype_field)
          ENDDO
        ENDIF
      ENDIF

      IF (IAND(dumpmask(22), code) .NE. 0) THEN
        CALL calc_charge_density(data, 0)
        CALL cfd_write_3d_cartesian_variable_parallel("Charge_Density", &
            "Derived", dims, stagger, "Grid", "Grid", data(1:nx, 1:ny, 1:nz), &
            subtype_field)
        IF (IAND(dumpmask(22), c_io_species) .NE. 0) THEN
          DO ispecies = 1, n_species
            CALL calc_charge_density(data, ispecies)
            WRITE(temp_name, '("Charge_Density_", a)') &
                TRIM(particle_species(ispecies)%name)
            CALL cfd_write_3d_cartesian_variable_parallel(&
                TRIM(ADJUSTL(temp_name)), "Derived", dims, stagger, "Grid", &
                "Grid", data(1:nx, 1:ny, 1:nz), subtype_field)
          ENDDO
        ENDIF
      ENDIF

      IF (IAND(dumpmask(23), code) .NE. 0) THEN
        CALL calc_number_density(data, 0)
        CALL cfd_write_3d_cartesian_variable_parallel("Number_Density", &
            "Derived", dims, stagger, "Grid", "Grid", data(1:nx, 1:ny, 1:nz), &
            subtype_field)
        IF (IAND(dumpmask(23), c_io_species) .NE. 0) THEN
          DO ispecies = 1, n_species
            CALL calc_number_density(data, ispecies)
            WRITE(temp_name, '("Number_Density_", a)') &
                TRIM(particle_species(ispecies)%name)
            CALL cfd_write_3d_cartesian_variable_parallel(&
                TRIM(ADJUSTL(temp_name)), "Derived", dims, stagger, "Grid", &
                "Grid", data(1:nx, 1:ny, 1:nz), subtype_field)
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
      CALL cfd_write_3d_cartesian_variable_parallel("Rank", "Processor", &
          dims, stagger, "Grid", "Grid", data(1:nx, 1:ny, 1:nz), subtype_field)
#endif

      IF (IAND(dumpmask(26), code) .NE. 0) THEN
        CALL write_dist_fns(code)
      ENDIF

#ifdef PARTICLE_PROBES
      IF (IAND(dumpmask(27), code) .NE. 0) THEN
        CALL write_probes(code)
      ENDIF
#endif

!!$      IF (IAND(dumpmask(28), code) .NE. 0) THEN
!!$        CALL calc_temperature(data, 0)
!!$        CALL cfd_write_2d_cartesian_variable_parallel("Temperature", &
!!$            "Derived", dims, stagger, "Grid", "Grid", data(1:nx, 1:ny), &
!!$            subtype_field)
!!$        IF (IAND(dumpmask(28), c_io_species) .NE. 0) THEN
!!$          DO ispecies = 1, n_species
!!$            CALL calc_temperature(data, ispecies)
!!$            WRITE(temp_name, '("Temperature_", a)') &
!!$                TRIM(particle_species(ispecies)%name)
!!$            CALL cfd_write_2d_cartesian_variable_parallel(&
!!$                TRIM(ADJUSTL(temp_name)), "Derived", dims, stagger, &
!!$                "Grid", "Grid", data(1:nx, 1:ny), subtype_field)
!!$          ENDDO
!!$        ENDIF
!!$      ENDIF

      ! CLOSE the file
      CALL cfd_close()

      output_file = output_file + 1
      IF (rank .EQ. 0) THEN
        WRITE(20, *) "Dumped data at", time, "at iteration", i, &
            "for dump", output_file-1
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

    REAL(num) :: dtx, dty, dtz

    dtx = dx/c
    dty = dy/c
    dtz = dz/c
    dt = MIN(dtx**2, dty**2, dtz**2)/SQRT(dtx**2+dty**2+dtz**2)
    IF (dt_plasma .NE. 0.0_num) dt = MIN(dt, dt_plasma)
    IF (dt_laser .NE. 0.0_num) dt = MIN(dt, dt_laser)
    dt = dt_multiplier * dt

  END SUBROUTINE set_dt



  SUBROUTINE energy_account()

  END SUBROUTINE energy_account



  SUBROUTINE write_species(filehandle, current_displacement)

    INTEGER, INTENT(IN) :: filehandle
    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: current_displacement

!!$    CALL MPI_FILE_SET_VIEW(filehandle, current_displacement, MPI_INTEGER, &
!!$        subtype_particle_int, "native", MPI_INFO_NULL, cfd_errcode)
!!$    CALL MPI_FILE_WRITE_ALL(filehandle, Part_Species, npart, MPI_INTEGER, &
!!$        status, errcode)

  END SUBROUTINE write_species



  ! iterator for particle positions
  SUBROUTINE iterate_particles(data, n_points, direction, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN) :: direction
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF
    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = cur%part_pos(direction)-window_shift(direction)
          ! IF (cur%part_pos(1) .EQ. cur%part_pos(2)) PRINT *, "PATBAD"
          cur=>cur%next
        ENDDO
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO
      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_particles



  ! iterator for particle charge
  SUBROUTINE iterate_charge(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF
    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGEMASS
          data(part_count) = cur%charge
#else
          data(part_count) = current_family%charge
#endif
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO
      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_charge



#ifdef PER_PARTICLE_WEIGHT

  ! iterator for particle weight
  ! Only present if you are using the PER_PARTICLE_WEIGHT
  ! Precompiler option
  SUBROUTINE iterate_weight(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF
    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = cur%weight
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO
      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_weight
#endif



  ! iterator for particle mass
  SUBROUTINE iterate_mass(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF
    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGEMASS
          data(part_count) = cur%mass
#else
          data(part_count) = current_family%mass
#endif
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO
      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_mass



#ifdef PART_DEBUG
  ! iterator for particle processor
  SUBROUTINE iterate_processor(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count
    REAL(num) :: root

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF
    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = REAL(cur%processor, num)
          IF (cur%processor .GE. nproc) PRINT *, "Bad Processor"
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO
      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_processor



  SUBROUTINE iterate_processor0(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count
    REAL(num) :: root

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF
    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = REAL(cur%processor_at_t0, num)
          IF (cur%processor .GE. nproc) PRINT *, "Bad Processor"
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO
      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_processor0
#endif



  ! iterator for particle processor
  SUBROUTINE iterate_species(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF
    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = REAL(current_family%id, num)
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO
      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_species



  ! iterator for particle velocities
  SUBROUTINE iterate_vx(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count
    REAL(num) :: root, part_m

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF
    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGEMASS
          part_m = cur%mass
#else
          part_m = current_family%mass
#endif
          root = SQRT(part_m**2 + &
              (cur%part_p(1)**2 + cur%part_p(2)**2 + cur%part_p(3)**2)/c**2)
          IF (root .NE. 0.0_num) root = 1.0_num/root
          data(part_count) = cur%part_p(1) * root
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO
      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_vx



  SUBROUTINE iterate_vy(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count
    REAL(num) :: root, part_m

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF
    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGEMASS
          part_m = cur%mass
#else
          part_m = current_family%mass
#endif
          root = SQRT(part_m**2 + &
              (cur%part_p(1)**2 + cur%part_p(2)**2 + cur%part_p(3)**2)/c**2)
          IF (root .NE. 0.0_num) root = 1.0_num/root
          data(part_count) = cur%part_p(2) * root
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO
      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_vy



  SUBROUTINE iterate_vz(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count
    REAL(num) :: root, part_m
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF
    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGEMASS
          part_m = cur%mass
#else
          part_m = current_family%mass
#endif
          root = SQRT(part_m**2 + &
              (cur%part_p(1)**2 + cur%part_p(2)**2 + cur%part_p(3)**2)/c**2)
          IF (root .NE. 0.0_num) root = 1.0_num/root
          data(part_count) = cur%part_p(3) * root
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO
      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_vz



  ! iterator for particle momenta
  SUBROUTINE iterate_px(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF
    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = cur%part_p(1)
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO
      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_px



  SUBROUTINE iterate_py(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF
    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = cur%part_p(2)
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO
      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_py



  SUBROUTINE iterate_pz(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF
    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = cur%part_p(3)
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO
      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_pz

END MODULE diagnostics
