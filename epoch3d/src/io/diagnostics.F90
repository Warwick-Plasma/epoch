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
    CHARACTER(LEN=8) :: dump_type
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: data
    REAL(num), DIMENSION(c_ndims) :: stagger = 0.0_num
    INTEGER(8) :: n_part_per_it = 100000, npart_local, npart_dump_global
    INTEGER :: ispecies, code
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER :: restart_flag

    dims = (/nx_global, ny_global, nz_global/)

    CALL io_test(i, print_arrays, last_call)
    ! Allows a maximum of 10^999 output dumps, should be enough for anyone
    ! (feel free to laugh when this isn't the case)
    WRITE(filename_desc, '("(a, ''/'', i", i3.3, ".", i3.3, ", ''.cfd'')")') &
        n_zeros, n_zeros
    WRITE(filename, filename_desc) TRIM(data_dir), output_file

    IF (print_arrays) THEN
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

      ! is a restart dump
      IF (IAND(code, c_io_restartable) .NE. 0) THEN
        restart_flag = 1
      ELSE
        restart_flag = 0
      ENDIF

      ALLOCATE(data(-2:nx+3, -2:ny+3, -2:nz+3))

      ! open the file
      ! (filename, rank of current process, MPI communicator (can be
      ! MPI_COMM_WORLD), file mode (c_cfd_read or c_cfd_write),
      ! cycle number, simulation time, job id)
      CALL cfd_open(filename, rank, comm, c_cfd_write, i, time, jobid)
      CALL cfd_write_job_info(restart_flag, sha1sum, 0)

      ! Write the snapshot information
      ! If you prefer the VisIt cycles to display the dump number, change i
      ! for output_file (code_time, n_iterations, rank used for writing)
      CALL cfd_write_snapshot_data(time, i, 0)

      IF (IAND(dumpmask(c_dump_part_grid), code) .NE. 0) &
          CALL cfd_write_nd_particle_grid_with_iterator_all("Particles", &
              "Part_Grid", iterate_particles, c_dimension_3d, npart_local, &
              npart_dump_global, npart_per_it, c_particle_cartesian, &
              subtype_particle_var)

      ! Write the cartesian mesh
      ! (mesh name, mesh class, x_array, y_array, rank used for writing)
      IF (IAND(dumpmask(c_dump_grid), code) .NE. 0) THEN
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

      ! (variable name, variable class, iterator function,
      ! global number of particles, number of particles to write per iteration,
      ! mesh name, mesh class, mpi type describing data distribution)
      IF (IAND(dumpmask(c_dump_part_species), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Species", &
              "Particles", iterate_species, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(c_dump_part_weight), code) .NE. 0) &
#ifdef PER_PARTICLE_WEIGHT
          CALL cfd_write_nd_particle_variable_with_iterator_all("Weight", &
              "Particles", iterate_weight, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
#else
          CALL cfd_write_real_constant("Weight", "Particles", weight, 0)
#endif
      IF (IAND(dumpmask(c_dump_part_px), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Px", &
              "Particles", iterate_px, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(c_dump_part_py), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Py", &
              "Particles", iterate_py, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(c_dump_part_pz), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Pz", &
              "Particles", iterate_pz, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(c_dump_part_vx), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Vx", &
              "Particles", iterate_vx, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(c_dump_part_vy), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Vy", &
              "Particles", iterate_vy, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(c_dump_part_vz), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Vz", &
              "Particles", iterate_vz, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(c_dump_part_charge), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Q", &
              "Particles", iterate_charge, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
      IF (IAND(dumpmask(c_dump_part_mass), code) .NE. 0) &
          CALL cfd_write_nd_particle_variable_with_iterator_all("Mass", &
              "Particles", iterate_mass, npart_dump_global, n_part_per_it, &
              "Particles", "Part_Grid", subtype_particle_var)
#ifdef PARTICLE_DEBUG
      CALL cfd_write_nd_particle_variable_with_iterator_all("Processor", &
          "Particles", iterate_processor, npart_dump_global, n_part_per_it, &
          "Particles", "Part_Grid", subtype_particle_var)
      CALL cfd_write_nd_particle_variable_with_iterator_all("Processor_at_t0", &
          "Particles", iterate_processor0, npart_dump_global, n_part_per_it, &
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
      CALL cfd_write_3d_cartesian_variable_parallel("Rank", "Processor", &
          dims, stagger, "Grid", "Grid", data(1:nx, 1:ny, 1:nz), subtype_field)
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
    ENDIF

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

    REAL(num) :: dtx, dty, dtz

    dtx = dx/c
    dty = dy/c
    dtz = dz/c
    dt = MIN(dtx**2, dty**2, dtz**2)/SQRT(dtx**2+dty**2+dtz**2)
    IF (dt_plasma_frequency .NE. 0.0_num) dt = MIN(dt, dt_plasma_frequency)
    IF (dt_laser .NE. 0.0_num) dt = MIN(dt, dt_laser)
    dt = dt_multiplier * dt

  END SUBROUTINE set_dt



  SUBROUTINE energy_account()

  END SUBROUTINE energy_account



  SUBROUTINE write_species(filehandle, current_displacement)

    INTEGER, INTENT(IN) :: filehandle
    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: current_displacement

!!$    CALL MPI_FILE_SET_VIEW(filehandle, current_displacement, MPI_INTEGER4, &
!!$        subtype_particle_int, "native", MPI_INFO_NULL, cfd_errcode)
!!$    CALL MPI_FILE_WRITE_ALL(filehandle, Part_Species, npart, MPI_INTEGER4, &
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



#ifdef PARTICLE_DEBUG
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
          root = SQRT(part_m**2 &
              + (cur%part_p(1)**2 + cur%part_p(2)**2 + cur%part_p(3)**2)/c**2)
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
          root = SQRT(part_m**2 &
              + (cur%part_p(1)**2 + cur%part_p(2)**2 + cur%part_p(3)**2)/c**2)
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
          root = SQRT(part_m**2 &
              + (cur%part_p(1)**2 + cur%part_p(2)**2 + cur%part_p(3)**2)/c**2)
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
 
 
 
  SUBROUTINE write_field(id, code, name, class, array)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), DIMENSION(:,:,:), INTENT(IN) :: array
    REAL(num), DIMENSION(c_ndims) :: stagger = 0.0_num
    INTEGER, DIMENSION(c_ndims) :: dims

    IF (IAND(dumpmask(id), code) .EQ. 0) RETURN

    dims = (/nx_global, ny_global, nz_global/)

    ! (variable name, variable class, global grid dimensions,
    ! grid stagger, mesh name, mesh class, variable,
    ! mpi type describing data distribution)
    CALL cfd_write_3d_cartesian_variable_parallel( &
        TRIM(name), TRIM(class), dims, stagger, &
        'Grid', 'Grid', array(1:nx,1:ny,1:nz), subtype_field)

  END SUBROUTINE write_field



  SUBROUTINE write_nspecies_field(id, code, name, class, func, data)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), DIMENSION(:,:,:), INTENT(INOUT) :: data
    REAL(num), DIMENSION(c_ndims) :: stagger = 0.0_num
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER :: ispecies
    CHARACTER(LEN=50) :: temp_name

    INTERFACE
      SUBROUTINE func(data_array, current_species)
        USE shared_data
        REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
        INTEGER, INTENT(IN) :: current_species
      END SUBROUTINE func
    END INTERFACE

    IF (IAND(dumpmask(id), code) .EQ. 0) RETURN

    dims = (/nx_global, ny_global, nz_global/)

    IF (IAND(dumpmask(id), c_io_no_intrinsic) .EQ. 0) THEN
      CALL func(data, 0)
      CALL cfd_write_3d_cartesian_variable_parallel( &
          TRIM(name), TRIM(class), dims, stagger, &
          'Grid', 'Grid', data(1:nx, 1:ny, 1:nz), subtype_field)
    ENDIF

    IF (IAND(dumpmask(id), c_io_species) .NE. 0) THEN
      DO ispecies = 1, n_species
        CALL func(data, ispecies)
        WRITE(temp_name, '(a, "_", a)') TRIM(name), &
            TRIM(particle_species(ispecies)%name)
        CALL cfd_write_3d_cartesian_variable_parallel( &
            TRIM(ADJUSTL(temp_name)), TRIM(class), dims, stagger, &
            'Grid', 'Grid', data(1:nx, 1:ny, 1:nz), subtype_field)
      ENDDO
    ENDIF

  END SUBROUTINE write_nspecies_field

END MODULE diagnostics
