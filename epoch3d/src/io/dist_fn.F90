MODULE dist_fn

  USE mpi
  USE sdf
  USE mpi_subtype_control
  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE attach_dist_fn(block)

    TYPE(distribution_function_block), POINTER :: block
    TYPE(distribution_function_block), POINTER :: current

    current=>dist_fns
    IF (.NOT. ASSOCIATED(current)) THEN
      ! This is the first distribution function to add
      dist_fns=>block
      RETURN
    ENDIF
    DO WHILE(ASSOCIATED(current%next))
      current=>current%next
    ENDDO
    current%next=>block

  END SUBROUTINE attach_dist_fn



  SUBROUTINE init_dist_fn(block)

    TYPE(distribution_function_block), POINTER :: block

    block%name = blank
    block%ndims = -1
    block%dumpmask = c_io_always
    block%directions = 0
    block%ranges = 1.0_num
    block%resolution = 1
    block%restrictions = 0.0_num
    block%use_restrictions = .FALSE.
    NULLIFY(block%next)
    ALLOCATE(block%use_species(n_species))
    block%use_species = .FALSE.

  END SUBROUTINE init_dist_fn



  SUBROUTINE clean_dist_fns()

    TYPE(distribution_function_block), POINTER :: current, next

    current=>dist_fns
    DO WHILE(ASSOCIATED(current))
      next=>current%next
      DEALLOCATE(current)
      current=>next
    ENDDO

  END SUBROUTINE clean_dist_fns



  SUBROUTINE write_dist_fns(sdf_handle, code)

    TYPE(sdf_file_handle) :: sdf_handle
    INTEGER, INTENT(IN) :: code

    INTEGER :: ispecies, errcode
    TYPE(distribution_function_block), POINTER :: current

    ! Write the distribution functions
    current=>dist_fns
    DO WHILE(ASSOCIATED(current))
      IF (IAND(current%dumpmask, code) .NE. 0) THEN
        DO ispecies = 1, n_species
          IF (.NOT. current%use_species(ispecies)) CYCLE

          CALL general_dist_fn(sdf_handle, current%name, current%directions, &
              current%ranges, current%resolution, ispecies, &
              current%restrictions, current%use_restrictions, current%ndims, &
              errcode)

          ! If there was an error writing the dist_fn then ignore it in future
          IF (errcode .NE. 0) current%dumpmask = c_io_never
        ENDDO
      ENDIF
      current=>current%next
    ENDDO

  END SUBROUTINE write_dist_fns



  SUBROUTINE general_dist_fn(sdf_handle, name, direction, ranges_in, &
      resolution_in, species, restrictions, use_restrictions, curdims, errcode)

    TYPE(sdf_file_handle) :: sdf_handle
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, DIMENSION(c_df_maxdims), INTENT(IN) :: direction
    REAL(num), DIMENSION(2,c_df_maxdims), INTENT(IN) :: ranges_in
    INTEGER, DIMENSION(c_df_maxdims), INTENT(IN) :: resolution_in
    INTEGER, INTENT(IN) :: species
    REAL(num), DIMENSION(2,c_df_maxdirs), INTENT(IN) :: restrictions
    LOGICAL, DIMENSION(c_df_maxdirs), INTENT(IN) :: use_restrictions
    INTEGER, INTENT(IN) :: curdims
    INTEGER, INTENT(OUT) :: errcode

    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: array, array_tmp
    REAL(num), DIMENSION(:), ALLOCATABLE :: grid1, grid2, grid3
    LOGICAL, DIMENSION(c_df_maxdims) :: parallel
    REAL(num), DIMENSION(c_df_maxdims) :: dgrid
    REAL(num) :: current_data, temp_data
    INTEGER :: idim, idir
    LOGICAL, DIMENSION(c_df_maxdims) :: calc_range
    LOGICAL :: calc_ranges

    LOGICAL :: use_x, use_y, use_z, need_reduce
    INTEGER, DIMENSION(c_df_maxdims) :: start_local, global_resolution
    INTEGER :: color, comm_new
    INTEGER :: new_type, array_type

    REAL(num), DIMENSION(2,c_df_maxdims) :: ranges
    INTEGER, DIMENSION(c_df_maxdims) :: resolution
    INTEGER, DIMENSION(c_df_maxdims) :: cell
    REAL(num) :: part_weight, part_mc, part_mc2, gamma_m1, start

    TYPE(particle), POINTER :: current
    CHARACTER(LEN=string_length) :: var_name
    CHARACTER(LEN=8), DIMENSION(c_df_maxdirs) :: labels, units
    REAL(num), DIMENSION(c_df_maxdirs) :: particle_data

    errcode = 0
    use_x = .FALSE.
    use_y = .FALSE.
    use_z = .FALSE.
    color = 0
    ranges = ranges_in
    resolution = resolution_in
    global_resolution = resolution
    parallel = .FALSE.
    start_local = 1
    calc_range = .FALSE.
    calc_ranges = .FALSE.

    current_data = 0.0_num
#ifndef PER_PARTICLE_CHARGE_MASS
    part_mc  = species_list(species)%mass * c
    part_mc2 = part_mc * c
#endif
#ifndef PER_PARTICLE_WEIGHT
    part_weight = species_list(species)%weight
#endif

    DO idim = 1, curdims
      IF (direction(idim) .EQ. c_dir_x) THEN
        use_x = .TRUE.
        resolution(idim) = nx
        ranges(1,idim) = x_min_local - 0.5_num * dx
        ranges(2,idim) = x_max_local + 0.5_num * dx
        start_local(idim) = cell_x_min(x_coords+1)
        global_resolution(idim) = nx_global
        dgrid(idim) = dx
        labels(idim) = 'X'
        units(idim)  = 'm'
        parallel(idim) = .TRUE.
        CYCLE

      ELSE IF (direction(idim) .EQ. c_dir_y) THEN
        use_y = .TRUE.
        resolution(idim) = ny
        ranges(1,idim) = y_min_local - 0.5_num * dy
        ranges(2,idim) = y_max_local + 0.5_num * dy
        start_local(idim) = cell_y_min(y_coords+1)
        global_resolution(idim) = ny_global
        dgrid(idim) = dy
        labels(idim) = 'Y'
        units(idim)  = 'm'
        parallel(idim) = .TRUE.
        CYCLE

      ELSE IF (direction(idim) .EQ. c_dir_z) THEN
        use_z = .TRUE.
        resolution(idim) = nz
        ranges(1,idim) = z_min_local - 0.5_num * dz
        ranges(2,idim) = z_max_local + 0.5_num * dz
        start_local(idim) = cell_z_min(z_coords+1)
        global_resolution(idim) = nz_global
        dgrid(idim) = dz
        labels(idim) = 'Z'
        units(idim)  = 'm'
        parallel(idim) = .TRUE.
        CYCLE

      ENDIF

      ! If we're here then this must be a momentum space direction
      ! So determine which momentum space directions are needed
      IF (ranges(1,idim) .EQ. ranges(2,idim)) THEN
        calc_range(idim) = .TRUE.
        calc_ranges = .TRUE.
      ENDIF

      IF (direction(idim) .EQ. c_dir_px) THEN
        labels(idim) = 'Px'
        units(idim)  = 'm/s'

      ELSE IF (direction(idim) .EQ. c_dir_py) THEN
        labels(idim) = 'Py'
        units(idim)  = 'm/s'

      ELSE IF (direction(idim) .EQ. c_dir_pz) THEN
        labels(idim) = 'Pz'
        units(idim)  = 'm/s'

      ELSE IF (direction(idim) .EQ. c_dir_en) THEN
        labels(idim) = 'en'
        units(idim)  = 'J'

      ELSE IF (direction(idim) .EQ. c_dir_gamma_m1) THEN
        labels(idim) = 'gamma-1'
        units(idim)  = ''

      ELSE
        IF (rank .EQ. 0) THEN
          WRITE(*,*) '*** WARNING ***'
          WRITE(*,*) 'Unable to write dist_fn. Ignoring.'
        ENDIF
        errcode = 1
        RETURN

      ENDIF
    ENDDO

    ! Calculate range for directions where needed
    IF (calc_ranges) THEN
      DO idim = 1, curdims
        IF (calc_range(idim)) THEN
          ranges(1,idim) = 1.0e6_num
          ranges(2,idim) = -1.0e6_num
        ENDIF
      ENDDO
      current=>species_list(species)%attached_list%head

      out1: DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc  = current%mass * c
        part_mc2 = part_mc * c
#endif
        gamma_m1 = SQRT(SUM((current%part_p / part_mc)**2) + 1.0_num) - 1.0_num

        particle_data(1:c_ndims) = current%part_pos
        particle_data(c_dir_px:c_dir_pz) = current%part_p
        particle_data(c_dir_en) = gamma_m1 * part_mc2
        particle_data(c_dir_gamma_m1) = gamma_m1

        current=>current%next

        DO idim = 1, curdims
          IF (calc_range(idim)) THEN
            DO idir = 1, c_df_maxdirs
              IF (use_restrictions(idir) &
                  .AND. (particle_data(idir) .LT. restrictions(1,idir) &
                  .OR. particle_data(idir) .GT. restrictions(2,idir))) &
                      CYCLE out1
            ENDDO

            current_data = particle_data(direction(idim))
            IF (current_data .LT. ranges(1,idim)) ranges(1,idim) = current_data
            IF (current_data .GT. ranges(2,idim)) ranges(2,idim) = current_data
          ENDIF
        ENDDO
      ENDDO out1

      DO idim = 1, curdims
        IF (.NOT. parallel(idim)) THEN
          ! If not parallel then this is a momentum dimension
          CALL MPI_ALLREDUCE(ranges(1,idim), temp_data, 1, mpireal, MPI_MIN, &
              comm, errcode)
          ranges(1,idim) = temp_data
          CALL MPI_ALLREDUCE(ranges(2,idim), temp_data, 1, mpireal, MPI_MAX, &
              comm, errcode)
          ranges(2,idim) = temp_data
        ENDIF
      ENDDO
    ENDIF

    DO idim = 1, curdims
      ! Fix so that if distribution function is zero then it picks an arbitrary
      ! scale in that direction
      IF (ranges(1,idim) .EQ. ranges(2,idim)) THEN
        ranges(1,idim) = -1.0_num
        ranges(2,idim) = 1.0_num
      ENDIF

      ! Calculate grid spacing
      IF (.NOT. parallel(idim)) dgrid(idim) = &
          (ranges(2,idim) - ranges(1,idim)) / REAL(resolution(idim), num)
    ENDDO

    ALLOCATE(array(resolution(1), resolution(2), resolution(3)))
    array = 0.0_num

    current=>species_list(species)%attached_list%head
    out2: DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
      part_mc  = current%mass * c
      part_mc2 = part_mc * c
#endif
#ifdef PER_PARTICLE_WEIGHT
      part_weight = current%weight
#endif
      gamma_m1 = SQRT(SUM((current%part_p / part_mc)**2) + 1.0_num) - 1.0_num

      particle_data(1:c_ndims) = current%part_pos
      particle_data(c_dir_px:c_dir_pz) = current%part_p
      particle_data(c_dir_en) = gamma_m1 * part_mc2
      particle_data(c_dir_gamma_m1) = gamma_m1

      current=>current%next

      DO idir = 1, c_df_maxdirs
        IF (use_restrictions(idir) &
            .AND. (particle_data(idir) .LT. restrictions(1,idir) &
            .OR. particle_data(idir) .GT. restrictions(2,idir))) &
                CYCLE out2
      ENDDO

      cell = 1
      DO idim = 1, curdims
        current_data = particle_data(direction(idim))
        cell(idim) = FLOOR((current_data - ranges(1,idim)) / dgrid(idim)) + 1
        IF (cell(idim) .LT. 1 .OR. cell(idim) .GT. resolution(idim)) &
            CYCLE out2
      ENDDO

      array(cell(1), cell(2), cell(3)) = &
          array(cell(1), cell(2), cell(3)) + part_weight ! * real_space_area
    ENDDO out2

    need_reduce = .TRUE.
    IF (use_x .AND. use_y .AND. use_z) need_reduce = .FALSE.

    IF (need_reduce) THEN
      ! If using x direction need to reduce across y, z
      IF (use_x) color = color + x_coords
      ! If using y direction need to reduce across x, z
      IF (use_y) color = color + nprocx * y_coords
      ! If using z direction need to reduce across x, y
      IF (use_z) color = color + nprocx * nprocy * z_coords

      CALL MPI_COMM_SPLIT(comm, color, rank, comm_new, errcode)
      ALLOCATE(array_tmp(resolution(1), resolution(2), resolution(3)))
      array_tmp = 0.0_num
      CALL MPI_ALLREDUCE(array, array_tmp, &
          resolution(1)*resolution(2)*resolution(3), mpireal, MPI_SUM, &
          comm_new, errcode)
      array = array_tmp
      DEALLOCATE(array_tmp)
      CALL MPI_COMM_FREE(comm_new, errcode)
    ENDIF

    ! Create grids
    ALLOCATE(grid1(global_resolution(1)))
    start = ranges(1,1) + 0.5_num * dgrid(1)
    DO idir = 1, global_resolution(1)
      grid1(idir) = start + (idir - 1) * dgrid(1)
    ENDDO

    IF (curdims .GE. 2) THEN
      ALLOCATE(grid2(global_resolution(2)))
      start = ranges(1,2) + 0.5_num * dgrid(2)
      DO idir = 1, global_resolution(2)
        grid2(idir) = start + (idir - 1) * dgrid(2)
      ENDDO
    ENDIF

    IF (curdims .GE. 3) THEN
      ALLOCATE(grid3(global_resolution(3)))
      start = ranges(1,3) + 0.5_num * dgrid(3)
      DO idir = 1, global_resolution(3)
        grid3(idir) = start + (idir - 1) * dgrid(3)
      ENDDO
    ENDIF

    var_name = TRIM(name) // '/' // TRIM(species_list(species)%name)

    IF (use_offset_grid) THEN
      IF (curdims .EQ. 1) THEN
        CALL sdf_write_srl_plain_mesh(sdf_handle, &
            'grid_full/' // TRIM(var_name), 'Grid_Full/' // TRIM(var_name), &
            grid1, labels, units)
      ELSE IF (curdims .EQ. 2) THEN
        CALL sdf_write_srl_plain_mesh(sdf_handle, &
            'grid_full/' // TRIM(var_name), 'Grid_Full/' // TRIM(var_name), &
            grid1, grid2, labels, units)
      ELSE IF (curdims .EQ. 3) THEN
        CALL sdf_write_srl_plain_mesh(sdf_handle, &
            'grid_full/' // TRIM(var_name), 'Grid_Full/' // TRIM(var_name), &
            grid1, grid2, grid3, labels, units)
      ENDIF
      IF (parallel(1)) grid1 = grid1 - ranges(1,1)
      IF (parallel(2)) grid2 = grid2 - ranges(1,2)
      IF (parallel(3)) grid3 = grid3 - ranges(1,3)
    ENDIF

    IF (curdims .EQ. 1) THEN
      CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid/' // TRIM(var_name), &
          'Grid/' // TRIM(var_name), grid1, labels, units)
      DEALLOCATE(grid1)
      new_type = &
          create_1d_array_subtype(resolution, global_resolution, start_local)
    ELSE IF (curdims .EQ. 2) THEN
      CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid/' // TRIM(var_name), &
          'Grid/' // TRIM(var_name), grid1, grid2, labels, units)
      DEALLOCATE(grid1, grid2)
      new_type = &
          create_2d_array_subtype(resolution, global_resolution, start_local)
    ELSE IF (curdims .EQ. 3) THEN
      CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid/' // TRIM(var_name), &
          'Grid/' // TRIM(var_name), grid1, grid2, grid3, labels, units)
      DEALLOCATE(grid1, grid2, grid3)
      new_type = &
          create_3d_array_subtype(resolution, global_resolution, start_local)
    ENDIF

    CALL MPI_TYPE_CONTIGUOUS(resolution(1) * resolution(2) * resolution(3), &
        mpireal, array_type, errcode)
    CALL MPI_TYPE_COMMIT(array_type, errcode)

    CALL sdf_write_plain_variable(sdf_handle, TRIM(var_name), &
        'dist_fn/' // TRIM(var_name), 'npart/cell', global_resolution, &
        c_stagger_vertex, 'grid/' // TRIM(var_name), array, new_type, &
        array_type)

    CALL MPI_TYPE_FREE(new_type, errcode)
    CALL MPI_TYPE_FREE(array_type, errcode)

    DEALLOCATE(array)

  END SUBROUTINE general_dist_fn

END MODULE dist_fn
