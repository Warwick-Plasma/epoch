MODULE dist_fn

  USE mpi_subtype_control

  IMPLICIT NONE

CONTAINS

  SUBROUTINE attach_dist_fn(block)

    TYPE(distribution_function_block), POINTER :: block
    TYPE(distribution_function_block), POINTER :: current

    current => dist_fns
    IF (.NOT. ASSOCIATED(current)) THEN
      ! This is the first distribution function to add
      dist_fns => block
      RETURN
    ENDIF
    DO WHILE(ASSOCIATED(current%next))
      current => current%next
    ENDDO
    current%next => block

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

    current => dist_fns
    DO WHILE(ASSOCIATED(current))
      next => current%next
      DEALLOCATE(current)
      current => next
    ENDDO

  END SUBROUTINE clean_dist_fns



  SUBROUTINE write_dist_fns(sdf_handle, code)

    TYPE(sdf_file_handle) :: sdf_handle
    INTEGER, INTENT(IN) :: code

    INTEGER :: ispecies, errcode
    TYPE(distribution_function_block), POINTER :: current
    LOGICAL :: convert

    ! Write the distribution functions
    current => dist_fns
    DO WHILE(ASSOCIATED(current))
      IF (IAND(current%dumpmask, code) .NE. 0) THEN
        DO ispecies = 1, n_species
          IF (.NOT. current%use_species(ispecies)) CYCLE

          convert = (IAND(IOR(dumpmask(c_dump_dist_fns),current%dumpmask), &
                          c_io_dump_single) .NE. 0)

          CALL general_dist_fn(sdf_handle, current%name, current%directions, &
              current%ranges, current%resolution, ispecies, &
              current%restrictions, current%use_restrictions, current%ndims, &
              convert, errcode)

          ! If there was an error writing the dist_fn then ignore it in future
          IF (errcode .NE. 0) current%dumpmask = c_io_never
        ENDDO
      ENDIF
      current => current%next
    ENDDO

  END SUBROUTINE write_dist_fns



  SUBROUTINE general_dist_fn(sdf_handle, name, direction, ranges_in, &
      resolution_in, species, restrictions, use_restrictions, curdims, &
      convert, errcode)

    TYPE(sdf_file_handle) :: sdf_handle
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, DIMENSION(c_df_maxdims), INTENT(IN) :: direction
    REAL(num), DIMENSION(2,c_df_maxdims), INTENT(IN) :: ranges_in
    INTEGER, DIMENSION(c_df_maxdims), INTENT(IN) :: resolution_in
    INTEGER, INTENT(IN) :: species
    REAL(num), DIMENSION(2,c_df_maxdirs), INTENT(IN) :: restrictions
    LOGICAL, DIMENSION(c_df_maxdirs), INTENT(IN) :: use_restrictions
    INTEGER, INTENT(IN) :: curdims
    LOGICAL, INTENT(IN) :: convert
    INTEGER, INTENT(OUT) :: errcode

    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: array, array_tmp
    REAL(num), DIMENSION(:), ALLOCATABLE :: grid1, grid2, grid3
    LOGICAL, DIMENSION(c_df_maxdims) :: parallel
    REAL(num), DIMENSION(c_df_maxdims) :: dgrid
    REAL(num) :: current_data, temp_data, theta, p, px, py, pz
    INTEGER :: idim, idir
    LOGICAL, DIMENSION(c_df_maxdims) :: calc_range
    LOGICAL :: calc_ranges

    LOGICAL :: use_x, need_reduce
    LOGICAL :: use_xy_angle, use_yz_angle, use_zx_angle
    INTEGER, DIMENSION(c_df_maxdims) :: start_local, global_resolution
    INTEGER :: new_type, array_type

    REAL(num), DIMENSION(2,c_df_maxdims) :: ranges
    INTEGER, DIMENSION(c_df_maxdims) :: resolution
    INTEGER, DIMENSION(c_df_maxdims) :: cell
    REAL(num) :: part_weight, part_mc, part_mc2, gamma_m1, start
    REAL(num) :: xy_max, yz_max, zx_max
    REAL(num), PARAMETER :: pi2 = 2.0_num * pi

    TYPE(particle), POINTER :: current, next
    CHARACTER(LEN=string_length) :: var_name
    CHARACTER(LEN=8), DIMENSION(c_df_maxdirs) :: labels, units
    REAL(num), DIMENSION(c_df_maxdirs) :: particle_data

    errcode = 0
    ! Update species count if necessary
    IF (io_list(species)%count_update_step .LT. step) THEN
      CALL MPI_ALLREDUCE(io_list(species)%attached_list%count, &
          io_list(species)%count, 1, MPI_INTEGER8, MPI_SUM, &
          comm, errcode)
      io_list(species)%count_update_step = step
    ENDIF

    IF (io_list(species)%count .LT. 1) RETURN

    use_x = .FALSE.
    ranges = ranges_in
    resolution = resolution_in
    global_resolution = resolution
    parallel = .FALSE.
    start_local = 1
    calc_range = .FALSE.
    calc_ranges = .FALSE.
    use_xy_angle = .FALSE.
    use_yz_angle = .FALSE.
    use_zx_angle = .FALSE.

    current_data = 0.0_num
#ifndef PER_PARTICLE_CHARGE_MASS
    part_mc  = io_list(species)%mass * c
    part_mc2 = part_mc * c
#endif
#ifndef PER_PARTICLE_WEIGHT
    part_weight = io_list(species)%weight
#endif

    DO idim = 1, curdims
      IF (direction(idim) .EQ. c_dir_x) THEN
        use_x = .TRUE.
        resolution(idim) = nx
        ranges(1,idim) = x_grid_min_local - 0.5_num * dx
        ranges(2,idim) = x_grid_max_local + 0.5_num * dx
        start_local(idim) = nx_global_min
        global_resolution(idim) = nx_global
        dgrid(idim) = dx
        labels(idim) = 'X'
        units(idim)  = 'm'
        parallel(idim) = .TRUE.
        CYCLE

      ENDIF

      ! If we're here then this must be a momentum space direction
      ! So determine which momentum space directions are needed
      IF (ABS(ranges(1,idim) - ranges(2,idim)) .LE. c_tiny) THEN
        calc_range(idim) = .TRUE.
        calc_ranges = .TRUE.
      ENDIF

      IF (direction(idim) .EQ. c_dir_px) THEN
        labels(idim) = 'Px'
        units(idim)  = 'kg.m/s'

      ELSE IF (direction(idim) .EQ. c_dir_py) THEN
        labels(idim) = 'Py'
        units(idim)  = 'kg.m/s'

      ELSE IF (direction(idim) .EQ. c_dir_pz) THEN
        labels(idim) = 'Pz'
        units(idim)  = 'kg.m/s'

      ELSE IF (direction(idim) .EQ. c_dir_en) THEN
        labels(idim) = 'en'
        units(idim)  = 'J'

      ELSE IF (direction(idim) .EQ. c_dir_gamma_m1) THEN
        labels(idim) = 'gamma-1'
        units(idim)  = ''

      ELSE IF (direction(idim) .EQ. c_dir_xy_angle) THEN
        use_xy_angle = .TRUE.
        labels(idim) = 'xy_angle'
        units(idim)  = 'radians'

      ELSE IF (direction(idim) .EQ. c_dir_yz_angle) THEN
        use_yz_angle = .TRUE.
        labels(idim) = 'yz_angle'
        units(idim)  = 'radians'

      ELSE IF (direction(idim) .EQ. c_dir_zx_angle) THEN
        use_zx_angle = .TRUE.
        labels(idim) = 'zx_angle'
        units(idim)  = 'radians'

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
          ranges(1,idim) =  HUGE(1.0_num) * 0.4_num
          ranges(2,idim) = -HUGE(1.0_num) * 0.4_num
        ENDIF
      ENDDO
      next => io_list(species)%attached_list%head

      out1: DO WHILE(ASSOCIATED(next))
        current => next
        next => current%next

#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc  = current%mass * c
        part_mc2 = part_mc * c
#endif
        gamma_m1 = SQRT(SUM((current%part_p / part_mc)**2) + 1.0_num) - 1.0_num
        px = current%part_p(1)
        py = current%part_p(2)
        pz = current%part_p(3)

        particle_data(1:c_ndims) = current%part_pos
        IF (io_list(species)%species_type .EQ. c_species_id_photon) THEN
          particle_data(c_dir_gamma_m1) = 0.0_num
#ifdef PHOTONS
          particle_data(c_dir_px) = px
          particle_data(c_dir_py) = py
          particle_data(c_dir_pz) = pz
          particle_data(c_dir_en) = current%particle_energy
#else
          particle_data(c_dir_px:c_dir_pz) = 0.0_num
          particle_data(c_dir_en) = 0.0_num
#endif
          ! Can't define gamma for photon so one is as good as anything
          particle_data(c_dir_gamma_m1) = 1.0_num
        ELSE
          particle_data(c_dir_px) = px
          particle_data(c_dir_py) = py
          particle_data(c_dir_pz) = pz
          particle_data(c_dir_en) = gamma_m1 * part_mc2
          particle_data(c_dir_gamma_m1) = gamma_m1
        ENDIF

        IF (use_xy_angle) THEN
          p = SQRT(px**2 + py**2)
          IF (p .LT. c_tiny) CYCLE

          theta = ASIN(py / p)
          IF (px .LT. 0.0_num) theta = SIGN(1.0_num,py) * pi - theta
          particle_data(c_dir_xy_angle) = theta
        ENDIF

        IF (use_yz_angle) THEN
          p = SQRT(py**2 + pz**2)
          IF (p .LT. c_tiny) CYCLE

          theta = ASIN(pz / p)
          IF (py .LT. 0.0_num) theta = SIGN(1.0_num,pz) * pi - theta
          particle_data(c_dir_yz_angle) = theta
        ENDIF

        IF (use_zx_angle) THEN
          p = SQRT(pz**2 + px**2)
          IF (p .LT. c_tiny) CYCLE

          theta = ASIN(px / p)
          IF (pz .LT. 0.0_num) theta = SIGN(1.0_num,px) * pi - theta
          particle_data(c_dir_zx_angle) = theta
        ENDIF

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
      IF (ABS(ranges(1,idim) - ranges(2,idim)) .LE. c_tiny) THEN
        ranges(1,idim) = -1.0_num
        ranges(2,idim) = 1.0_num
      ENDIF

      IF (direction(idim) .EQ. c_dir_xy_angle) THEN
        xy_max = ranges(2,idim)
      ELSE IF (direction(idim) .EQ. c_dir_yz_angle) THEN
        yz_max = ranges(2,idim)
      ELSE IF (direction(idim) .EQ. c_dir_zx_angle) THEN
        zx_max = ranges(2,idim)
      ENDIF

      ! Calculate grid spacing
      IF (.NOT. parallel(idim)) dgrid(idim) = &
          (ranges(2,idim) - ranges(1,idim)) / REAL(resolution(idim), num)
    ENDDO

    ALLOCATE(array(resolution(1), resolution(2), resolution(3)))
    array = 0.0_num

    next => io_list(species)%attached_list%head

    out2: DO WHILE(ASSOCIATED(next))
      current => next
      next => current%next

#ifdef PER_PARTICLE_CHARGE_MASS
      part_mc  = current%mass * c
      part_mc2 = part_mc * c
#endif
#ifdef PER_PARTICLE_WEIGHT
      part_weight = current%weight
#endif
      gamma_m1 = SQRT(SUM((current%part_p / part_mc)**2) + 1.0_num) - 1.0_num
      px = current%part_p(1)
      py = current%part_p(2)
      pz = current%part_p(3)

      particle_data(1:c_ndims) = current%part_pos
      IF (io_list(species)%species_type .EQ. c_species_id_photon) THEN
        particle_data(c_dir_gamma_m1) = 0.0_num
#ifdef PHOTONS
        particle_data(c_dir_px) = px
        particle_data(c_dir_py) = py
        particle_data(c_dir_pz) = pz
        particle_data(c_dir_en) = current%particle_energy
#else
        particle_data(c_dir_px:c_dir_pz) = 0.0_num
        particle_data(c_dir_en) = 0.0_num
#endif
        ! Can't define gamma for photon so one is as good as anything
        particle_data(c_dir_gamma_m1) = 1.0_num
      ELSE
        particle_data(c_dir_px) = px
        particle_data(c_dir_py) = py
        particle_data(c_dir_pz) = pz
        particle_data(c_dir_en) = gamma_m1 * part_mc2
        particle_data(c_dir_gamma_m1) = gamma_m1
      ENDIF

      IF (use_xy_angle) THEN
        p = SQRT(px**2 + py**2)
        IF (p .LT. c_tiny) CYCLE

        theta = ASIN(py / p)
        IF (px .LT. 0.0_num) theta = SIGN(1.0_num,py) * pi - theta
        IF ((theta + pi2) .LT. xy_max) theta = theta + pi2
        particle_data(c_dir_xy_angle) = theta
      ENDIF

      IF (use_yz_angle) THEN
        p = SQRT(py**2 + pz**2)
        IF (p .LT. c_tiny) CYCLE

        theta = ASIN(pz / p)
        IF (py .LT. 0.0_num) theta = SIGN(1.0_num,pz) * pi - theta
        IF ((theta + pi2) .LT. yz_max) theta = theta + pi2
        particle_data(c_dir_yz_angle) = theta
      ENDIF

      IF (use_zx_angle) THEN
        p = SQRT(pz**2 + px**2)
        IF (p .LT. c_tiny) CYCLE

        theta = ASIN(px / p)
        IF (pz .LT. 0.0_num) theta = SIGN(1.0_num,px) * pi - theta
        IF ((theta + pi2) .LT. zx_max) theta = theta + pi2
        particle_data(c_dir_zx_angle) = theta
      ENDIF

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
    IF (use_x) need_reduce = .FALSE.

    IF (need_reduce) THEN
      ALLOCATE(array_tmp(resolution(1), resolution(2), resolution(3)))
      array_tmp = 0.0_num
      CALL MPI_ALLREDUCE(array, array_tmp, &
          resolution(1)*resolution(2)*resolution(3), mpireal, MPI_SUM, &
          comm, errcode)
      array = array_tmp
      DEALLOCATE(array_tmp)
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

    var_name = TRIM(name) // '/' // TRIM(io_list(species)%name)

    IF (use_offset_grid) THEN
      IF (curdims .EQ. 1) THEN
        CALL sdf_write_srl_plain_mesh(sdf_handle, &
            'grid_full/' // TRIM(var_name), 'Grid_Full/' // TRIM(var_name), &
            grid1, convert, labels, units)
      ELSE IF (curdims .EQ. 2) THEN
        CALL sdf_write_srl_plain_mesh(sdf_handle, &
            'grid_full/' // TRIM(var_name), 'Grid_Full/' // TRIM(var_name), &
            grid1, grid2, convert, labels, units)
      ELSE IF (curdims .EQ. 3) THEN
        CALL sdf_write_srl_plain_mesh(sdf_handle, &
            'grid_full/' // TRIM(var_name), 'Grid_Full/' // TRIM(var_name), &
            grid1, grid2, grid3, convert, labels, units)
      ENDIF
      IF (parallel(1)) grid1 = grid1 - ranges(1,1)
      IF (parallel(2)) grid2 = grid2 - ranges(1,2)
      IF (parallel(3)) grid3 = grid3 - ranges(1,3)
    ENDIF

    IF (curdims .EQ. 1) THEN
      CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid/' // TRIM(var_name), &
          'Grid/' // TRIM(var_name), grid1, convert, labels, units)
      DEALLOCATE(grid1)
      new_type = create_1d_array_subtype(mpireal, resolution, &
          global_resolution, start_local)
    ELSE IF (curdims .EQ. 2) THEN
      CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid/' // TRIM(var_name), &
          'Grid/' // TRIM(var_name), grid1, grid2, convert, labels, units)
      DEALLOCATE(grid1, grid2)
      new_type = create_2d_array_subtype(mpireal, resolution, &
          global_resolution, start_local)
    ELSE IF (curdims .EQ. 3) THEN
      CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid/' // TRIM(var_name), &
          'Grid/' // TRIM(var_name), grid1, grid2, grid3, convert, labels, &
          units)
      DEALLOCATE(grid1, grid2, grid3)
      new_type = create_3d_array_subtype(mpireal, resolution, &
          global_resolution, start_local)
    ENDIF

    CALL MPI_TYPE_CONTIGUOUS(resolution(1) * resolution(2) * resolution(3), &
        mpireal, array_type, errcode)
    CALL MPI_TYPE_COMMIT(array_type, errcode)

    CALL sdf_write_plain_variable(sdf_handle, TRIM(var_name), &
        'dist_fn/' // TRIM(var_name), 'npart/cell', global_resolution, &
        c_stagger_vertex, 'grid/' // TRIM(var_name), array, new_type, &
        array_type, convert)

    CALL MPI_TYPE_FREE(new_type, errcode)
    CALL MPI_TYPE_FREE(array_type, errcode)

    DEALLOCATE(array)

  END SUBROUTINE general_dist_fn

END MODULE dist_fn
