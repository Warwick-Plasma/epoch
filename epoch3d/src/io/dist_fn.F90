MODULE dist_fn

  USE cfd
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



  SUBROUTINE setup_dist_fn(block)

    TYPE(distribution_function_block), POINTER :: block

    NULLIFY(block%next)
    ALLOCATE(block%use_species(n_species))
    block%use_species = .FALSE.
    block%ranges = 1.0_num
    block%use_restrictions = .FALSE.
    block%ndims = -1

  END SUBROUTINE setup_dist_fn



  SUBROUTINE clean_dist_fns()

    TYPE(distribution_function_block), POINTER :: current, next

    current=>dist_fns
    DO WHILE(ASSOCIATED(current))
      next=>current%next
      DEALLOCATE(current)
      current=>next
    ENDDO

  END SUBROUTINE clean_dist_fns



  SUBROUTINE write_dist_fns(code)

    INTEGER, INTENT(IN) :: code

    INTEGER :: ispecies
    REAL(num), DIMENSION(2,c_df_maxdims) :: ranges
    LOGICAL, DIMENSION(c_df_maxdirs) :: use_restrictions
    REAL(num), DIMENSION(2,c_df_maxdirs) :: restrictions
    INTEGER, DIMENSION(c_df_maxdims) :: resolution = (/100, 100, 100/)
    TYPE(distribution_function_block), POINTER :: current

    ! Write the distribution functions
    current=>dist_fns
    DO WHILE(ASSOCIATED(current))
      IF (IAND(current%dumpmask, code) .NE. 0) THEN
        DO ispecies = 1, n_species
          IF (.NOT. current%use_species(ispecies)) CYCLE
          ranges = current%ranges
          resolution = current%resolution
          restrictions = current%restrictions
          use_restrictions = current%use_restrictions

          IF (current%ndims .EQ. 1) THEN
            CALL general_1d_dist_fn(current%name, current%directions, &
                ranges, resolution, ispecies, restrictions, use_restrictions)
          ELSE IF (current%ndims .EQ. 2) THEN
            CALL general_2d_dist_fn(current%name, current%directions, &
                ranges, resolution, ispecies, restrictions, use_restrictions)
          ELSE IF (current%ndims .EQ. 3) THEN
            CALL general_3d_dist_fn(current%name, current%directions, &
                ranges, resolution, ispecies, restrictions, use_restrictions)
          ENDIF
        ENDDO
      ENDIF
      current=>current%next
    ENDDO

  END SUBROUTINE write_dist_fns



  SUBROUTINE general_1d_dist_fn(name, direction, ranges, resolution, species, &
      restrictions, use_restrictions)

    INTEGER, PARAMETER :: c_df_curdims = 1
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, DIMENSION(c_df_curdims), INTENT(IN) :: direction
    REAL(num), DIMENSION(2,c_df_curdims), INTENT(INOUT) :: ranges
    INTEGER, DIMENSION(c_df_curdims), INTENT(INOUT) :: resolution
    INTEGER, INTENT(IN) :: species
    REAL(num), DIMENSION(2,c_df_maxdirs), INTENT(IN) :: restrictions
    LOGICAL, DIMENSION(c_df_maxdirs), INTENT(IN) :: use_restrictions

    REAL(num), DIMENSION(:), ALLOCATABLE :: data, data2
    REAL(num), DIMENSION(:), ALLOCATABLE :: grid1
    LOGICAL, DIMENSION(c_df_curdims) :: parallel
    REAL(num), DIMENSION(c_df_curdims) :: dgrid
    REAL(num) :: current_data, temp_data
    INTEGER :: idim, idir
    LOGICAL, DIMENSION(c_df_curdims) :: calc_range
    LOGICAL :: calc_ranges

    LOGICAL :: use_x, use_y, use_z, need_reduce
    INTEGER, DIMENSION(c_df_curdims) :: start_local, global_resolution
    INTEGER :: color, comm_new
    INTEGER :: new_type, array_type

    LOGICAL, DIMENSION(c_df_curdims, c_df_maxdirs) :: use_direction
    LOGICAL, DIMENSION(c_df_curdims) :: calc_mod
    INTEGER, DIMENSION(c_df_curdims) :: p_count
    INTEGER, DIMENSION(c_df_curdims) :: l_direction
    REAL(num), DIMENSION(c_df_curdims) :: conv
    INTEGER, DIMENSION(c_df_curdims) :: cell
    LOGICAL :: use_this
    REAL(num) :: real_space_area, part_weight
    REAL(num) :: part_mass, part_mass_c, gamma_mass_c

    TYPE(particle), POINTER :: current
    CHARACTER(LEN=string_length) :: grid_name, norm_grid_name, var_name
    REAL(num), DIMENSION(c_df_curdims) :: stagger = 0.0_num
    REAL(num), DIMENSION(c_df_maxdirs) :: particle_data

    REAL(num) :: max_p_conv

    INTEGER :: ind

    use_x = .FALSE.
    use_y = .FALSE.
    color = 0
    global_resolution = resolution
    parallel = .FALSE.
    start_local = 1
    calc_range = .FALSE.
    calc_ranges = .FALSE.
    p_count = 0
    use_direction = .FALSE.
    l_direction = 0

    real_space_area = 1.0_num
    current_data = 0.0_num
    part_mass = particle_species(species)%mass
    part_mass_c = part_mass * c

    DO idim = 1, c_df_curdims
      IF (direction(idim) .EQ. c_dir_x) THEN
        use_x = .TRUE.
        resolution(idim) = nx
        ranges(1,idim) = x_min_local
        ranges(2,idim) = x_max_local
        start_local(idim) = cell_x_min(coordinates(c_ndims)+1)
        global_resolution(idim) = nx_global
        dgrid(idim) = dx
        parallel(idim) = .TRUE.
        use_direction(idim,1) = .TRUE.
        l_direction(idim) = 1
        conv(idim) = MAX(length_x, length_y, length_z)
        CYCLE

      ELSE IF (direction(idim) .EQ. c_dir_y) THEN
        use_y = .TRUE.
        resolution(idim) = ny
        ranges(1,idim) = y_min_local
        ranges(2,idim) = y_max_local
        start_local(idim) = cell_y_min(coordinates(c_ndims-1)+1)
        global_resolution(idim) = ny_global
        dgrid(idim) = dy
        parallel(idim) = .TRUE.
        use_direction(idim,2) = .TRUE.
        l_direction(idim) = 2
        conv(idim) = MAX(length_x, length_y, length_z)
        CYCLE

      ELSE IF (direction(idim) .EQ. c_dir_z) THEN
        use_z = .TRUE.
        resolution(idim) = nz
        ranges(1,idim) = z_min_local
        ranges(2,idim) = z_max_local
        start_local(idim) = cell_z_min(coordinates(c_ndims-2)+1)
        global_resolution(idim) = nz_global
        dgrid(idim) = dz
        parallel(idim) = .TRUE.
        use_direction(idim,3) = .TRUE.
        l_direction(idim) = 3
        conv(idim) = MAX(length_x, length_y, length_z)
        CYCLE

      ENDIF

      conv(idim) = c*m0
      ! If we're here then this must be a momentum space direction
      ! So determine which momentum space directions are needed
      IF (ranges(1,idim) .EQ. ranges(2,idim)) THEN
        calc_range(idim) = .TRUE.
        calc_ranges = .TRUE.
      ENDIF

      IF (IAND(direction(idim), c_dir_px) .NE. 0) THEN
        use_direction(idim, c_ndims+1) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+1
      ENDIF

      IF (IAND(direction(idim), c_dir_py) .NE. 0) THEN
        use_direction(idim, c_ndims+2) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+2
      ENDIF

      IF (IAND(direction(idim), c_dir_pz) .NE. 0) THEN
        use_direction(idim, c_ndims+3) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+3
      ENDIF

      IF (IAND(direction(idim), c_dir_en) .NE. 0) THEN
        use_direction(idim, c_ndims+4) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+4
      ENDIF

      IF (IAND(direction(idim), c_dir_gamma_m1) .NE. 0) THEN
        use_direction(idim, c_ndims+5) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+5
      ENDIF
    ENDDO

    IF (.NOT. use_x) real_space_area = real_space_area * dx
    IF (.NOT. use_y) real_space_area = real_space_area * dy
    IF (.NOT. use_z) real_space_area = real_space_area * dz

    DO idim = 1, c_df_curdims
      calc_mod(idim) = .FALSE.
      IF (p_count(idim) .GT. 1) calc_mod(idim) = .TRUE.
    ENDDO

    ! Calculate range for directions where needed
    IF (calc_ranges) THEN
      DO idim = 1, c_df_curdims
        IF (calc_range(idim)) THEN
          ranges(1,idim) = 1.0e6_num
          ranges(2,idim) = -1.0e6_num
        ENDIF
      ENDDO
      current=>particle_species(species)%attached_list%head
      ind = 0

      DO WHILE(ASSOCIATED(current))
        particle_data(1:c_ndims) = current%part_pos
        particle_data(c_ndims+1:c_ndims+3) = current%part_p
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mass = current%mass
        part_mass_c = part_mass * c
#endif
        gamma_mass_c = SQRT(SUM(current%part_p**2) + part_mass_c**2)
        particle_data(c_ndims+4) = c * (gamma_mass_c - part_mass_c)
        particle_data(c_ndims+5) = gamma_mass_c / part_mass_c - 1.0_num
        DO idim = 1, c_df_curdims
          IF (calc_range(idim)) THEN
            IF (calc_mod(idim)) THEN
              DO idir = 1, c_df_maxdirs
                IF (use_direction(idim, idir)) &
                    current_data = current_data + particle_data(idir)**2
              ENDDO
              current_data = SQRT(current_data)
            ELSE
              current_data = particle_data(l_direction(idim))
            ENDIF
            use_this = .TRUE.
            DO idir = 1, c_df_maxdirs
              IF (use_restrictions(idir) &
                  .AND. (particle_data(idir) .LT. restrictions(1,idir) &
                  .OR. particle_data(idir) .GT. restrictions(2,idir))) &
                      use_this = .FALSE.
            ENDDO
            IF (use_this) THEN
              IF (current_data .LT. ranges(1,idim)) &
                  ranges(1,idim) = current_data
              IF (current_data .GT. ranges(2,idim)) &
                  ranges(2,idim) = current_data
            ENDIF
          ENDIF
        ENDDO
        ind = ind + 1
        current=>current%next
      ENDDO
    ENDIF

    max_p_conv = -10.0_num
    DO idim = 1, c_df_curdims
      IF (.NOT. parallel(idim)) THEN
        ! If not parallel then this is a momentum dimension
        CALL MPI_ALLREDUCE(ranges(1,idim), temp_data, 1, mpireal, MPI_MIN, &
            comm, errcode)
        ranges(1,idim) = temp_data
        CALL MPI_ALLREDUCE(ranges(2,idim), temp_data, 1, mpireal, MPI_MAX, &
            comm, errcode)
        ranges(2,idim) = temp_data
      ENDIF
      ! Fix so that if distribution function is zero then it picks an arbitrary
      ! scale in that direction
      IF (ranges(1,idim) .EQ. ranges(2,idim)) THEN
        ranges(1,idim) = -1.0_num
        ranges(2,idim) = 1.0_num
      ENDIF
      ! Calculate the maximum range of a momentum direction
      IF (ranges(2,idim) - ranges(1,idim) .GT. max_p_conv &
          .AND. .NOT. parallel(idim)) &
              max_p_conv = ranges(2,idim) - ranges(1,idim)
    ENDDO

    DO idim = 1, c_df_curdims
      IF (.NOT. parallel(idim)) conv(idim) = max_p_conv
    ENDDO

    ! Calculate grid spacing
    DO idim = 1, c_df_curdims
      IF (.NOT. parallel(idim)) dgrid(idim) = &
          (ranges(2,idim) - ranges(1,idim)) / REAL(resolution(idim)-1, num)
    ENDDO

    ALLOCATE(data(resolution(1)))
    data = 0.0_num

    current=>particle_species(species)%attached_list%head
    part_weight = weight
    DO WHILE(ASSOCIATED(current))
      particle_data(1:c_ndims) = current%part_pos
      particle_data(c_ndims+1:c_ndims+3) = current%part_p
#ifdef PER_PARTICLE_CHARGE_MASS
      part_mass = current%mass
      part_mass_c = part_mass * c
#endif
      gamma_mass_c = SQRT(SUM(current%part_p**2) + part_mass_c**2)
      particle_data(c_ndims+4) = c * (gamma_mass_c - part_mass_c)
      particle_data(c_ndims+5) = gamma_mass_c / part_mass_c - 1.0_num
      use_this = .TRUE.
      DO idir = 1, c_df_maxdirs
        IF (use_restrictions(idir) &
            .AND. (particle_data(idir) .LT. restrictions(1,idir) &
            .OR. particle_data(idir) .GT. restrictions(2,idir))) &
                use_this = .FALSE.
      ENDDO
      IF (use_this) THEN
        DO idim = 1, c_df_curdims
          IF (calc_mod(idim)) THEN
            DO idir = 1, c_df_maxdirs
              IF (use_direction(idim, idir)) &
                  current_data = current_data + particle_data(idir)**2
            ENDDO
            current_data = SQRT(current_data)
          ELSE
            current_data = particle_data(l_direction(idim))
          ENDIF
          cell(idim) = NINT((current_data - ranges(1,idim)) / dgrid(idim)) + 1
          IF (cell(idim) .LT. 1 .OR. cell(idim) .GT. resolution(idim)) &
              use_this = .FALSE.
        ENDDO
#ifdef PER_PARTICLE_WEIGHT
        part_weight = current%weight
#endif
        IF (use_this) data(cell(1)) = &
            data(cell(1)) + part_weight ! * real_space_area
      ENDIF
      current=>current%next
    ENDDO

    need_reduce = .TRUE.
    IF (use_x .AND. use_y .AND. use_z) need_reduce = .FALSE.

    IF (need_reduce) THEN
      ! If using x direction need to reduce across y, z
      IF (use_x) color = color + coordinates(c_ndims)
      ! If using y direction need to reduce across x, z
      IF (use_y) color = color + nprocx * coordinates(c_ndims-1)
      ! If using z direction need to reduce across x, y
      IF (use_z) color = color + nprocx * nprocy * coordinates(c_ndims-2)

      CALL MPI_COMM_SPLIT(comm, color, rank, comm_new, errcode)
      ALLOCATE(data2(resolution(1)))
      data2 = 0.0_num
      CALL MPI_ALLREDUCE(data, data2, &
          resolution(1), mpireal, MPI_SUM, &
          comm_new, errcode)
      data = data2
      DEALLOCATE(data2)
      CALL MPI_COMM_FREE(comm_new, errcode)
    ENDIF

    ! Create grids
    ALLOCATE(grid1(global_resolution(1)))

    DO idir = 1, global_resolution(1)
      grid1(idir) = ranges(1,1) + (idir - 1) * dgrid(1)
    ENDDO

    grid_name = "Grid_" // TRIM(name) // "_" &
        // TRIM(particle_species(species)%name)
    norm_grid_name = "Norm_Grid_" // TRIM(name) // "_" &
        // TRIM(particle_species(species)%name)
    var_name = TRIM(name) // "_" // TRIM(particle_species(species)%name)

    CALL cfd_write_1d_cartesian_grid(TRIM(grid_name), "Grid", &
        grid1, 0)
    IF (use_offset_grid) THEN
      IF (parallel(1)) grid1 = grid1 - ranges(1,1)
    ENDIF
    CALL cfd_write_1d_cartesian_grid(TRIM(norm_grid_name), "Grid", &
        grid1/conv(1), 0)

    DEALLOCATE(grid1)

    new_type = &
        create_1d_array_subtype(resolution, global_resolution, start_local)
    CALL MPI_TYPE_CONTIGUOUS(resolution(1), &
        mpireal, array_type, errcode)
    CALL MPI_TYPE_COMMIT(array_type, errcode)

    CALL cfd_write_1d_cartesian_variable_parallel(TRIM(var_name), "dist_fn", &
        global_resolution, stagger, TRIM(norm_grid_name), "Grid", data, &
        new_type, array_type)

    CALL MPI_TYPE_FREE(new_type, errcode)
    CALL MPI_TYPE_FREE(array_type, errcode)

    DEALLOCATE(data)

  END SUBROUTINE general_1d_dist_fn



  SUBROUTINE general_2d_dist_fn(name, direction, ranges, resolution, species, &
      restrictions, use_restrictions)

    INTEGER, PARAMETER :: c_df_curdims = 2
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, DIMENSION(c_df_curdims), INTENT(IN) :: direction
    REAL(num), DIMENSION(2,c_df_curdims), INTENT(INOUT) :: ranges
    INTEGER, DIMENSION(c_df_curdims), INTENT(INOUT) :: resolution
    INTEGER, INTENT(IN) :: species
    REAL(num), DIMENSION(2,c_df_maxdirs), INTENT(IN) :: restrictions
    LOGICAL, DIMENSION(c_df_maxdirs), INTENT(IN) :: use_restrictions

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: data, data2
    REAL(num), DIMENSION(:), ALLOCATABLE :: grid1, grid2
    LOGICAL, DIMENSION(c_df_curdims) :: parallel
    REAL(num), DIMENSION(c_df_curdims) :: dgrid
    REAL(num) :: current_data, temp_data
    INTEGER :: idim, idir
    LOGICAL, DIMENSION(c_df_curdims) :: calc_range
    LOGICAL :: calc_ranges

    LOGICAL :: use_x, use_y, use_z, need_reduce
    INTEGER, DIMENSION(c_df_curdims) :: start_local, global_resolution
    INTEGER :: color, comm_new
    INTEGER :: new_type, array_type

    LOGICAL, DIMENSION(c_df_curdims, c_df_maxdirs) :: use_direction
    LOGICAL, DIMENSION(c_df_curdims) :: calc_mod
    INTEGER, DIMENSION(c_df_curdims) :: p_count
    INTEGER, DIMENSION(c_df_curdims) :: l_direction
    REAL(num), DIMENSION(c_df_curdims) :: conv
    INTEGER, DIMENSION(c_df_curdims) :: cell
    LOGICAL :: use_this
    REAL(num) :: real_space_area, part_weight
    REAL(num) :: part_mass, part_mass_c, gamma_mass_c

    TYPE(particle), POINTER :: current
    CHARACTER(LEN=string_length) :: grid_name, norm_grid_name, var_name
    REAL(num), DIMENSION(c_df_curdims) :: stagger = 0.0_num
    REAL(num), DIMENSION(c_df_maxdirs) :: particle_data

    REAL(num) :: max_p_conv

    INTEGER :: ind

    use_x = .FALSE.
    use_y = .FALSE.
    color = 0
    global_resolution = resolution
    parallel = .FALSE.
    start_local = 1
    calc_range = .FALSE.
    calc_ranges = .FALSE.
    p_count = 0
    use_direction = .FALSE.
    l_direction = 0

    real_space_area = 1.0_num
    current_data = 0.0_num
    part_mass = particle_species(species)%mass
    part_mass_c = part_mass * c

    DO idim = 1, c_df_curdims
      IF (direction(idim) .EQ. c_dir_x) THEN
        use_x = .TRUE.
        resolution(idim) = nx
        ranges(1,idim) = x_min_local
        ranges(2,idim) = x_max_local
        start_local(idim) = cell_x_min(coordinates(c_ndims)+1)
        global_resolution(idim) = nx_global
        dgrid(idim) = dx
        parallel(idim) = .TRUE.
        use_direction(idim,1) = .TRUE.
        l_direction(idim) = 1
        conv(idim) = MAX(length_x, length_y, length_z)
        CYCLE

      ELSE IF (direction(idim) .EQ. c_dir_y) THEN
        use_y = .TRUE.
        resolution(idim) = ny
        ranges(1,idim) = y_min_local
        ranges(2,idim) = y_max_local
        start_local(idim) = cell_y_min(coordinates(c_ndims-1)+1)
        global_resolution(idim) = ny_global
        dgrid(idim) = dy
        parallel(idim) = .TRUE.
        use_direction(idim,2) = .TRUE.
        l_direction(idim) = 2
        conv(idim) = MAX(length_x, length_y, length_z)
        CYCLE

      ELSE IF (direction(idim) .EQ. c_dir_z) THEN
        use_z = .TRUE.
        resolution(idim) = nz
        ranges(1,idim) = z_min_local
        ranges(2,idim) = z_max_local
        start_local(idim) = cell_z_min(coordinates(c_ndims-2)+1)
        global_resolution(idim) = nz_global
        dgrid(idim) = dz
        parallel(idim) = .TRUE.
        use_direction(idim,3) = .TRUE.
        l_direction(idim) = 3
        conv(idim) = MAX(length_x, length_y, length_z)
        CYCLE

      ENDIF

      conv(idim) = c*m0
      ! If we're here then this must be a momentum space direction
      ! So determine which momentum space directions are needed
      IF (ranges(1,idim) .EQ. ranges(2,idim)) THEN
        calc_range(idim) = .TRUE.
        calc_ranges = .TRUE.
      ENDIF

      IF (IAND(direction(idim), c_dir_px) .NE. 0) THEN
        use_direction(idim, c_ndims+1) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+1
      ENDIF

      IF (IAND(direction(idim), c_dir_py) .NE. 0) THEN
        use_direction(idim, c_ndims+2) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+2
      ENDIF

      IF (IAND(direction(idim), c_dir_pz) .NE. 0) THEN
        use_direction(idim, c_ndims+3) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+3
      ENDIF

      IF (IAND(direction(idim), c_dir_en) .NE. 0) THEN
        use_direction(idim, c_ndims+4) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+4
      ENDIF

      IF (IAND(direction(idim), c_dir_gamma_m1) .NE. 0) THEN
        use_direction(idim, c_ndims+5) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+5
      ENDIF
    ENDDO

    IF (.NOT. use_x) real_space_area = real_space_area * dx
    IF (.NOT. use_y) real_space_area = real_space_area * dy
    IF (.NOT. use_z) real_space_area = real_space_area * dz

    DO idim = 1, c_df_curdims
      calc_mod(idim) = .FALSE.
      IF (p_count(idim) .GT. 1) calc_mod(idim) = .TRUE.
    ENDDO

    ! Calculate range for directions where needed
    IF (calc_ranges) THEN
      DO idim = 1, c_df_curdims
        IF (calc_range(idim)) THEN
          ranges(1,idim) = 1.0e6_num
          ranges(2,idim) = -1.0e6_num
        ENDIF
      ENDDO
      current=>particle_species(species)%attached_list%head
      ind = 0

      DO WHILE(ASSOCIATED(current))
        particle_data(1:c_ndims) = current%part_pos
        particle_data(c_ndims+1:c_ndims+3) = current%part_p
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mass = current%mass
        part_mass_c = part_mass * c
#endif
        gamma_mass_c = SQRT(SUM(current%part_p**2) + part_mass_c**2)
        particle_data(c_ndims+4) = c * (gamma_mass_c - part_mass_c)
        particle_data(c_ndims+5) = gamma_mass_c / part_mass_c - 1.0_num
        DO idim = 1, c_df_curdims
          IF (calc_range(idim)) THEN
            IF (calc_mod(idim)) THEN
              DO idir = 1, c_df_maxdirs
                IF (use_direction(idim, idir)) &
                    current_data = current_data + particle_data(idir)**2
              ENDDO
              current_data = SQRT(current_data)
            ELSE
              current_data = particle_data(l_direction(idim))
            ENDIF
            use_this = .TRUE.
            DO idir = 1, c_df_maxdirs
              IF (use_restrictions(idir) &
                  .AND. (particle_data(idir) .LT. restrictions(1,idir) &
                  .OR. particle_data(idir) .GT. restrictions(2,idir))) &
                      use_this = .FALSE.
            ENDDO
            IF (use_this) THEN
              IF (current_data .LT. ranges(1,idim)) &
                  ranges(1,idim) = current_data
              IF (current_data .GT. ranges(2,idim)) &
                  ranges(2,idim) = current_data
            ENDIF
          ENDIF
        ENDDO
        ind = ind + 1
        current=>current%next
      ENDDO
    ENDIF

    max_p_conv = -10.0_num
    DO idim = 1, c_df_curdims
      IF (.NOT. parallel(idim)) THEN
        ! If not parallel then this is a momentum dimension
        CALL MPI_ALLREDUCE(ranges(1,idim), temp_data, 1, mpireal, MPI_MIN, &
            comm, errcode)
        ranges(1,idim) = temp_data
        CALL MPI_ALLREDUCE(ranges(2,idim), temp_data, 1, mpireal, MPI_MAX, &
            comm, errcode)
        ranges(2,idim) = temp_data
      ENDIF
      ! Fix so that if distribution function is zero then it picks an arbitrary
      ! scale in that direction
      IF (ranges(1,idim) .EQ. ranges(2,idim)) THEN
        ranges(1,idim) = -1.0_num
        ranges(2,idim) = 1.0_num
      ENDIF
      ! Calculate the maximum range of a momentum direction
      IF (ranges(2,idim) - ranges(1,idim) .GT. max_p_conv &
          .AND. .NOT. parallel(idim)) &
              max_p_conv = ranges(2,idim) - ranges(1,idim)
    ENDDO

    DO idim = 1, c_df_curdims
      IF (.NOT. parallel(idim)) conv(idim) = max_p_conv
    ENDDO

    ! Calculate grid spacing
    DO idim = 1, c_df_curdims
      IF (.NOT. parallel(idim)) dgrid(idim) = &
          (ranges(2,idim) - ranges(1,idim)) / REAL(resolution(idim)-1, num)
    ENDDO

    ALLOCATE(data(resolution(1), resolution(2)))
    data = 0.0_num

    current=>particle_species(species)%attached_list%head
    part_weight = weight
    DO WHILE(ASSOCIATED(current))
      particle_data(1:c_ndims) = current%part_pos
      particle_data(c_ndims+1:c_ndims+3) = current%part_p
#ifdef PER_PARTICLE_CHARGE_MASS
      part_mass = current%mass
      part_mass_c = part_mass * c
#endif
      gamma_mass_c = SQRT(SUM(current%part_p**2) + part_mass_c**2)
      particle_data(c_ndims+4) = c * (gamma_mass_c - part_mass_c)
      particle_data(c_ndims+5) = gamma_mass_c / part_mass_c - 1.0_num
      use_this = .TRUE.
      DO idir = 1, c_df_maxdirs
        IF (use_restrictions(idir) &
            .AND. (particle_data(idir) .LT. restrictions(1,idir) &
            .OR. particle_data(idir) .GT. restrictions(2,idir))) &
                use_this = .FALSE.
      ENDDO
      IF (use_this) THEN
        DO idim = 1, c_df_curdims
          IF (calc_mod(idim)) THEN
            DO idir = 1, c_df_maxdirs
              IF (use_direction(idim, idir)) &
                  current_data = current_data + particle_data(idir)**2
            ENDDO
            current_data = SQRT(current_data)
          ELSE
            current_data = particle_data(l_direction(idim))
          ENDIF
          cell(idim) = NINT((current_data - ranges(1,idim)) / dgrid(idim)) + 1
          IF (cell(idim) .LT. 1 .OR. cell(idim) .GT. resolution(idim)) &
              use_this = .FALSE.
        ENDDO
#ifdef PER_PARTICLE_WEIGHT
        part_weight = current%weight
#endif
        IF (use_this) data(cell(1), cell(2)) = &
            data(cell(1), cell(2)) + part_weight ! * real_space_area
      ENDIF
      current=>current%next
    ENDDO

    need_reduce = .TRUE.
    IF (use_x .AND. use_y .AND. use_z) need_reduce = .FALSE.

    IF (need_reduce) THEN
      ! If using x direction need to reduce across y, z
      IF (use_x) color = color + coordinates(c_ndims)
      ! If using y direction need to reduce across x, z
      IF (use_y) color = color + nprocx * coordinates(c_ndims-1)
      ! If using z direction need to reduce across x, y
      IF (use_z) color = color + nprocx * nprocy * coordinates(c_ndims-2)

      CALL MPI_COMM_SPLIT(comm, color, rank, comm_new, errcode)
      ALLOCATE(data2(resolution(1), resolution(2)))
      data2 = 0.0_num
      CALL MPI_ALLREDUCE(data, data2, &
          resolution(1)*resolution(2), mpireal, MPI_SUM, &
          comm_new, errcode)
      data = data2
      DEALLOCATE(data2)
      CALL MPI_COMM_FREE(comm_new, errcode)
    ENDIF

    ! Create grids
    ALLOCATE(grid1(global_resolution(1)))
    ALLOCATE(grid2(global_resolution(2)))

    DO idir = 1, global_resolution(1)
      grid1(idir) = ranges(1,1) + (idir - 1) * dgrid(1)
    ENDDO

    DO idir = 1, global_resolution(2)
      grid2(idir) = ranges(1,2) + (idir - 1) * dgrid(2)
    ENDDO

    grid_name = "Grid_" // TRIM(name) // "_" &
        // TRIM(particle_species(species)%name)
    norm_grid_name = "Norm_Grid_" // TRIM(name) // "_" &
        // TRIM(particle_species(species)%name)
    var_name = TRIM(name) // "_" // TRIM(particle_species(species)%name)

    CALL cfd_write_2d_cartesian_grid(TRIM(grid_name), "Grid", &
        grid1, grid2, 0)
    IF (use_offset_grid) THEN
      IF (parallel(1)) grid1 = grid1 - ranges(1,1)
      IF (parallel(2)) grid2 = grid2 - ranges(1,2)
    ENDIF
    CALL cfd_write_2d_cartesian_grid(TRIM(norm_grid_name), "Grid", &
        grid1/conv(1), grid2/conv(2), 0)

    DEALLOCATE(grid1, grid2)

    new_type = &
        create_2d_array_subtype(resolution, global_resolution, start_local)
    CALL MPI_TYPE_CONTIGUOUS(resolution(1) * resolution(2), &
        mpireal, array_type, errcode)
    CALL MPI_TYPE_COMMIT(array_type, errcode)

    CALL cfd_write_2d_cartesian_variable_parallel(TRIM(var_name), "dist_fn", &
        global_resolution, stagger, TRIM(norm_grid_name), "Grid", data, &
        new_type, array_type)

    CALL MPI_TYPE_FREE(new_type, errcode)
    CALL MPI_TYPE_FREE(array_type, errcode)

    DEALLOCATE(data)

  END SUBROUTINE general_2d_dist_fn



  SUBROUTINE general_3d_dist_fn(name, direction, ranges, resolution, species, &
      restrictions, use_restrictions)

    INTEGER, PARAMETER :: c_df_curdims = 3
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, DIMENSION(c_df_curdims), INTENT(IN) :: direction
    REAL(num), DIMENSION(2,c_df_curdims), INTENT(INOUT) :: ranges
    INTEGER, DIMENSION(c_df_curdims), INTENT(INOUT) :: resolution
    INTEGER, INTENT(IN) :: species
    REAL(num), DIMENSION(2,c_df_maxdirs), INTENT(IN) :: restrictions
    LOGICAL, DIMENSION(c_df_maxdirs), INTENT(IN) :: use_restrictions

    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: data, data2
    REAL(num), DIMENSION(:), ALLOCATABLE :: grid1, grid2, grid3
    LOGICAL, DIMENSION(c_df_curdims) :: parallel
    REAL(num), DIMENSION(c_df_curdims) :: dgrid
    REAL(num) :: current_data, temp_data
    INTEGER :: idim, idir
    LOGICAL, DIMENSION(c_df_curdims) :: calc_range
    LOGICAL :: calc_ranges

    LOGICAL :: use_x, use_y, use_z, need_reduce
    INTEGER, DIMENSION(c_df_curdims) :: start_local, global_resolution
    INTEGER :: color, comm_new
    INTEGER :: new_type, array_type

    LOGICAL, DIMENSION(c_df_curdims, c_df_maxdirs) :: use_direction
    LOGICAL, DIMENSION(c_df_curdims) :: calc_mod
    INTEGER, DIMENSION(c_df_curdims) :: p_count
    INTEGER, DIMENSION(c_df_curdims) :: l_direction
    REAL(num), DIMENSION(c_df_curdims) :: conv
    INTEGER, DIMENSION(c_df_curdims) :: cell
    LOGICAL :: use_this
    REAL(num) :: real_space_area, part_weight
    REAL(num) :: part_mass, part_mass_c, gamma_mass_c

    TYPE(particle), POINTER :: current
    CHARACTER(LEN=string_length) :: grid_name, norm_grid_name, var_name
    REAL(num), DIMENSION(c_df_curdims) :: stagger = 0.0_num
    REAL(num), DIMENSION(c_df_maxdirs) :: particle_data

    REAL(num) :: max_p_conv

    INTEGER :: ind

    use_x = .FALSE.
    use_y = .FALSE.
    color = 0
    global_resolution = resolution
    parallel = .FALSE.
    start_local = 1
    calc_range = .FALSE.
    calc_ranges = .FALSE.
    p_count = 0
    use_direction = .FALSE.
    l_direction = 0

    real_space_area = 1.0_num
    current_data = 0.0_num
    part_mass = particle_species(species)%mass
    part_mass_c = part_mass * c

    DO idim = 1, c_df_curdims
      IF (direction(idim) .EQ. c_dir_x) THEN
        use_x = .TRUE.
        resolution(idim) = nx
        ranges(1,idim) = x_min_local
        ranges(2,idim) = x_max_local
        start_local(idim) = cell_x_min(coordinates(c_ndims)+1)
        global_resolution(idim) = nx_global
        dgrid(idim) = dx
        parallel(idim) = .TRUE.
        use_direction(idim,1) = .TRUE.
        l_direction(idim) = 1
        conv(idim) = MAX(length_x, length_y, length_z)
        CYCLE

      ELSE IF (direction(idim) .EQ. c_dir_y) THEN
        use_y = .TRUE.
        resolution(idim) = ny
        ranges(1,idim) = y_min_local
        ranges(2,idim) = y_max_local
        start_local(idim) = cell_y_min(coordinates(c_ndims-1)+1)
        global_resolution(idim) = ny_global
        dgrid(idim) = dy
        parallel(idim) = .TRUE.
        use_direction(idim,2) = .TRUE.
        l_direction(idim) = 2
        conv(idim) = MAX(length_x, length_y, length_z)
        CYCLE

      ELSE IF (direction(idim) .EQ. c_dir_z) THEN
        use_z = .TRUE.
        resolution(idim) = nz
        ranges(1,idim) = z_min_local
        ranges(2,idim) = z_max_local
        start_local(idim) = cell_z_min(coordinates(c_ndims-2)+1)
        global_resolution(idim) = nz_global
        dgrid(idim) = dz
        parallel(idim) = .TRUE.
        use_direction(idim,3) = .TRUE.
        l_direction(idim) = 3
        conv(idim) = MAX(length_x, length_y, length_z)
        CYCLE

      ENDIF

      conv(idim) = c*m0
      ! If we're here then this must be a momentum space direction
      ! So determine which momentum space directions are needed
      IF (ranges(1,idim) .EQ. ranges(2,idim)) THEN
        calc_range(idim) = .TRUE.
        calc_ranges = .TRUE.
      ENDIF

      IF (IAND(direction(idim), c_dir_px) .NE. 0) THEN
        use_direction(idim, c_ndims+1) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+1
      ENDIF

      IF (IAND(direction(idim), c_dir_py) .NE. 0) THEN
        use_direction(idim, c_ndims+2) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+2
      ENDIF

      IF (IAND(direction(idim), c_dir_pz) .NE. 0) THEN
        use_direction(idim, c_ndims+3) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+3
      ENDIF

      IF (IAND(direction(idim), c_dir_en) .NE. 0) THEN
        use_direction(idim, c_ndims+4) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+4
      ENDIF

      IF (IAND(direction(idim), c_dir_gamma_m1) .NE. 0) THEN
        use_direction(idim, c_ndims+5) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = c_ndims+5
      ENDIF
    ENDDO

    IF (.NOT. use_x) real_space_area = real_space_area * dx
    IF (.NOT. use_y) real_space_area = real_space_area * dy
    IF (.NOT. use_z) real_space_area = real_space_area * dz

    DO idim = 1, c_df_curdims
      calc_mod(idim) = .FALSE.
      IF (p_count(idim) .GT. 1) calc_mod(idim) = .TRUE.
    ENDDO

    ! Calculate range for directions where needed
    IF (calc_ranges) THEN
      DO idim = 1, c_df_curdims
        IF (calc_range(idim)) THEN
          ranges(1,idim) = 1.0e6_num
          ranges(2,idim) = -1.0e6_num
        ENDIF
      ENDDO
      current=>particle_species(species)%attached_list%head
      ind = 0

      DO WHILE(ASSOCIATED(current))
        particle_data(1:c_ndims) = current%part_pos
        particle_data(c_ndims+1:c_ndims+3) = current%part_p
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mass = current%mass
        part_mass_c = part_mass * c
#endif
        gamma_mass_c = SQRT(SUM(current%part_p**2) + part_mass_c**2)
        particle_data(c_ndims+4) = c * (gamma_mass_c - part_mass_c)
        particle_data(c_ndims+5) = gamma_mass_c / part_mass_c - 1.0_num
        DO idim = 1, c_df_curdims
          IF (calc_range(idim)) THEN
            IF (calc_mod(idim)) THEN
              DO idir = 1, c_df_maxdirs
                IF (use_direction(idim, idir)) &
                    current_data = current_data + particle_data(idir)**2
              ENDDO
              current_data = SQRT(current_data)
            ELSE
              current_data = particle_data(l_direction(idim))
            ENDIF
            use_this = .TRUE.
            DO idir = 1, c_df_maxdirs
              IF (use_restrictions(idir) &
                  .AND. (particle_data(idir) .LT. restrictions(1,idir) &
                  .OR. particle_data(idir) .GT. restrictions(2,idir))) &
                      use_this = .FALSE.
            ENDDO
            IF (use_this) THEN
              IF (current_data .LT. ranges(1,idim)) &
                  ranges(1,idim) = current_data
              IF (current_data .GT. ranges(2,idim)) &
                  ranges(2,idim) = current_data
            ENDIF
          ENDIF
        ENDDO
        ind = ind + 1
        current=>current%next
      ENDDO
    ENDIF

    max_p_conv = -10.0_num
    DO idim = 1, c_df_curdims
      IF (.NOT. parallel(idim)) THEN
        ! If not parallel then this is a momentum dimension
        CALL MPI_ALLREDUCE(ranges(1,idim), temp_data, 1, mpireal, MPI_MIN, &
            comm, errcode)
        ranges(1,idim) = temp_data
        CALL MPI_ALLREDUCE(ranges(2,idim), temp_data, 1, mpireal, MPI_MAX, &
            comm, errcode)
        ranges(2,idim) = temp_data
      ENDIF
      ! Fix so that if distribution function is zero then it picks an arbitrary
      ! scale in that direction
      IF (ranges(1,idim) .EQ. ranges(2,idim)) THEN
        ranges(1,idim) = -1.0_num
        ranges(2,idim) = 1.0_num
      ENDIF
      ! Calculate the maximum range of a momentum direction
      IF (ranges(2,idim) - ranges(1,idim) .GT. max_p_conv &
          .AND. .NOT. parallel(idim)) &
              max_p_conv = ranges(2,idim) - ranges(1,idim)
    ENDDO

    DO idim = 1, c_df_curdims
      IF (.NOT. parallel(idim)) conv(idim) = max_p_conv
    ENDDO

    ! Calculate grid spacing
    DO idim = 1, c_df_curdims
      IF (.NOT. parallel(idim)) dgrid(idim) = &
          (ranges(2,idim) - ranges(1,idim)) / REAL(resolution(idim)-1, num)
    ENDDO

    ALLOCATE(data(resolution(1), resolution(2), resolution(3)))
    data = 0.0_num

    current=>particle_species(species)%attached_list%head
    part_weight = weight
    DO WHILE(ASSOCIATED(current))
      particle_data(1:c_ndims) = current%part_pos
      particle_data(c_ndims+1:c_ndims+3) = current%part_p
#ifdef PER_PARTICLE_CHARGE_MASS
      part_mass = current%mass
      part_mass_c = part_mass * c
#endif
      gamma_mass_c = SQRT(SUM(current%part_p**2) + part_mass_c**2)
      particle_data(c_ndims+4) = c * (gamma_mass_c - part_mass_c)
      particle_data(c_ndims+5) = gamma_mass_c / part_mass_c - 1.0_num
      use_this = .TRUE.
      DO idir = 1, c_df_maxdirs
        IF (use_restrictions(idir) &
            .AND. (particle_data(idir) .LT. restrictions(1,idir) &
            .OR. particle_data(idir) .GT. restrictions(2,idir))) &
                use_this = .FALSE.
      ENDDO
      IF (use_this) THEN
        DO idim = 1, c_df_curdims
          IF (calc_mod(idim)) THEN
            DO idir = 1, c_df_maxdirs
              IF (use_direction(idim, idir)) &
                  current_data = current_data + particle_data(idir)**2
            ENDDO
            current_data = SQRT(current_data)
          ELSE
            current_data = particle_data(l_direction(idim))
          ENDIF
          cell(idim) = NINT((current_data - ranges(1,idim)) / dgrid(idim)) + 1
          IF (cell(idim) .LT. 1 .OR. cell(idim) .GT. resolution(idim)) &
              use_this = .FALSE.
        ENDDO
#ifdef PER_PARTICLE_WEIGHT
        part_weight = current%weight
#endif
        IF (use_this) data(cell(1), cell(2), cell(3)) = &
            data(cell(1), cell(2), cell(3)) + part_weight ! * real_space_area
      ENDIF
      current=>current%next
    ENDDO

    need_reduce = .TRUE.
    IF (use_x .AND. use_y .AND. use_z) need_reduce = .FALSE.

    IF (need_reduce) THEN
      ! If using x direction need to reduce across y, z
      IF (use_x) color = color + coordinates(c_ndims)
      ! If using y direction need to reduce across x, z
      IF (use_y) color = color + nprocx * coordinates(c_ndims-1)
      ! If using z direction need to reduce across x, y
      IF (use_z) color = color + nprocx * nprocy * coordinates(c_ndims-2)

      CALL MPI_COMM_SPLIT(comm, color, rank, comm_new, errcode)
      ALLOCATE(data2(resolution(1), resolution(2), resolution(3)))
      data2 = 0.0_num
      CALL MPI_ALLREDUCE(data, data2, &
          resolution(1)*resolution(2)*resolution(3), mpireal, MPI_SUM, &
          comm_new, errcode)
      data = data2
      DEALLOCATE(data2)
      CALL MPI_COMM_FREE(comm_new, errcode)
    ENDIF

    ! Create grids
    ALLOCATE(grid1(global_resolution(1)))
    ALLOCATE(grid2(global_resolution(2)))
    ALLOCATE(grid3(global_resolution(3)))

    DO idir = 1, global_resolution(1)
      grid1(idir) = ranges(1,1) + (idir - 1) * dgrid(1)
    ENDDO

    DO idir = 1, global_resolution(2)
      grid2(idir) = ranges(1,2) + (idir - 1) * dgrid(2)
    ENDDO

    DO idir = 1, global_resolution(3)
      grid3(idir) = ranges(1,3) + (idir - 1) * dgrid(3)
    ENDDO

    grid_name = "Grid_" // TRIM(name) // "_" &
        // TRIM(particle_species(species)%name)
    norm_grid_name = "Norm_Grid_" // TRIM(name) // "_" &
        // TRIM(particle_species(species)%name)
    var_name = TRIM(name) // "_" // TRIM(particle_species(species)%name)

    CALL cfd_write_3d_cartesian_grid(TRIM(grid_name), "Grid", &
        grid1, grid2, grid3, 0)
    IF (use_offset_grid) THEN
      IF (parallel(1)) grid1 = grid1 - ranges(1,1)
      IF (parallel(2)) grid2 = grid2 - ranges(1,2)
      IF (parallel(3)) grid3 = grid3 - ranges(1,3)
    ENDIF
    CALL cfd_write_3d_cartesian_grid(TRIM(norm_grid_name), "Grid", &
        grid1/conv(1), grid2/conv(2), grid3/conv(3), 0)

    DEALLOCATE(grid1, grid2, grid3)

    new_type = &
        create_3d_array_subtype(resolution, global_resolution, start_local)
    CALL MPI_TYPE_CONTIGUOUS(resolution(1) * resolution(2) * resolution(3), &
        mpireal, array_type, errcode)
    CALL MPI_TYPE_COMMIT(array_type, errcode)

    CALL cfd_write_3d_cartesian_variable_parallel(TRIM(var_name), "dist_fn", &
        global_resolution, stagger, TRIM(norm_grid_name), "Grid", data, &
        new_type, array_type)

    CALL MPI_TYPE_FREE(new_type, errcode)
    CALL MPI_TYPE_FREE(array_type, errcode)

    DEALLOCATE(data)

  END SUBROUTINE general_3d_dist_fn

END MODULE dist_fn
