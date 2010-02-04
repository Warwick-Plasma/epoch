MODULE dist_fn

  USE shared_data
  USE partlist
  USE output_cartesian
  USE output
  USE iocontrol
  USE mpi_subtype_control
  IMPLICIT NONE

CONTAINS

  SUBROUTINE attach_dist_fn(BLOCK)

    TYPE(distribution_function_block), POINTER :: BLOCK
    TYPE(distribution_function_block), POINTER :: current

    current=>dist_fns
    IF (.NOT. ASSOCIATED(current)) THEN
      ! This is the first distribution function to add
      dist_fns=>BLOCK
      RETURN
    ENDIF
    DO WHILE(ASSOCIATED(current%next))
      current=>current%next
    ENDDO
    current%next=>BLOCK

  END SUBROUTINE attach_dist_fn



  SUBROUTINE setup_dist_fn(BLOCK)

    TYPE(distribution_function_block), POINTER :: BLOCK

    NULLIFY(block%next)
    ALLOCATE(block%use_species(1:n_species))
    block%use_species = .FALSE.
    block%ranges = 1.0_num
    block%use_restrictions = .FALSE.
    block%ndims = -1

  END SUBROUTINE setup_dist_fn



  SUBROUTINE clean_dist_fns

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
    REAL(num), DIMENSION(3, 2) :: ranges
    LOGICAL, DIMENSION(6) :: use_restrictions
    REAL(num), DIMENSION(6, 2) :: restrictions
    INTEGER, DIMENSION(3) :: resolution = (/100, 100, 100/)
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

          IF (current%ndims .EQ. 2) THEN
            CALL general_2d_dist_fn(current%name, current%directions(1:2), &
                ranges(1:2,:), resolution(1:2), ispecies, restrictions, &
                use_restrictions)
          ELSE
            CALL general_3d_dist_fn(current%name, current%directions, &
                ranges, resolution, ispecies, restrictions, use_restrictions)
          ENDIF
          IF (current%store_ranges) current%ranges = ranges
        ENDDO
      ENDIF
      current=>current%next
    ENDDO

  END SUBROUTINE write_dist_fns



  SUBROUTINE general_3d_dist_fn(name, direction, range, resolution, species, &
      restrictions, use_restrictions)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, DIMENSION(3), INTENT(IN) :: direction
    REAL(num), DIMENSION(3, 2), INTENT(INOUT) :: range
    INTEGER, DIMENSION(3), INTENT(INOUT) :: resolution
    INTEGER, INTENT(IN) :: species
    REAL(num), DIMENSION(6, 2), INTENT(IN) :: restrictions
    LOGICAL, DIMENSION(6), INTENT(IN) :: use_restrictions

    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: DATA, data2
    REAL(num), DIMENSION(:), ALLOCATABLE :: grid1, grid2, grid3
    LOGICAL, DIMENSION(3) :: parallel
    REAL(num), DIMENSION(3) :: dgrid
    REAL(num) :: current_data, temp_data
    INTEGER :: idim, idir
    LOGICAL, DIMENSION(3) :: calc_range
    LOGICAL :: calc_ranges

    LOGICAL :: use_x, use_y, use_z, need_reduce
    INTEGER, DIMENSION(3) :: start_local, global_resolution
    INTEGER :: color
    INTEGER :: comm_new, new_type

    LOGICAL, DIMENSION(3, 6) :: use_direction
    LOGICAL, DIMENSION(3) :: calc_mod
    INTEGER, DIMENSION(3) :: p_count
    INTEGER, DIMENSION(3) :: l_direction
    REAL(num), DIMENSION(3) :: conv
    INTEGER, DIMENSION(3) :: cell
    LOGICAL :: use_this

    TYPE(particle), POINTER :: current
    CHARACTER(LEN=string_length) :: grid_name, norm_grid_name, var_name
    REAL(num), DIMENSION(3) :: stagger = 0.0_num
    REAL(num), DIMENSION(6) :: particle_data

    REAL(num) :: max_p_conv

    INTEGER :: ind

    use_x = .FALSE.
    use_y = .FALSE.
    need_reduce = .TRUE.
    color = 0
    global_resolution = resolution
    parallel = .FALSE.
    start_local = 1
    calc_range = .FALSE.
    calc_ranges = .FALSE.
    p_count = 0
    use_direction = .FALSE.
    l_direction = 0

    DO idim = 1, 3
      IF (direction(idim) .EQ. c_dir_x) THEN
        use_x = .TRUE.
        resolution(idim) = nx
        RANGE(idim, 1) = x_start_local
        RANGE(idim, 2) = x_end_local
        start_local(idim) = cell_x_start(coordinates(3)+1)
        global_resolution(idim) = nx_global
        dgrid(idim) = dx
        parallel(idim) = .TRUE.
        use_direction(idim, 1) = .TRUE.
        l_direction(idim) = 1
        conv(idim) = MAX(length_x, length_y, length_z)
        CYCLE
      ENDIF
      IF (direction(idim) .EQ. c_dir_y)  THEN
        use_y = .TRUE.
        resolution(idim) = ny
        RANGE(idim, 1) = y_start_local
        RANGE(idim, 2) = y_end_local
        start_local(idim) = cell_y_start(coordinates(2)+1)
        global_resolution(idim) = ny_global
        dgrid(idim) = dy
        parallel(idim) = .TRUE.
        use_direction(idim, 2) = .TRUE.
        l_direction(idim) = 2
        conv(idim) = MAX(length_x, length_y, length_z)
        CYCLE
      ENDIF
      IF (direction(idim) .EQ. c_dir_z)  THEN
        use_z = .TRUE.
        resolution(idim) = ny
        RANGE(idim, 1) = z_start_local
        RANGE(idim, 2) = z_end_local
        start_local(idim) = cell_z_start(coordinates(1)+1)
        global_resolution(idim) = nz_global
        dgrid(idim) = dz
        parallel(idim) = .TRUE.
        use_direction(idim, 3) = .TRUE.
        l_direction(idim) = 3
        conv(idim) = MAX(length_x, length_y, length_z)
        CYCLE
      ENDIF
      conv(idim) = c*m0
      ! If we're here then this must be a momentum space direction
      ! So determine which momentum space directions are needed
      IF (RANGE(idim, 1) .EQ. RANGE(idim, 2)) THEN
        calc_range(idim) = .TRUE.
        calc_ranges = .TRUE.
      ENDIF
      IF (IAND(direction(idim), c_dir_px) .NE. 0) THEN
        use_direction(idim, 4) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = 4
      ENDIF
      IF (IAND(direction(idim), c_dir_py) .NE. 0) THEN
        use_direction(idim, 5) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = 5
      ENDIF
      IF (IAND(direction(idim), c_dir_pz) .NE. 0) THEN
        use_direction(idim, 6) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = 6
      ENDIF
    ENDDO

    DO idim = 1, 3
      IF (p_count(idim) .GT. 1) calc_mod(idim) = .TRUE.
    ENDDO
    ! Calculate ranges for directions where needed
    IF (calc_ranges) THEN
      DO idim = 1, 3
        IF (calc_range(idim)) THEN
          RANGE(idim, 1) = 1.0e6_num
          RANGE(idim, 2) = -1.0e6_num
        ENDIF
      ENDDO
      current=>particle_species(species)%attached_list%head
      ind = 0
      DO WHILE(ASSOCIATED(current))
        particle_data(1:3) = current%part_pos
        particle_data(4:6) = current%part_p
        DO idim = 1, 3
          IF (calc_range(idim)) THEN
            IF (calc_mod(idim)) THEN
              DO idir = 1, 6
                IF (use_direction(idim, idir)) &
                    current_data = current_data+particle_data(idir)
              ENDDO
              current_data = SQRT(current_data)
            ELSE
              current_data = particle_data(l_direction(idim))
            ENDIF
            use_this = .TRUE.
            DO idir = 1, 6
              IF (use_restrictions(idir) .AND. &
                  (particle_data(idir) .LT. restrictions(idir, 1) .OR. &
                  particle_data(idir) .GT. restrictions(idir, 2))) &
                      use_this = .FALSE.
            ENDDO
            IF (use_this) THEN
              IF (current_data .LT. RANGE(idim, 1)) &
                  RANGE(idim, 1) = current_data
              IF (current_data .GT. RANGE(idim, 2)) &
                  RANGE(idim, 2) = current_data
            ENDIF
          ENDIF
        ENDDO
        ind = ind+1
        current=>current%next
      ENDDO
    ENDIF

    max_p_conv = -10.0_num
    DO idim = 1, 3
      IF (.NOT. parallel(idim)) THEN
        ! If not parallel then this is a momentum DIMENSION
        CALL MPI_ALLREDUCE(RANGE(idim, 1), temp_data, 1, mpireal, MPI_MIN, &
            comm, errcode)
        RANGE(idim, 1) = temp_data
        CALL MPI_ALLREDUCE(RANGE(idim, 2), temp_data, 1, mpireal, MPI_MAX, &
            comm, errcode)
        RANGE(idim, 2) = temp_data
      ENDIF
      ! Fix so that if distribution function is zero then it picks an arbitrary
      ! scale in that direction
      IF (RANGE(idim, 1) .EQ. RANGE(idim, 2)) THEN
        RANGE(idim, 1) = -1.0_num
        RANGE(idim, 2) = 1.0_num
      ENDIF
      ! Calculate the maxmium range of a momentum direction
      IF (RANGE(idim, 2)-RANGE(idim, 1) .GT. max_p_conv) &
          max_p_conv = RANGE(idim, 2)-RANGE(idim, 1)
    ENDDO

    DO idim = 1, 2
      IF (.NOT. parallel(idim)) conv(idim) = max_p_conv
    ENDDO

    ! Setup MPI
    IF (use_x .AND. use_y .AND. use_z) need_reduce = .FALSE.
    color = 1 ! coordinates(3) + nprocx * coordinates(2) + &
              ! nprocx*nprocy*coordinates(1)
    IF (use_z) color = color+nprocx*nprocy*coordinates(1)
    ! If using x direction need to reduce only across all y
    IF (use_y) color = color+nprocx * coordinates(2)
    ! If using y direction need to reduce only across all x
    IF (use_x) color = color+coordinates(3)

    IF (need_reduce) THEN
      CALL MPI_COMM_SPLIT(comm, color, rank, comm_new, errcode)
    ELSE
      comm_new = MPI_COMM_NULL
    ENDIF

    new_type = create_3d_array_subtype(resolution, global_resolution, &
        start_local)
    ! Create grids
    DO idim = 1, 3
      IF (.NOT. parallel(idim)) &
          dgrid(idim) = (RANGE(idim, 2) - &
              RANGE(idim, 1))/REAL(resolution(idim)-1, num)
    ENDDO
    ALLOCATE(grid1(0:global_resolution(1)), grid2(0:global_resolution(2)))
    ALLOCATE(grid3(0:global_resolution(3)))

    grid1(0) = -dgrid(1) + RANGE(1, 1)
    DO idir = 1, global_resolution(1)
      grid1(idir) = grid1(idir-1)+dgrid(1)
    ENDDO

    grid2(0) = -dgrid(2) + RANGE(2, 1)
    DO idir = 1, global_resolution(2)
      grid2(idir) = grid2(idir-1)+dgrid(2)
    ENDDO

    grid3(0) = -dgrid(3) + RANGE(3, 1)
    DO idir = 1, global_resolution(3)
      grid3(idir) = grid3(idir-1)+dgrid(3)
    ENDDO

    ALLOCATE(DATA(1:resolution(1), 1:resolution(2), 1:resolution(3)))
    DATA = 0.0_num

    current=>particle_species(species)%attached_list%head
    DO WHILE(ASSOCIATED(current))
      particle_data(1:3) = current%part_pos
      particle_data(4:6) = current%part_p
      use_this = .TRUE.
      DO idir = 1, 6
        IF (use_restrictions(idir) .AND. &
            (particle_data(idir) .LT. restrictions(idir, 1) .OR. &
            particle_data(idir) .GT. restrictions(idir, 2))) use_this = .FALSE.
      ENDDO
      IF (use_this) THEN
        DO idim = 1, 3
          IF (calc_mod(idim)) THEN
            DO idir = 1, 6
              IF (use_direction(idim, idir)) &
                  current_data = current_data+particle_data(idir)
            ENDDO
            current_data = SQRT(current_data)
          ELSE
            current_data = particle_data(l_direction(idim))
          ENDIF
          cell(idim) = NINT((current_data-RANGE(idim, 1))/dgrid(idim))+1
          IF (cell(idim) .LT. 1 .OR. cell(idim) .GT. resolution(idim)) &
              use_this = .FALSE.
        ENDDO
        IF (use_this) &
            DATA(cell(1), cell(2), cell(3)) = &
                DATA(cell(1), cell(2), cell(3))+current%weight
      ENDIF
      current=>current%next
    ENDDO

    IF (need_reduce) THEN
      ALLOCATE(data2(1:resolution(1), 1:resolution(2), 1:resolution(3)))
      data2 = 0.0_num
      CALL MPI_ALLREDUCE(DATA, data2, &
          resolution(1)*resolution(2)*resolution(3), mpireal, &
          MPI_SUM, comm_new, errcode)
      DATA = data2
      DEALLOCATE(data2)
    ENDIF

    grid_name = "Grid_" // TRIM(name) // "_" // &
        TRIM(particle_species(species)%name)
    norm_grid_name = "Norm_Grid_" // TRIM(name) // "_" // &
        TRIM(particle_species(species)%name)
    var_name = TRIM(name) // "_" // TRIM(particle_species(species)%name)

    CALL cfd_write_3d_cartesian_grid(TRIM(grid_name), "Grid", &
        grid1(1:global_resolution(1)), grid2(1:global_resolution(2)), &
        grid3(1:global_resolution(3)), 0)
    CALL cfd_write_3d_cartesian_grid(TRIM(norm_grid_name), "Grid", &
        grid1(1:global_resolution(1))/conv(1), &
        grid2(1:global_resolution(2))/conv(2), &
        grid3(1:global_resolution(3))/conv(3), 0)

    CALL cfd_write_3d_cartesian_variable_parallel(TRIM(var_name), "dist_fn", &
        global_resolution, stagger, TRIM(norm_grid_name), "Grid", &
        DATA, new_type)
    CALL MPI_TYPE_FREE(new_type, errcode)
    IF (need_reduce) CALL MPI_COMM_FREE(comm_new, errcode)

    DEALLOCATE(DATA)
    DEALLOCATE(grid1, grid2, grid3)

  END SUBROUTINE general_3d_dist_fn



  SUBROUTINE general_2d_dist_fn(name, direction, range, resolution, species, &
      restrictions, use_restrictions)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, DIMENSION(2), INTENT(IN) :: direction
    REAL(num), DIMENSION(2, 2), INTENT(INOUT) :: range
    INTEGER, DIMENSION(2), INTENT(INOUT) :: resolution
    INTEGER, INTENT(IN) :: species
    REAL(num), DIMENSION(6, 2), INTENT(IN) :: restrictions
    LOGICAL, DIMENSION(6), INTENT(IN) :: use_restrictions

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: DATA, data2
    REAL(num), DIMENSION(:), ALLOCATABLE :: grid1, grid2
    LOGICAL, DIMENSION(2) :: parallel
    REAL(num), DIMENSION(2) :: dgrid
    REAL(num) :: current_data, temp_data
    INTEGER :: idim, idir
    LOGICAL, DIMENSION(2) :: calc_range
    LOGICAL :: calc_ranges

    LOGICAL :: use_x, use_y, use_z, need_reduce
    INTEGER, DIMENSION(2) :: start_local, global_resolution
    INTEGER :: color
    INTEGER :: comm_new, new_type

    LOGICAL, DIMENSION(2, 6) :: use_direction
    LOGICAL, DIMENSION(2) :: calc_mod
    INTEGER, DIMENSION(2) :: p_count
    INTEGER, DIMENSION(2) :: l_direction
    REAL(num), DIMENSION(2) :: conv
    INTEGER, DIMENSION(2) :: cell
    LOGICAL :: use_this

    TYPE(particle), POINTER :: current
    CHARACTER(LEN=string_length) :: grid_name, norm_grid_name, var_name
    REAL(num), DIMENSION(2) :: stagger = 0.0_num
    REAL(num), DIMENSION(6) :: particle_data
    REAL(num) :: max_p_conv

    INTEGER :: ind

    use_x = .FALSE.
    use_y = .FALSE.
    use_z = .FALSE.

    need_reduce = .TRUE.
    color = 0
    global_resolution = resolution
    parallel = .FALSE.
    start_local = 1
    calc_range = .FALSE.
    calc_ranges = .FALSE.
    p_count = 0
    use_direction = .FALSE.
    l_direction = 0

    DO idim = 1, 2

      IF (direction(idim) .EQ. c_dir_x) THEN
        use_x = .TRUE.
        resolution(idim) = nx
        RANGE(idim, 1) = x_start_local
        RANGE(idim, 2) = x_end_local
        ! IF (use_offset_grid) range(idim,:) = range(idim,:)-x_start
        start_local(idim) = cell_x_start(coordinates(2)+1)
        global_resolution(idim) = nx_global
        dgrid(idim) = dx
        parallel(idim) = .TRUE.
        use_direction(idim, 1) = .TRUE.
        l_direction(idim) = 1
        conv(idim) = MAX(length_x, length_y)
        CYCLE
      ENDIF
      IF (direction(idim) .EQ. c_dir_y)  THEN
        use_y = .TRUE.
        resolution(idim) = ny
        RANGE(idim, 1) = y_start_local
        RANGE(idim, 2) = y_end_local
        ! IF (use_offset_grid) range(idim,:) = range(idim,:)-y_start
        start_local(idim) = cell_y_start(coordinates(1)+1)
        global_resolution(idim) = ny_global
        dgrid(idim) = dy
        parallel(idim) = .TRUE.
        use_direction(idim, 2) = .TRUE.
        l_direction(idim) = 2
        conv(idim) = MAX(length_x, length_y)
        CYCLE
      ENDIF
      IF (direction(idim) .EQ. c_dir_z)  THEN
        use_z = .TRUE.
        resolution(idim) = ny
        RANGE(idim, 1) = z_start_local
        RANGE(idim, 2) = z_end_local
        start_local(idim) = cell_z_start(coordinates(1)+1)
        global_resolution(idim) = nz_global
        dgrid(idim) = dz
        parallel(idim) = .TRUE.
        use_direction(idim, 3) = .TRUE.
        l_direction(idim) = 3
        conv(idim) = MAX(length_x, length_y, length_z)
        CYCLE
      ENDIF
      conv(idim) = c*m0
      ! If we're here then this must be a momentum space direction
      ! So determine which momentum space directions are needed
      IF (RANGE(idim, 1) .EQ. RANGE(idim, 2)) THEN
        calc_range(idim) = .TRUE.
        calc_ranges = .TRUE.
      ENDIF
      IF (IAND(direction(idim), c_dir_px) .NE. 0) THEN
        use_direction(idim, 4) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = 4
      ENDIF
      IF (IAND(direction(idim), c_dir_py) .NE. 0) THEN
        use_direction(idim, 5) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = 5
      ENDIF
      IF (IAND(direction(idim), c_dir_pz) .NE. 0) THEN
        use_direction(idim, 6) = .TRUE.
        p_count(idim) = p_count(idim)+1
        l_direction(idim) = 6
      ENDIF
    ENDDO

    DO idim = 1, 2
      IF (p_count(idim) .GT. 1) calc_mod(idim) = .TRUE.
    ENDDO
    ! Calculate ranges for directions where needed
    IF (calc_ranges) THEN
      DO idim = 1, 2
        IF (calc_range(idim)) THEN
          RANGE(idim, 1) = 1.0e6_num
          RANGE(idim, 2) = -1.0e6_num
        ENDIF
      ENDDO
      current=>particle_species(species)%attached_list%head
      ind = 0
      DO WHILE(ASSOCIATED(current))
        particle_data(1:3) = current%part_pos
        particle_data(4:6) = current%part_p
        DO idim = 1, 2
          IF (calc_range(idim)) THEN
            IF (calc_mod(idim)) THEN
              DO idir = 1, 6
                IF (use_direction(idim, idir)) &
                    current_data = current_data+particle_data(idir)
              ENDDO
              current_data = SQRT(current_data)
            ELSE
              current_data = particle_data(l_direction(idim))
            ENDIF
            use_this = .TRUE.
            DO idir = 1, 6
              IF (use_restrictions(idir) .AND. &
                  (particle_data(idir) .LT. restrictions(idir, 1) .OR. &
                  particle_data(idir) .GT. restrictions(idir, 2))) &
                      use_this = .FALSE.
            ENDDO
            IF (use_this) THEN
              IF (current_data .LT. RANGE(idim, 1)) &
                  RANGE(idim, 1) = current_data
              IF (current_data .GT. RANGE(idim, 2)) &
                  RANGE(idim, 2) = current_data
            ENDIF
          ENDIF
        ENDDO
        ind = ind+1
        current=>current%next
      ENDDO
    ENDIF

    max_p_conv = -10.0_num
    DO idim = 1, 2
      IF (.NOT. parallel(idim)) THEN
        ! If not parallel then this is a momentum DIMENSION
        CALL MPI_ALLREDUCE(RANGE(idim, 1), temp_data, 1, mpireal, &
            MPI_MIN, comm, errcode)
        RANGE(idim, 1) = temp_data
        CALL MPI_ALLREDUCE(RANGE(idim, 2), temp_data, 1, mpireal, &
            MPI_MAX, comm, errcode)
        RANGE(idim, 2) = temp_data

        ! Calculate the maxmium range of a momentum direction
        IF (RANGE(idim, 2)-RANGE(idim, 1) .GT. max_p_conv) &
            max_p_conv = RANGE(idim, 2)-RANGE(idim, 1)
      ENDIF
      ! Fix so that if distribution function is zero then it picks an arbitrary
      ! scale in that direction
      IF (RANGE(idim, 1) .EQ. RANGE(idim, 2)) THEN
        RANGE(idim, 1) = -1.0_num
        RANGE(idim, 2) = 1.0_num
      ENDIF
    ENDDO

    DO idim = 1, 2
      IF (.NOT. parallel(idim)) conv(idim) = max_p_conv
    ENDDO

    ! Setup MPI
    IF (use_x .AND. use_y .AND. use_z) need_reduce = .FALSE.
    color = 1 ! coordinates(3) + nprocx * coordinates(2) + &
              ! nprocx*nprocy*coordinates(1)
    IF (use_z) color = color+nprocx*nprocy*coordinates(1)
    ! If using x direction need to reduce only across all y
    IF (use_y) color = color+nprocx * coordinates(2)
    ! If using y direction need to reduce only across all x
    IF (use_x) color = color+coordinates(3)

!!$    ! If using x direction need to reduce only across all y
!!$    IF (.NOT. use_y) color = color+coordinates(2)
!!$    ! If using y direction need to reduce only across all x
!!$    IF (.NOT. use_x) color = color+nprocx*coordinates(1)

    IF (need_reduce) THEN
      CALL MPI_COMM_SPLIT(comm, color, rank, comm_new, errcode)
    ELSE
      comm_new = MPI_COMM_NULL
    ENDIF

    new_type = create_2d_array_subtype(resolution, global_resolution, &
        start_local)
    ! Create grids
    DO idim = 1, 2
      IF (.NOT. parallel(idim)) &
          dgrid(idim) = (RANGE(idim, 2)-RANGE(idim, 1)) / &
              REAL(resolution(idim)-1, num)
    ENDDO
    ALLOCATE(grid1(0:global_resolution(1)), grid2(0:global_resolution(2)))

    grid1(0) = -dgrid(1) + RANGE(1, 1)
    DO idir = 1, global_resolution(1)
      grid1(idir) = grid1(idir-1)+dgrid(1)
    ENDDO

    grid2(0) = -dgrid(2) + RANGE(2, 1)
    DO idir = 1, global_resolution(2)
      grid2(idir) = grid2(idir-1)+dgrid(2)
    ENDDO

    ALLOCATE(DATA(1:resolution(1), 1:resolution(2)))
    DATA = 0.0_num

    current=>particle_species(species)%attached_list%head
    DO WHILE(ASSOCIATED(current))
      particle_data(1:3) = current%part_pos
      particle_data(4:6) = current%part_p
      use_this = .TRUE.
      DO idir = 1, 6
        IF (use_restrictions(idir) .AND. &
            (particle_data(idir) .LT. restrictions(idir, 1) .OR. &
            particle_data(idir) .GT. restrictions(idir, 2))) use_this = .FALSE.
      ENDDO
      IF (use_this) THEN
        DO idim = 1, 2
          IF (calc_mod(idim)) THEN
            DO idir = 1, 5
              IF (use_direction(idim, idir)) &
                  current_data = current_data+particle_data(idir)
            ENDDO
            current_data = SQRT(current_data)
          ELSE
            current_data = particle_data(l_direction(idim))
          ENDIF
          cell(idim) = NINT((current_data-RANGE(idim, 1))/dgrid(idim))+1
          IF (cell(idim) .LT. 1 .OR. cell(idim) .GT. resolution(idim)) &
              use_this = .FALSE.
        ENDDO

        IF (use_this) &
            DATA(cell(1), cell(2)) = DATA(cell(1), cell(2))+current%weight
      ENDIF
      current=>current%next
    ENDDO

    IF (need_reduce) THEN
      ALLOCATE(data2(1:resolution(1), 1:resolution(2)))
      data2 = 0.0_num
      CALL MPI_ALLREDUCE(DATA, data2, resolution(1)*resolution(2), &
          mpireal, MPI_SUM, comm_new, errcode)
      DATA = data2
      DEALLOCATE(data2)
    ENDIF

    grid_name = "Grid_" // TRIM(name) // "_" // &
        TRIM(particle_species(species)%name)
    norm_grid_name = "Norm_Grid_" // TRIM(name) // "_" // &
        TRIM(particle_species(species)%name)
    var_name = TRIM(name) // "_" // TRIM(particle_species(species)%name)

    CALL cfd_write_2d_cartesian_grid(TRIM(grid_name), "Grid", &
        grid1(1:global_resolution(1)), grid2(1:global_resolution(2)), 0)
    CALL cfd_write_2d_cartesian_grid(TRIM(norm_grid_name), "Grid", &
        grid1(1:global_resolution(1))/conv(1), &
        grid2(1:global_resolution(2))/conv(2), 0)

    CALL cfd_write_2d_cartesian_variable_parallel(TRIM(var_name), "dist_fn", &
        global_resolution, stagger, TRIM(norm_grid_name), "Grid", &
        DATA, new_type)
    CALL MPI_TYPE_FREE(new_type, errcode)
    IF (need_reduce) CALL MPI_COMM_FREE(comm_new, errcode)

    DEALLOCATE(DATA)
    DEALLOCATE(grid1, grid2)

  END SUBROUTINE general_2d_dist_fn

END MODULE dist_fn
