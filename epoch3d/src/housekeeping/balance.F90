MODULE balance

  USE partlist
  USE boundary
  USE mpi_subtype_control

  IMPLICIT NONE

CONTAINS

  SUBROUTINE balance_workload(over_ride)

    ! This subroutine determines whether or not the code needs rebalancing,
    ! calculates where to split the domain and calls other subroutines to
    ! actually rearrange the fields and particles onto the new processors

    ! This is really, really hard to do properly
    ! So cheat

    LOGICAL, INTENT(IN) :: over_ride
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: density_x
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: density_y
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: density_z
    INTEGER, DIMENSION(:), ALLOCATABLE :: starts_x, ends_x
    INTEGER, DIMENSION(:), ALLOCATABLE :: starts_y, ends_y
    INTEGER, DIMENSION(:), ALLOCATABLE :: starts_z, ends_z
    INTEGER :: new_cell_x_min, new_cell_x_max
    INTEGER :: new_cell_y_min, new_cell_y_max
    INTEGER :: new_cell_z_min, new_cell_z_max
    REAL(num) :: balance_frac, npart_av
    REAL(num) :: balance_frac_x, balance_frac_y, balance_frac_z
    INTEGER(KIND=8) :: min_x, max_x, min_y, max_y, min_z, max_z
    INTEGER(KIND=8) :: npart_local, sum_npart, max_npart, wk
    INTEGER :: iproc
    INTEGER, DIMENSION(c_ndims,2) :: domain
#ifdef PARTICLE_DEBUG
    TYPE(particle), POINTER :: current
    INTEGER :: ispecies
#endif

    ! On one processor do nothing to save time
    IF (nproc .EQ. 1) RETURN

    ! This parameter allows selecting the mode of the autobalancing between
    ! leftsweep, rightsweep, auto(best of leftsweep and rightsweep) or both
    balance_mode = c_lb_all

    ! count particles
    npart_local = get_total_local_particles()

    ! The over_ride flag allows the code to force a load balancing sweep
    ! at t = 0
    IF (.NOT. over_ride) THEN
      CALL MPI_ALLREDUCE(npart_local, max_npart, 1, MPI_INTEGER8, MPI_MAX, &
          comm, errcode)
      CALL MPI_ALLREDUCE(npart_local, sum_npart, 1, MPI_INTEGER8, MPI_SUM, &
          comm, errcode)
      npart_av = REAL(sum_npart, num) / nproc
      balance_frac = (npart_av + SQRT(npart_av)) / REAL(max_npart, num)
      IF (balance_frac .GT. dlb_threshold) RETURN
      IF (rank .EQ. 0) PRINT *, "Load balancing with fraction", balance_frac
    ENDIF

    ALLOCATE(starts_x(1:nprocx), ends_x(1:nprocx))
    ALLOCATE(starts_y(1:nprocy), ends_y(1:nprocy))
    ALLOCATE(starts_z(1:nprocz), ends_z(1:nprocz))

    ! Sweep in X
    IF (IAND(balance_mode, c_lb_x) .NE. 0 &
        .OR. IAND(balance_mode, c_lb_auto) .NE. 0) THEN
      ! Rebalancing in X
      ALLOCATE(density_x(0:nx_global+1))
      CALL get_density_in_x(density_x)
      CALL calculate_breaks(density_x, nprocx, starts_x, ends_x)
    ELSE
      ! Just keep the original lengths
      starts_x = cell_x_min
      ends_x = cell_x_max
    ENDIF

    ! Sweep in Y
    IF (IAND(balance_mode, c_lb_y) .NE. 0 &
        .OR. IAND(balance_mode, c_lb_auto) .NE. 0) THEN
      ! Rebalancing in Y
      ALLOCATE(density_y(0:ny_global+1))
      CALL get_density_in_y(density_y)
      CALL calculate_breaks(density_y, nprocy, starts_y, ends_y)
    ELSE
      ! Just keep the original lengths
      starts_y = cell_y_min
      ends_y = cell_y_max
    ENDIF

    ! Sweep in Z
    IF (IAND(balance_mode, c_lb_z) .NE. 0 &
        .OR. IAND(balance_mode, c_lb_auto) .NE. 0) THEN
      ! Rebalancing in Z
      ALLOCATE(density_z(0:nz_global+1))
      CALL get_density_in_z(density_z)
      CALL calculate_breaks(density_z, nprocz, starts_z, ends_z)
    ELSE
      ! Just keep the original lengths
      starts_z = cell_z_min
      ends_z = cell_z_max
    ENDIF

    ! In the autobalancer then determine whether to balance in X, Y or Z
    ! Is this worth keeping?
    IF (IAND(balance_mode, c_lb_auto) .NE. 0 ) THEN

      ! Code is auto load balancing
      max_x = 0
      min_x = npart_global
      DO iproc = 1, nprocx
        wk = SUM(density_x(starts_x(iproc):ends_x(iproc)))
        IF (wk .GT. max_x) max_x = wk
        IF (wk .LT. min_x) min_x = wk
      ENDDO

      max_y = 0
      min_y = npart_global
      DO iproc = 1, nprocy
        wk = SUM(density_y(starts_y(iproc):ends_y(iproc)))
        IF (wk .GT. max_y) max_y = wk
        IF (wk .LT. min_y) min_y = wk
      ENDDO

      max_z = 0
      min_z = npart_global
      DO iproc = 1, nprocz
        wk = SUM(density_z(starts_z(iproc):ends_z(iproc)))
        IF (wk .GT. max_z) max_z = wk
        IF (wk .LT. min_z) min_z = wk
      ENDDO

      balance_frac_x = REAL(min_x, num) / REAL(max_x, num)
      balance_frac_y = REAL(min_y, num) / REAL(max_y, num)
      balance_frac_z = REAL(min_z, num) / REAL(max_z, num)

      IF (balance_frac_x .LT. balance_frac_y) THEN
        IF (balance_frac_x .LT. balance_frac_z) THEN
          starts_x = cell_x_min
          ends_x = cell_x_max
        ELSE
          starts_z = cell_z_min
          ends_z = cell_z_max
        ENDIF
      ELSE
        IF (balance_frac_y .LT. balance_frac_z) THEN
          starts_y = cell_y_min
          ends_y = cell_y_max
        ELSE
          starts_z = cell_z_min
          ends_z = cell_z_max
        ENDIF
      ENDIF

    ENDIF

    ! Now need to calculate the start and end points for the new domain on
    ! the current processor
    new_cell_x_min = starts_x(coordinates(c_ndims)+1)
    new_cell_x_max = ends_x(coordinates(c_ndims)+1)

    new_cell_y_min = starts_y(coordinates(c_ndims-1)+1)
    new_cell_y_max = ends_y(coordinates(c_ndims-1)+1)

    new_cell_z_min = starts_z(coordinates(c_ndims-2)+1)
    new_cell_z_max = ends_z(coordinates(c_ndims-2)+1)

    domain(1,:) = (/new_cell_x_min, new_cell_x_max/)
    domain(2,:) = (/new_cell_y_min, new_cell_y_max/)
    domain(3,:) = (/new_cell_z_min, new_cell_z_max/)

    ! Redistribute the field variables
    CALL redistribute_fields(domain)

    ! Copy the new lengths into the permanent variables
    cell_x_min = starts_x
    cell_y_min = starts_y
    cell_z_min = starts_z
    cell_x_max = ends_x
    cell_y_max = ends_y
    cell_z_max = ends_z

    ! Set the new nx, ny, nz
    nx = new_cell_x_max - new_cell_x_min + 1
    ny = new_cell_y_max - new_cell_y_min + 1
    nz = new_cell_z_max - new_cell_z_min + 1

    ! Do X, Y, Z arrays separately because we already have global copies
    DEALLOCATE(x, y, z)
    ALLOCATE(x(-2:nx+3), y(-2:ny+3), z(-2:nz+3))
    x(0:nx+1) = x_global(new_cell_x_min-1:new_cell_x_max+1)
    y(0:ny+1) = y_global(new_cell_y_min-1:new_cell_y_max+1)
    z(0:nz+1) = z_global(new_cell_z_min-1:new_cell_z_max+1)

    ! Recalculate x_mins and x_maxs so that rebalancing works next time
    DO iproc = 0, nprocx - 1
      x_mins(iproc) = x_global(cell_x_min(iproc+1))
      x_maxs(iproc) = x_global(cell_x_max(iproc+1))
    ENDDO
    ! Same for y
    DO iproc = 0, nprocy - 1
      y_mins(iproc) = y_global(cell_y_min(iproc+1))
      y_maxs(iproc) = y_global(cell_y_max(iproc+1))
    ENDDO
    ! Same for z
    DO iproc = 0, nprocz - 1
      z_mins(iproc) = z_global(cell_z_min(iproc+1))
      z_maxs(iproc) = z_global(cell_z_max(iproc+1))
    ENDDO

    ! Set the lengths of the current domain so that the particle balancer
    ! works properly
    x_min_local = x_mins(coordinates(c_ndims))
    x_max_local = x_maxs(coordinates(c_ndims))
    y_min_local = y_mins(coordinates(c_ndims-1))
    y_max_local = y_maxs(coordinates(c_ndims-1))
    z_min_local = z_mins(coordinates(c_ndims-2))
    z_max_local = z_maxs(coordinates(c_ndims-2))

    ! Redistribute the particles onto their new processors
    CALL distribute_particles

    ! If running with particle debugging then set the t = 0 processor if
    ! over_ride = true
#ifdef PARTICLE_DEBUG
    IF (over_ride) THEN
      DO ispecies = 1, n_species
        current=>particle_species(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          current%processor_at_t0 = rank
          current=>current%next
        ENDDO
      ENDDO
    ENDIF
#endif

  END SUBROUTINE balance_workload



  SUBROUTINE redistribute_fields(new_domain)

    ! This subroutine redistributes the 2D field variables over the new
    ! processor layout. If using a 2D field of your own then set the
    ! redistribute_field subroutine to implement it. 1D fields, you're on
    ! your own (have global copies and use those to repopulate?)

    INTEGER :: nx_new, ny_new, nz_new
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: new_domain
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temp2d
    TYPE(laser_block), POINTER :: current

    nx_new = new_domain(1,2) - new_domain(1,1) + 1
    ny_new = new_domain(2,2) - new_domain(2,1) + 1
    nz_new = new_domain(3,2) - new_domain(3,1) + 1

    ALLOCATE(temp(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))

    temp = 0.0_num
    CALL redistribute_field(new_domain, ex, temp)
    DEALLOCATE(ex)
    ALLOCATE(ex(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    ex = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, ey, temp)
    DEALLOCATE(ey)
    ALLOCATE(ey(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    ey = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, ez, temp)
    DEALLOCATE(ez)
    ALLOCATE(ez(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    ez = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, bx, temp)
    DEALLOCATE(bx)
    ALLOCATE(bx(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    bx = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, by, temp)
    DEALLOCATE(by)
    ALLOCATE(by(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    by = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, bz, temp)
    DEALLOCATE(bz)
    ALLOCATE(bz(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    bz = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, jx, temp)
    DEALLOCATE(jx)
    ALLOCATE(jx(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    jx = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, jy, temp)
    DEALLOCATE(jy)
    ALLOCATE(jy(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    jy = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, jz, temp)
    DEALLOCATE(jz)
    ALLOCATE(jz(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    jz = temp

    DEALLOCATE(temp)

    ALLOCATE(temp2d(-2:ny_new+3, -2:nz_new+3))

    current=>laser_x_min
    DO WHILE(ASSOCIATED(current))
      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_x, current%profile, temp2d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:ny_new+3, -2:nz_new+3))
      current%profile = temp2d

      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_x, current%phase, temp2d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:ny_new+3, -2:nz_new+3))
      current%phase = temp2d

      current=>current%next
    ENDDO

    current=>laser_x_max
    DO WHILE(ASSOCIATED(current))
      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_x, current%profile, temp2d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:ny_new+3, -2:nz_new+3))
      current%profile = temp2d

      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_x, current%phase, temp2d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:ny_new+3, -2:nz_new+3))
      current%phase = temp2d

      current=>current%next
    ENDDO

    DEALLOCATE(temp2d)
    ALLOCATE(temp2d(-2:nx_new+3, -2:nz_new+3))

    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_y, current%profile, temp2d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:nx_new+3, -2:nz_new+3))
      current%profile = temp2d

      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_y, current%phase, temp2d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:nx_new+3, -2:nz_new+3))
      current%phase = temp2d

      current=>current%next
    ENDDO

    current=>laser_y_max
    DO WHILE(ASSOCIATED(current))
      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_y, current%profile, temp2d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:nx_new+3, -2:nz_new+3))
      current%profile = temp2d

      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_y, current%phase, temp2d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:nx_new+3, -2:nz_new+3))
      current%phase = temp2d

      current=>current%next
    ENDDO

    DEALLOCATE(temp2d)
    ALLOCATE(temp2d(-2:nx_new+3, -2:ny_new+3))

    current=>laser_z_min
    DO WHILE(ASSOCIATED(current))
      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_z, current%profile, temp2d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:nx_new+3, -2:ny_new+3))
      current%profile = temp2d

      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_z, current%phase, temp2d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:nx_new+3, -2:ny_new+3))
      current%phase = temp2d

      current=>current%next
    ENDDO

    current=>laser_z_max
    DO WHILE(ASSOCIATED(current))
      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_z, current%profile, temp2d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:nx_new+3, -2:ny_new+3))
      current%profile = temp2d

      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_z, current%phase, temp2d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:nx_new+3, -2:ny_new+3))
      current%phase = temp2d

      current=>current%next
    ENDDO

    DEALLOCATE(temp2d)

  END SUBROUTINE redistribute_fields



  SUBROUTINE redistribute_field(domain, field_in, field_out)

    ! This subroutine redistributes the fields over the new processor layout
    ! The current version works by writing the field to a file and then each
    ! processor loads back in it's own part. This is better than the previous
    ! version where each processor produced it's own copy of the global array
    ! and then took its own subsection
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: domain
    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: field_out
    INTEGER :: nx_new, ny_new, nz_new
    INTEGER :: subarray_write, subarray_read
    INTEGER :: subtype_write, subtype_read, fh
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0
    CHARACTER(LEN=9+data_dir_max_length+n_zeros) :: filename

    WRITE(filename, '(a, "/balance.dat")') TRIM(data_dir)

    nx_new = domain(1,2) - domain(1,1) + 1
    ny_new = domain(2,2) - domain(2,1) + 1
    nz_new = domain(3,2) - domain(3,1) + 1

    subarray_write = create_current_field_subarray()
    subtype_write  = create_current_field_subtype()

    subarray_read = create_field_subarray(nx_new, ny_new, nz_new)
    subtype_read  = create_field_subtype(nx_new, ny_new, nz_new, &
        domain(1,1), domain(2,1), domain(3,1))

    CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_RDWR+MPI_MODE_CREATE, &
        MPI_INFO_NULL, fh, errcode)

    CALL MPI_FILE_SET_VIEW(fh, offset, subarray_write, subtype_write, &
        "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(fh, field_in, 1, subarray_write, status, errcode)

    CALL MPI_FILE_SET_VIEW(fh, offset, subarray_read, subtype_read, &
        "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_READ_ALL(fh, field_out, 1, subarray_read, status, errcode)

    CALL MPI_FILE_CLOSE(fh, errcode)

    CALL MPI_TYPE_FREE(subarray_write, errcode)
    CALL MPI_TYPE_FREE(subtype_write, errcode)
    CALL MPI_TYPE_FREE(subarray_read, errcode)
    CALL MPI_TYPE_FREE(subtype_read, errcode)

    CALL do_field_mpi_with_lengths(field_out, nx_new, ny_new, nz_new)

  END SUBROUTINE redistribute_field



  SUBROUTINE redistribute_field_2d(domain, direction, field_in, field_out)

    ! This subroutine redistributes the fields over the new processor layout
    ! The current version works by writing the field to a file and then each
    ! processor loads back in it's own part. This is better than the previous
    ! version where each processor produced it's own copy of the global array
    ! and then took its own subsection
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: domain
    INTEGER, INTENT(IN) :: direction
    REAL(num), DIMENSION(-2:,-2:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(-2:,-2:), INTENT(OUT) :: field_out
    INTEGER :: n1, n2, n1_new, n2_new, n1_global, n2_global
    INTEGER :: n1_old_start, n2_old_start, n1_dir, n2_dir
    INTEGER :: subtype_write, subtype_read, fh
    INTEGER :: subarray_write, subarray_read
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0
    CHARACTER(LEN=9+data_dir_max_length+n_zeros) :: filename

    IF (direction .EQ. c_dir_x) THEN
      n1_dir = 2
      n2_dir = 3
      n1 = ny
      n2 = nz
      n1_global = ny_global
      n2_global = nz_global
      n1_old_start = cell_y_min(coordinates(c_ndims-1)+1)
      n2_old_start = cell_z_min(coordinates(c_ndims-2)+1)
    ELSE IF (direction .EQ. c_dir_y) THEN
      n1_dir = 1
      n2_dir = 3
      n1 = nx
      n2 = nz
      n1_global = nx_global
      n2_global = nz_global
      n1_old_start = cell_x_min(coordinates(c_ndims)+1)
      n2_old_start = cell_z_min(coordinates(c_ndims-2)+1)
    ELSE
      n1_dir = 1
      n2_dir = 2
      n1 = nx
      n2 = ny
      n1_global = nx_global
      n2_global = ny_global
      n1_old_start = cell_x_min(coordinates(c_ndims)+1)
      n2_old_start = cell_y_min(coordinates(c_ndims-1)+1)
    ENDIF

    WRITE(filename, '(a, "/balance.dat")') TRIM(data_dir)

    n1_new = domain(n1_dir,2) - domain(n1_dir,1) + 1
    n2_new = domain(n2_dir,2) - domain(n2_dir,1) + 1

    subarray_write = create_field_subarray(n1, n2)
    subtype_write = create_2d_array_subtype((/n1, n2/), &
        (/n1_global, n2_global/), (/n1_old_start, n2_old_start/))

    subarray_read = create_field_subarray(n1_new, n2_new)
    subtype_read  = create_2d_array_subtype((/n1_new, n2_new/), &
        (/n1_global, n2_global/), (/domain(n1_dir,1), domain(n2_dir,1)/))

    CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_RDWR+MPI_MODE_CREATE, &
        MPI_INFO_NULL, fh, errcode)

    CALL MPI_FILE_SET_VIEW(fh, offset, mpireal, subtype_write, "native", &
        MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(fh, field_in, 1, subarray_write, status, errcode)

    CALL MPI_FILE_SET_VIEW(fh, offset, mpireal, subtype_read, "native", &
        MPI_INFO_NULL, errcode)
    CALL MPI_FILE_READ_ALL(fh, field_out, 1, subarray_read, status, errcode)

    CALL MPI_FILE_CLOSE(fh, errcode)

    CALL MPI_TYPE_FREE(subtype_write, errcode)
    CALL MPI_TYPE_FREE(subtype_read, errcode)
    CALL MPI_TYPE_FREE(subarray_write, errcode)
    CALL MPI_TYPE_FREE(subarray_read, errcode)

  END SUBROUTINE redistribute_field_2d



  SUBROUTINE get_density_in_x(density)

    ! Calculate total particle density across the X direction
    ! Summed in the Y and Z directions

    INTEGER(KIND=8), DIMENSION(:), INTENT(INOUT) :: density
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: temp
    TYPE(particle), POINTER :: current
    REAL(num) :: part_x
    INTEGER :: cell_x1, ispecies

    density = 0.0_num

    DO ispecies = 1, n_species
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so x_min, NOT x_min_local
        part_x = current%part_pos(1) - x_min
        cell_x1 = NINT(part_x / dx) + 1
        density(cell_x1) = density(cell_x1) + 1
        current=>current%next
      ENDDO
    ENDDO

    ! Now have local densities, so add using MPI
    ALLOCATE(temp(0:nx_global+1))
    CALL MPI_ALLREDUCE(density, temp, nx_global+2, MPI_INTEGER8, &
        MPI_SUM, comm, errcode)
    density = temp

    DEALLOCATE(temp)

  END SUBROUTINE get_density_in_x



  SUBROUTINE get_density_in_y(density)

    ! Calculate total particle density across the Y direction
    ! Summed in the X and Z directions

    INTEGER(KIND=8), DIMENSION(:), INTENT(INOUT) :: density
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: temp
    TYPE(particle), POINTER :: current
    REAL(num) :: part_y
    INTEGER :: cell_y1, ispecies

    density = 0.0_num

    DO ispecies = 1, n_species
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so y_min, NOT y_min_local
        part_y = current%part_pos(2) - y_min
        cell_y1 = NINT(part_y / dy) + 1
        density(cell_y1) = density(cell_y1) + 1
        current=>current%next
      ENDDO
    ENDDO

    ! Now have local densities, so add using MPI
    ALLOCATE(temp(0:ny_global+1))
    CALL MPI_ALLREDUCE(density, temp, ny_global+2, MPI_INTEGER8, &
        MPI_SUM, comm, errcode)
    density = temp

    DEALLOCATE(temp)

  END SUBROUTINE get_density_in_y



  SUBROUTINE get_density_in_z(density)

    ! Calculate total particle density across the Z direction
    ! Summed in the X and Y directions

    INTEGER(KIND=8), DIMENSION(:), INTENT(INOUT) :: density
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: temp
    TYPE(particle), POINTER :: current
    REAL(num) :: part_z
    INTEGER :: cell_z1, ispecies

    density = 0.0_num

    DO ispecies = 1, n_species
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so z_min, NOT z_min_local
        part_z = current%part_pos(3) - z_min
        cell_z1 = NINT(part_z / dz) + 1
        density(cell_z1) = density(cell_z1) + 1
        current=>current%next
      ENDDO
    ENDDO

    ! Now have local densities, so add using MPI
    ALLOCATE(temp(0:nz_global+1))
    CALL MPI_ALLREDUCE(density, temp, nz_global+2, MPI_INTEGER8, &
        MPI_SUM, comm, errcode)
    density = temp

    DEALLOCATE(temp)

  END SUBROUTINE get_density_in_z



  SUBROUTINE calculate_breaks(density, nproc, starts, ends)

    ! This subroutine calculates the places in a given density profile to split
    ! The domain to give the most even subdivision possible

    INTEGER(KIND=8), INTENT(IN), DIMENSION(:) :: density
    INTEGER, INTENT(IN) :: nproc
    INTEGER, DIMENSION(:), INTENT(OUT) :: starts, ends
    INTEGER :: sz, idim, partition
    INTEGER(KIND=8) :: total
    INTEGER :: npart_per_proc_ideal

    ! -2 because of ghost cells at each end
    sz = SIZE(density) - 2
    IF (nproc .EQ. 1) THEN
      starts = 1
      ends = sz
    ENDIF

    npart_per_proc_ideal = INT(SUM(density) / nproc)
    partition = 2
    starts(1) = 1
    total = 0
    DO idim = 1, sz
      IF (partition .GT. nproc) EXIT
      IF (total .GE. npart_per_proc_ideal &
          .OR. ABS(total + density(idim) - npart_per_proc_ideal) &
          .GT. ABS(total - npart_per_proc_ideal)  .OR. idim .EQ. sz) THEN
        total = density(idim)
        starts(partition) = idim - 1
        partition = partition + 1
        ! If you've reached the last processor, have already done the best
        ! you can, so just leave
      ELSE
        total = total + density(idim)
      ENDIF
    ENDDO

    ends(nproc) = sz
    DO idim = 1, nproc - 1
      ends(idim) = starts(idim+1) - 1
    ENDDO

  END SUBROUTINE calculate_breaks



  FUNCTION get_particle_processor(a_particle)

    ! This subroutine calculates which processor a given particles resides on

    TYPE(particle), INTENT(IN) :: a_particle
    INTEGER :: get_particle_processor
    INTEGER :: iproc, coords(c_ndims)

    get_particle_processor = -1
    coords = -1

    ! This could be replaced by a bisection method, but for the moment I
    ! just don't care

    DO iproc = 0, nprocx - 1
      IF (a_particle%part_pos(1) .GE. x_mins(iproc) - dx / 2.0_num &
          .AND. a_particle%part_pos(1) .LT. x_maxs(iproc) + dx / 2.0_num) THEN
        coords(c_ndims) = iproc
        EXIT
      ENDIF
    ENDDO

    DO iproc = 0, nprocy - 1
      IF (a_particle%part_pos(2) .GE. y_mins(iproc) - dy / 2.0_num &
          .AND. a_particle%part_pos(2) .LT. y_maxs(iproc) + dy / 2.0_num) THEN
        coords(c_ndims-1) = iproc
        EXIT
      ENDIF
    ENDDO

    DO iproc = 0, nprocz - 1
      IF (a_particle%part_pos(3) .GE. z_mins(iproc) - dz / 2.0_num &
          .AND. a_particle%part_pos(3) .LT. z_maxs(iproc) + dz / 2.0_num) THEN
        coords(c_ndims-2) = iproc
        EXIT
      ENDIF
    ENDDO

    IF (MINVAL(coords) .LT. 0) PRINT *, "UNLOCATABLE PARTICLE", coords
    IF (MINVAL(coords) .LT. 0) RETURN
    CALL MPI_CART_RANK(comm, coords, get_particle_processor, errcode)
    ! IF (get_particle_processor .NE. rank) PRINT *,

  END FUNCTION get_particle_processor



  ! This subroutine is used to rearrange particles over processors
  SUBROUTINE distribute_particles

    ! This subroutine actually moves particles which are on the wrong processor
    ! And moves then to the correct processor.

    TYPE(particle_list), DIMENSION(:), ALLOCATABLE :: pointers_send
    TYPE(particle_list), DIMENSION(:), ALLOCATABLE :: pointers_recv
    TYPE(particle), POINTER :: current, next
    INTEGER :: part_proc, iproc_recv, iproc_send, ispecies, ierr

    ALLOCATE(pointers_send(0:nproc-1), pointers_recv(0:nproc-1))
    DO ispecies = 1, n_species
      current=>particle_species(ispecies)%attached_list%head
      DO iproc_send = 0, nproc - 1
        CALL create_empty_partlist(pointers_send(iproc_send))
        CALL create_empty_partlist(pointers_recv(iproc_send))
      ENDDO

      DO WHILE(ASSOCIATED(current))
        next=>current%next
        part_proc = get_particle_processor(current)
        IF (part_proc .LT. 0) THEN
          PRINT *, "Unlocatable particle on processor", rank, current%part_pos
          CALL MPI_ABORT(comm, errcode, ierr)
          STOP
        ENDIF
#ifdef PARTICLE_DEBUG
        current%processor = part_proc
#endif
        IF (part_proc .NE. rank) THEN
          CALL remove_particle_from_partlist(&
              particle_species(ispecies)%attached_list, current)
          CALL add_particle_to_partlist(pointers_send(part_proc), current)
        ENDIF
        current=>next
      ENDDO

      DO iproc_send = 0, nproc - 1
        DO iproc_recv = 0, nproc - 1
          IF (iproc_send .NE. iproc_recv) THEN
            IF (rank .EQ. iproc_send) THEN
              CALL partlist_send(pointers_send(iproc_recv), iproc_recv)
              CALL destroy_partlist(pointers_send(iproc_recv))
            ENDIF
            IF (rank .EQ. iproc_recv) &
                CALL partlist_recv(pointers_recv(iproc_send), iproc_send)
          ENDIF
        ENDDO
      ENDDO

      DO iproc_recv = 0, nproc - 1
        CALL append_partlist(particle_species(ispecies)%attached_list, &
            pointers_recv(iproc_recv))
      ENDDO
    ENDDO

    DEALLOCATE(pointers_send, pointers_recv)

  END SUBROUTINE distribute_particles

END MODULE balance
