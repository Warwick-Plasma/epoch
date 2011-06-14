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
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: load_x
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: load_y
    INTEGER, DIMENSION(:), ALLOCATABLE :: starts_x, ends_x
    INTEGER, DIMENSION(:), ALLOCATABLE :: starts_y, ends_y
    REAL(num) :: balance_frac, npart_av
    REAL(num) :: balance_frac_x, balance_frac_y
    INTEGER(KIND=8) :: min_x, max_x, min_y, max_y
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
      IF (max_npart .LE. 0) RETURN
      CALL MPI_ALLREDUCE(npart_local, sum_npart, 1, MPI_INTEGER8, MPI_SUM, &
          comm, errcode)
      npart_av = REAL(sum_npart, num) / nproc
      balance_frac = (npart_av + SQRT(npart_av)) / REAL(max_npart, num)
      IF (balance_frac .GT. dlb_threshold) RETURN
      IF (rank .EQ. 0) PRINT *, "Load balancing with fraction", balance_frac
    ENDIF

    ALLOCATE(starts_x(1:nprocx), ends_x(1:nprocx))
    ALLOCATE(starts_y(1:nprocy), ends_y(1:nprocy))

    ! Sweep in X
    IF (IAND(balance_mode, c_lb_x) .NE. 0 &
        .OR. IAND(balance_mode, c_lb_auto) .NE. 0) THEN
      ! Rebalancing in X
      ALLOCATE(load_x(nx_global))
      CALL get_load_in_x(load_x)
      CALL calculate_breaks(load_x, nprocx, starts_x, ends_x)
    ELSE
      ! Just keep the original lengths
      starts_x = cell_x_min
      ends_x = cell_x_max
    ENDIF

    ! Sweep in Y
    IF (IAND(balance_mode, c_lb_y) .NE. 0 &
        .OR. IAND(balance_mode, c_lb_auto) .NE. 0) THEN
      ! Rebalancing in Y
      ALLOCATE(load_y(ny_global))
      CALL get_load_in_y(load_y)
      CALL calculate_breaks(load_y, nprocy, starts_y, ends_y)
    ELSE
      ! Just keep the original lengths
      starts_y = cell_y_min
      ends_y = cell_y_max
    ENDIF

    ! In the autobalancer then determine whether to balance in X or Y
    ! Is this worth keeping?
    IF (IAND(balance_mode, c_lb_auto) .NE. 0 ) THEN

      ! Code is auto load balancing
      max_x = 0
      min_x = npart_global
      DO iproc = 1, nprocx
        wk = SUM(load_x(starts_x(iproc):ends_x(iproc)))
        IF (wk .GT. max_x) max_x = wk
        IF (wk .LT. min_x) min_x = wk
      ENDDO

      max_y = 0
      min_y = npart_global
      DO iproc = 1, nprocy
        wk = SUM(load_y(starts_y(iproc):ends_y(iproc)))
        IF (wk .GT. max_y) max_y = wk
        IF (wk .LT. min_y) min_y = wk
      ENDDO

      balance_frac_x = REAL(min_x, num) / REAL(max_x, num)
      balance_frac_y = REAL(min_y, num) / REAL(max_y, num)

      IF (balance_frac_x .LT. balance_frac_y) THEN
        starts_x = cell_x_min
        ends_x = cell_x_max
      ELSE
        starts_y = cell_y_min
        ends_y = cell_y_max
      ENDIF

    ENDIF

    IF (ALLOCATED(load_x)) DEALLOCATE(load_x)
    IF (ALLOCATED(load_y)) DEALLOCATE(load_y)

    ! Now need to calculate the start and end points for the new domain on
    ! the current processor
    nx_global_min = starts_x(x_coords+1)
    nx_global_max = ends_x(x_coords+1)

    ny_global_min = starts_y(y_coords+1)
    ny_global_max = ends_y(y_coords+1)

    domain(1,:) = (/nx_global_min, nx_global_max/)
    domain(2,:) = (/ny_global_min, ny_global_max/)

    ! Redistribute the field variables
    CALL redistribute_fields(domain)

    ! Copy the new lengths into the permanent variables
    cell_x_min = starts_x
    cell_y_min = starts_y
    cell_x_max = ends_x
    cell_y_max = ends_y

    ! Set the new nx, ny
    nx = nx_global_max - nx_global_min + 1
    ny = ny_global_max - ny_global_min + 1

    ! Do X, Y arrays separately because we already have global copies
    DEALLOCATE(x, y)
    ALLOCATE(x(-2:nx+3), y(-2:ny+3))
    x(-2:nx+3) = x_global(nx_global_min-3:nx_global_max+3)
    y(-2:ny+3) = y_global(ny_global_min-3:ny_global_max+3)

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

    ! Set the lengths of the current domain so that the particle balancer
    ! works properly
    x_min_local = x_mins(x_coords)
    x_max_local = x_maxs(x_coords)
    y_min_local = y_mins(y_coords)
    y_max_local = y_maxs(y_coords)

    ! Redistribute the particles onto their new processors
    CALL distribute_particles

    ! If running with particle debugging then set the t = 0 processor if
    ! over_ride = true
#ifdef PARTICLE_DEBUG
    IF (over_ride) THEN
      DO ispecies = 1, n_species
        current=>species_list(ispecies)%attached_list%head
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

    INTEGER :: nx_new, ny_new
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: new_domain
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp3d
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temp
    REAL(num), DIMENSION(:), ALLOCATABLE :: temp1d
    TYPE(laser_block), POINTER :: current
    INTEGER :: ispecies, index, n_species_local

    nx_new = new_domain(1,2) - new_domain(1,1) + 1
    ny_new = new_domain(2,2) - new_domain(2,1) + 1

    ALLOCATE(temp(-2:nx_new+3, -2:ny_new+3))

    temp = 0.0_num
    CALL redistribute_field(new_domain, ex, temp)
    DEALLOCATE(ex)
    ALLOCATE(ex(-2:nx_new+3, -2:ny_new+3))
    ex = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, ey, temp)
    DEALLOCATE(ey)
    ALLOCATE(ey(-2:nx_new+3, -2:ny_new+3))
    ey = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, ez, temp)
    DEALLOCATE(ez)
    ALLOCATE(ez(-2:nx_new+3, -2:ny_new+3))
    ez = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, bx, temp)
    DEALLOCATE(bx)
    ALLOCATE(bx(-2:nx_new+3, -2:ny_new+3))
    bx = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, by, temp)
    DEALLOCATE(by)
    ALLOCATE(by(-2:nx_new+3, -2:ny_new+3))
    by = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, bz, temp)
    DEALLOCATE(bz)
    ALLOCATE(bz(-2:nx_new+3, -2:ny_new+3))
    bz = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, jx, temp)
    DEALLOCATE(jx)
    ALLOCATE(jx(-2:nx_new+3, -2:ny_new+3))
    jx = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, jy, temp)
    DEALLOCATE(jy)
    ALLOCATE(jy(-2:nx_new+3, -2:ny_new+3))
    jy = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, jz, temp)
    DEALLOCATE(jz)
    ALLOCATE(jz(-2:nx_new+3, -2:ny_new+3))
    jz = temp

    DEALLOCATE(temp)

    DO index = 1, num_vars_to_dump
      IF (.NOT. ASSOCIATED(averaged_data(index)%array)) CYCLE

      n_species_local = 0
      IF (IAND(dumpmask(index), c_io_no_sum) .EQ. 0) &
          n_species_local = 1
      IF (IAND(dumpmask(index), c_io_species) .NE. 0) &
          n_species_local = n_species_local + n_species

      ALLOCATE(temp3d(-2:nx_new+3, -2:ny_new+3, 1:n_species_local))

      DO ispecies = 1, n_species_local
        CALL redistribute_field(new_domain, &
            averaged_data(index)%array(:,:,ispecies), temp3d(:,:,ispecies))
      ENDDO

      DEALLOCATE(averaged_data(index)%array)
      ALLOCATE(averaged_data(index)%array(-2:nx_new+3, -2:ny_new+3, &
          n_species_local))

      averaged_data(index)%array = temp3d

      DEALLOCATE(temp3d)
    ENDDO

    ALLOCATE(temp1d(-2:ny_new+3))

    current=>laser_x_min
    DO WHILE(ASSOCIATED(current))
      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain, c_dir_x, current%profile, temp1d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:ny_new+3))
      current%profile = temp1d

      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain, c_dir_x, current%phase, temp1d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:ny_new+3))
      current%phase = temp1d

      current=>current%next
    ENDDO

    current=>laser_x_max
    DO WHILE(ASSOCIATED(current))
      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain, c_dir_x, current%profile, temp1d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:ny_new+3))
      current%profile = temp1d

      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain, c_dir_x, current%phase, temp1d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:ny_new+3))
      current%phase = temp1d

      current=>current%next
    ENDDO

    DEALLOCATE(temp1d)
    ALLOCATE(temp1d(-2:nx_new+3))

    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain, c_dir_y, current%profile, temp1d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:nx_new+3))
      current%profile = temp1d

      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain, c_dir_y, current%phase, temp1d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:nx_new+3))
      current%phase = temp1d

      current=>current%next
    ENDDO

    current=>laser_y_max
    DO WHILE(ASSOCIATED(current))
      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain, c_dir_y, current%profile, temp1d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:nx_new+3))
      current%profile = temp1d

      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain, c_dir_y, current%phase, temp1d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:nx_new+3))
      current%phase = temp1d

      current=>current%next
    ENDDO

    DEALLOCATE(temp1d)

  END SUBROUTINE redistribute_fields



  SUBROUTINE redistribute_field(domain, field_in, field_out)

    ! This subroutine redistributes the fields over the new processor layout
    ! The current version works by writing the field to a file and then each
    ! processor loads back in it's own part. This is better than the previous
    ! version where each processor produced it's own copy of the global array
    ! and then took its own subsection
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: domain
    REAL(num), DIMENSION(-2:,-2:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(-2:,-2:), INTENT(OUT) :: field_out
    INTEGER :: nx_new, ny_new
    INTEGER :: subarray_write, subarray_read
    INTEGER :: subtype_write, subtype_read, fh
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0
    CHARACTER(LEN=9+data_dir_max_length+n_zeros) :: filename

    WRITE(filename, '(a, "/balance.dat")') TRIM(data_dir)

    nx_new = domain(1,2) - domain(1,1) + 1
    ny_new = domain(2,2) - domain(2,1) + 1

    subarray_write = create_current_field_subarray()
    subtype_write  = create_current_field_subtype()

    subarray_read = create_field_subarray(nx_new, ny_new)
    subtype_read  = create_field_subtype(nx_new, ny_new, &
        domain(1,1), domain(2,1))

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

    CALL do_field_mpi_with_lengths(field_out, nx_new, ny_new)

  END SUBROUTINE redistribute_field



  SUBROUTINE redistribute_field_1d(domain, direction, field_in, field_out)

    ! This subroutine redistributes a 1D field over the new processor layout
    ! The current version works by producing a global copy on each processor
    ! And then extracting the required part for the local processor.
    ! in 1D, this is probably OK
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: domain
    INTEGER, INTENT(IN) :: direction
    REAL(num), DIMENSION(-2:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(-2:), INTENT(OUT) :: field_out
    REAL(num), DIMENSION(:), ALLOCATABLE :: field_new, field_temp
    INTEGER :: old_start, new_start
    INTEGER :: old_pts, new_pts, dir
    INTEGER :: npts_global, new_comm, color, coord, i
    INTEGER :: ghost_start, ghost_end

    ghost_start = 0
    ghost_end = 0

    IF (direction .EQ. c_dir_x) THEN
      dir = 2
      color = x_coords
      coord = y_coords
      old_start = cell_y_min(y_coords+1)
      old_pts = ny
      npts_global = ny_global
      IF (y_max_boundary) ghost_end = 3
    ELSE
      dir = 1
      color = y_coords
      coord = x_coords
      old_start = cell_x_min(x_coords+1)
      old_pts = nx
      npts_global = nx_global
      IF (x_max_boundary) ghost_end = 3
    ENDIF

    IF (coord .EQ. 0) ghost_start = -3

    new_start = domain(dir,1)
    new_pts = domain(dir,2) - new_start + 1

    CALL MPI_COMM_SPLIT(comm, color, rank, new_comm, errcode)

    ! Create a global copy of the whole array
    ALLOCATE(field_new(-2:npts_global+3), field_temp(-2:npts_global+3))

    field_new = 0.0_num
    DO i = 1 + ghost_start, old_pts + ghost_end
      field_new(i+old_start-1) = field_in(i)
    ENDDO

    CALL MPI_ALLREDUCE(field_new, field_temp, npts_global+6, &
        mpireal, MPI_SUM, new_comm, errcode)
    DEALLOCATE(field_new)

    DO i = -2, new_pts + 3
      field_out(i) = field_temp(i+new_start-1)
    ENDDO

    DEALLOCATE(field_temp)
    CALL MPI_COMM_FREE(new_comm, errcode)

  END SUBROUTINE redistribute_field_1d



  SUBROUTINE get_load_in_x(load)

    ! Calculate total load across the X direction
    ! Summed in the Y direction

    INTEGER(KIND=8), DIMENSION(:), INTENT(OUT) :: load
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: temp
    TYPE(particle), POINTER :: current
    INTEGER :: cell, ispecies, sz

    load = 0

    DO ispecies = 1, n_species
      current=>species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so x_min, NOT x_min_local
#ifdef PARTICLE_SHAPE_TOPHAT
        cell = FLOOR((current%part_pos(1) - x_min) / dx) + 1
#else
        cell = FLOOR((current%part_pos(1) - x_min) / dx + 1.5_num)
#endif
        load(cell) = load(cell) + 1
        current=>current%next
      ENDDO
    ENDDO

    ! Now have local densities, so add using MPI
    sz = SIZE(load)
    ALLOCATE(temp(sz))
    CALL MPI_ALLREDUCE(load, temp, sz, MPI_INTEGER8, MPI_SUM, comm, errcode)

    ! Adjust the load of pushing one particle relative to the load
    ! of updating one field cell, then add on the field load.
    ! The push_per_field factor will be updated automatically in future.
    load = push_per_field * temp + ny

    DEALLOCATE(temp)

  END SUBROUTINE get_load_in_x



  SUBROUTINE get_load_in_y(load)

    ! Calculate total load across the Y direction
    ! Summed in the X direction

    INTEGER(KIND=8), DIMENSION(:), INTENT(OUT) :: load
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: temp
    TYPE(particle), POINTER :: current
    INTEGER :: cell, ispecies, sz

    load = 0

    DO ispecies = 1, n_species
      current=>species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so y_min, NOT y_min_local
#ifdef PARTICLE_SHAPE_TOPHAT
        cell = FLOOR((current%part_pos(2) - y_min) / dy) + 1
#else
        cell = FLOOR((current%part_pos(2) - y_min) / dy + 1.5_num)
#endif
        load(cell) = load(cell) + 1
        current=>current%next
      ENDDO
    ENDDO

    ! Now have local densities, so add using MPI
    sz = SIZE(load)
    ALLOCATE(temp(sz))
    CALL MPI_ALLREDUCE(load, temp, sz, MPI_INTEGER8, MPI_SUM, comm, errcode)

    ! Adjust the load of pushing one particle relative to the load
    ! of updating one field cell, then add on the field load.
    ! The push_per_field factor will be updated automatically in future.
    load = push_per_field * temp + nx

    DEALLOCATE(temp)

  END SUBROUTINE get_load_in_y



  SUBROUTINE calculate_breaks(load, nproc, starts, ends)

    ! This subroutine calculates the places in a given load profile to split
    ! The domain to give the most even subdivision possible

    INTEGER(KIND=8), INTENT(IN), DIMENSION(:) :: load
    INTEGER, INTENT(IN) :: nproc
    INTEGER, DIMENSION(:), INTENT(OUT) :: starts, ends
    INTEGER :: sz, idim, partition, proc, old, shift
    INTEGER(KIND=8) :: total, total_old
    INTEGER :: load_per_proc_ideal

    sz = SIZE(load)
    ends = sz

    load_per_proc_ideal = FLOOR((SUM(load) + 0.5d0) / nproc)

    proc = 1
    total = 0
    DO idim = 1, sz
      total_old = total
      total = total + load(idim)
      IF (total .GE. load_per_proc_ideal) THEN
        ! Pick the split that most closely matches the load
        IF (load_per_proc_ideal - total_old &
            .LT. total - load_per_proc_ideal) THEN
          ends(proc) = idim - 1
        ELSE
          ends(proc) = idim
        ENDIF
        proc = proc + 1
        total = total - load_per_proc_ideal
        IF (proc .EQ. nproc) EXIT
      ENDIF
    ENDDO

    ! Sanity check. Must be one cell of separation between each endpoint.
    ! Forwards (unnecessary?)
    !old = 0
    !DO proc = 1, nproc
    !  IF (ends(proc) - old .LE. 0) THEN
    !    ends(proc) = old + 1
    !  ENDIF
    !  old = ends(proc)
    !ENDDO

    ! Backwards
    old = sz + 1
    DO proc = nproc, 1, -1
      IF (old - ends(proc) .LE. 0) THEN
        ends(proc) = old - 1
      ENDIF
      old = ends(proc)
    ENDDO

    ! Set starts
    starts(1) = 1
    DO proc = 2, nproc
      starts(proc) = ends(proc-1) + 1
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

    IF (MINVAL(coords) .LT. 0) THEN
      WRITE(*,*) 'UNLOCATABLE PARTICLE', coords
      RETURN
    ENDIF
    CALL MPI_CART_RANK(comm, coords, get_particle_processor, errcode)
    ! IF (get_particle_processor .NE. rank) PRINT *,

  END FUNCTION get_particle_processor



  ! This subroutine is used to rearrange particles over processors
  SUBROUTINE distribute_particles

    ! This subroutine moves particles which are on the wrong processor
    ! to the correct processor.

    TYPE(particle_list), DIMENSION(:), ALLOCATABLE :: pointers_send
    TYPE(particle_list), DIMENSION(:), ALLOCATABLE :: pointers_recv
    TYPE(particle), POINTER :: current, next
    INTEGER :: part_proc, iproc_recv, iproc_send, ispecies, ierr

    ALLOCATE(pointers_send(0:nproc-1), pointers_recv(0:nproc-1))
    DO ispecies = 1, n_species
      current=>species_list(ispecies)%attached_list%head
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
              species_list(ispecies)%attached_list, current)
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
        CALL append_partlist(species_list(ispecies)%attached_list, &
            pointers_recv(iproc_recv))
      ENDDO
    ENDDO

    DEALLOCATE(pointers_send, pointers_recv)

  END SUBROUTINE distribute_particles

END MODULE balance
