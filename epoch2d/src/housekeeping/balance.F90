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
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: npart_each_rank
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: density_x, density_y
    INTEGER, DIMENSION(:), ALLOCATABLE :: starts_x, ends_x
    INTEGER, DIMENSION(:), ALLOCATABLE :: starts_y, ends_y
    INTEGER :: new_cell_x_min, new_cell_x_max
    INTEGER :: new_cell_y_min, new_cell_y_max
    REAL(num) :: balance_frac, balance_frac_x, balance_frac_y
    INTEGER(KIND=8) :: max_x, max_y, wk, min_x, min_y, npart_local
    INTEGER :: iproc
    INTEGER, DIMENSION(2, 2) :: domain
#ifdef PARTICLE_DEBUG
    TYPE(particle), POINTER :: current
    INTEGER :: ispecies
#endif

    ! On one processor do nothing to SAVE time
    IF (nproc .EQ. 1) RETURN

    ! This parameter allows selecting the mode of the autobalancing between
    ! leftsweep, rightsweep, auto(best of leftsweep and rightsweep) or both
    balance_mode = c_lb_both

    ! count particles
    npart_local = get_total_local_particles()

    ! The over_ride flag allows the code to force a load balancing sweep
    ! at t = 0
    IF (.NOT. over_ride) THEN
      ALLOCATE(npart_each_rank(1:nproc))
      ! Get npart for each rank
      CALL MPI_ALLGATHER(npart_local, 1, MPI_INTEGER8, npart_each_rank, 1, &
          MPI_INTEGER8, comm, errcode)
      ! Determine ratio of npart on between most loaded and least loaded
      ! processor. Maybe this can be replaced by and MPI_ALLREDUCE to
      ! find min/max?
      balance_frac = REAL(MINVAL(npart_each_rank), num) &
          / REAL(MAXVAL(npart_each_rank), num)
      IF (balance_frac .GT. dlb_threshold) RETURN
      IF (rank .EQ. 0) PRINT *, "Load balancing with fraction", balance_frac
      DEALLOCATE(npart_each_rank)
    ENDIF

    ALLOCATE(starts_x(1:nprocx), ends_x(1:nprocx))
    ALLOCATE(starts_y(1:nprocy), ends_y(1:nprocy))

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

    ! In the autobalancer then determine whether to balance in X or Y
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

      balance_frac_x = REAL(min_x, num)/REAL(max_x, num)
      balance_frac_y = REAL(min_y, num)/REAL(max_y, num)

      IF (balance_frac_y .LT. balance_frac_x) THEN
        starts_y = cell_y_min
        ends_y = cell_y_max
      ELSE
        starts_x = cell_x_min
        ends_x = cell_x_max
      ENDIF

    ENDIF

    ! Now need to calculate the start and end points for the new domain on
    ! the current processor
    new_cell_x_min = starts_x(coordinates(2)+1)
    new_cell_x_max = ends_x(coordinates(2)+1)

    new_cell_y_min = starts_y(coordinates(1)+1)
    new_cell_y_max = ends_y(coordinates(1)+1)

    ! Redeistribute the field variables
    domain(1,:) = (/new_cell_x_min, new_cell_x_max/)
    domain(2,:) = (/new_cell_y_min, new_cell_y_max/)
    CALL redistribute_fields(domain)

    ! Copy the new lengths into the permanent variables
    cell_x_min = starts_x
    cell_y_min = starts_y
    cell_x_max = ends_x
    cell_y_max = ends_y

    ! Set the new nx and ny
    nx = new_cell_x_max-new_cell_x_min+1
    ny = new_cell_y_max-new_cell_y_min+1

    ! Do X and Y arrays separatly because we already have global copies
    ! of X and Y
    DEALLOCATE(x, y)
    ALLOCATE(x(-2:nx+3), y(-2:ny+3))
    x(0:nx+1) = x_global(new_cell_x_min-1:new_cell_x_max+1)
    y(0:ny+1) = y_global(new_cell_y_min-1:new_cell_y_max+1)

    ! Reallocate the kinetic energy calculation
    DEALLOCATE(ekbar, ekbar_sum, ct)
    ALLOCATE(ekbar(1:nx, 1:ny, 1:n_species))
    ALLOCATE(ekbar_sum(-2:nx+3, -2:ny+3, 1:n_species))
    ALLOCATE(ct(-2:nx+3, -2:ny+3, 1:n_species))

    ! Recalculate x_mins and y_mins so that rebalancing works next time
    DO iproc = 0, nprocx-1
      x_mins(iproc) = x_global(cell_x_min(iproc+1))
      x_maxs(iproc) = x_global(cell_x_max(iproc+1))
    ENDDO
    ! Same for y
    DO iproc = 0, nprocy-1
      y_mins(iproc) = y_global(cell_y_min(iproc+1))
      y_maxs(iproc) = y_global(cell_y_max(iproc+1))
    ENDDO

    ! Set the lengths of the current domain so that the particle balancer
    ! works properly
    x_min_local = x_mins(coordinates(2))
    x_max_local = x_maxs(coordinates(2))
    y_min_local = y_mins(coordinates(1))
    y_max_local = y_maxs(coordinates(1))

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
    INTEGER :: nx_new, ny_new
    INTEGER, DIMENSION(2, 2), INTENT(IN) :: new_domain
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temp
    REAL(num), DIMENSION(:), ALLOCATABLE :: temp1d
    TYPE(laser_block), POINTER :: current

    nx_new = new_domain(1, 2)-new_domain(1, 1)+1
    ny_new = new_domain(2, 2)-new_domain(2, 1)+1

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

#ifndef FULL_LASER_SETTINGS
    ALLOCATE(temp1d(-2:ny_new+3))
    current=>laser_x_min
    DO WHILE(ASSOCIATED(current))
      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain(2, 1), new_domain(2, 2), &
          cell_y_min(coordinates(1)+1), cell_y_max(coordinates(1)+1), &
          ny_global, current%profile, temp1d, c_dir_x)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:ny_new+3))
      current%profile = temp1d
      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain(2, 1), new_domain(2, 2), &
          cell_y_min(coordinates(1)+1), cell_y_max(coordinates(1)+1), &
          ny_global, current%phase, temp1d, c_dir_x)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:ny_new+3))
      current%phase = temp1d

      current=>current%next
    ENDDO

    current=>laser_x_max
    DO WHILE(ASSOCIATED(current))
      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain(2, 1), new_domain(2, 2), &
          cell_y_min(coordinates(1)+1), cell_y_max(coordinates(1)+1), &
          ny_global, current%profile, temp1d, c_dir_x)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:ny_new+3))
      current%profile = temp1d
      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain(2, 1), new_domain(2, 2), &
          cell_y_min(coordinates(1)+1), cell_y_max(coordinates(1)+1), &
          ny_global, current%phase, temp1d, c_dir_x)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:ny_new+3))
      current%phase = temp1d

      current=>current%next
    ENDDO

    ! 1D arrays
    DEALLOCATE(temp1d)
    ALLOCATE(temp1d(-2:nx_new+3))
    current=>laser_y_max
    DO WHILE(ASSOCIATED(current))
      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain(1, 1), new_domain(1, 2), &
          cell_x_min(coordinates(2)+1), cell_x_max(coordinates(2)+1), &
          nx_global, current%profile, temp1d, c_dir_y)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:nx_new+3))
      current%profile = temp1d
      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain(1, 1), new_domain(1, 2), &
          cell_x_min(coordinates(2)+1), cell_x_max(coordinates(2)+1), &
          nx_global, current%phase, temp1d, c_dir_y)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:nx_new+3))
      current%phase = temp1d
      current=>current%next
    ENDDO

    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain(1, 1), new_domain(1, 2), &
          cell_x_min(coordinates(2)+1), cell_x_max(coordinates(2)+1), &
          nx_global, current%profile, temp1d, c_dir_y)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:nx_new+3))
      current%profile = temp1d
      temp1d = 0.0_num
      CALL redistribute_field_1d(new_domain(1, 1), new_domain(1, 2), &
          cell_x_min(coordinates(2)+1), cell_x_max(coordinates(2)+1), &
          nx_global, current%phase, temp1d, c_dir_y)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:nx_new+3))
      current%phase = temp1d
      current=>current%next
    ENDDO
#endif

  END SUBROUTINE redistribute_fields



  SUBROUTINE redistribute_field(domain, field, new_field)

    ! This subroutine redistributes the fields over the new processor layout
    ! The current version works by writing the field to a file and then each
    ! processor loads back in it's own part. This is better than the previous
    ! version where each processor produced it's own copy of the global array
    ! and then took its own subsection
    INTEGER, DIMENSION(2, 2), INTENT(IN) :: domain
    REAL(num), DIMENSION(-2:,-2:), INTENT(IN) :: field
    REAL(num), DIMENSION(-2:,-2:), INTENT(OUT) :: new_field
    INTEGER :: nx_new, ny_new
    INTEGER :: subtype_write, subtype_read, fh
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0
    CHARACTER(LEN=9+data_dir_max_length+n_zeros) :: filename

    WRITE(filename, '(a, "/balance.dat")') TRIM(data_dir)

    nx_new = domain(1, 2)-domain(1, 1)+1
    ny_new = domain(2, 2)-domain(2, 1)+1

    CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_RDWR+MPI_MODE_CREATE, &
        MPI_INFO_NULL, fh, errcode)
    subtype_write = create_current_field_subtype()
    subtype_read  = create_field_subtype(nx_new, ny_new, domain(1, 1), &
        domain(2, 1))

    CALL MPI_FILE_SET_VIEW(fh, offset, mpireal, subtype_write, "native", &
        MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(fh, field(1:nx, 1:ny), nx*ny, mpireal, &
        status, errcode)
    CALL MPI_BARRIER(comm, errcode)
    CALL MPI_FILE_SEEK(fh, offset, MPI_SEEK_SET, errcode)
    CALL MPI_FILE_SET_VIEW(fh, offset, mpireal, subtype_read, "native", &
        MPI_INFO_NULL, errcode)
    CALL MPI_FILE_READ_ALL(fh, new_field(1:nx_new, 1:ny_new), nx_new*ny_new, &
        mpireal, status, errcode)
    CALL MPI_FILE_CLOSE(fh, errcode)
    CALL MPI_BARRIER(comm, errcode)

    CALL MPI_TYPE_FREE(subtype_write, errcode)
    CALL MPI_TYPE_FREE(subtype_read, errcode)

    CALL do_field_mpi_with_lengths(new_field, nx_new, ny_new)

  END SUBROUTINE redistribute_field



  SUBROUTINE redistribute_field_1d(new_start, new_end, old_start, old_end, &
      npts_global, field_in, field_out, direction)

    ! This subroutine redistributes a 1D field over the new processor layout
    ! The current version works by producing a global copy on each processor
    ! And then extracting the required part for the local processor.
    ! in 1D, this is probably OK
    INTEGER, INTENT(IN) :: new_start, new_end, npts_global
    INTEGER, INTENT(IN) :: old_start, old_end, direction
    REAL(num), DIMENSION(-2:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(-2:), INTENT(OUT) :: field_out
    REAL(num), DIMENSION(:), ALLOCATABLE :: field_new, field_temp
    INTEGER :: new_pts, old_pts
    INTEGER :: new_comm, color

    IF (IAND(direction, c_dir_x) .EQ. 0) color = coordinates(1)
    IF (IAND(direction, c_dir_y) .EQ. 0) color = coordinates(2)

    CALL MPI_COMM_SPLIT(comm, color, rank, new_comm, errcode)

    new_pts = new_end-new_start+1
    old_pts = old_end-old_start+1

    ! Create a global copy of the whole array
    ALLOCATE(field_new(1:npts_global), field_temp(0:npts_global+1))
    field_new = 0.0_num
    field_new(old_start:old_end) = field_in(1:old_pts)
    CALL MPI_ALLREDUCE(field_new, field_temp(1:npts_global), npts_global, &
        mpireal, MPI_SUM, new_comm, errcode)
    DEALLOCATE(field_new)

    field_out(0:new_pts+1) = field_temp(new_start-1:new_end+1)
    DEALLOCATE(field_temp)
    CALL MPI_COMM_FREE(new_comm, errcode)

  END SUBROUTINE redistribute_field_1d



  SUBROUTINE get_density_in_x(density)

    ! Calculate total particle density across the X direction
    ! Summed in the Y direction

    INTEGER(KIND=8), DIMENSION(:), INTENT(INOUT) :: density
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: temp
    TYPE(particle), POINTER :: current
    REAL(num) :: part_x
    INTEGER :: cell_x1, ispecies

    density = 0.0_num

    DO ispecies = 1, n_species
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so not x_min, NOT x_min_local
        part_x = current%part_pos(1)-x_min
        cell_x1 = NINT(part_x/dx)+1
        density(cell_x1) = density(cell_x1)+1
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
    ! Summed in the X direction

    INTEGER(KIND=8), DIMENSION(:), INTENT(INOUT) :: density
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: temp
    TYPE(particle), POINTER :: current
    REAL(num) :: part_y
    INTEGER :: cell_y1, ispecies

    density = 0.0_num

    DO ispecies = 1, n_species
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so not x_min, NOT x_min_local
        part_y = current%part_pos(2)-y_min
        cell_y1 = NINT(part_y/dy)+1
        density(cell_y1) = density(cell_y1)+1
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
    sz = SIZE(density)-2
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
          .OR. ABS(total + density(idim) -npart_per_proc_ideal) &
          .GT. ABS(total-npart_per_proc_ideal)  .OR. idim .EQ. sz) THEN
        total = density(idim)
        starts(partition) = idim-1
        partition = partition+1
        ! If you've reached the last processor, have already done the best
        ! you can, so just leave
      ELSE
        total = total+density(idim)
      ENDIF
    ENDDO

    ends(nproc) = sz
    DO idim = 1, nproc-1
      ends(idim) = starts(idim+1)-1
    ENDDO

  END SUBROUTINE calculate_breaks



  FUNCTION get_particle_processor(a_particle)

    ! This subroutine calculates which processor a given particles resides on

    TYPE(particle), INTENT(IN) :: a_particle
    INTEGER :: get_particle_processor
    INTEGER :: iproc, coords(2)

    get_particle_processor = -1
    coords = -1

    ! This could be replaced by a bisection method, but for the moment I
    ! just don't care

    DO iproc = 0, nprocx-1
      IF (a_particle%part_pos(1) .GE. x_mins(iproc) - dx/2.0_num &
          .AND. a_particle%part_pos(1) .LE. x_maxs(iproc) + dx/2.0_num) THEN
        coords(2) = iproc
        EXIT
      ENDIF
    ENDDO

    DO iproc = 0, nprocy-1
      IF (a_particle%part_pos(2) .GE. y_mins(iproc) -dy/2.0_num &
          .AND. a_particle%part_pos(2) .LE. y_maxs(iproc) + dy/2.0_num) THEN
        coords(1) = iproc
        EXIT
      ENDIF
    ENDDO

    IF (MINVAL(coords) .LT. 0) PRINT *, "UNLOCATABLE PARTICLE"
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
    INTEGER :: part_proc, iproc_recv, iproc_send, ispecies

    ALLOCATE(pointers_send(0:nproc-1), pointers_recv(0:nproc-1))
    DO ispecies = 1, n_species
      current=>particle_species(ispecies)%attached_list%head
      DO iproc_send = 0, nproc-1
        CALL create_empty_partlist(pointers_send(iproc_send))
        CALL create_empty_partlist(pointers_recv(iproc_send))
      ENDDO

      DO WHILE(ASSOCIATED(current))
        next=>current%next
        part_proc = get_particle_processor(current)
        IF (part_proc .LT. 0) THEN
          PRINT *, "Unlocatable particle on processor", rank, &
              current%part_pos, dx, dy
          CALL MPI_BARRIER(comm, errcode)
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

      DO iproc_send = 0, nproc-1
        DO iproc_recv = 0, nproc-1
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

      DO iproc_recv = 0, nproc-1
        CALL append_partlist(particle_species(ispecies)%attached_list, &
            pointers_recv(iproc_recv))
      ENDDO
    ENDDO

    DEALLOCATE(pointers_send, pointers_recv)

  END SUBROUTINE distribute_particles

END MODULE balance
