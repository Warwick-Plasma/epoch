MODULE balance

  USE mpi
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
    INTEGER, DIMENSION(:), ALLOCATABLE :: starts_x, ends_x
    REAL(num) :: balance_frac, npart_av
    INTEGER(KIND=8) :: npart_local, sum_npart, max_npart
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

    IF (ALLOCATED(load_x)) DEALLOCATE(load_x)

    ! Now need to calculate the start and end points for the new domain on
    ! the current processor
    nx_global_min = starts_x(x_coords+1)
    nx_global_max = ends_x(x_coords+1)

    domain(1,:) = (/nx_global_min, nx_global_max/)

    ! Redistribute the field variables
    CALL redistribute_fields(domain)

    ! Copy the new lengths into the permanent variables
    cell_x_min = starts_x
    cell_x_max = ends_x

    ! Set the new nx
    nx = nx_global_max - nx_global_min + 1

    ! Do X array separately because we already have global copies
    DEALLOCATE(x)
    ALLOCATE(x(-2:nx+3))
    x(-2:nx+3) = x_global(nx_global_min-3:nx_global_max+3)

    ! Recalculate x_mins and x_maxs so that rebalancing works next time
    DO iproc = 0, nprocx - 1
      x_mins(iproc) = x_global(cell_x_min(iproc+1))
      x_maxs(iproc) = x_global(cell_x_max(iproc+1))
    ENDDO

    ! Set the lengths of the current domain so that the particle balancer
    ! works properly
    x_min_local = x_mins(x_coords)
    x_max_local = x_maxs(x_coords)

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

    INTEGER :: nx_new
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: new_domain
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temp2d
    REAL(num), DIMENSION(:), ALLOCATABLE :: temp
    INTEGER :: ispecies, index, n_species_local

    nx_new = new_domain(1,2) - new_domain(1,1) + 1

    ALLOCATE(temp(-2:nx_new+3))

    temp = 0.0_num
    CALL redistribute_field(new_domain, ex, temp)
    DEALLOCATE(ex)
    ALLOCATE(ex(-2:nx_new+3))
    ex = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, ey, temp)
    DEALLOCATE(ey)
    ALLOCATE(ey(-2:nx_new+3))
    ey = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, ez, temp)
    DEALLOCATE(ez)
    ALLOCATE(ez(-2:nx_new+3))
    ez = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, bx, temp)
    DEALLOCATE(bx)
    ALLOCATE(bx(-2:nx_new+3))
    bx = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, by, temp)
    DEALLOCATE(by)
    ALLOCATE(by(-2:nx_new+3))
    by = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, bz, temp)
    DEALLOCATE(bz)
    ALLOCATE(bz(-2:nx_new+3))
    bz = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, jx, temp)
    DEALLOCATE(jx)
    ALLOCATE(jx(-2:nx_new+3))
    jx = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, jy, temp)
    DEALLOCATE(jy)
    ALLOCATE(jy(-2:nx_new+3))
    jy = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, jz, temp)
    DEALLOCATE(jz)
    ALLOCATE(jz(-2:nx_new+3))
    jz = temp

    IF (cpml_boundaries) THEN
      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_e_psiyx, temp)
      DEALLOCATE(cpml_e_psiyx)
      ALLOCATE(cpml_e_psiyx(-2:nx_new+3))
      cpml_e_psiyx = temp

      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_e_psizx, temp)
      DEALLOCATE(cpml_e_psizx)
      ALLOCATE(cpml_e_psizx(-2:nx_new+3))
      cpml_e_psizx = temp

      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_b_psiyx, temp)
      DEALLOCATE(cpml_b_psiyx)
      ALLOCATE(cpml_b_psiyx(-2:nx_new+3))
      cpml_b_psiyx = temp

      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_b_psizx, temp)
      DEALLOCATE(cpml_b_psizx)
      ALLOCATE(cpml_b_psizx(-2:nx_new+3))
      cpml_b_psizx = temp

      CALL deallocate_cpml_helpers
      CALL set_cpml_helpers
    ENDIF

    DEALLOCATE(temp)

    DO index = 1, num_vars_to_dump
      IF (.NOT. ASSOCIATED(averaged_data(index)%array)) CYCLE

      n_species_local = 0
      IF (IAND(dumpmask(index), c_io_no_sum) .EQ. 0) &
          n_species_local = 1
      IF (IAND(dumpmask(index), c_io_species) .NE. 0) &
          n_species_local = n_species_local + n_species

      IF (n_species_local .LE. 0) CYCLE

      ALLOCATE(temp2d(-2:nx_new+3, 1:n_species_local))

      DO ispecies = 1, n_species_local
        CALL redistribute_field(new_domain, &
            averaged_data(index)%array(:,ispecies), temp2d(:,ispecies))
      ENDDO

      DEALLOCATE(averaged_data(index)%array)
      ALLOCATE(averaged_data(index)%array(-2:nx_new+3,n_species_local))

      averaged_data(index)%array = temp2d

      DEALLOCATE(temp2d)
    ENDDO

    ! No need to rebalance lasers in 1D, lasers are just a point!
    ! Same goes for moving windows and thermal boundaries.

  END SUBROUTINE redistribute_fields



  SUBROUTINE redistribute_field(domain, field_in, field_out)

    ! This subroutine redistributes a 1D field over the new processor layout
    ! The current version works by producing a global copy on each processor
    ! And then extracting the required part for the local processor.
    ! in 1D, this is probably OK
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: domain
    REAL(num), DIMENSION(-2:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(-2:), INTENT(OUT) :: field_out
    REAL(num), DIMENSION(:), ALLOCATABLE :: field_new, field_temp
    INTEGER :: old_start, new_start
    INTEGER :: old_pts, new_pts, dir
    INTEGER :: npts_global, new_comm, color, coord, i
    INTEGER :: ghost_start, ghost_end

    ghost_start = 0
    ghost_end = 0

    dir = 1
    color = 1
    coord = x_coords
    old_start = cell_x_min(x_coords+1)
    old_pts = nx
    npts_global = nx_global
    IF (x_max_boundary) ghost_end = 3

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

  END SUBROUTINE redistribute_field



  SUBROUTINE get_load_in_x(load)

    ! Calculate total load across the X direction

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
        cell = FLOOR((current%part_pos - x_min) / dx) + 1
#else
        cell = FLOOR((current%part_pos - x_min) / dx + 1.5_num)
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
    load = push_per_field * temp + 1

    DEALLOCATE(temp)

  END SUBROUTINE get_load_in_x



  SUBROUTINE calculate_breaks(load, nproc, starts, ends)

    ! This subroutine calculates the places in a given load profile to split
    ! The domain to give the most even subdivision possible

    INTEGER(KIND=8), INTENT(IN), DIMENSION(:) :: load
    INTEGER, INTENT(IN) :: nproc
    INTEGER, DIMENSION(:), INTENT(OUT) :: starts, ends
    INTEGER :: sz, idim, proc, old
    INTEGER(KIND=8) :: total, total_old, load_per_proc_ideal

    sz = SIZE(load)
    ends = sz

    load_per_proc_ideal = FLOOR((SUM(load) + 0.5d0) / nproc, 8)

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
    old = 0
    DO proc = 1, nproc
      IF (ends(proc) - old .LE. 0) THEN
        ends(proc) = old + 1
      ENDIF
      old = ends(proc)
    ENDDO

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
      IF (a_particle%part_pos .GE. x_mins(iproc) - dx / 2.0_num &
          .AND. a_particle%part_pos .LT. x_maxs(iproc) + dx / 2.0_num) THEN
        coords(c_ndims) = iproc
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
    INTEGER :: part_proc, iproc, ispecies, i, oldrank, newrank, nodd, nevn, ierr
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: sendcounts, recvcounts
    INTEGER, DIMENSION(:), ALLOCATABLE :: oddlist, evnlist

    nodd = (nproc + 1) / 2
    nevn = nproc / 2

    ALLOCATE(pointers_send(0:nproc-1), pointers_recv(0:nproc-1))
    ALLOCATE(sendcounts(0:nproc-1), recvcounts(0:nproc-1))
    ALLOCATE(oddlist(0:nodd-1), evnlist(0:nevn-1))

    DO ispecies = 1, n_species
      current=>species_list(ispecies)%attached_list%head
      DO iproc = 0, nproc - 1
        CALL create_empty_partlist(pointers_send(iproc))
        CALL create_empty_partlist(pointers_recv(iproc))
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

      DO iproc = 0, nproc - 1
        sendcounts(iproc) = pointers_send(iproc)%count
      ENDDO

      ! Split the processors into odd and even lists.
      ! Even processors send then receive, Odd processors vice-versa.
      ! Next, on even processors split the even list in two and on
      ! odd processors split the odd list in two.
      ! Repeat until the list size is equal to one.

      ! If the number of processors is not divisible by two then the
      ! odd list has one extra entry.

      nodd = (nproc + 1) / 2
      nevn = nproc / 2
      DO i = 0, nevn - 1
        oddlist(i) = 2 * i
        evnlist(i) = 2 * i + 1
        IF (oddlist(i) .EQ. rank .OR. evnlist(i) .EQ. rank) THEN
          newrank = i
        ENDIF
      ENDDO

      IF (nodd .NE. nevn) oddlist(nodd-1) = nproc - 1

      CALL MPI_ALLTOALL(sendcounts, 1, MPI_INTEGER8, recvcounts, 1, &
          MPI_INTEGER8, comm, errcode)

      oldrank = rank
      DO WHILE(nevn .GT. 0)
        IF (MOD(oldrank,2) .EQ. 0) THEN
          oldrank = newrank
          DO i = 0, nevn - 1
            iproc = evnlist(i)
            IF (sendcounts(iproc) .GT. 0) THEN
              CALL partlist_send_nocount(pointers_send(iproc), iproc)
              CALL destroy_partlist(pointers_send(iproc))
            ENDIF
          ENDDO

          DO i = 0, nevn - 1
            iproc = evnlist(i)
            IF (iproc .GE. nproc) EXIT
            IF (recvcounts(iproc) .GT. 0) THEN
              CALL partlist_recv_nocount(pointers_recv(iproc), iproc, &
                  recvcounts(iproc))
            ENDIF
          ENDDO

          nevn = nodd / 2
          nodd = (nodd + 1) / 2
          DO i = 0, nevn - 1
            oddlist(i) = oddlist(2*i)
            evnlist(i) = oddlist(2*i+1)
            IF (oddlist(i) .EQ. rank .OR. evnlist(i) .EQ. rank) THEN
              newrank = i
            ENDIF
          ENDDO
          IF (nodd .NE. nevn) oddlist(nodd-1) = oddlist(2*nodd-2)

        ELSE
          oldrank = newrank
          DO i = 0, nodd - 1
            iproc = oddlist(i)
            IF (iproc .LT. 0) EXIT
            IF (recvcounts(iproc) .GT. 0) THEN
              CALL partlist_recv_nocount(pointers_recv(iproc), iproc, &
                  recvcounts(iproc))
            ENDIF
          ENDDO

          DO i = 0, nodd - 1
            iproc = oddlist(i)
            IF (iproc .LT. 0) EXIT
            IF (sendcounts(iproc) .GT. 0) THEN
              CALL partlist_send_nocount(pointers_send(iproc), iproc)
              CALL destroy_partlist(pointers_send(iproc))
            ENDIF
          ENDDO

          nodd = (nevn + 1) / 2
          nevn = nevn / 2
          DO i = 0, nevn - 1
            oddlist(i) = evnlist(2*i)
            evnlist(i) = evnlist(2*i+1)
            IF (oddlist(i) .EQ. rank .OR. evnlist(i) .EQ. rank) THEN
              newrank = i
            ENDIF
          ENDDO
          IF (nodd .NE. nevn) oddlist(nodd-1) = evnlist(2*nodd-2)
        ENDIF
      ENDDO

      DO iproc = 0, nproc - 1
        CALL append_partlist(species_list(ispecies)%attached_list, &
            pointers_recv(iproc))
      ENDDO
    ENDDO

    DEALLOCATE(sendcounts, recvcounts, oddlist, evnlist)
    DEALLOCATE(pointers_send, pointers_recv)

  END SUBROUTINE distribute_particles

END MODULE balance
