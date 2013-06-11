MODULE balance

  USE boundary
  USE mpi_subtype_control
  USE redblack_module

  IMPLICIT NONE

  INTEGER, DIMENSION(:), ALLOCATABLE :: new_cell_x_min, new_cell_x_max
  LOGICAL :: overriding

CONTAINS

  SUBROUTINE balance_workload(over_ride)

    ! This subroutine determines whether or not the code needs rebalancing,
    ! calculates where to split the domain and calls other subroutines to
    ! actually rearrange the fields and particles onto the new processors

    ! This is really, really hard to do properly
    ! So cheat

    LOGICAL, INTENT(IN) :: over_ride
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: load_x
    REAL(num) :: balance_frac, npart_av
    INTEGER(i8) :: npart_local, sum_npart, max_npart
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
      IF (rank .EQ. 0) PRINT *, 'Load balancing with fraction', balance_frac
    ENDIF

    IF (.NOT.use_exact_restart) THEN
      overriding = over_ride

      ALLOCATE(new_cell_x_min(nprocx), new_cell_x_max(nprocx))

      new_cell_x_min = cell_x_min
      new_cell_x_max = cell_x_max

      ! Sweep in X
      IF (nprocx .GT. 1) THEN
        IF (IAND(balance_mode, c_lb_x) .NE. 0 &
            .OR. IAND(balance_mode, c_lb_auto) .NE. 0) THEN
          ! Rebalancing in X
          ALLOCATE(load_x(nx_global))
          CALL get_load_in_x(load_x)
          CALL calculate_breaks(load_x, nprocx, new_cell_x_min, new_cell_x_max)
        ENDIF
      ENDIF

      IF (ALLOCATED(load_x)) DEALLOCATE(load_x)

      ! Now need to calculate the start and end points for the new domain on
      ! the current processor

      domain(1,:) = (/new_cell_x_min(x_coords+1), new_cell_x_max(x_coords+1)/)

      ! Redistribute the field variables
      CALL redistribute_fields(domain)

      ! Copy the new lengths into the permanent variables
      cell_x_min = new_cell_x_min
      cell_x_max = new_cell_x_max

      ! Set the new nx
      nx_global_min = cell_x_min(x_coords+1)
      nx_global_max = cell_x_max(x_coords+1)
      n_global_min(1) = nx_global_min
      n_global_max(1) = nx_global_max

      nx = nx_global_max - nx_global_min + 1

      DEALLOCATE(new_cell_x_min, new_cell_x_max)

      ! Do X array separately because we already have global copies
      DEALLOCATE(x)
      ALLOCATE(x(-2:nx+3))
      x(-2:nx+3) = x_global(nx_global_min-3:nx_global_max+3)

      ! Recalculate x_grid_mins/maxs so that rebalancing works next time
      DO iproc = 0, nprocx - 1
        x_grid_mins(iproc) = x_global(cell_x_min(iproc+1))
        x_grid_maxs(iproc) = x_global(cell_x_max(iproc+1))
      ENDDO

      ! Set the lengths of the current domain so that the particle balancer
      ! works properly
      x_grid_min_local = x_grid_mins(x_coords)
      x_grid_max_local = x_grid_maxs(x_coords)

      x_min_local = x_grid_min_local + (cpml_x_min_offset - 0.5_num) * dx
      x_max_local = x_grid_max_local - (cpml_x_max_offset - 0.5_num) * dx
    ENDIF

    ! Redistribute the particles onto their new processors
    CALL distribute_particles

    ! If running with particle debugging then set the t = 0 processor if
    ! over_ride = true
#ifdef PARTICLE_DEBUG
    IF (over_ride) THEN
      DO ispecies = 1, n_species
        current => species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          current%processor_at_t0 = rank
          current => current%next
        ENDDO
      ENDDO
    ENDIF
#endif

    use_exact_restart = .FALSE.

  END SUBROUTINE balance_workload



  SUBROUTINE redistribute_fields(new_domain)

    ! This subroutine redistributes the field variables over the new
    ! processor layout. If using a field of your own then set the
    ! redistribute_field subroutine to implement it.

    INTEGER :: nx_new
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: new_domain
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temp_sum
    REAL(r4), DIMENSION(:,:), ALLOCATABLE :: r4temp_sum
    REAL(num), DIMENSION(:), ALLOCATABLE :: temp, temp2
    INTEGER :: i, ispecies, io, id, nspec_local, mask

    nx_new = new_domain(1,2) - new_domain(1,1) + 1

    ! The following code is quite messy and repetitive. Unfortunately, the
    ! F90 standard does not allow the ALLOCATABLE attribute for subroutine
    ! arguments and POINTER arrays are not as fast.

    ! Full domain arrays

    ALLOCATE(temp(-2:nx_new+3))

    ! Current will be recalculated during the particle push, so there
    ! is no need to copy the contents of the old arrays.
    ! If overriding, then we may not be doing a particle push next
    ! so we still have to balance the arrays.
    ! It is done slightly differently since the arrays may be
    ! a different size.

    IF (overriding) THEN
      ALLOCATE(temp2(-2:nx+3))

      temp2(0:nx+1) = jx(0:nx+1)
      CALL remap_field(temp2, temp)
      DEALLOCATE(jx)
      ALLOCATE(jx(1-jng:nx_new+jng))
      jx(0:nx_new+1) = temp(0:nx_new+1)

      temp2(0:nx+1) = jy(0:nx+1)
      CALL remap_field(temp2, temp)
      DEALLOCATE(jy)
      ALLOCATE(jy(1-jng:nx_new+jng))
      jy(0:nx_new+1) = temp(0:nx_new+1)

      temp2(0:nx+1) = jz(0:nx+1)
      CALL remap_field(temp2, temp)
      DEALLOCATE(jz)
      ALLOCATE(jz(1-jng:nx_new+jng))
      jz(0:nx_new+1) = temp(0:nx_new+1)

      DEALLOCATE(temp2)
    ELSE
      DEALLOCATE(jx)
      DEALLOCATE(jy)
      DEALLOCATE(jz)
      ALLOCATE(jx(1-jng:nx_new+jng))
      ALLOCATE(jy(1-jng:nx_new+jng))
      ALLOCATE(jz(1-jng:nx_new+jng))
    ENDIF

    CALL remap_field(ex, temp)
    DEALLOCATE(ex)
    ALLOCATE(ex(-2:nx_new+3))
    ex = temp

    CALL remap_field(ey, temp)
    DEALLOCATE(ey)
    ALLOCATE(ey(-2:nx_new+3))
    ey = temp

    CALL remap_field(ez, temp)
    DEALLOCATE(ez)
    ALLOCATE(ez(-2:nx_new+3))
    ez = temp

    CALL remap_field(bx, temp)
    DEALLOCATE(bx)
    ALLOCATE(bx(-2:nx_new+3))
    bx = temp

    CALL remap_field(by, temp)
    DEALLOCATE(by)
    ALLOCATE(by(-2:nx_new+3))
    by = temp

    CALL remap_field(bz, temp)
    DEALLOCATE(bz)
    ALLOCATE(bz(-2:nx_new+3))
    bz = temp

    DO ispecies = 1, n_species
      IF (species_list(ispecies)%migrate%fluid) THEN
        CALL remap_field(species_list(ispecies)%migrate%fluid_energy, temp)
        DEALLOCATE(species_list(ispecies)%migrate%fluid_energy)
        ALLOCATE(species_list(ispecies)%migrate%fluid_energy(-2:nx_new+3))
        species_list(ispecies)%migrate%fluid_energy = temp

        CALL remap_field(species_list(ispecies)%migrate%fluid_density, temp)
        DEALLOCATE(species_list(ispecies)%migrate%fluid_density)
        ALLOCATE(species_list(ispecies)%migrate%fluid_density(-2:nx_new+3))
        species_list(ispecies)%migrate%fluid_density = temp
      ENDIF
    ENDDO

    IF (cpml_boundaries) THEN
      CALL remap_field(cpml_psi_eyx, temp)
      DEALLOCATE(cpml_psi_eyx)
      ALLOCATE(cpml_psi_eyx(-2:nx_new+3))
      cpml_psi_eyx = temp

      CALL remap_field(cpml_psi_byx, temp)
      DEALLOCATE(cpml_psi_byx)
      ALLOCATE(cpml_psi_byx(-2:nx_new+3))
      cpml_psi_byx = temp

      CALL remap_field(cpml_psi_ezx, temp)
      DEALLOCATE(cpml_psi_ezx)
      ALLOCATE(cpml_psi_ezx(-2:nx_new+3))
      cpml_psi_ezx = temp

      CALL remap_field(cpml_psi_bzx, temp)
      DEALLOCATE(cpml_psi_bzx)
      ALLOCATE(cpml_psi_bzx(-2:nx_new+3))
      cpml_psi_bzx = temp

      CALL deallocate_cpml_helpers
      CALL set_cpml_helpers(nx_new, new_domain(1,1), new_domain(1,2))
    ENDIF

    DEALLOCATE(temp)

    ! Full domain arrays with an additional index

    DO id = 1, num_vars_to_dump
      io = averaged_var_block(id)
      IF (io .EQ. 0) CYCLE

      mask = io_block_list(io)%dumpmask(id)
      nspec_local = 0
      IF (IAND(mask, c_io_no_sum) .EQ. 0) &
          nspec_local = 1
      IF (IAND(mask, c_io_species) .NE. 0) &
          nspec_local = nspec_local + n_species

      IF (nspec_local .LE. 0) CYCLE

      IF (io_block_list(io)%averaged_data(id)%dump_single) THEN
        IF (.NOT. ASSOCIATED(io_block_list(io)%averaged_data(id)%r4array)) CYCLE

        ALLOCATE(r4temp_sum(-2:nx_new+3, nspec_local))

        DO i = 1, nspec_local
          CALL remap_field_r4(&
              io_block_list(io)%averaged_data(id)%r4array(:,i), &
              r4temp_sum(:,i))
        ENDDO

        DEALLOCATE(io_block_list(io)%averaged_data(id)%r4array)
        ALLOCATE(io_block_list(io)%averaged_data(id)&
            %r4array(-2:nx_new+3, nspec_local))

        io_block_list(io)%averaged_data(id)%r4array = r4temp_sum

        DEALLOCATE(r4temp_sum)
      ELSE
        IF (.NOT. ASSOCIATED(io_block_list(io)%averaged_data(id)%array)) CYCLE

        ALLOCATE(temp_sum(-2:nx_new+3, nspec_local))

        DO i = 1, nspec_local
          CALL remap_field(&
              io_block_list(io)%averaged_data(id)%array(:,i), &
              temp_sum(:,i))
        ENDDO

        DEALLOCATE(io_block_list(io)%averaged_data(id)%array)
        ALLOCATE(io_block_list(io)%averaged_data(id)&
            %array(-2:nx_new+3, nspec_local))

        io_block_list(io)%averaged_data(id)%array = temp_sum

        DEALLOCATE(temp_sum)
      ENDIF
    ENDDO

  END SUBROUTINE redistribute_fields



  SUBROUTINE remap_field(field_in, field_out)

    ! This is a wrapper for the field redistribution routine
    REAL(num), DIMENSION(:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(c_ndims) :: n_new, cdim
    INTEGER :: i

    n_new = SHAPE(field_out) - 2 * 3

    DO i = 1, c_ndims
      cdim(i) = c_ndims + 1 - i
    ENDDO

    CALL redistribute_field_1d(field_in, field_out, cdim, &
        cell_x_min, cell_x_max, new_cell_x_min, new_cell_x_max)

    CALL do_field_mpi_with_lengths(field_out, ng, n_new(1))

  END SUBROUTINE remap_field



  SUBROUTINE remap_field_r4(field_in, field_out)

    ! This is a wrapper for the field redistribution routine
    REAL(r4), DIMENSION(:), INTENT(IN) :: field_in
    REAL(r4), DIMENSION(:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(c_ndims) :: n_new, cdim
    INTEGER :: i

    n_new = SHAPE(field_out) - 2 * 3

    DO i = 1, c_ndims
      cdim(i) = c_ndims + 1 - i
    ENDDO

    CALL redistribute_field_1d_r4(field_in, field_out, cdim, &
        cell_x_min, cell_x_max, new_cell_x_min, new_cell_x_max)

    CALL do_field_mpi_with_lengths_r4(field_out, ng, n_new(1))

  END SUBROUTINE remap_field_r4



  SUBROUTINE redistribute_field_1d(field_in, field_out, cdim, &
      old_cell_min1, old_cell_max1, new_cell_min1, new_cell_max1)

    ! This subroutine redistributes the fields over the new processor layout
    INTEGER, PARAMETER :: nd = 1
    REAL(num), DIMENSION(1-ng:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(1-ng:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(nd), INTENT(IN) :: cdim
    INTEGER, DIMENSION(:), INTENT(IN) :: old_cell_min1, old_cell_max1
    INTEGER, DIMENSION(:), INTENT(IN) :: new_cell_min1, new_cell_max1
    INTEGER :: irank, basetype, n, ng0, ng1
    INTEGER :: i, iproc, inew
    INTEGER, DIMENSION(nd) :: type_min, type_max, old_0, old_1, new_0
    INTEGER, DIMENSION(nd) :: n_global, n_local, start, nprocs
    INTEGER, DIMENSION(nd) :: old_min, old_max, new_min, new_max
    INTEGER, DIMENSION(c_ndims) :: coord
    INTEGER, DIMENSION(nd) :: our_coords, nmin, nmax
    INTEGER, DIMENSION(:), ALLOCATABLE :: sendtypes, recvtypes

    basetype = mpireal

    ALLOCATE(sendtypes(0:nproc-1))
    ALLOCATE(recvtypes(0:nproc-1))

    DO i = 1, nd
      our_coords(i) = coordinates(cdim(i))
    ENDDO

    nprocs(1) = SIZE(old_cell_min1)

    old_min(1) = old_cell_min1(our_coords(1)+1)
    old_max(1) = old_cell_max1(our_coords(1)+1)
    new_min(1) = new_cell_min1(our_coords(1)+1)
    new_max(1) = new_cell_max1(our_coords(1)+1)

    nmin(1) = new_cell_min1(1)
    nmax(1) = new_cell_max1(nprocs(1))

    tag = 0
    sendtypes = 0
    recvtypes = 0
    coord = coordinates

    ! Create array of sendtypes

    DO i = 1,nd
      n_global(i) = old_max(i) - old_min(i) + 2 * ng + 1
    ENDDO

    n = 1
    type_min(n) = old_min(n)
    type_max(n) = old_min(n)

    ! Find the new processor on which the old x_min resides
    ! This could be sped up by using bisection.
    DO iproc = 1, nprocs(n)-1
      IF (new_cell_min1(iproc) .LE. old_min(n) &
          .AND. new_cell_max1(iproc) .GE. old_min(n)) EXIT
    ENDDO

    DO WHILE(type_max(n) .LE. old_max(n))
      coord(cdim(n)) = iproc - 1
      type_max(n) = new_cell_max1(iproc)
      IF (type_max(n) .GT. old_max(n)) type_max(n) = old_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) .EQ. nmin(n)) ng0 = ng
      IF (type_max(n) .EQ. nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - old_min(n) + ng - ng0 + 1

      CALL MPI_CART_RANK(comm, coord, irank, errcode)

      IF (rank .NE. irank) THEN
        sendtypes(irank) = create_1d_array_subtype(basetype, n_local, &
            n_global, start)
      ELSE
        ! New domain is on the same processor as the old domain.
        ! Just copy the region rather than using MPI.
        DO i = 1,nd
          old_0(i) = start(i) - ng
          old_1(i) = old_0(i) + n_local(i) - 1
        ENDDO
      ENDIF

      n = 1
      IF (type_max(n) .EQ. old_max(n)) EXIT
      iproc = iproc + 1
      type_min(n) = new_cell_min1(iproc)
    ENDDO

    ! Create array of recvtypes

    DO i = 1,nd
      n_global(i) = new_max(i) - new_min(i) + 2 * ng + 1
    ENDDO

    n = 1
    type_min(n) = new_min(n)
    type_max(n) = new_min(n)

    ! Find the old processor on which the new x_min resides
    ! This could be sped up by using bisection.
    DO iproc = 1, nprocs(n)-1
      IF (old_cell_min1(iproc) .LE. new_min(n) &
          .AND. old_cell_max1(iproc) .GE. new_min(n)) EXIT
    ENDDO

    DO WHILE(type_max(n) .LE. new_max(n))
      coord(cdim(n)) = iproc - 1
      type_max(n) = old_cell_max1(iproc)
      IF (type_max(n) .GT. new_max(n)) type_max(n) = new_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) .EQ. nmin(n)) ng0 = ng
      IF (type_max(n) .EQ. nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - new_min(n) + ng - ng0 + 1

      CALL MPI_CART_RANK(comm, coord, irank, errcode)

      IF (rank .NE. irank) THEN
        recvtypes(irank) = create_1d_array_subtype(basetype, n_local, &
            n_global, start)
      ELSE
        ! New domain is on the same processor as the old domain.
        ! Just copy the region rather than using MPI.
        DO i = 1,nd
          new_0(i) = start(i) - ng
        ENDDO
        DO i = old_0(1),old_1(1)
          inew = new_0(1) + i - old_0(1)
          field_out(inew) = field_in(i)
        ENDDO
      ENDIF

      n = 1
      IF (type_max(n) .EQ. new_max(n)) EXIT
      iproc = iproc + 1
      type_min(n) = old_cell_min1(iproc)
    ENDDO

    CALL redblack(field_in, field_out, sendtypes, recvtypes)

    DO i = 0,nproc-1
      IF (sendtypes(i) .NE. 0) CALL MPI_TYPE_FREE(sendtypes(i), errcode)
      IF (recvtypes(i) .NE. 0) CALL MPI_TYPE_FREE(recvtypes(i), errcode)
    ENDDO

    DEALLOCATE(sendtypes)
    DEALLOCATE(recvtypes)

  END SUBROUTINE redistribute_field_1d



  SUBROUTINE redistribute_field_1d_r4(field_in, field_out, cdim, &
      old_cell_min1, old_cell_max1, new_cell_min1, new_cell_max1)

    ! This subroutine redistributes the fields over the new processor layout
    INTEGER, PARAMETER :: nd = 1
    REAL(r4), DIMENSION(1-ng:), INTENT(IN) :: field_in
    REAL(r4), DIMENSION(1-ng:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(nd), INTENT(IN) :: cdim
    INTEGER, DIMENSION(:), INTENT(IN) :: old_cell_min1, old_cell_max1
    INTEGER, DIMENSION(:), INTENT(IN) :: new_cell_min1, new_cell_max1
    INTEGER :: irank, basetype, n, ng0, ng1
    INTEGER :: i, iproc, inew
    INTEGER, DIMENSION(nd) :: type_min, type_max, old_0, old_1, new_0
    INTEGER, DIMENSION(nd) :: n_global, n_local, start, nprocs
    INTEGER, DIMENSION(nd) :: old_min, old_max, new_min, new_max
    INTEGER, DIMENSION(c_ndims) :: coord
    INTEGER, DIMENSION(nd) :: our_coords, nmin, nmax
    INTEGER, DIMENSION(:), ALLOCATABLE :: sendtypes, recvtypes

    basetype = MPI_REAL4

    ALLOCATE(sendtypes(0:nproc-1))
    ALLOCATE(recvtypes(0:nproc-1))

    DO i = 1, nd
      our_coords(i) = coordinates(cdim(i))
    ENDDO

    nprocs(1) = SIZE(old_cell_min1)

    old_min(1) = old_cell_min1(our_coords(1)+1)
    old_max(1) = old_cell_max1(our_coords(1)+1)
    new_min(1) = new_cell_min1(our_coords(1)+1)
    new_max(1) = new_cell_max1(our_coords(1)+1)

    nmin(1) = new_cell_min1(1)
    nmax(1) = new_cell_max1(nprocs(1))

    tag = 0
    sendtypes = 0
    recvtypes = 0
    coord = coordinates

    ! Create array of sendtypes

    DO i = 1,nd
      n_global(i) = old_max(i) - old_min(i) + 2 * ng + 1
    ENDDO

    n = 1
    type_min(n) = old_min(n)
    type_max(n) = old_min(n)

    ! Find the new processor on which the old x_min resides
    ! This could be sped up by using bisection.
    DO iproc = 1, nprocs(n)-1
      IF (new_cell_min1(iproc) .LE. old_min(n) &
          .AND. new_cell_max1(iproc) .GE. old_min(n)) EXIT
    ENDDO

    DO WHILE(type_max(n) .LE. old_max(n))
      coord(cdim(n)) = iproc - 1
      type_max(n) = new_cell_max1(iproc)
      IF (type_max(n) .GT. old_max(n)) type_max(n) = old_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) .EQ. nmin(n)) ng0 = ng
      IF (type_max(n) .EQ. nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - old_min(n) + ng - ng0 + 1

      CALL MPI_CART_RANK(comm, coord, irank, errcode)

      IF (rank .NE. irank) THEN
        sendtypes(irank) = create_1d_array_subtype(basetype, n_local, &
            n_global, start)
      ELSE
        ! New domain is on the same processor as the old domain.
        ! Just copy the region rather than using MPI.
        DO i = 1,nd
          old_0(i) = start(i) - ng
          old_1(i) = old_0(i) + n_local(i) - 1
        ENDDO
      ENDIF

      n = 1
      IF (type_max(n) .EQ. old_max(n)) EXIT
      iproc = iproc + 1
      type_min(n) = new_cell_min1(iproc)
    ENDDO

    ! Create array of recvtypes

    DO i = 1,nd
      n_global(i) = new_max(i) - new_min(i) + 2 * ng + 1
    ENDDO

    n = 1
    type_min(n) = new_min(n)
    type_max(n) = new_min(n)

    ! Find the old processor on which the new x_min resides
    ! This could be sped up by using bisection.
    DO iproc = 1, nprocs(n)-1
      IF (old_cell_min1(iproc) .LE. new_min(n) &
          .AND. old_cell_max1(iproc) .GE. new_min(n)) EXIT
    ENDDO

    DO WHILE(type_max(n) .LE. new_max(n))
      coord(cdim(n)) = iproc - 1
      type_max(n) = old_cell_max1(iproc)
      IF (type_max(n) .GT. new_max(n)) type_max(n) = new_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) .EQ. nmin(n)) ng0 = ng
      IF (type_max(n) .EQ. nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - new_min(n) + ng - ng0 + 1

      CALL MPI_CART_RANK(comm, coord, irank, errcode)

      IF (rank .NE. irank) THEN
        recvtypes(irank) = create_1d_array_subtype(basetype, n_local, &
            n_global, start)
      ELSE
        ! New domain is on the same processor as the old domain.
        ! Just copy the region rather than using MPI.
        DO i = 1,nd
          new_0(i) = start(i) - ng
        ENDDO
        DO i = old_0(1),old_1(1)
          inew = new_0(1) + i - old_0(1)
          field_out(inew) = field_in(i)
        ENDDO
      ENDIF

      n = 1
      IF (type_max(n) .EQ. new_max(n)) EXIT
      iproc = iproc + 1
      type_min(n) = old_cell_min1(iproc)
    ENDDO

    CALL redblack(field_in, field_out, sendtypes, recvtypes)

    DO i = 0,nproc-1
      IF (sendtypes(i) .NE. 0) CALL MPI_TYPE_FREE(sendtypes(i), errcode)
      IF (recvtypes(i) .NE. 0) CALL MPI_TYPE_FREE(recvtypes(i), errcode)
    ENDDO

    DEALLOCATE(sendtypes)
    DEALLOCATE(recvtypes)

  END SUBROUTINE redistribute_field_1d_r4



  SUBROUTINE get_load_in_x(load)

    ! Calculate total load across the X direction

    INTEGER(i8), DIMENSION(:), INTENT(OUT) :: load
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: temp
    TYPE(particle), POINTER :: current
    INTEGER :: cell, ispecies, sz

    load = 0

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so x_grid_min, NOT x_grid_min_local
        cell = FLOOR((current%part_pos - x_grid_min) / dx + 1.5_num)

        load(cell) = load(cell) + 1
        current => current%next
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



  SUBROUTINE calculate_breaks(load, nproc, mins, maxs)

    ! This subroutine calculates the places in a given load profile to split
    ! The domain to give the most even subdivision possible

    INTEGER(i8), INTENT(IN), DIMENSION(:) :: load
    INTEGER, INTENT(IN) :: nproc
    INTEGER, DIMENSION(:), INTENT(OUT) :: mins, maxs
    INTEGER :: sz, idim, proc, old, nextra
    INTEGER(i8) :: total, total_old, load_per_proc_ideal

    sz = SIZE(load)
    maxs = sz

    load_per_proc_ideal = FLOOR((SUM(load) + 0.5d0) / nproc, i8)

    proc = 1
    total = 0
    old = 1
    nextra = 0
    DO idim = 1, sz
      IF (nextra .GT. 0) THEN
        nextra = nextra - 1
        CYCLE
      ENDIF
      total_old = total
      total = total + load(idim)
      IF (total .GE. load_per_proc_ideal) THEN
        ! Pick the split that most closely matches the load
        IF (load_per_proc_ideal - total_old &
            .LT. total - load_per_proc_ideal) THEN
          maxs(proc) = idim - 1
        ELSE
          maxs(proc) = idim
        ENDIF
        ! To communicate ghost cell information correctly, each domain must
        ! contain at least ng cells.
        nextra = old - maxs(proc) + ng
        IF (nextra .GT. 0) THEN
          maxs(proc) = maxs(proc) + nextra
        ENDIF
        proc = proc + 1
        IF (proc .EQ. nproc) EXIT
        total = 0
        old = maxs(proc-1)
      ENDIF
    ENDDO
    maxs(nproc) = sz

    ! Sanity check. Must be one cell of separation between each endpoint.
    ! Backwards
    old = sz
    DO proc = nproc-1, 1, -1
      IF (old - maxs(proc) .LT. ng) THEN
        maxs(proc) = old - ng
      ENDIF
      old = maxs(proc)
    ENDDO

    ! Forwards (unnecessary?)
    old = 0
    DO proc = 1, nproc-1
      IF (maxs(proc) - old .LT. ng) THEN
        maxs(proc) = old + ng
      ENDIF
      old = maxs(proc)
    ENDDO

    ! Set mins
    mins(1) = 1
    DO proc = 2, nproc
      mins(proc) = maxs(proc-1) + 1
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
      IF (a_particle%part_pos .GE. x_grid_mins(iproc) - dx / 2.0_num &
          .AND. a_particle%part_pos .LT. x_grid_maxs(iproc) + dx / 2.0_num) THEN
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
    INTEGER :: part_proc, iproc, ispecies, ierr
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: sendcounts, recvcounts

    ALLOCATE(pointers_send(0:nproc-1), pointers_recv(0:nproc-1))
    ALLOCATE(sendcounts(0:nproc-1), recvcounts(0:nproc-1))

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO iproc = 0, nproc - 1
        CALL create_empty_partlist(pointers_send(iproc))
        CALL create_empty_partlist(pointers_recv(iproc))
      ENDDO

      DO WHILE(ASSOCIATED(current))
        next => current%next
        part_proc = get_particle_processor(current)
        IF (part_proc .LT. 0) THEN
          PRINT *, 'Unlocatable particle on processor', rank, current%part_pos
          CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
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
        current => next
      ENDDO

      DO iproc = 0, nproc - 1
        sendcounts(iproc) = pointers_send(iproc)%count
      ENDDO

      CALL MPI_ALLTOALL(sendcounts, 1, MPI_INTEGER8, recvcounts, 1, &
          MPI_INTEGER8, comm, errcode)

      CALL redblack(pointers_send, pointers_recv, sendcounts, recvcounts)

      DO iproc = 0, nproc - 1
        CALL append_partlist(species_list(ispecies)%attached_list, &
            pointers_recv(iproc))
      ENDDO
    ENDDO

    DEALLOCATE(sendcounts, recvcounts)
    DEALLOCATE(pointers_send, pointers_recv)

  END SUBROUTINE distribute_particles

END MODULE balance
