MODULE mpi_routines

  USE mpi
  USE partlist
  USE helper

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpi_initialise, mpi_close, mpi_minimal_init, setup_communicator

  REAL(dbl) :: start_time, end_time

CONTAINS

  SUBROUTINE mpi_minimal_init

    CALL MPI_INIT(errcode)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, errcode)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, errcode)

  END SUBROUTINE mpi_minimal_init



  SUBROUTINE setup_communicator

    INTEGER, PARAMETER :: ndims = 1
    INTEGER :: dims(ndims), idim
    LOGICAL :: periods(ndims), reorder, op
    INTEGER :: test_coords(ndims)
    INTEGER :: ix

    IF (comm .NE. MPI_COMM_NULL) CALL MPI_COMM_FREE(comm, errcode)

    dims = (/nprocx/)
    CALL MPI_DIMS_CREATE(nproc, ndims, dims, errcode)

    periods = .FALSE.
    reorder = .TRUE.

    ! Set boundary to be periodic if *any* boundary condition requires it.
    ! Once there are per-species boundary conditions then this will be true
    ! if any of the species are periodic

    IF (bc_field(c_bd_x_min) .EQ. c_bc_periodic &
        .OR. bc_x_min_after_move .EQ. c_bc_periodic &
        .OR. bc_particle(c_bd_x_min) .EQ. c_bc_periodic) &
            periods(c_ndims) = .TRUE.

    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, reorder, &
        comm, errcode)
    CALL MPI_COMM_RANK(comm, rank, errcode)
    CALL MPI_CART_COORDS(comm, rank, ndims, coordinates, errcode)
    CALL MPI_CART_SHIFT(comm, 0, 1, proc_x_min, proc_x_max, errcode)

    nprocx = dims(1)
    nprocdir = dims

    x_coords = coordinates(c_ndims)
    x_min_boundary = .FALSE.
    x_max_boundary = .FALSE.
    IF (x_coords .EQ. 0) x_min_boundary = .TRUE.
    IF (x_coords .EQ. nprocx - 1) x_max_boundary = .TRUE.

    neighbour = MPI_PROC_NULL
    DO ix = -1, 1
      test_coords = coordinates
      test_coords(1) = test_coords(1)+ix
      op = .TRUE.

      ! For some stupid reason MPI_CART_RANK returns an error rather than
      ! MPI_PROC_NULL if the coords are out of range.
      DO idim = 1, ndims
        IF ((test_coords(idim) .LT. 0 &
            .OR. test_coords(idim) .GE. dims(idim)) &
            .AND. .NOT. periods(idim)) op = .FALSE.
      ENDDO
      IF (op) CALL MPI_CART_RANK(comm, test_coords, neighbour(ix), errcode)
    ENDDO

  END SUBROUTINE setup_communicator



  SUBROUTINE mpi_initialise

    INTEGER :: ispecies, idim
    INTEGER :: nx0, nxp

    IF (.NOT.cpml_boundaries) cpml_thickness = 0

    CALL setup_communicator

    nx_global = nx_global + 2 * cpml_thickness

    nx0 = nx_global / nprocx

    nx  = nx0

    ! If the number of gridpoints cannot be exactly subdivided then fix
    ! The first nxp processors have nx0 grid points
    ! The remaining processors have nx0+1 grid points
    IF (nx * nprocx .NE. nx_global) THEN
      nxp = (nx + 1) * nprocx - nx_global
      IF (x_coords .GE. nxp) nx = nx + 1
    ELSE
      nxp = nprocx
    ENDIF

    ALLOCATE(npart_each_rank(nproc))
    ALLOCATE(x_mins(0:nprocx-1), x_maxs(0:nprocx-1))
    ALLOCATE(cell_x_min(nprocx), cell_x_max(nprocx))

    DO idim = 1, nxp
      cell_x_min(idim) = (idim - 1) * nx0 + 1
      cell_x_max(idim) = idim * nx0
    ENDDO
    DO idim = nxp + 1, nprocx
      cell_x_min(idim) = nxp * nx0 + (idim - nxp - 1) * (nx0 + 1) + 1
      cell_x_max(idim) = nxp * nx0 + (idim - nxp) * (nx0 + 1)
    ENDDO

    subtype_field = 0

    ALLOCATE(x(1-ng:nx+ng))
    ALLOCATE(x_global(1-ng:nx_global+ng))
    ALLOCATE(xb_global(nx_global+1))
    ALLOCATE(xb_offset_global(nx_global+1))
    ALLOCATE(ex(1-ng:nx+ng))
    ALLOCATE(ey(1-ng:nx+ng))
    ALLOCATE(ez(1-ng:nx+ng))
    ALLOCATE(bx(1-ng:nx+ng))
    ALLOCATE(by(1-ng:nx+ng))
    ALLOCATE(bz(1-ng:nx+ng))
    ! Current may need an extra layer of ghostcells.
    ALLOCATE(jx(1-jng:nx+jng))
    ALLOCATE(jy(1-jng:nx+jng))
    ALLOCATE(jz(1-jng:nx+jng))

    ! Setup the particle lists
    IF (n_species .GT. 0) &
        NULLIFY(species_list(1)%prev, species_list(n_species)%next)
    DO ispecies = 1, n_species-1
      species_list(ispecies)%next => species_list(ispecies+1)
    ENDDO
    DO ispecies = 2, n_species
      species_list(ispecies)%prev => species_list(ispecies-1)
    ENDDO
    DO ispecies = 1, n_species
      species_list(ispecies)%id = ispecies
#ifdef PARTICLE_PROBES
      NULLIFY(species_list(ispecies)%attached_probes)
#endif
      NULLIFY(species_list(ispecies)%attached_list%next)
      NULLIFY(species_list(ispecies)%attached_list%prev)
      CALL create_empty_partlist(species_list(ispecies)%attached_list)

      IF (bc_particle(c_bd_x_min) .EQ. c_bc_thermal) THEN
        ALLOCATE(species_list(ispecies)%ext_temp_x_min(1:3))
      ENDIF
      IF (bc_particle(c_bd_x_max) .EQ. c_bc_thermal) THEN
        ALLOCATE(species_list(ispecies)%ext_temp_x_max(1:3))
      ENDIF
    ENDDO

    CALL allocate_ic

    start_time = MPI_WTIME()

  END SUBROUTINE mpi_initialise



  SUBROUTINE mpi_close

    INTEGER :: seconds, minutes, hours, total

    IF (rank .EQ. 0) THEN
      end_time = MPI_WTIME()
      total = INT(end_time - start_time)
      seconds = MOD(total, 60)
      minutes = MOD(total / 60, 60)
      hours = total / 3600
      WRITE(stat_unit, *)
      WRITE(stat_unit, '(''runtime = '', i4, ''h '', i2, ''m '', i2, &
          & ''s on '', i4, '' process elements.'')') hours, minutes, seconds, &
          nproc
    ENDIF

    CALL MPI_BARRIER(comm, errcode)

  END SUBROUTINE mpi_close

END MODULE mpi_routines
