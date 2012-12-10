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

    INTEGER, PARAMETER :: ndims = 2
    INTEGER :: dims(ndims), idim
    LOGICAL :: periods(ndims), reorder, op
    INTEGER :: test_coords(ndims)
    INTEGER :: ix, iy
    INTEGER :: nxsplit, nysplit
    INTEGER :: area, minarea

    IF (comm .NE. MPI_COMM_NULL) CALL MPI_COMM_FREE(comm, errcode)

    IF (MAX(nprocx,1) * MAX(nprocy,1) .GT. nproc) THEN
      IF (rank .EQ. 0) THEN
        PRINT *, 'Unable to use requested processor subdivision. Using ' &
            // 'default division.'
      ENDIF
      nprocx = 0
      nprocy = 0
    ENDIF

    IF (nprocx * nprocy .EQ. 0) THEN
      ! Find the processor split which minimizes surface area of
      ! the resulting domain

      minarea = nx_global + ny_global

      DO ix = 1, nproc
        iy = nproc / ix
        IF (ix * iy .NE. nproc) CYCLE

        nxsplit = nx_global / ix
        nysplit = ny_global / iy

        area = nxsplit + nysplit
        IF (area .LT. minarea) THEN
          nprocx = ix
          nprocy = iy
          minarea = area
        ENDIF
      ENDDO
    ENDIF

    dims = (/nprocy, nprocx/)
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

    IF (bc_field(c_bd_y_min) .EQ. c_bc_periodic &
        .OR. bc_particle(c_bd_y_min) .EQ. c_bc_periodic) &
            periods(c_ndims-1) = .TRUE.

    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, reorder, &
        comm, errcode)
    CALL MPI_COMM_RANK(comm, rank, errcode)
    CALL MPI_CART_COORDS(comm, rank, ndims, coordinates, errcode)
    CALL MPI_CART_SHIFT(comm, 1, 1, proc_x_min, proc_x_max, errcode)
    CALL MPI_CART_SHIFT(comm, 0, 1, proc_y_min, proc_y_max, errcode)

    nprocx = dims(2)
    nprocy = dims(1)
    nprocdir = dims

    IF (rank .EQ. 0) THEN
      PRINT *, 'Processor subdivision is ', (/nprocx, nprocy/)
    ENDIF

    x_coords = coordinates(c_ndims)
    x_min_boundary = .FALSE.
    x_max_boundary = .FALSE.
    IF (x_coords .EQ. 0) x_min_boundary = .TRUE.
    IF (x_coords .EQ. nprocx - 1) x_max_boundary = .TRUE.

    y_coords = coordinates(c_ndims-1)
    y_min_boundary = .FALSE.
    y_max_boundary = .FALSE.
    IF (y_coords .EQ. 0) y_min_boundary = .TRUE.
    IF (y_coords .EQ. nprocy - 1) y_max_boundary = .TRUE.

    neighbour = MPI_PROC_NULL
    DO iy = -1, 1
      DO ix = -1, 1
        test_coords = coordinates
        test_coords(1) = test_coords(1)+iy
        test_coords(2) = test_coords(2)+ix
        op = .TRUE.

        ! For some stupid reason MPI_CART_RANK returns an error rather than
        ! MPI_PROC_NULL if the coords are out of range.
        DO idim = 1, ndims
          IF ((test_coords(idim) .LT. 0 &
              .OR. test_coords(idim) .GE. dims(idim)) &
              .AND. .NOT. periods(idim)) op = .FALSE.
        ENDDO
        IF (op) CALL MPI_CART_RANK(comm, test_coords, neighbour(ix,iy), errcode)
      ENDDO
    ENDDO

  END SUBROUTINE setup_communicator



  SUBROUTINE mpi_initialise

    INTEGER :: ispecies, idim
    INTEGER :: nx0, nxp
    INTEGER :: ny0, nyp

    IF (.NOT.cpml_boundaries) cpml_thickness = 0

    CALL setup_communicator

    nx_global = nx_global + 2 * cpml_thickness
    ny_global = ny_global + 2 * cpml_thickness

    nx0 = nx_global / nprocx
    ny0 = ny_global / nprocy

    nx  = nx0
    ny  = ny0

    ! If the number of gridpoints cannot be exactly subdivided then fix
    ! The first nxp processors have nx0 grid points
    ! The remaining processors have nx0+1 grid points
    IF (nx * nprocx .NE. nx_global) THEN
      nxp = (nx + 1) * nprocx - nx_global
      IF (x_coords .GE. nxp) nx = nx + 1
    ELSE
      nxp = nprocx
    ENDIF

    IF (ny * nprocy .NE. ny_global) THEN
      nyp = (ny + 1) * nprocy - ny_global
      IF (y_coords .GE. nyp) ny = ny + 1
    ELSE
      nyp = nprocy
    ENDIF

    ALLOCATE(npart_each_rank(nproc))
    ALLOCATE(x_mins(0:nprocx-1), x_maxs(0:nprocx-1))
    ALLOCATE(y_mins(0:nprocy-1), y_maxs(0:nprocy-1))
    ALLOCATE(cell_x_min(nprocx), cell_x_max(nprocx))
    ALLOCATE(cell_y_min(nprocy), cell_y_max(nprocy))

    DO idim = 1, nxp
      cell_x_min(idim) = (idim - 1) * nx0 + 1
      cell_x_max(idim) = idim * nx0
    ENDDO
    DO idim = nxp + 1, nprocx
      cell_x_min(idim) = nxp * nx0 + (idim - nxp - 1) * (nx0 + 1) + 1
      cell_x_max(idim) = nxp * nx0 + (idim - nxp) * (nx0 + 1)
    ENDDO

    DO idim = 1, nyp
      cell_y_min(idim) = (idim - 1) * ny0 + 1
      cell_y_max(idim) = idim * ny0
    ENDDO
    DO idim = nyp + 1, nprocy
      cell_y_min(idim) = nyp * ny0 + (idim - nyp - 1) * (ny0 + 1) + 1
      cell_y_max(idim) = nyp * ny0 + (idim - nyp) * (ny0 + 1)
    ENDDO

    subtype_field = 0

    ALLOCATE(x(1-ng:nx+ng), y(1-ng:ny+ng))
    ALLOCATE(x_global(1-ng:nx_global+ng))
    ALLOCATE(y_global(1-ng:ny_global+ng))
    ALLOCATE(xb_global(nx_global+1))
    ALLOCATE(yb_global(ny_global+1))
    ALLOCATE(xb_offset_global(nx_global+1))
    ALLOCATE(yb_offset_global(ny_global+1))
    ALLOCATE(ex(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(ey(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(ez(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(bx(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(by(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(bz(1-ng:nx+ng, 1-ng:ny+ng))
    ! Current may need an extra layer of ghostcells.
    ALLOCATE(jx(1-jng:nx+jng, 1-jng:ny+jng))
    ALLOCATE(jy(1-jng:nx+jng, 1-jng:ny+jng))
    ALLOCATE(jz(1-jng:nx+jng, 1-jng:ny+jng))

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
        ALLOCATE(species_list(ispecies)%ext_temp_x_min(1-ng:ny+ng,1:3))
      ENDIF
      IF (bc_particle(c_bd_x_max) .EQ. c_bc_thermal) THEN
        ALLOCATE(species_list(ispecies)%ext_temp_x_max(1-ng:ny+ng,1:3))
      ENDIF
      IF (bc_particle(c_bd_y_min) .EQ. c_bc_thermal) THEN
        ALLOCATE(species_list(ispecies)%ext_temp_y_min(1-ng:nx+ng,1:3))
      ENDIF
      IF (bc_particle(c_bd_y_max) .EQ. c_bc_thermal) THEN
        ALLOCATE(species_list(ispecies)%ext_temp_y_max(1-ng:nx+ng,1:3))
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
