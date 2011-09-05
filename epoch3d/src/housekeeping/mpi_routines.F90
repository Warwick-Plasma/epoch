MODULE mpi_routines

  USE mpi
  USE partlist

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

    INTEGER, PARAMETER :: ndims = 3
    INTEGER :: dims(ndims), idim
    LOGICAL :: periods(ndims), reorder, op
    INTEGER :: test_coords(ndims)
    INTEGER :: ix, iy, iz
    INTEGER :: nxsplit, nysplit, nzsplit
    INTEGER :: area, minarea, nprocyz

    IF (comm .NE. MPI_COMM_NULL) CALL MPI_COMM_FREE(comm, errcode)

    IF (MAX(nprocx,1) * MAX(nprocy,1) * MAX(nprocz,1) .GT. nproc) THEN
      IF (rank .EQ. 0) THEN
        PRINT *, 'Unable to use requested processor subdivision. Using ' &
            // 'default division.'
      ENDIF
      nprocx = 0
      nprocy = 0
      nprocz = 0
    ENDIF

    IF (nprocx * nprocy * nprocz .EQ. 0) THEN
      ! Find the processor split which minimizes surface area of
      ! the resulting domain

      minarea = nx_global * ny_global + ny_global * nz_global &
          + nz_global * nx_global

      DO ix = 1, nproc
        nprocyz = nproc / ix
        IF (ix * nprocyz .NE. nproc) CYCLE

        nxsplit = nx_global / ix

        DO iy = 1, nprocyz
          iz = nprocyz / iy
          IF (iy * iz .NE. nprocyz) CYCLE

          nysplit = ny_global / iy
          nzsplit = nz_global / iz

          area = nxsplit * nysplit + nysplit * nzsplit + nzsplit * nxsplit
          IF (area .LT. minarea) THEN
            nprocx = ix
            nprocy = iy
            nprocz = iz
            minarea = area
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    dims = (/nprocz, nprocy, nprocx/)
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

    IF (bc_field(c_bd_z_min) .EQ. c_bc_periodic &
        .OR. bc_particle(c_bd_z_min) .EQ. c_bc_periodic) &
            periods(c_ndims-2) = .TRUE.

    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, reorder, &
        comm, errcode)
    CALL MPI_COMM_RANK(comm, rank, errcode)
    CALL MPI_CART_COORDS(comm, rank, ndims, coordinates, errcode)
    CALL MPI_CART_SHIFT(comm, 2, 1, proc_x_min, proc_x_max, errcode)
    CALL MPI_CART_SHIFT(comm, 1, 1, proc_y_min, proc_y_max, errcode)
    CALL MPI_CART_SHIFT(comm, 0, 1, proc_z_min, proc_z_max, errcode)

    nprocx = dims(3)
    nprocy = dims(2)
    nprocz = dims(1)

    IF (rank .EQ. 0) THEN
      PRINT *, 'Processor subdivision is ', (/nprocx, nprocy, nprocz/)
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

    z_coords = coordinates(c_ndims-2)
    z_min_boundary = .FALSE.
    z_max_boundary = .FALSE.
    IF (z_coords .EQ. 0) z_min_boundary = .TRUE.
    IF (z_coords .EQ. nprocz - 1) z_max_boundary = .TRUE.

    neighbour = MPI_PROC_NULL
    DO iz = -1, 1
      DO iy = -1, 1
        DO ix = -1, 1
          test_coords = coordinates
          test_coords(1) = test_coords(1)+iz
          test_coords(2) = test_coords(2)+iy
          test_coords(3) = test_coords(3)+ix
          op = .TRUE.
          ! For some stupid reason MPI_CART_RANK returns an error rather than
          ! MPI_PROC_NULL if the coords are out of range.
          DO idim = 1, ndims
            IF ((test_coords(idim) .LT. 0 &
                .OR. test_coords(idim) .GE. dims(idim)) &
                .AND. .NOT. periods(idim)) op = .FALSE.
          ENDDO
          IF (op) THEN
            CALL MPI_CART_RANK(comm, test_coords, neighbour(ix,iy,iz), errcode)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE setup_communicator



  SUBROUTINE mpi_initialise

    INTEGER :: ispecies, idim
    INTEGER :: nx0, nxp
    INTEGER :: ny0, nyp
    INTEGER :: nz0, nzp

    IF (.NOT.cpml_boundaries) cpml_thickness = 0

    CALL setup_communicator

    nx_global = nx_global + 2 * cpml_thickness
    ny_global = ny_global + 2 * cpml_thickness
    nz_global = nz_global + 2 * cpml_thickness

    nx0 = nx_global / nprocx
    ny0 = ny_global / nprocy
    nz0 = nz_global / nprocz

    nx  = nx0
    ny  = ny0
    nz  = nz0

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

    IF (nz * nprocz .NE. nz_global) THEN
      nzp = (nz + 1) * nprocz - nz_global
      IF (z_coords .GE. nzp) nz = nz + 1
    ELSE
      nzp = nprocz
    ENDIF

    ALLOCATE(npart_each_rank(1:nproc))
    ALLOCATE(x_mins(0:nprocx-1), x_maxs(0:nprocx-1))
    ALLOCATE(y_mins(0:nprocy-1), y_maxs(0:nprocy-1))
    ALLOCATE(z_mins(0:nprocz-1), z_maxs(0:nprocz-1))
    ALLOCATE(cell_x_min(1:nprocx), cell_x_max(1:nprocx))
    ALLOCATE(cell_y_min(1:nprocy), cell_y_max(1:nprocy))
    ALLOCATE(cell_z_min(1:nprocz), cell_z_max(1:nprocz))

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

    DO idim = 1, nzp
      cell_z_min(idim) = (idim - 1) * nz0 + 1
      cell_z_max(idim) = idim * nz0
    ENDDO
    DO idim = nzp + 1, nprocz
      cell_z_min(idim) = nzp * nz0 + (idim - nzp - 1) * (nz0 + 1) + 1
      cell_z_max(idim) = nzp * nz0 + (idim - nzp) * (nz0 + 1)
    ENDDO

    subtype_field = 0

    ALLOCATE(x(-2:nx+3), y(-2:ny+3), z(-2:nz+3))
    ALLOCATE(x_global(-2:nx_global+3))
    ALLOCATE(y_global(-2:ny_global+3))
    ALLOCATE(z_global(-2:nz_global+3))
    ALLOCATE(xb_global(nx_global+1))
    ALLOCATE(yb_global(ny_global+1))
    ALLOCATE(zb_global(nz_global+1))
    ALLOCATE(xb_offset_global(nx_global+1))
    ALLOCATE(yb_offset_global(ny_global+1))
    ALLOCATE(zb_offset_global(nz_global+1))
    ALLOCATE(ex(-2:nx+3, -2:ny+3, -2:nz+3))
    ALLOCATE(ey(-2:nx+3, -2:ny+3, -2:nz+3))
    ALLOCATE(ez(-2:nx+3, -2:ny+3, -2:nz+3))
    ALLOCATE(bx(-2:nx+3, -2:ny+3, -2:nz+3))
    ALLOCATE(by(-2:nx+3, -2:ny+3, -2:nz+3))
    ALLOCATE(bz(-2:nx+3, -2:ny+3, -2:nz+3))
    ALLOCATE(jx(-2:nx+3, -2:ny+3, -2:nz+3))
    ALLOCATE(jy(-2:nx+3, -2:ny+3, -2:nz+3))
    ALLOCATE(jz(-2:nx+3, -2:ny+3, -2:nz+3))

    ! Setup the particle lists
    IF (n_species .GT. 0) &
        NULLIFY(species_list(1)%prev, species_list(n_species)%next)
    DO ispecies = 1, n_species-1
      species_list(ispecies)%next=>species_list(ispecies+1)
    ENDDO
    DO ispecies = 2, n_species
      species_list(ispecies)%prev=>species_list(ispecies-1)
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
        ALLOCATE(species_list(ispecies)%ext_temp_x_min(-2:ny+3,-2:nz+3,1:3))
      ENDIF
      IF (bc_particle(c_bd_x_max) .EQ. c_bc_thermal) THEN
        ALLOCATE(species_list(ispecies)%ext_temp_x_max(-2:ny+3,-2:nz+3,1:3))
      ENDIF
      IF (bc_particle(c_bd_y_min) .EQ. c_bc_thermal) THEN
        ALLOCATE(species_list(ispecies)%ext_temp_y_min(-2:nx+3,-2:nz+3,1:3))
      ENDIF
      IF (bc_particle(c_bd_y_max) .EQ. c_bc_thermal) THEN
        ALLOCATE(species_list(ispecies)%ext_temp_y_max(-2:nx+3,-2:nz+3,1:3))
      ENDIF
      IF (bc_particle(c_bd_z_min) .EQ. c_bc_thermal) THEN
        ALLOCATE(species_list(ispecies)%ext_temp_z_min(-2:nx+3,-2:ny+3,1:3))
      ENDIF
      IF (bc_particle(c_bd_z_max) .EQ. c_bc_thermal) THEN
        ALLOCATE(species_list(ispecies)%ext_temp_z_max(-2:nx+3,-2:ny+3,1:3))
      ENDIF
    ENDDO

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
          ''s on '', i4, '' process elements.'')') hours, minutes, seconds, &
          nproc
    ENDIF

    CALL MPI_BARRIER(comm, errcode)

  END SUBROUTINE mpi_close

END MODULE mpi_routines
