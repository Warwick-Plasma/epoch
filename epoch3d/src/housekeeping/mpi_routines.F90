MODULE mpi_routines

  USE partlist

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpi_initialise, mpi_close, mpi_minimal_init, setup_communicator

  REAL(dbl) :: start_time, end_time

CONTAINS

  SUBROUTINE mpi_minimal_init

    INTEGER :: s

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

    dims = (/nprocz, nprocy, nprocx/)
    CALL MPI_DIMS_CREATE(nproc, ndims, dims, errcode)

    periods = .TRUE.
    reorder = .TRUE.

    IF (bc_field(c_bd_x_min)  .NE. c_bc_periodic &
        .OR. bc_field(c_bd_x_max) .NE. c_bc_periodic &
        .OR. bc_particle(c_bd_x_min) .NE. c_bc_periodic &
        .OR. bc_particle(c_bd_x_max) .NE. c_bc_periodic) &
            periods(c_ndims) = .FALSE.

    IF (bc_field(c_bd_y_min)  .NE. c_bc_periodic &
        .OR. bc_field(c_bd_y_max) .NE. c_bc_periodic &
        .OR. bc_particle(c_bd_y_min) .NE. c_bc_periodic &
        .OR. bc_particle(c_bd_y_max) .NE. c_bc_periodic) &
            periods(c_ndims-1) = .FALSE.

    IF (bc_field(c_bd_z_min)  .NE. c_bc_periodic &
        .OR. bc_field(c_bd_z_max) .NE. c_bc_periodic &
        .OR. bc_particle(c_bd_z_min) .NE. c_bc_periodic &
        .OR. bc_particle(c_bd_z_max) .NE. c_bc_periodic) &
            periods(c_ndims-2) = .FALSE.

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
      PRINT *, "Processor subdivision is ", (/nprocx, nprocy, nprocz/)
    ENDIF
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
          IF (op) &
              CALL MPI_CART_RANK(comm, test_coords, neighbour(ix, iy, iz), &
                  errcode)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE setup_communicator



  SUBROUTINE mpi_initialise

    INTEGER :: ispecies, idim
    INTEGER :: nx_big, nx_little
    INTEGER :: ny_big, ny_little
    INTEGER :: nz_big, nz_little

    CALL setup_communicator

    nx = nx_global / nprocx
    ny = ny_global / nprocy
    nz = nz_global / nprocz

    ! If the number of gridpoints cannot be exactly subdivided then fix
    IF (nx * nprocx .NE. nx_global) THEN
      nx_big = nx + 1
      nx_little = nx_global - (nx + 1) * (nprocx - 1)
      IF (coordinates(3) .NE. nprocx-1) THEN
        nx = nx_big
      ELSE
        nx = nx_little
      ENDIF
    ELSE
      nx_big = nx
      nx_little = nx
    ENDIF

    IF (ny * nprocy .NE. ny_global) THEN
      ny_big = ny + 1
      ny_little = ny_global - (ny + 1) * (nprocy - 1)
      IF (coordinates(2) .NE. nprocy-1) THEN
        ny = ny_big
      ELSE
        ny = ny_little
      ENDIF
    ELSE
      ny_big = ny
      ny_little = ny
    ENDIF

    IF (nz * nprocz .NE. nz_global) THEN
      nz_big = nz + 1
      nz_little = nz_global - (nz + 1) * (nprocz - 1)
      IF (coordinates(1) .NE. nprocz-1) THEN
        nz = nz_big
      ELSE
        nz = nz_little
      ENDIF
    ELSE
      nz_big = nz
      nz_little = nz
    ENDIF

    ALLOCATE(npart_each_rank(1:nproc))
    ALLOCATE(x_mins(0:nprocx-1), x_maxs(0:nprocx-1))
    ALLOCATE(y_mins(0:nprocy-1), y_maxs(0:nprocy-1))
    ALLOCATE(z_mins(0:nprocz-1), z_maxs(0:nprocz-1))
    ALLOCATE(cell_x_min(1:nprocx), cell_x_max(1:nprocx))
    ALLOCATE(cell_y_min(1:nprocy), cell_y_max(1:nprocy))
    ALLOCATE(cell_z_min(1:nprocz), cell_z_max(1:nprocz))

    DO idim = 1, nprocx-1
      cell_x_min(idim) = nx_big * (idim - 1) + 1
      cell_x_max(idim) = nx_big * idim
    ENDDO
    cell_x_min(nprocx) = (nprocx - 1) * nx_big + 1
    cell_x_max(nprocx) = (nprocx - 1) * nx_big + nx_little

    DO idim = 1, nprocy-1
      cell_y_min(idim) = ny_big * (idim - 1) + 1
      cell_y_max(idim) = ny_big * idim
    ENDDO
    cell_y_min(nprocy) = (nprocy - 1) * ny_big + 1
    cell_y_max(nprocy) = (nprocy - 1) * ny_big + ny_little

    DO idim = 1, nprocz-1
      cell_z_min(idim) = nz_big * (idim - 1) + 1
      cell_z_max(idim) = nz_big * idim
    ENDDO
    cell_z_min(nprocz) = (nprocz - 1) * nz_big + 1
    cell_z_max(nprocz) = (nprocz - 1) * nz_big + nz_little

    subtype_field = 0

    ALLOCATE(x(-2:nx+3), y(-2:ny+3), z(-2:nz+3))
    ALLOCATE(x_global(-2:nx_global+3))
    ALLOCATE(y_global(-2:ny_global+3))
    ALLOCATE(z_global(-2:nz_global+3))
    ALLOCATE(x_offset_global(-2:nx_global+3))
    ALLOCATE(y_offset_global(-2:ny_global+3))
    ALLOCATE(z_offset_global(-2:nz_global+3))
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
        NULLIFY(particle_species(1)%prev, particle_species(n_species)%next)
    DO ispecies = 1, n_species-1
      particle_species(ispecies)%next=>particle_species(ispecies+1)
    ENDDO
    DO ispecies = 2, n_species
      particle_species(ispecies)%prev=>particle_species(ispecies-1)
    ENDDO
    DO ispecies = 1, n_species
      particle_species(ispecies)%id = ispecies
#ifdef PARTICLE_PROBES
      NULLIFY(particle_species(ispecies)%attached_probes)
#endif
      NULLIFY(particle_species(ispecies)%attached_list%next)
      NULLIFY(particle_species(ispecies)%attached_list%prev)
      CALL create_empty_partlist(particle_species(ispecies)%attached_list)
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
      WRITE(20, *)
      WRITE(20, '("runtime = ", i4, "h ", i2, "m ", i2, "s on ", i4, &
          &" process elements.")') hours, minutes, seconds, nproc
    ENDIF

    CALL MPI_BARRIER(comm, errcode)

  END SUBROUTINE mpi_close

END MODULE mpi_routines
