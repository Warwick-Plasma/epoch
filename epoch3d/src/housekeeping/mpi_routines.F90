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

    IF (rank .EQ. 0) THEN
      OPEN(unit=10, file="./input.deck", iostat=s, status='OLD')
      deckfile = (s .EQ. 0)
      CLOSE(10)
    ENDIF

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

    IF (bc_x_min  .NE. c_bc_periodic &
        .OR. bc_x_max .NE. c_bc_periodic) periods(3) = .FALSE.
    IF (bc_y_min  .NE. c_bc_periodic &
        .OR. bc_y_max    .NE. c_bc_periodic) periods(2) = .FALSE.
    IF (bc_z_min .NE. c_bc_periodic &
        .OR. bc_z_max .NE. c_bc_periodic) periods(1) = .FALSE.

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

    CALL setup_communicator

    nx = nx_global/nprocx
    ny = ny_global/nprocy
    nz = nz_global/nprocz

    ALLOCATE(nx_each_rank(1:nproc))
    ALLOCATE(ny_each_rank(1:nproc))
    ALLOCATE(nz_each_rank(1:nproc))
    ALLOCATE(npart_each_rank(1:nproc))
    subtype_field = 0
    subtype_particle = 0

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
    ALLOCATE(ekbar(1:nx, 1:ny, 1:nz, 1:n_species))
    ALLOCATE(ekbar_sum(-2:nx+3, -2:ny+3, -2:nz+3, 1:n_species))
    ALLOCATE(ct(-2:nx+3, -2:ny+3, -2:nz+3, 1:n_species))
    ALLOCATE(starts_x(0:nprocx-1), ends_x(0:nprocx-1))
    ALLOCATE(starts_y(0:nprocy-1), ends_y(0:nprocy-1))
    ALLOCATE(starts_z(0:nprocz-1), ends_z(0:nprocz-1))
    ALLOCATE(cell_x_min(1:nprocx), cell_x_max(1:nprocx))
    ALLOCATE(cell_y_min(1:nprocy), cell_y_max(1:nprocy))
    ALLOCATE(cell_z_min(1:nprocz), cell_z_max(1:nprocz))

    DO idim = 0, nprocx-1
      cell_x_min(idim+1) = nx*idim+1
      cell_x_max(idim+1) = nx*(idim+1)
    ENDDO

    DO idim = 0, nprocy-1
      cell_y_min(idim+1) = ny*idim+1
      cell_y_max(idim+1) = ny*(idim+1)
    ENDDO

    DO idim = 0, nprocz-1
      cell_z_min(idim+1) = nz*idim+1
      cell_z_max(idim+1) = nz*(idim+1)
    ENDDO

    ! Setup the particle lists
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
