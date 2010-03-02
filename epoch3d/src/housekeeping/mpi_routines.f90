MODULE mpi_routines

  USE shared_data
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

    IF (comm .NE. MPI_COMM_NULL) CALL MPI_COMM_FREE(comm, errcode)

    IF (MAX(nprocx, 1)*MAX(nprocy, 1)*MAX(nprocz, 1) .GT. nproc) THEN
      IF (rank .EQ. 0) &
          PRINT *, "Unable to use requested processor subdivision. Using &
              &default division."
      nprocx = 0
      nprocy = 0
      nprocz = 0
    ENDIF

    dims = (/nprocz, nprocy, nprocx/)
    CALL MPI_DIMS_CREATE(nproc, ndims, dims, errcode)

    periods = .TRUE.
    reorder = .TRUE.

    IF (xbc_left  .NE. c_bc_periodic .OR. &
        xbc_right .NE. c_bc_periodic) periods(3) = .FALSE.
    IF (ybc_down  .NE. c_bc_periodic .OR. &
        ybc_up    .NE. c_bc_periodic) periods(2) = .FALSE.
    IF (zbc_front .NE. c_bc_periodic .OR. &
        zbc_back  .NE. c_bc_periodic) periods(1) = .FALSE.

    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, reorder, &
        comm, errcode)
    CALL MPI_COMM_RANK(comm, rank, errcode)
    CALL MPI_CART_COORDS(comm, rank, ndims, coordinates, errcode)
    CALL MPI_CART_SHIFT(comm, 2, 1, left, right, errcode)
    CALL MPI_CART_SHIFT(comm, 1, 1, down, up, errcode)
    CALL MPI_CART_SHIFT(comm, 0, 1, back, front, errcode)

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
            IF ((test_coords(idim) .LT. 0 .OR. &
                test_coords(idim) .GE. dims(idim)) .AND. &
                .NOT. periods(idim)) op = .FALSE.
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
    INTEGER(KIND=8) :: npart_this_species, npart

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

    npart = npart_global/nproc
    IF (npart*nproc .NE. npart_global) THEN
      IF (rank .EQ. 0) PRINT *, "Unable to divide particles at t = 0. Quitting."
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode)
    ENDIF

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
!!$    ALLOCATE(start_each_rank(0:nproc-1, 1:3), end_each_rank(0:nproc-1, 1:3))
    ALLOCATE(x_starts(0:nprocx-1), x_ends(0:nprocx-1))
    ALLOCATE(y_starts(0:nprocy-1), y_ends(0:nprocy-1))
    ALLOCATE(z_starts(0:nprocz-1), z_ends(0:nprocz-1))
    ALLOCATE(cell_x_start(1:nprocx), cell_x_end(1:nprocx))
    ALLOCATE(cell_y_start(1:nprocy), cell_y_end(1:nprocy))
    ALLOCATE(cell_z_start(1:nprocz), cell_z_end(1:nprocz))

    DO idim = 0, nprocx-1
      cell_x_start(idim+1) = nx*idim+1
      cell_x_end(idim+1) = nx*(idim+1)
    ENDDO

    DO idim = 0, nprocy-1
      cell_y_start(idim+1) = ny*idim+1
      cell_y_end(idim+1) = ny*(idim+1)
    ENDDO

    DO idim = 0, nprocz-1
      cell_z_start(idim+1) = nz*idim+1
      cell_z_end(idim+1) = nz*(idim+1)
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
      npart_this_species = particle_species(ispecies)%count
      NULLIFY(particle_species(ispecies)%attached_list%next)
      NULLIFY(particle_species(ispecies)%attached_list%prev)
      IF (restart .OR. IOR(ictype, c_ic_autoload) .NE. 0) THEN
        CALL create_empty_partlist(particle_species(ispecies)%attached_list)
      ELSE
        CALL create_allocated_partlist(&
            particle_species(ispecies)%attached_list, npart_this_species)
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
      WRITE(20, *)
      WRITE(20, '("runtime = ", i4, "h ", i2, "m ", i2, "s on ", i4, &
          &" process elements.")') hours, minutes, seconds, nproc
    ENDIF

    CALL MPI_BARRIER(comm, errcode)

  END SUBROUTINE mpi_close

END MODULE mpi_routines
