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

    INTEGER, PARAMETER :: ndims = 1
    INTEGER :: dims(ndims), idim
    LOGICAL :: periods(ndims), reorder, op
    INTEGER :: test_coords(ndims)
    INTEGER :: ix

    IF (comm .NE. MPI_COMM_NULL) CALL MPI_COMM_FREE(comm, errcode)

    dims = (/nprocx/)
    CALL MPI_DIMS_CREATE(nproc, ndims, dims, errcode)

    periods = .TRUE.
    reorder = .TRUE.

    IF (bc_x_min .NE. c_bc_periodic &
        .OR. bc_x_max .NE. c_bc_periodic) periods(1) = .FALSE.

    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, reorder, &
        comm, errcode)
    CALL MPI_COMM_RANK(comm, rank, errcode)
    CALL MPI_CART_COORDS(comm, rank, 1, coordinates, errcode)
    CALL MPI_CART_SHIFT(comm, 0, 1, proc_x_min, proc_x_max, errcode)

    nprocx = dims(1)

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

    CALL setup_communicator
    nx = nx_global/nprocx

    ALLOCATE(nx_each_rank(1:nproc), npart_each_rank(1:nproc))
    subtype_field = 0
    subtype_particle = 0

    ALLOCATE(x(-2:nx+3))
    ALLOCATE(x_global(-2:nx_global+3))
    ALLOCATE(x_offset_global(-2:nx_global+3))
    ALLOCATE(ex(-2:nx+3), ey(-2:nx+3), ez(-2:nx+3))
    ALLOCATE(bx(-2:nx+3), by(-2:nx+3), bz(-2:nx+3))
    ALLOCATE(jx(-2:nx+3), jy(-2:nx+3), jz(-2:nx+3))
    ALLOCATE(ekbar(1:nx, 1:n_species), ekbar_sum(-2:nx+3, 1:n_species))
    ALLOCATE(ct(-2:nx+3, 1:n_species))
    ALLOCATE(starts_x(0:nprocx-1), ends_x(0:nprocx-1))
    ALLOCATE(cell_x_min(1:nprocx), cell_x_max(1:nprocx))

    DO idim = 0, nprocx-1
      cell_x_min(idim+1) = nx*idim+1
      cell_x_max(idim+1) = nx*(idim+1)
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
