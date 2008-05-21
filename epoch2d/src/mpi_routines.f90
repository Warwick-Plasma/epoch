MODULE mpi_routines

  USE shared_data
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpi_initialise, mpi_close,mpi_minimal_init

  REAL(dbl) :: start_time, end_time

CONTAINS

  SUBROUTINE mpi_minimal_init
    INTEGER s

    CALL MPI_INIT(errcode)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, errcode)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, errcode)

    IF (rank .EQ. 0) THEN
      open(unit=10,file="./input.deck",iostat=s,status='OLD')
      deckfile=(s .EQ. 0)
      close(10)
    ENDIF

  END SUBROUTINE mpi_minimal_init


  SUBROUTINE mpi_initialise

    INTEGER :: ndims, dims(3)
    LOGICAL :: periods(3), reorder
    INTEGER :: starts(2), sizes(2), subsizes(2)

!    PRINT *,"IN MPI"
    ndims=1

    npart=npart_global/nproc
    IF (npart*nproc /= npart_global) THEN
       IF (rank .EQ.0) PRINT *,"Unable to divide particles at t=0. Quitting."
       CALL MPI_ABORT(MPI_COMM_WORLD,errcode)
    ENDIF
    comm=MPI_COMM_WORLD

    ! Create the subarray for this problem: subtype decribes where this
    ! process's data fits into the global picture.

    ! set up the starting point for my subgrid (assumes arrays start at 0)
    starts(1) = rank * npart

    ! the grid sizes
    subsizes = (/npart,0/)
    sizes = (/npart_global,0/)

    ! set up and commit the subarray type
    CALL MPI_TYPE_CREATE_SUBARRAY(ndims, sizes, subsizes, starts, &
         MPI_ORDER_FORTRAN, mpireal, subtype, errcode)
    CALL MPI_TYPE_COMMIT(subtype,errcode)

    CALL MPI_TYPE_CREATE_SUBARRAY(ndims, sizes, subsizes, starts, &
         MPI_ORDER_FORTRAN, MPI_INTEGER, subtype_int, errcode)
    CALL MPI_TYPE_COMMIT(subtype_int,errcode)


    !Now for the particle mesh
    starts = (/rank * npart,1/)
    sizes=(/npart_global,3/)
    subsizes=(/npart,3/)

    CALL MPI_TYPE_CREATE_SUBARRAY(ndims, sizes, subsizes, starts, &
         MPI_ORDER_FORTRAN, mpireal, subtype_particle_mesh, errcode)

    CALL MPI_TYPE_COMMIT(subtype_particle_mesh,errcode)



    ALLOCATE(x(-2:nx+3),y(-2:ny+3))
    ALLOCATE(Ex(-2:nx+3,-2:ny+3),Ey(-2:nx+3,-2:ny+3),Ez(-2:nx+3,-2:ny+3))
    ALLOCATE(Bx(-2:nx+3,-2:ny+3),By(-2:nx+3,-2:ny+3),Bz(-2:nx+3,-2:ny+3))
    ALLOCATE(Jx(-2:nx+3,-2:ny+3),Jy(-2:nx+3,-2:ny+3),Jz(-2:nx+3,-2:ny+3))
    ALLOCATE(temperature(1:nx,1:ny,1:nspecies),sum(-2:nx+3,-2:ny+3,1:nspecies))
    ALLOCATE(ct(-2:nx+3,-2:ny+3,1:nspecies))

    ALLOCATE(Part_pos(1:npart,1:2),Part_P(1:npart,1:3),Part_species(1:npart))

    start_time=MPI_WTIME()

  END SUBROUTINE mpi_initialise



  SUBROUTINE mpi_close
    INTEGER :: seconds, minutes, hours, total

    IF(rank == 0) THEN
       end_time = MPI_WTIME()
       total = INT(end_time - start_time)
       seconds = MOD(total,60)
       minutes = MOD(total / 60, 60)
       hours = total / 3600
       WRITE(20,*)
       WRITE(20,'("runtime = ",i4,"h ",i2,"m ",i2,"s on ",i4," process elements.")') &
            hours,minutes,seconds,nproc
    END IF

    CALL MPI_BARRIER(comm, errcode)


  END SUBROUTINE mpi_close



END MODULE mpi_routines
