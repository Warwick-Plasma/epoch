MODULE mpi_routines

  USE shared_data
  USE partlist
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

  END SUBROUTINE mpi_minimal_init


  SUBROUTINE mpi_initialise

    INTEGER :: ndims, dims(1)
    LOGICAL :: periods(1), reorder
    INTEGER :: starts(1), sizes(1), subsizes(1), types(1)
    TYPE(PARTICLE), POINTER :: cur

    !PRINT *,"IN"

    IF (domain == DO_DECOMPOSED) THEN
       ndims=1
       dims=nproc
       periods=.TRUE.
       IF (xbc_left .NE.BC_PERIODIC .OR. xbc_right .NE. BC_PERIODIC) periods=.FALSE.
       CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods,  &
            reorder, comm, errcode)
       CALL MPI_COMM_RANK(comm, rank, errcode)
       CALL MPI_CART_COORDS(comm, rank, 1, coordinates, errcode)
       CALL MPI_CART_SHIFT(comm, 0, 1, left, right, errcode)
       nprocx=nproc
       nx=nx_global/nprocx
    ELSE
       nx=nx_global
       ndims=1
       comm=MPI_COMM_WORLD
       IF (xbc_left .EQ. BC_PERIODIC .OR. xbc_right .EQ. BC_PERIODIC) THEN
          left=rank
          right=rank
       ELSE
          left=MPI_PROC_NULL
          right=MPI_PROC_NULL
       ENDIF
    ENDIF

    ALLOCATE(nx_each_rank(1:nproc),npart_each_rank(1:nproc))
    subtype_field=0
    subtype_particle=0

    !Initially setup npart = npart_global/nproc then rebalance after initial conditions
    npart=npart_global/nproc
    IF (npart*nproc /= npart_global) THEN
       IF (rank .EQ.0) PRINT *,"Unable to divide particles at t=0. Quitting."
       CALL MPI_ABORT(MPI_COMM_WORLD,errcode)
    ENDIF

    ALLOCATE(x(-2:nx+3),x_global(-2:nx_global+3))
    ALLOCATE(Ex(-2:nx+3),Ey(-2:nx+3),Ez(-2:nx+3))
    ALLOCATE(Bx(-2:nx+3),By(-2:nx+3),Bz(-2:nx+3))
    ALLOCATE(Jx(-2:nx+3),Jy(-2:nx+3),Jz(-2:nx+3))
    ALLOCATE(ek_bar(1:nx,1:nspecies),ekbar_sum(-2:nx+3,1:nspecies))
    ALLOCATE(ct(-2:nx+3,1:nspecies))
    ALLOCATE(x_starts(0:nproc-1), x_ends(0:nproc-1))

    jx_sum=0.0_num
    jx_sum2=0.0_num

    CALL Create_Allocated_PartList(MainRoot, npart)

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
       !WRITE(20,'("runtime = ",i4,"h ",i2,"m ",i2,"s on ",i4," process elements.")') &
      !      hours,minutes,seconds,nproc
    END IF

    CALL MPI_BARRIER(comm, errcode)


  END SUBROUTINE mpi_close



END MODULE mpi_routines
