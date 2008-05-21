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

    INTEGER :: ndims, dims(3),idim
    LOGICAL :: periods(3), reorder
    INTEGER :: starts(3), sizes(3), subsizes(3),test_coords(3)
    LOGICAL :: op
    TYPE(particle),POINTER :: Cur


    ndims=3
    dims=0
    periods=.FALSE.
    reorder=.TRUE.
    IF (xbc_left .EQ. periodic .OR. xbc_right .EQ. periodic) periods(3)=.TRUE.
    IF (ybc_up .EQ. periodic .OR. ybc_down .EQ. periodic) periods(2)=.TRUE.
    IF (zbc_front .EQ. periodic .OR. zbc_back .EQ. periodic) periods(1)=.TRUE.

    dims = (/nprocz,nprocy,nprocx/)
    IF (MAX(dims(1),1)*MAX(dims(2),1)*MAX(dims(3),1) .GT. nproc) THEN
       dims=0
       IF (rank .EQ. 0) THEN
          PRINT *,"Too many processors requested in override."
          PRINT *,"Reverting to automatic decomposition."
          PRINT *,"******************************************"
          PRINT *,""
       ENDIF
    ENDIF

    CALL MPI_DIMS_CREATE(nproc,ndims,dims,errcode)

    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods,  &
         reorder, comm, errcode)
    CALL MPI_COMM_RANK(comm, rank, errcode)
    CALL MPI_CART_COORDS(comm, rank, ndims, coordinates, errcode)
    CALL MPI_CART_SHIFT(comm, 2, 1, left, right, errcode)
    CALL MPI_CART_SHIFT(comm, 1, 1, down, up, errcode)
    CALL MPI_CART_SHIFT(comm, 0, 1, front, back, errcode)

    neighbour=MPI_PROC_NULL
    DO iz=-1,1
       DO iy=-1,1
          DO ix=-1,1
             test_coords=coordinates
             test_coords(1)=test_coords(1)+iz
             test_coords(2)=test_coords(2)+iy
             test_coords(3)=test_coords(3)+ix
             op=.TRUE.
             !For some stupid reason MPI_CART_RANK returns an error rather than MPI_PROC_NULL
             !If the coords are out of range.
             DO idim=1,ndims
                IF ((test_coords(idim) .LT. 0 .OR. test_coords(idim) .GE. dims(idim)) .AND. .NOT. periods(idim)) op=.FALSE.
             ENDDO
             IF (op) CALL MPI_CART_RANK(comm,test_coords,neighbour(ix,iy,iz),errcode)
          ENDDO
       ENDDO
    ENDDO

!!$    WRITE(rank+10,*),left,right,up,down,front,back
!!$    DO iz=-1,1
!!$       DO iy=-1,1
!!$          DO ix=-1,1
!!$             IF (neighbour(ix,iy,iz) .NE. MPI_PROC_NULL) WRITE(rank+10,'("(",i3,i3,i3,")",i3)'),ix,iy,iz,neighbour(ix,iy,iz)
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO

!!$    OPEN(unit=rank+10,form="formatted")
!!$    WRITE(rank+10,*) neighbour(0,0,:)
!!$    CLOSE(rank+10)

    nprocx=dims(3)
    nprocy=dims(2)
    nprocz=dims(1)

    nx=nx_global/nprocx
    ny=ny_global/nprocy
    nz=nz_global/nprocz

    IF (nx*nprocx /= nx_global .OR. ny*nprocy /= ny_global .OR. nz*nprocz /= nz_global) THEN
       IF (rank .EQ.0) PRINT *,"Unable to divide grid at t=0. Quitting."
       CALL MPI_ABORT(MPI_COMM_WORLD,errcode)
    ENDIF

    starts=(/coordinates(3)*nx,coordinates(2)*ny,coordinates(1)*nz/)
    sizes=(/nx_global,ny_global,nz_global/)
    subsizes=(/nx,ny,nz/)

    CALL MPI_TYPE_CREATE_SUBARRAY(ndims,sizes,subsizes,starts,&
         MPI_ORDER_FORTRAN, mpireal, subtype_field, errcode)
    CALL MPI_TYPE_COMMIT(subtype_field,errcode)

    npart=npart_global/nproc
    IF (npart*nproc /= npart_global) THEN
       IF (rank .EQ.0) PRINT *,"Unable to divide particles at t=0. Quitting."
       CALL MPI_ABORT(MPI_COMM_WORLD,errcode)
    ENDIF

    !Subtypes for particles have been moved to balance.f90 since they are recalculated every time

    ALLOCATE(x_global(-2:nx_global+3),y_global(-2:ny_global+3),z_global(-2:nz_global+3))
    ALLOCATE(x(-2:nx+3),y(-2:ny+3),z(-2:nz+3))

    ALLOCATE(Ex(-2:nx+3,-2:ny+3,-2:nz+3),Ey(-2:nx+3,-2:ny+3,-2:nz+3),Ez(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(Bx(-2:nx+3,-2:ny+3,-2:nz+3),By(-2:nx+3,-2:ny+3,-2:nz+3),Bz(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(Jx(-2:nx+3,-2:ny+3,-2:nz+3),Jy(-2:nx+3,-2:ny+3,-2:nz+3),Jz(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(temperature(1:nx,1:ny,1:nz,1:nspecies),sum(-2:nx+3,-2:ny+3,-2:nz+3,1:nspecies))
    ALLOCATE(ct(-2:nx+3,-2:ny+3,-2:nz+3,1:nspecies))

    ALLOCATE(Part_species(1:npart))
    subtype_particle=0
    subtype_particle_var=0
    subtype_particle_int=0


    !Allocate particle list
    ALLOCATE(Head)
    NULLIFY(Head%Next,Head%Prev)
    Cur=>Head
    DO ix=2,npart
       ALLOCATE(Cur%Next)
       Cur%Part_P=0.0_num
       Cur%Part_Pos=0.0_num
       !       Cur%Next%Prev=>Cur
       NULLIFY(Cur%Next%Next)
       Cur=>Cur%Next
       !This means that when finished, Tail will point to last element
       Tail=>Cur
    ENDDO

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
