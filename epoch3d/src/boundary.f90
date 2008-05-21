MODULE boundary

  USE shared_data
  USE particle_swap

  IMPLICIT NONE

CONTAINS

  SUBROUTINE DoField_MPI(Array)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: Array

    !Send left
    CALL MPI_SENDRECV(Array(1:2,1:ny,1:nz),2*ny*nz,mpireal,left,tag,Array(nx+1:nx+2,1:ny,1:nz),2*ny*nz&
         ,mpireal,right,tag,comm,status,errcode) 
    !Send right
    CALL MPI_SENDRECV(Array(nx-1:nx,1:ny,1:nz),2*ny*nz,mpireal,right,tag,Array(-1:0,1:ny,1:nz),2*ny*nz&
         ,mpireal,left,tag,comm,status,errcode) 

    !Send up
    CALL MPI_SENDRECV(Array(1:nx,1:2,1:nz),2*nx*nz,mpireal,up,tag,Array(1:nx,ny+1:ny+2,1:nz),2*nx*nz&
         ,mpireal,down,tag,comm,status,errcode) 
    !Send down
    CALL MPI_SENDRECV(Array(1:nx,ny-1:ny,1:nz),2*nx*nz,mpireal,down,tag,Array(1:nx,-1:0,1:nz),2*nx*nz&
         ,mpireal,up,tag,comm,status,errcode) 

    !Send front
    CALL MPI_SENDRECV(Array(1:nx,1:ny,1:2),2*nx*ny,mpireal,front,tag,Array(1:nx,1:ny,nz+1:nz+2),2*nx*ny&
         ,mpireal,back,tag,comm,status,errcode) 
    !Send back
    CALL MPI_SENDRECV(Array(1:nx,1:ny,nz-1:nz),2*nx*ny,mpireal,back,tag,Array(1:nx,1:ny,-1:0),2*nx*ny&
         ,mpireal,front,tag,comm,status,errcode) 

  END SUBROUTINE DoField_MPI

  SUBROUTINE Efield_bcs

    CALL DoField_MPI(Ex)
    CALL DoField_MPI(Ey)
    CALL DoField_MPI(Ez)

    IF (xbc_left == other .AND. left == MPI_PROC_NULL) THEN
       Ex(-1:0,:,:)=0.0_num
       DO iz=1,nz
          Ez(1,1:ny,iz)=las_amp*sin(las_omega * time) * exp(-((y(1:ny)-(y_start+y_end)/2.0_num)/(0.1_num * length_y))**2)&
               * exp(-((z(iz)-(z_start+z_end)/2.0_num)/(0.1_num * length_z))**2)
       ENDDO
       Ey(0,:,:)=Ey(1,:,:)

       Ey(-1:0,:,:)=0.0_num
    ENDIF

    IF (xbc_right == other .AND. right == MPI_PROC_NULL) THEN
       Ex(nx+1:nx+2,:,:)=0.0_num
       Ey(nx+1:nx+2,:,:)=0.0_num
       Ez(nx+1:nx+2,:,:)=0.0_num
    ENDIF


    IF (ybc_down == other .AND. down == MPI_PROC_NULL) THEN
       Ex(:,-1:0,:)=0.0_num
       Ey(:,-1:0,:)=0.0_num
       Ez(:,-1:0,:)=0.0_num
    ENDIF

    IF (ybc_up == other .AND. up == MPI_PROC_NULL) THEN
       Ex(:,ny+1:ny+2,:)=0.0_num
       Ey(:,ny+1:ny+2,:)=0.0_num
       Ez(:,ny+1:ny+2,:)=0.0_num
    ENDIF


    IF (zbc_back == other .AND. back == MPI_PROC_NULL) THEN
       Ex(:,-1:0,:)=0.0_num
       Ey(:,-1:0,:)=0.0_num
       Ez(:,-1:0,:)=0.0_num
    ENDIF

    IF (zbc_front == other .AND. front == MPI_PROC_NULL) THEN
       Ex(:,ny+1:ny+2,:)=0.0_num
       Ey(:,ny+1:ny+2,:)=0.0_num
       Ez(:,ny+1:ny+2,:)=0.0_num
    ENDIF


  END SUBROUTINE Efield_bcs

  SUBROUTINE Bfield_bcs

    CALL DoField_MPI(Bx)
    CALL DoField_MPI(By)
    CALL DoField_MPI(Bz)    

    IF (xbc_left == other .AND. left == MPI_PROC_NULL) THEN
       Bx(-1:0,:,:)=0.0_num
       By(-1:0,:,:)=0.0_num
       Bz(-1:0,:,:)=0.0_num
    ENDIF

    IF (xbc_right == other .AND. right == MPI_PROC_NULL) THEN
       Bx(nx+1:nx+2,:,:)=0.0_num
       By(nx+1:nx+2,:,:)=0.0_num
       Bz(nx+1:nx+2,:,:)=0.0_num
    ENDIF


    IF (ybc_down == other .AND. down == MPI_PROC_NULL) THEN
       Bx(:,-1:0,:)=0.0_num
       By(:,-1:0,:)=0.0_num
       Bz(:,-1:0,:)=0.0_num
    ENDIF

    IF (ybc_up == other .AND. up == MPI_PROC_NULL) THEN
       Bx(:,ny+1:ny+2,:)=0.0_num
       By(:,ny+1:ny+2,:)=0.0_num
       Bz(:,ny+1:ny+2,:)=0.0_num
    ENDIF


    IF (zbc_back == other .AND. back == MPI_PROC_NULL) THEN
       Bx(:,-1:0,:)=0.0_num
       By(:,-1:0,:)=0.0_num
       Bz(:,-1:0,:)=0.0_num
    ENDIF

    IF (zbc_front == other .AND. front == MPI_PROC_NULL) THEN
       Bx(:,ny+1:ny+2,:)=0.0_num
       By(:,ny+1:ny+2,:)=0.0_num
       Bz(:,ny+1:ny+2,:)=0.0_num
    ENDIF


  END SUBROUTINE Bfield_bcs

  SUBROUTINE Particle_bcs

    TYPE:: particlepointer

       TYPE(particle),POINTER :: Pointer

    END TYPE particlepointer

    TYPE(particle),POINTER :: Cur,last_good,Next
    TYPE(particlepointer),DIMENSION(-1:1,-1:1,-1:1) :: heads
    INTEGER(8),DIMENSION(-1:1,-1:1,-1:1) :: counts
    INTEGER :: xbd,ybd,zbd,ixp,iyp,izp,ct
    LOGICAL :: Out_Of_Bounds
    LOGICAL,DIMENSION(-1:1,-1:1,-1:1) :: done

!    PRINT *,"entering boundary"

    Cur=>Head
    NULLIFY(last_good)

    DO iz=-1,1
       DO iy=-1,1
          DO ix=-1,1
             NULLIFY(heads(ix,iy,iz)%Pointer)
             counts(ix,iy,iz)=0
          ENDDO
       ENDDO
    ENDDO

    DO WHILE (ASSOCIATED(Cur))
       Out_Of_Bounds=.FALSE.
       Next=>Cur%Next

       xbd=0
       ybd=0
       zbd=0

       !These conditions apply if a particle has passed a physical boundary
       !Not a processor boundary or a periodic boundary
       IF (Cur%Part_Pos(1) .LT. x_start .AND. left == MPI_PROC_NULL .AND. xbc_left == other) THEN
          !Particle has crossed left boundary
          cur%part_pos(1) =  x_start - (cur%part_pos(1)-x_start)
          IF (cur%part_pos(1) .LT. x_start) WRITE(10+rank,*) "BAD PARTICLE LOW X"
          Cur%part_p(1) = - Cur%part_p(1)
       ENDIF
       IF (Cur%Part_Pos(1) .GT. x_end .AND. right == MPI_PROC_NULL .AND. xbc_right == other) THEN
          !Particle has crossed right boundary
          Cur%part_pos(1) =  x_end - (Cur%part_pos(1)-x_end)
          IF (cur%part_pos(1) .GT. x_end) WRITE(10+rank,*) "BAD PARTICLE HIGH X"
          Cur%part_p(1) = - Cur%part_p(1)
       ENDIF
       IF (Cur%Part_Pos(2) .LT. y_start .AND. down == MPI_PROC_NULL .AND. ybc_down == other) THEN
          !Particle has crossed bottom boundary
          WRITE(10+rank,*)Cur%Part_Pos(2)
          cur%part_pos(2) =  y_start - (cur%part_pos(2)-y_start)
          IF (cur%part_pos(2) .LT. y_start) WRITE(10+rank,*) "BAD PARTICLE LOW Y",cur%part_pos(2),y_start,y(0)
          Cur%part_p(2) = - Cur%part_p(2)
       ENDIF
       IF (Cur%Part_Pos(2) .GT. y_end .AND. up == MPI_PROC_NULL .AND. ybc_up == other) THEN
          !Particle has crossed top boundary
          Cur%part_pos(2) =  y_end - (Cur%part_pos(2)-y_end)
          IF (cur%part_pos(2) .GT. y_end) WRITE(10+rank,*) "BAD PARTICLE HIGH Y"
          Cur%part_p(2) = - Cur%part_p(2)
       ENDIF
       IF (Cur%Part_Pos(3) .LT. z_start .AND. back == MPI_PROC_NULL .AND. zbc_back == other) THEN
          !Particle has crossed back boundary
          cur%part_pos(3) =  z_start - (cur%part_pos(3)-z_start)
          IF (cur%part_pos(3) .LT. z_start) WRITE(10+rank,*) "BAD PARTICLE LOW Z"
          Cur%part_p(3) = - Cur%part_p(3)
       ENDIF
       IF (Cur%Part_Pos(3) .GT. z_end .AND. front == MPI_PROC_NULL .AND. zbc_front == other) THEN
          Cur%part_pos(3) =  z_end - (Cur%part_pos(3)-z_end)
          IF (cur%part_pos(3) .GT. z_end) WRITE(10+rank,*) "BAD PARTICLE HIGH Y"
          Cur%part_p(3) = - Cur%part_p(3)
       ENDIF

       IF (Cur%Part_Pos(1) .LT. x_start_local) xbd=-1
       IF (Cur%Part_Pos(1) .GT. x_end_local)   xbd=1
       IF (Cur%Part_Pos(2) .LT. y_start_local) ybd=-1
       IF (Cur%Part_Pos(2) .GT. y_end_local)   ybd=1
       IF (Cur%Part_Pos(3) .LT. z_start_local) zbd=-1
       IF (Cur%Part_Pos(3) .GT. z_end_local)   zbd=1

       IF (xbd .GT. 0 .AND. right == MPI_PROC_NULL) PRINT *,"Problem with particle"

       IF (ABS(xbd) + ABS(ybd) + ABS(zbd) .GT. 0) THEN
          !Particle has left box
          counts(xbd,ybd,zbd) = counts(xbd,ybd,zbd)+1
          Cur%Next=>heads(xbd,ybd,zbd)%Pointer
          heads(xbd,ybd,zbd)%Pointer=>Cur
          IF (ASSOCIATED(last_good)) THEN
             last_good%Next=>Next
          ELSE
             Head=>Next
          ENDIF
          !If this is the last particle then set the tail marker to the last known good particle
          IF (ASSOCIATED(Tail,TARGET=Cur)) Tail=>Last_good
       ELSE
          last_good=>Cur
       ENDIF

       !Move to next particle
       Cur=>Next
    ENDDO


    !swap Particles
    Done=.FALSE.
    DO iz=-1,1
       DO iy=-1,1
          DO ix=-1,1
             IF (ABS(ix)+ABS(iy)+ABS(iz) .EQ. 0) CYCLE
             ixp=-ix
             iyp=-iy
             izp=-iz
!!$             IF (counts(ix,iy,iz) .NE. 0) WRITE(rank+10,*) "Sending",counts(ix,iy,iz),"to rank",neighbour(ix,iy,iz),"coords",coordinates+(/ix,iy,iz/)
!!$             IF (counts(ixp,iyp,izp) .NE. 0) WRITE(rank+10,*) "Sending",counts(ixp,iyp,izp),"to rank",neighbour(ixp,iyp,izp),"coords",coordinates+(/ixp,iyp,izp/)
             IF (counts(ix,iy,iz) .NE. 0 .AND. neighbour(ix,iy,iz) .EQ. MPI_PROC_NULL) WRITE(rank+10,*) "***ERROR*** Sending particle to nowhere"

             CALL Swap_Particles(Heads(ixp,iyp,izp)%Pointer,Heads(ix,iy,iz)%Pointer,counts(ixp,iyp,izp),counts(ix,iy,iz),neighbour(ixp,iyp,izp),neighbour(ix,iy,iz))
          ENDDO
       ENDDO
    ENDDO

!    PRINT *,"All swaps done",rank



    !Particles should only lie outside boundaries if the periodic boundaries are turned on
    !This now moves them to within the boundaries
    Cur=>Head
    ct=0
    DO WHILE(ASSOCIATED(Cur))
       !       IF (rank .EQ. 0) PRINT *,Cur%label
       IF(Cur%Part_Pos(1) .GT. x_end .AND. xbc_left == periodic) Cur%Part_Pos(1)=Cur%Part_Pos(1)-length_x
       IF(Cur%Part_Pos(1) .LT. x_start .AND. xbc_right == periodic) Cur%Part_Pos(1)=Cur%Part_Pos(1)+length_x
       IF(Cur%Part_Pos(2) .GT. y_end) Cur%Part_Pos(2)=Cur%Part_Pos(2)-length_y
       IF(Cur%Part_Pos(2) .LT. y_start) Cur%Part_Pos(2)=Cur%Part_Pos(2)+length_y
       IF(Cur%Part_Pos(3) .GT. z_end) Cur%Part_Pos(3)=Cur%Part_Pos(3)-length_z
       IF(Cur%Part_Pos(3) .LT. z_start) Cur%Part_Pos(3)=Cur%Part_Pos(3)+length_z

       Cur=>Cur%Next
    ENDDO

    !PRINT *,"Particle bcs_done",rank

  END SUBROUTINE Particle_bcs

END MODULE boundary
