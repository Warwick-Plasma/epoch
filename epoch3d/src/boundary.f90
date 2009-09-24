MODULE boundary

  USE shared_data
  USE partlist
  USE shared_parser_data
  IMPLICIT NONE

CONTAINS

  SUBROUTINE Setup_Particle_Boundaries

    AnyOpen=.FALSE.

    !For some types of boundary, fields and particles are treated in different ways, deal with that here
    IF (xbc_right == BC_PERIODIC) THEN
       xbc_right_particle=BC_PERIODIC
       xbc_right_field=BC_PERIODIC
    ENDIF
    IF (xbc_left == BC_PERIODIC) THEN
       xbc_left_particle=BC_PERIODIC
       xbc_left_field=BC_PERIODIC
    ENDIF

    IF (ybc_up == BC_PERIODIC) THEN
       ybc_up_particle=BC_PERIODIC
       ybc_up_field=BC_PERIODIC
    ENDIF
    IF (ybc_down == BC_PERIODIC) THEN
       ybc_down_particle=BC_PERIODIC
       ybc_down_field=BC_PERIODIC
    ENDIF

    IF (zbc_front == BC_PERIODIC) THEN
       zbc_front_particle=BC_PERIODIC
       zbc_front_field=BC_PERIODIC
    ENDIF
    IF (zbc_back == BC_PERIODIC) THEN
       zbc_back_particle=BC_PERIODIC
       zbc_back_field=BC_PERIODIC
    ENDIF

    !For some types of boundary, fields and particles are treated in different ways, deal with that here
    IF (xbc_right == BC_OTHER) THEN
       xbc_right_particle=BC_REFLECT
       xbc_right_field=BC_CLAMP
    ENDIF
    IF (xbc_left == BC_OTHER) THEN
       xbc_left_particle=BC_REFLECT
       xbc_left_field=BC_CLAMP
    ENDIF

    IF (ybc_up == BC_OTHER) THEN
       ybc_up_particle=BC_REFLECT
       ybc_up_field=BC_CLAMP
    ENDIF
    IF (ybc_down == BC_OTHER) THEN
       ybc_down_particle=BC_REFLECT
       ybc_down_field=BC_CLAMP
    ENDIF

    IF (zbc_front == BC_OTHER) THEN
       zbc_front_particle=BC_REFLECT
       zbc_front_field=BC_CLAMP
    ENDIF
    IF (zbc_back == BC_OTHER) THEN
       zbc_back_particle=BC_REFLECT
       zbc_back_field=BC_CLAMP
    ENDIF

    !Laser boundaries reflect particles off a hard wall
    IF (xbc_left == BC_SIMPLE_LASER .OR. xbc_left == BC_SIMPLE_OUTFLOW) THEN
       xbc_left_particle=BC_OPEN
       xbc_left_field=BC_ZERO_GRADIENT
       AnyOpen=.TRUE.
    ENDIF

    IF (xbc_right == BC_SIMPLE_LASER .OR. xbc_right == BC_SIMPLE_OUTFLOW) THEN
       xbc_right_particle=BC_OPEN
       xbc_left_field=BC_ZERO_GRADIENT
       AnyOpen=.TRUE.
    ENDIF

    IF (ybc_up == BC_SIMPLE_LASER .OR. ybc_up == BC_SIMPLE_OUTFLOW) THEN
       ybc_up_particle=BC_OPEN
       ybc_up_field=BC_ZERO_GRADIENT
       AnyOpen=.TRUE.
    ENDIF

    IF (ybc_down == BC_SIMPLE_LASER .OR. ybc_down == BC_SIMPLE_OUTFLOW) THEN
       ybc_down_particle=BC_OPEN
       ybc_down_field=BC_ZERO_GRADIENT
       AnyOpen=.TRUE.
    ENDIF

    IF (zbc_front == BC_SIMPLE_LASER .OR. zbc_front == BC_SIMPLE_OUTFLOW) THEN
       zbc_front_particle=BC_OPEN
       zbc_front_field=BC_ZERO_GRADIENT
       AnyOpen=.TRUE.
    ENDIF

    IF (zbc_back == BC_SIMPLE_LASER .OR. zbc_back == BC_SIMPLE_OUTFLOW) THEN
       zbc_back_particle=BC_OPEN
       zbc_back_field=BC_ZERO_GRADIENT
       AnyOpen=.TRUE.
    ENDIF

    IF (AnyOpen) CALL Create_Empty_Partlist(Ejected_Particles)

  END SUBROUTINE Setup_Particle_Boundaries

  !Exchanges field values at processor boundaries and applies field boundary conditions
  SUBROUTINE Field_BC(Field)

    REAL(num),DIMENSION(-2:,-2:,-2:),INTENT(INOUT) :: Field

    CALL Do_Field_MPI_With_Lengths(Field,nx,ny,nz)

  END SUBROUTINE Field_BC

  SUBROUTINE Do_Field_MPI_With_Lengths(Field,nx_local,ny_local,nz_local)

    REAL(num),DIMENSION(-2:,-2:,-2:),INTENT(INOUT) :: Field
    INTEGER,INTENT(IN) :: nx_local, ny_local, nz_local

    CALL MPI_SENDRECV(Field(1:3,:,:),3*(ny_local+6)*(nz_local+6),mpireal,left,tag,&
         Field(nx_local+1:nx_local+3,:,:),3*(ny_local+6)*(nz_local+6),mpireal,&
         right,tag,comm,status,errcode)
    CALL MPI_SENDRECV(Field(nx_local-2:nx_local,:,:),3*(ny_local+6)*(nz_local+6),&
         mpireal,right,tag,Field(-2:0,:,:),3*(ny_local+6)*(nz_local+6),mpireal,&
         left,tag,comm,status,errcode)

    CALL MPI_SENDRECV(Field(:,1:3,:),3*(nx_local+6)*(nz_local+6),mpireal,down,tag,&
         Field(:,ny_local+1:ny_local+3,:),3*(nx_local+6)*(nz_local+6),mpireal,up,&
         tag,comm,status,errcode)
    CALL MPI_SENDRECV(Field(:,ny_local-2:ny_local,:),3*(nx_local+6)*(nz_local+6),&
         mpireal,up,tag,Field(:,-2:0,:),3*(nx_local+6)*(nz_local+6),mpireal,down&
         ,tag,comm,status,errcode)

    CALL MPI_SENDRECV(Field(:,:,1:3),3*(nx_local+6)*(ny_local+6),mpireal,back,tag,&
         Field(:,:,nz_local+1:nz_local+3),3*(nx_local+6)*(ny_local+6),mpireal,front,&
         tag,comm,status,errcode)
    CALL MPI_SENDRECV(Field(:,:,nz_local-2:nz_local),3*(nx_local+6)*(ny_local+6),&
         mpireal,front,tag,Field(:,:,-2:0),3*(nx_local+6)*(ny_local+6),mpireal,back&
         ,tag,comm,status,errcode)

  END SUBROUTINE Do_Field_MPI_With_Lengths

  SUBROUTINE Field_Zero_Gradient(Field,Force)

    REAL(num),DIMENSION(-2:,-2:,-2:),INTENT(INOUT) :: Field
    LOGICAL,INTENT(IN) :: Force


    IF ((xbc_left_field == BC_ZERO_GRADIENT .OR. Force) .AND. left == MPI_PROC_NULL) THEN
       Field(0,:,:)=Field(1,:,:)
    ENDIF

    IF ((xbc_right_field == BC_ZERO_GRADIENT .OR. Force)  .AND. right == MPI_PROC_NULL) THEN
       Field(nx+1,:,:)=Field(nx,:,:)
    ENDIF

    IF ((ybc_down_field == BC_ZERO_GRADIENT .OR. Force) .AND. down == MPI_PROC_NULL) THEN
       Field(:,0,:)=Field(:,1,:)
    ENDIF

    IF ((ybc_up_field == BC_ZERO_GRADIENT .OR. Force) .AND. up == MPI_PROC_NULL) THEN
       Field(:,ny+1,:)=Field(:,ny,:)
    ENDIF

    IF ((zbc_back_field == BC_ZERO_GRADIENT .OR. Force) .AND. back == MPI_PROC_NULL) THEN
       Field(:,:,0)=Field(:,:,1)
    ENDIF

    IF ((zbc_front_field == BC_ZERO_GRADIENT .OR. Force) .AND. front == MPI_PROC_NULL) THEN
       Field(:,:,nz+1)=Field(:,:,nz)
    ENDIF


  END SUBROUTINE Field_Zero_Gradient

  SUBROUTINE Field_Clamp_Zero(Field)

    REAL(num),DIMENSION(-2:,-2:,-2:),INTENT(INOUT) :: Field


    IF (xbc_left_field == BC_CLAMP .AND. left == MPI_PROC_NULL) THEN
       Field(0,:,:)=0.0_num
    ENDIF

    IF (xbc_right_field == BC_CLAMP .AND. right == MPI_PROC_NULL) THEN
       Field(nx+1,:,:)=0.0_num
    ENDIF

    IF (ybc_down_field == BC_CLAMP .AND. down == MPI_PROC_NULL) THEN
       Field(:,0,:)=0.0_num
    ENDIF

    IF (ybc_up_field == BC_CLAMP .AND. up == MPI_PROC_NULL) THEN
       Field(:,ny+1,:)=0.0_num
    ENDIF

    IF (zbc_back_field == BC_CLAMP .AND. back == MPI_PROC_NULL) THEN
       Field(:,:,0)=0.0_num
    ENDIF

    IF (zbc_front_field == BC_CLAMP .AND. front == MPI_PROC_NULL) THEN
       Field(:,:,nz+1)=0.0_num
    ENDIF


  END SUBROUTINE Field_Clamp_Zero

  SUBROUTINE Processor_Summation_BCS(array)

    REAL(num),DIMENSION(-2:,-2:,-2:),INTENT(INOUT) :: array
    REAL(num),DIMENSION(:,:,:),ALLOCATABLE :: temp
    INTEGER :: nxp,nyp,nzp

    INTEGER,DIMENSION(-1:1,-1:1,-1:1) :: Sizes,xStart,xEnd,yStart,yEnd,zStart,zEnd,xshift,yshift,zshift
    INTEGER :: xs,xe,ys,ye,zs,ze,xf,yf,zf

    nxp=nx+1
    nyp=ny+1
    nzp=nz+1

    Sizes=0
    xStart=0
    yStart=0
    zStart=0
    xEnd=0
    yEnd=0
    zEnd=0
    xshift=0
    yshift=0
    zshift=0

    DO iz=-1,1
       DO iy=-1,1
          DO ix=-1,1
             Sizes(ix,iy,iz)=1
             IF (ix == 0) THEN
                Sizes(ix,iy,iz) = Sizes(ix,iy,iz) * (nx+6)
                xStart(ix,iy,iz) = -2
                xEnd(ix,iy,iz) = nx+3
             ELSE IF (ix==1) THEN
                Sizes(ix,iy,iz)=Sizes(ix,iy,iz) * 3
                xStart(ix,iy,iz) = nx+1
                xEnd(ix,iy,iz) = nx+3
                xShift(ix,iy,iz) = -nx
             ELSE
                Sizes(ix,iy,iz)=Sizes(ix,iy,iz) * 3
                xStart(ix,iy,iz) = -2
                xEnd(ix,iy,iz) = 0
                xShift(ix,iy,iz) = nx
             ENDIF
             IF (iy == 0) THEN
                Sizes(ix,iy,iz) = Sizes(ix,iy,iz) * (ny+6)
                yStart(ix,iy,iz) = -2
                yEnd(ix,iy,iz) = ny+3
             ELSE IF (iy==1) THEN
                Sizes(ix,iy,iz)=Sizes(ix,iy,iz) * 3
                yStart(ix,iy,iz) = ny+1
                yEnd(ix,iy,iz) = ny+3
                yShift(ix,iy,iz) = -ny
             ELSE
                Sizes(ix,iy,iz)=Sizes(ix,iy,iz) * 3
                yStart(ix,iy,iz) = -2
                yEnd(ix,iy,iz) = 0
                yShift(ix,iy,iz) = ny
             ENDIF
             IF (iz == 0) THEN
                Sizes(ix,iy,iz) = Sizes(ix,iy,iz) * (nz+6)
                zStart(ix,iy,iz) = -2
                zEnd(ix,iy,iz) = nz+3
             ELSE IF (iz==1) THEN
                Sizes(ix,iy,iz)=Sizes(ix,iy,iz) * 3
                zStart(ix,iy,iz) = nz+1
                zEnd(ix,iy,iz) = nz+3
                zShift(ix,iy,iz) = -nz
             ELSE
                Sizes(ix,iy,iz)=Sizes(ix,iy,iz) * 3
                zStart(ix,iy,iz) = -2
                zEnd(ix,iy,iz) = 0
                zShift(ix,iy,iz) = nz
             ENDIF

          ENDDO
       ENDDO
    ENDDO


    DO iz=-1,1
       DO iy=-1,1
          DO ix=-1,1
             IF (ix == 0 .AND. iy == 0 .AND. iz ==0) CYCLE
             IF (ABS(ix)+ABS(iy)+ABS(iz) .GT. 1) CYCLE
!!$             IF (sizes(ix,iy,iz) .EQ. 0) THEN
!!$                WRITE(rank+10,*) "Zero size",ix,iy,iz
!!$                CYCLE
!!$             ENDIF
             !Copy the starts into variables with shorter names, or this is HORRIFIC to read
             xs=xstart(ix,iy,iz)
             ys=ystart(ix,iy,iz)
             zs=zstart(ix,iy,iz)
             xe=xend(ix,iy,iz)
             ye=yend(ix,iy,iz)
             ze=zend(ix,iy,iz)
             xf=xshift(ix,iy,iz)
             yf=yshift(ix,iy,iz)
             zf=zshift(ix,iy,iz)
             ALLOCATE(temp(xs:xe,ys:ye,zs:ze))
             temp=0.0_num
             CALL MPI_SENDRECV(Array(xs:xe,ys:ye,zs:ze),sizes(ix,iy,iz),mpireal,neighbour(ix,iy,iz),tag,temp,sizes(-ix,-iy,-iz),&
                  mpireal,neighbour(-ix,-iy,-iz),tag,comm,status,errcode)
             Array(xs+xf:xe+xf,ys+yf:ye+yf,zs+zf:ze+zf) = Array(xs+xf:xe+xf,ys+yf:ye+yf,zs+zf:ze+zf) + temp
             DEALLOCATE(temp)
          ENDDO
       ENDDO
    ENDDO
    CALL Field_BC(Array)

  END SUBROUTINE Processor_Summation_BCS

  SUBROUTINE Efield_bcs

    !These are the MPI boundaries
    CALL Field_BC(Ex)
    CALL Field_BC(Ey)
    CALL Field_BC(Ez)

    !These apply zero field boundary conditions on the edges
    CALL Field_Clamp_Zero(Ex)
    CALL Field_Clamp_Zero(Ey)
    CALL Field_Clamp_Zero(Ez)
    !These apply zero field gradient boundary conditions on the edges
    CALL Field_Zero_Gradient(Ex,.FALSE.)
    CALL Field_Zero_Gradient(Ey,.FALSE.)
    CALL Field_Zero_Gradient(Ez,.FALSE.)

  END SUBROUTINE Efield_bcs

  SUBROUTINE Bfield_bcs(MPI_Only)

    LOGICAL, INTENT(IN) :: MPI_Only

    !These are the MPI boundaries
    CALL Field_BC(Bx)
    CALL Field_BC(By)
    CALL Field_BC(Bz)

    IF (.NOT. MPI_Only) THEN
       !These apply zero field boundary conditions on the edges
       CALL Field_Clamp_Zero(Bx)
       CALL Field_Clamp_Zero(By)
       CALL Field_Clamp_Zero(Bz)
       !These apply zero field boundary conditions on the edges
       CALL Field_Zero_Gradient(Bx,.FALSE.)
       CALL Field_Zero_Gradient(By,.FALSE.)
       CALL Field_Zero_Gradient(Bz,.FALSE.)
    ENDIF

  END SUBROUTINE Bfield_bcs

  SUBROUTINE Particle_bcs

    TYPE(particle),POINTER :: Cur,last_good,Next
    TYPE(particlelist),DIMENSION(-1:1,-1:1,-1:1) :: send,recv
    INTEGER :: xbd,ybd,zbd,ixp,iyp,izp
    LOGICAL :: Out_Of_Bounds
    LOGICAL,DIMENSION(-1:1,-1:1,-1:1) :: done
    INTEGER :: iSpecies

    DO iSpecies=1,nspecies
       Cur=>ParticleSpecies(iSpecies)%AttachedList%Head
       NULLIFY(last_good)

       DO iz=-1,1
          DO iy=-1,1
             DO ix=-1,1
                CALL Create_Empty_Partlist(send(ix,iy,iz))
                CALL Create_Empty_Partlist(recv(ix,iy,iz))
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
          IF (Cur%Part_Pos(1) .LE. x_start-dx/2.0_num .AND. left == MPI_PROC_NULL .AND. xbc_left_particle == BC_REFLECT) THEN
             !Particle has crossed left boundary
             cur%part_pos(1) =  2.0_num * (x_start-dx/2.0_num) - cur%part_pos(1)
             Cur%part_p(1) = - Cur%part_p(1)
          ENDIF
          IF (Cur%Part_Pos(1) .GE. x_end+dx/2.0_num .AND. right == MPI_PROC_NULL .AND. xbc_right_particle == BC_REFLECT) THEN
             !Particle has crossed right boundary
             Cur%part_pos(1) =  2.0_num *(x_end+dx/2.0_num) - Cur%part_pos(1)
             Cur%part_p(1) = - Cur%part_p(1)
          ENDIF
          IF (Cur%Part_Pos(2) .LE. y_start-dy/2.0_num .AND. down == MPI_PROC_NULL .AND. ybc_down_particle == BC_REFLECT) THEN
             !Particle has crossed bottom boundary
             cur%part_pos(2) =  2.0_num * (y_start-dy/2.0_num) - cur%part_pos(2)
             Cur%part_p(2) = - Cur%part_p(2)
          ENDIF
          IF (Cur%Part_Pos(2) .GE. y_end+dy/2.0_num .AND. up == MPI_PROC_NULL .AND. ybc_up_particle == BC_REFLECT) THEN
             !          PRINT *,"Reflecting"
             !Particle has crossed top boundary
             Cur%part_pos(2) =  2.0_num * (y_end + dy/2.0_num) - Cur%part_pos(2)
             !IF (cur%part_pos(2) .GT. y_end) WRITE(10+rank,*) "BAD PARTICLE HIGH Y"
             Cur%part_p(2) = - Cur%part_p(2)
          ENDIF
          IF (Cur%Part_Pos(3) .LT. z_start+dz/2.0_num .AND. back == MPI_PROC_NULL .AND. zbc_back_particle == BC_OTHER) THEN
             !Particle has crossed back boundary
             Cur%part_pos(3) =  2.0_num * (z_start-dz/2.0_num) - cur%part_pos(3)
             Cur%part_p(3) = - Cur%part_p(3)
          ENDIF
          IF (Cur%Part_Pos(3) .GT. z_end+dz/2.0_num .AND. front == MPI_PROC_NULL .AND. zbc_front_particle == BC_OTHER) THEN
             !Particle has crossed front boundary
             Cur%part_pos(3) =  2.0_num * (z_end + dz/2.0_num) - Cur%part_pos(2)
             Cur%part_p(3) = - Cur%part_p(3)
          ENDIF

          IF (Cur%Part_Pos(1) .LT. x_start_local - dx/2.0_num) xbd=-1
          IF (Cur%Part_Pos(1) .GT. x_end_local + dx/2.0_num )  xbd=1
          IF (Cur%Part_Pos(2) .LT. y_start_local - dy/2.0_num) ybd=-1
          IF (Cur%Part_Pos(2) .GT. y_end_local + dy/2.0_num)   ybd=1
          IF (Cur%Part_Pos(3) .LT. z_start_local - dz/2.0_num) zbd=-1
          IF (Cur%Part_Pos(3) .GT. z_end_local + dz/2.0_num)   zbd=1

          IF ((Cur%Part_Pos(1) .LT. x_start - dx/2.0_num) .AND. (xbc_left_particle == BC_OPEN)) out_of_bounds = .TRUE.
          IF ((Cur%Part_Pos(1) .GT. x_end + dx/2.0_num) .AND. (xbc_right_particle == BC_OPEN)) out_of_bounds = .TRUE.
          IF ((Cur%Part_Pos(2) .LT. y_start - dy/2.0_num) .AND. (ybc_down_particle == BC_OPEN)) out_of_bounds = .TRUE.
          IF ((Cur%Part_Pos(2) .GT. y_end +dy/2.0_num) .AND. (ybc_up_particle == BC_OPEN)) out_of_bounds = .TRUE.
          IF ((Cur%Part_Pos(3) .LT. z_start - dz/2.0_num) .AND. (zbc_back_particle == BC_OPEN)) out_of_bounds = .TRUE.
          IF ((Cur%Part_Pos(3) .GT. z_end +dz/2.0_num) .AND. (zbc_front_particle == BC_OPEN)) out_of_bounds = .TRUE.

          IF (ABS(xbd) + ABS(ybd) + ABS(zbd) .GT. 0) THEN
             !Particle has left box
             CALL Remove_Particle_From_PartList(ParticleSpecies(iSpecies)%AttachedList, Cur)
             IF (.NOT. out_of_bounds) THEN
                CALL Add_Particle_To_PartList(send(xbd,ybd,zbd),Cur)
             ELSE
                !CALL Add_Particle_To_PartList(ejected_particles,Cur)
                DEALLOCATE(Cur)
             ENDIF
          ENDIF

          !Move to next particle
          Cur=>Next
       ENDDO



       !swap Particles
       Done=.FALSE.
       DO iz=-1,1
          DO iy=-1,1
             DO ix=-1,1
                IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
                ixp=-ix
                iyp=-iy
                izp=-iz
                CALL PartList_SendRecv(send(ix,iy,iz),recv(ixp,iyp,izp),Neighbour(ix,iy,iz),Neighbour(ixp,iyp,izp))
                CALL Append_PartList(ParticleSpecies(ispecies)%AttachedList,recv(ixp,iyp,izp))
             ENDDO
          ENDDO
       ENDDO


       !Particles should only lie outside boundaries if the periodic boundaries are turned on
       !This now moves them to within the boundaries
       Cur=>ParticleSpecies(iSpecies)%AttachedList%Head
       ct=0
       DO WHILE(ASSOCIATED(Cur))
          IF(Cur%Part_Pos(1) .GT. x_end+dx/2.0_num .AND. xbc_left_particle == BC_PERIODIC) &
               Cur%Part_Pos(1)=Cur%Part_Pos(1)-length_x - dx
          IF(Cur%Part_Pos(1) .LT. x_start-dx/2.0_num .AND. xbc_right_particle == BC_PERIODIC) &
               Cur%Part_Pos(1)=Cur%Part_Pos(1)+length_x + dx
          IF(Cur%Part_Pos(2) .GT. y_end+dy/2.0_num .AND. ybc_up_particle == BC_PERIODIC) &
               Cur%Part_Pos(2)=Cur%Part_Pos(2)-length_y - dy
          IF(Cur%Part_Pos(2) .LT. y_start-dy/2.0_num .AND. ybc_down_particle == BC_PERIODIC) &
               Cur%Part_Pos(2)=Cur%Part_Pos(2)+length_y + dy
          IF(Cur%Part_Pos(3) .GT. z_end+dz/2.0_num .AND. zbc_front_particle == BC_PERIODIC) &
               Cur%Part_Pos(3)=Cur%Part_Pos(3)-length_z - dz
          IF(Cur%Part_Pos(3) .LT. z_start-dz/2.0_num .AND. zbc_back_particle == BC_PERIODIC) &
               Cur%Part_Pos(3)=Cur%Part_Pos(3)+length_z + dz
          Cur=>Cur%Next
       ENDDO
    ENDDO
    !PRINT *,"Particle bcs_done",rank

  END SUBROUTINE Particle_bcs


END MODULE boundary
