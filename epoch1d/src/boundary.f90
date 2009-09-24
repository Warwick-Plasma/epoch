
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

    !For some types of boundary, fields and particles are treated in different ways, deal with that here
    IF (xbc_right == BC_OTHER) THEN
       xbc_right_particle=BC_REFLECT
       xbc_right_field=BC_CLAMP
    ENDIF
    IF (xbc_left == BC_OTHER) THEN
       xbc_left_particle=BC_REFLECT
       xbc_left_field=BC_CLAMP
    ENDIF

    !Laser boundaries reflect particles off a hard wall
    IF (xbc_left == BC_SIMPLE_LASER .OR. xbc_left == BC_SIMPLE_OUTFLOW) THEN
       xbc_left_particle=BC_OPEN
       xbc_left_field=BC_CLAMP
       AnyOpen=.TRUE.
    ENDIF

    IF (xbc_right == BC_SIMPLE_LASER .OR. xbc_right == BC_SIMPLE_OUTFLOW) THEN
       xbc_right_particle=BC_OPEN
       xbc_right_field=BC_CLAMP
       AnyOpen=.TRUE.
    ENDIF

    IF (AnyOpen) CALL Create_Empty_Partlist(Ejected_Particles)

  END SUBROUTINE Setup_Particle_Boundaries

  !Exchanges field values at processor boundaries and applies field boundary conditions
  SUBROUTINE Field_BC(Field)

    REAL(num),DIMENSION(-2:),INTENT(INOUT) :: Field

    CALL Do_Field_MPI_With_Lengths(Field,nx)

  END SUBROUTINE Field_BC

  SUBROUTINE Do_Field_MPI_With_Lengths(Field,nx_local)

    REAL(num),DIMENSION(-2:),INTENT(INOUT) :: Field
    INTEGER,INTENT(IN) :: nx_local

    CALL MPI_SENDRECV(Field(1:3),3,mpireal,left,tag,&
         Field(nx_local+1:nx_local+3),3,mpireal,&
         right,tag,comm,status,errcode)
    CALL MPI_SENDRECV(Field(nx_local-2:nx_local),3,&
         mpireal,right,tag,Field(-2:0),3,mpireal,&
         left,tag,comm,status,errcode)

  END SUBROUTINE Do_Field_MPI_With_Lengths

  SUBROUTINE Field_Zero_Gradient(Field,Force)

    REAL(num),DIMENSION(-2:),INTENT(INOUT) :: Field
    LOGICAL,INTENT(IN) :: Force


    IF ((xbc_left_field == BC_ZERO_GRADIENT .OR. Force) .AND. left == MPI_PROC_NULL) THEN
       Field(0)=Field(1)
       Field(-1)=Field(2)
    ENDIF

    IF ((xbc_right_field == BC_ZERO_GRADIENT .OR. Force)  .AND. right == MPI_PROC_NULL) THEN
       Field(nx+1)=Field(nx)
       Field(nx+2)=Field(nx-1)
    ENDIF

  END SUBROUTINE Field_Zero_Gradient

  SUBROUTINE Field_Clamp_Zero(Field,Stagger)

    REAL(num),DIMENSION(-2:),INTENT(INOUT) :: Field
    INTEGER,DIMENSION(1),INTENT(IN) :: Stagger


    IF (xbc_left_field == BC_CLAMP .AND. left == MPI_PROC_NULL) THEN
       IF (stagger(1) .EQ. 1) THEN
          Field(0)=0.0_num
          Field(-1)=-Field(1)
       ELSE
          Field(0)=-Field(1)
          Field(-1)=-Field(2)
       ENDIF
    ENDIF

    IF (xbc_right_field == BC_CLAMP .AND. right == MPI_PROC_NULL) THEN
       IF (stagger(1) .EQ. 1) THEN
          Field(nx)=0.0_num
       ELSE
          Field(nx+1)=-Field(nx)
          Field(nx+2)=-Field(nx-1)
       ENDIF
    ENDIF

  END SUBROUTINE Field_Clamp_Zero

  SUBROUTINE Processor_Summation_BCS(array)

    REAL(num),DIMENSION(-2:),INTENT(INOUT) :: array
    REAL(num),DIMENSION(:),ALLOCATABLE :: temp
    INTEGER :: nxp,nyp,nzp

    INTEGER,DIMENSION(-1:1) :: Sizes,xStart,xEnd,xshift
    INTEGER :: xs,xe,ys,ye,xf,yf

    nxp=nx+3

    Sizes=0
    xStart=0
    xEnd=0
    xshift=0

       DO ix=-1,1
          Sizes(ix)=1
          IF (ix == 0) THEN
             Sizes(ix) = Sizes(ix) * (nx+6)
             xStart(ix) = -2
             xEnd(ix) = nx+3
          ELSE IF (ix==1) THEN
             Sizes(ix) = Sizes(ix) *3
             xStart(ix) = nx+1
             xEnd(ix) = nx+3
             xShift(ix) = -nx
          ELSE
             Sizes(ix) = Sizes(ix) *3
             xStart(ix) = -2
             xEnd(ix) = 0
             xShift(ix) = nx
          ENDIF
       ENDDO

!!$    DO iy=-1,1
!!$       DO ix=-1,1
!!$          sizes(ix)=(xEnd(ix)-xStart(ix)+1)*(yEnd(ix)-yStart(ix)+1)
!!$       ENDDO
!!$    ENDDO

    !PRINT *,rank,ystart(0,1),yend(0,1),yshift(0,1)

       DO ix=-1,1
          IF (ix == 0 ) CYCLE
          !Copy the starts into variables with shorter names, or this is HORRIFIC to read
          xs=xstart(ix)
          xe=xend(ix)
          xf=xshift(ix)
          ALLOCATE(temp(xs:xe))
          temp=0.0_num
          !IF (neighbour(ix) .EQ. MPI_PROC_NULL) PRINT *,"BAD NEIGHBOUR",ix,iy
          CALL MPI_SENDRECV(Array(xs:xe),sizes(ix),mpireal,neighbour(ix),tag,temp,sizes(-ix),&
               mpireal,neighbour(-ix),tag,comm,status,errcode)
          Array(xs+xf:xe+xf) = Array(xs+xf:xe+xf) + temp
          DEALLOCATE(temp)
       ENDDO

    CALL Field_BC(Array)

  END SUBROUTINE Processor_Summation_BCS

  SUBROUTINE Efield_bcs

    !These are the MPI boundaries
    CALL Field_BC(Ex)
    CALL Field_BC(Ey)
    CALL Field_BC(Ez)

    !These apply zero field boundary conditions on the edges
    CALL Field_Clamp_Zero(Ex,(/1/))
    CALL Field_Clamp_Zero(Ey,(/0/))
    CALL Field_Clamp_Zero(Ez,(/0/))
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
       CALL Field_Clamp_Zero(Bx,(/0/))
       CALL Field_Clamp_Zero(By,(/1/))
       CALL Field_Clamp_Zero(Bz,(/1/))
       !These apply zero field boundary conditions on the edges
       CALL Field_Zero_Gradient(Bx,.FALSE.)
       CALL Field_Zero_Gradient(By,.FALSE.)
       CALL Field_Zero_Gradient(Bz,.FALSE.)
    ENDIF

  END SUBROUTINE Bfield_bcs

  SUBROUTINE Particle_bcs

    TYPE(particle),POINTER :: Cur,last_good,Next
    TYPE(particlelist),DIMENSION(-1:1) :: send,recv
    INTEGER :: xbd,ixp
    LOGICAL :: Out_Of_Bounds
    LOGICAL,DIMENSION(-1:1) :: done
    INTEGER :: iSpecies

    DO iSpecies=1,nspecies
       Cur=>ParticleSpecies(iSpecies)%AttachedList%Head
       NULLIFY(last_good)


          DO ix=-1,1
             CALL Create_Empty_Partlist(send(ix))
             CALL Create_Empty_Partlist(recv(ix))
          ENDDO

       DO WHILE (ASSOCIATED(Cur))
          Out_Of_Bounds=.FALSE.
          Next=>Cur%Next

          xbd=0

          !These conditions apply if a particle has passed a physical boundary
          !Not a processor boundary or a periodic boundary
          IF (Cur%Part_Pos .LE. x_start-dx/2.0_num .AND. left == MPI_PROC_NULL .AND. xbc_left_particle == BC_REFLECT) THEN
             !Particle has crossed left boundary
             cur%part_pos =  2.0_num * (x_start-dx/2.0_num) - cur%part_pos
             Cur%part_p(1) = - Cur%part_p(1)
          ENDIF
          IF (Cur%Part_Pos .GE. x_end+dx/2.0_num .AND. right == MPI_PROC_NULL .AND. xbc_right_particle == BC_REFLECT) THEN
             !Particle has crossed right boundary
             Cur%part_pos =  2.0_num *(x_end+dx/2.0_num) - Cur%part_pos
             Cur%part_p(1) = - Cur%part_p(1)
          ENDIF

          IF (Cur%part_pos .LT. x_start_local - dx/2.0_num) xbd=-1
          IF (Cur%part_pos .GT. x_end_local + dx/2.0_num )   xbd=1

          IF ((Cur%part_pos .LT. x_start - dx/2.0_num) .AND. (xbc_left_particle == BC_OPEN)) out_of_bounds = .TRUE.
          IF ((Cur%part_pos .GT. x_end + dx/2.0_num) .AND. (xbc_right_particle == BC_OPEN)) out_of_bounds = .TRUE.

          IF (ABS(xbd) .GT. 0) THEN
             !Particle has left box
             CALL Remove_Particle_From_PartList(ParticleSpecies(iSpecies)%AttachedList, Cur)
             IF (.NOT. out_of_bounds) THEN
                CALL Add_Particle_To_PartList(send(xbd),Cur)
             ELSE
                IF (DumpMask(29) .NE. IO_NEVER) THEN
                   CALL Add_Particle_To_PartList(ejected_particles,Cur)
                ELSE
                   DEALLOCATE(Cur)
                ENDIF
             ENDIF
          ENDIF

          !Move to next particle
          Cur=>Next
       ENDDO



       !swap Particles
       Done=.FALSE.
          DO ix=-1,1
             IF (ABS(ix) .EQ. 0) CYCLE
             ixp=-ix
             CALL PartList_SendRecv(send(ix),recv(ixp),Neighbour(ix),Neighbour(ixp))
             CALL Append_PartList(ParticleSpecies(ispecies)%AttachedList,recv(ixp))
       ENDDO


       !Particles should only lie outside boundaries if the periodic boundaries are turned on
       !This now moves them to within the boundaries
       Cur=>ParticleSpecies(iSpecies)%AttachedList%Head
       ct=0
       DO WHILE(ASSOCIATED(Cur))
          IF(Cur%part_pos .GT. x_end+dx/2.0_num .AND. xbc_left_particle == BC_PERIODIC) &
               Cur%part_pos=Cur%part_pos-length_x - dx
          IF(Cur%part_pos .LT. x_start-dx/2.0_num .AND. xbc_right_particle == BC_PERIODIC) &
               Cur%part_pos=Cur%part_pos+length_x + dx
          Cur=>Cur%Next
       ENDDO
    ENDDO
    !PRINT *,"Particle bcs_done",rank

  END SUBROUTINE Particle_bcs


END MODULE boundary
