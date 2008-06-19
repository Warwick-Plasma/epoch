MODULE boundary

  USE shared_data
  USE partlist

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Setup_Particle_Boundaries

    !For some types of boundary, fields and particles are treated in different ways, deal with that here
    xbc_right_particle=xbc_right
    xbc_left_particle=xbc_left

    IF (xbc_left == BC_SIMPLE_LASER) THEN
       xbc_left_particle=BC_OTHER
    ENDIF

    IF (xbc_right == BC_SIMPLE_LASER) THEN
       xbc_right_particle=BC_OTHER
    ENDIF

  END SUBROUTINE Setup_Particle_Boundaries

  SUBROUTINE Field_Reduction(Val)

    REAL(num),DIMENSION(-2:),INTENT(INOUT) :: Val
    REAL(num),DIMENSION(:),ALLOCATABLE ::Temp

    ALLOCATE(Temp(-2:nx_global+3))
    CALL MPI_ALLREDUCE(Val,Temp,nx_global+6,mpireal,MPI_SUM,comm,errcode)
    Val=Temp
    DEALLOCATE(Temp)

  END SUBROUTINE Field_Reduction

  !Exchanges field values at processor boundaries and applies field boundary conditions
  SUBROUTINE Field_BC(Field)

    REAL(num),DIMENSION(-2:),INTENT(INOUT) :: Field

    CALL MPI_SENDRECV(Field(1:2),2,mpireal,left,tag,Field(nx+1:nx+2),2,mpireal,right,tag,comm,status,errcode)
    CALL MPI_SENDRECV(Field(nx-1:nx),2,mpireal,right,tag,Field(-1:0),2,mpireal,left,tag,comm,status,errcode)

    IF (xbc_left == BC_OTHER .AND. left == MPI_PROC_NULL) THEN
       Field(-1:0)=0.0_num
       Field(-1:1)=0.0_num 
       Field(-1:0)=0.0_num
    ENDIF

    IF (xbc_right == BC_OTHER .AND. right == MPI_PROC_NULL) THEN
       Field(nx+1:nx+2)=0.0_num
       Field(nx+1:nx+2)=0.0_num
       Field(nx+1:nx+2)=0.0_num
    ENDIF

  END SUBROUTINE Field_BC

  !Deals with quantities such as current which are summed at processor boundaries
  SUBROUTINE Processor_Summation_bcs(Val)
    REAL(num),DIMENSION(-2:), INTENT(INOUT) :: Val
    REAL(num),DIMENSION(3) :: Data

    Data=0.0_num
    CALL MPI_SENDRECV(Val(nx+1:nx+3),3,mpireal,right,tag,Data,3,mpireal,left,tag,comm,status,errcode)
    Val(1:3)=Val(1:3) + Data

    Data=0.0_num
    CALL MPI_SENDRECV(Val(-2:0),3,mpireal,left,tag,Data,3,mpireal,right,tag,comm,status,errcode)
    Val(nx-2:nx)=Val(nx-2:nx) + Data

  END SUBROUTINE Processor_Summation_bcs


  !Summation on periodic boundaries on shared domains
  !This routine DOES NOT WORK on decomposed domains
  SUBROUTINE Periodic_Summation_bcs(Val)

    REAL(num),DIMENSION(-2:), INTENT(INOUT) :: Val

    IF (xbc_left .EQ. BC_PERIODIC .OR. xbc_right == BC_PERIODIC) THEN
       Val(1:3)=Val(1:3)+Val(nx+1:nx+3)
       Val(nx-2:nx)=Val(nx-2:nx)+Val(-2:0)
    ENDIF

  END SUBROUTINE Periodic_Summation_Bcs

  !E and B fields have different boundary conditions to allow for drivers
  SUBROUTINE Efield_bcs

    CALL Field_BC(Ex)
    CALL Field_BC(Ey)
    CALL Field_BC(Ez)

    IF (xbc_left == BC_OTHER .AND. left == MPI_PROC_NULL) THEN
       Ex(-1:0)=0.0_num
       Ey(-1:1)=0.0_num 
       Ez(-1:0)=0.0_num
    ENDIF

    IF (xbc_left==BC_SIMPLE_LASER .AND. left == MPI_PROC_NULL) THEN
       Ex(-1:0)=0.0_num
       Ey(-1:1)=Laser_Amp * SIN(Laser_Omega * time)
       Ez(-1:0)=0.0_num
    ENDIF

    IF (xbc_right == BC_OTHER .AND. right == MPI_PROC_NULL) THEN
       Ex(nx+1:nx+2)=0.0_num
       Ey(nx+1:nx+2)=0.0_num
       Ez(nx+1:nx+2)=0.0_num
    ENDIF

    IF (xbc_right == BC_SIMPLE_LASER .AND. right == MPI_PROC_NULL) THEN
       Ex(nx+1:nx+2)=0.0_num
       Ey(nx+1:nx+2)=Laser_Amp * SIN(Laser_Omega * time)
       Ez(nx+1:nx+2)=0.0_num
    ENDIF

  END SUBROUTINE Efield_bcs

  SUBROUTINE Bfield_bcs

    CALL Field_BC(Bx)
    CALL Field_BC(By)
    CALL Field_BC(Bz)

    IF ((xbc_left == BC_OTHER .OR. xbc_left == BC_SIMPLE_LASER) .AND. left == MPI_PROC_NULL) THEN
       Bx(-1:0)=0.0_num
       By(-1:0)=0.0_num
       Bz(-1:0)=0.0_num
    ENDIF

    IF ((xbc_right == BC_OTHER .OR. xbc_right == BC_SIMPLE_LASER) .AND. right == MPI_PROC_NULL) THEN
       Bx(nx+1:nx+2)=0.0_num
       By(nx+1:nx+2)=0.0_num
       Bz(nx+1:nx+2)=0.0_num
    ENDIF

  END SUBROUTINE Bfield_bcs

  SUBROUTINE Particle_bcs

    TYPE(Particle), POINTER :: Current, Next
    TYPE(ParticlePointer) :: Left_Send_List,Right_Send_List
    TYPE(ParticlePointer) :: Left_Recv_List, Right_Recv_List

    CALL Create_Empty_Partlist(Left_Send_List)
    CALL Create_Empty_Partlist(Left_Recv_List)
    CALL Create_Empty_Partlist(Right_Send_List)
    CALL Create_Empty_Partlist(Right_Recv_List)

    Current=>MainRoot%Head
    DO WHILE(ASSOCIATED(Current))
       Next=>Current%Next
       IF (Current%Part_Pos .GT. x_end_local) THEN
          !PRINT *,"Moving",rank,Current%Part_Pos
          IF (right .EQ. MPI_PROC_NULL) THEN
             IF (xbc_right_particle .EQ. BC_OTHER) THEN
                !Reached right boundary with reflecting boundaries
                Current%part_pos =  x_end - (Current%part_pos-x_end)
                Current%part_p(1) = - Current%part_p(1)
             ELSE IF (xbc_right_particle .EQ. BC_PERIODIC) THEN
                !Reached right boundary in periodic domain without domain decomposition
                Current%Part_Pos=Current%Part_Pos-length_x
             ENDIF
          ELSE
             !Flip if gone over periodic boundary
             IF (Current%Part_Pos .GT. x_end) THEN
                Current%Part_Pos=Current%Part_Pos-length_x
             ENDIF
             !Reached right processor boundary, put on list to go right
             CALL Remove_Particle_From_PartList(MainRoot, Current)
             CALL Add_Particle_To_PartList(Right_Send_List,Current)
          ENDIF
       ELSE IF (Current%Part_Pos .LT. x_start_local) THEN
          IF (left .EQ. MPI_PROC_NULL) THEN
             IF (xbc_left_particle .EQ. BC_OTHER) THEN
                !Reached left boundary with reflecting boundaries
                Current%part_pos  = x_start - (Current%part_pos-x_start)
                Current%part_p(1) = - Current%part_p(1)
             ELSE IF (xbc_left_particle .EQ. BC_PERIODIC) THEN
                !Reached left boundary in periodic domain without domain decomposition
                Current%Part_Pos=Current%Part_Pos+length_x
             ENDIF
          ELSE
             !Flip if gone over periodic boundary
             IF (Current%Part_Pos .LT. x_start) THEN
                Current%Part_Pos=Current%Part_Pos+length_x
             ENDIF
             !Reached left processor boundary, put on list to go left
             CALL Remove_Particle_From_PartList(MainRoot, Current)
             CALL Add_Particle_To_PartList(Left_Send_List,Current)
          ENDIF
       ELSE

       ENDIF
       Current=>Next
    ENDDO

    !Now have all the particles, so do the actual swapping
    !The particles to be sent are destroyed by this operation
    CALL PartList_SendRecv(Left_Send_List,Right_Recv_List,left,right)
    CALL PartList_SendRecv(Right_Send_List,Left_Recv_List,right,left)
    
    !Append the newly received particles to the main list of particles
    CALL Append_PartList(MainRoot,Right_Recv_List)
    CALL Append_PartList(MainRoot,Left_Recv_List)

    !Calculate new npart
    npart=MainRoot%Count

  END SUBROUTINE Particle_bcs

END MODULE boundary
