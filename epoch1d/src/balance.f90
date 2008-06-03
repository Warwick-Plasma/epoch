MODULE balance

  USE shared_data
  USE partlist

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Balance_Workload(OverRide)

    !In theory, this is a complete, working load balancing system for PIC1D
    !But it is totally untested


    LOGICAL,INTENT(IN) :: OverRide
    REAL(num) :: balance_frac
    INTEGER(8),DIMENSION(:),ALLOCATABLE :: npart_per_cell,npart_per_cell_global
    INTEGER(8) :: npart_per_rank_ideal,count
    INTEGER,DIMENSION(:),ALLOCATABLE :: displs,nx_each_rank_new,startpoint,endpoint
    INTEGER :: cell_x
    REAL(num),DIMENSION(:),ALLOCATABLE :: Temp_field
    TYPE(PARTICLE),POINTER :: Cur,NewHead,NewCur,Next,Temp,Tail
    INTEGER, DIMENSION(3) :: length,disp,type
    INTEGER :: subtype_int_field,iproc,partition,isrc,idest
    REAL(num) :: start_x,end_x
    REAL(num), DIMENSION(:), ALLOCATABLE :: DataIn,DataOut
    INTEGER,PARAMETER :: nvar=5

    !If running on one processor, you can do nothing to balance the load, so just return
    IF (nproc .EQ. 0) RETURN

    !If not overriding then check npart on each processor
    !If the ratio of the largest number of particles to the smallest number is less than
    !dlb_threshold then just leave
    IF (.NOT. OverRide) THEN
       !Get npart for each rank
       CALL MPI_ALLGATHER(npart,1,MPI_INTEGER8,npart_each_rank,1,MPI_INTEGER8,comm,errcode)
       balance_frac=REAL(MINVAL(npart_each_rank),num)/REAL(MAXVAL(npart_each_rank),num)
       IF (balance_frac .GT. dlb_threshold) RETURN
    ENDIF

    IF (rank .EQ. 0) THEN
       IF (.NOT. OverRide) THEN
          PRINT *,"Code rebalancing with balance fraction of",balance_frac
       ENDIF
    ENDIF

    CALL MPI_ALLGATHER(nx,1,MPI_INTEGER,nx_each_rank,1,MPI_INTEGER,comm,errcode)

    !First count number of particles per cell
    ALLOCATE(npart_per_cell(1:nx))
    npart_per_cell=0
    Cur=>MainRoot%Head
    DO WHILE(ASSOCIATED(Cur))
       cell_x = NINT((Cur%Part_Pos-x_start_local)/dx)+1
       IF (cell_x .GE. 1 .AND. cell_x .LE. nx) &
            npart_per_cell(cell_x)=npart_per_cell(cell_x)+1
       Cur=>Cur%Next
    ENDDO

    ALLOCATE(npart_per_cell_global(1:nx_global))
    npart_per_cell_global = 0

    !gather all the nparts arrays
    ALLOCATE(displs(1:nproc))
    displs=0
    DO iproc=2,nproc
       displs(iproc)=displs(iproc-1) + nx_each_rank(iproc-1)
    ENDDO
    CALL MPI_ALLGATHERV(npart_per_cell,nx,MPI_INTEGER8,npart_per_cell_global,nx_each_rank,displs,MPI_INTEGER8,comm,errcode)

    DEALLOCATE(npart_per_cell)


    !Now calculate the balanced number of particles per processor
    npart_per_rank_ideal=npart_global/nproc


    !This just counts through the number of particles on a cell by cell basis and 
    !Chops the current processor off when it has enough particles
    ALLOCATE(nx_each_rank_new(-1:nproc-1),startpoint(0:nproc),endpoint(-1:nproc-1))
    partition=0
    startpoint=1
    endpoint=0
    count=0
    nx_each_rank_new(-1)=0
    DO ix=1,nx_global
       count=count+npart_per_cell_global(ix)
       IF (count .GE. npart_per_rank_ideal) THEN
          count=0
          nx_each_rank_new(partition)=ix-SUM(nx_each_rank_new(0:partition-1))
          partition=partition+1
          IF (partition .EQ. nproc-1) EXIT
       ENDIF
    ENDDO
    nx_each_rank_new(nproc-1)=nx_global - SUM(nx_each_rank_new(0:nproc-2))

!!$    nx_each_rank_new=nx

    startpoint(0)=1

    DO iproc=0,nproc-1
       endpoint(iproc)=startpoint(iproc) + nx_each_rank_new(iproc) - 1
       startpoint(iproc+1)=endpoint(iproc) + 1
    ENDDO

    !No longer need npart_per_cell
    DEALLOCATE(npart_per_cell_global)
    ALLOCATE(Temp_Field(1:nx_global))


    !This isn't the most elegant way of doing this, but in 1D it suffices.
    !Just gather the global field variables and take your own local part out

    !Rebalance Ex
    CALL MPI_ALLGATHERV(Ex(1:nx),nx,mpireal,Temp_Field,nx_each_rank,displs,mpireal,comm,errcode)
    DEALLOCATE(Ex)
    ALLOCATE(Ex(-2:nx_each_rank_new(rank)+2))
    Ex(1:nx_each_rank_new(rank)) = Temp_Field(startpoint(rank):endpoint(rank))

    !Rebalance Ey
    CALL MPI_ALLGATHERV(Ey(1:nx),nx,mpireal,Temp_Field,nx_each_rank,displs,mpireal,comm,errcode)
    DEALLOCATE(Ey)
    ALLOCATE(Ey(-2:nx_each_rank_new(rank)+2))
    Ey(1:nx_each_rank_new(rank))=Temp_Field(startpoint(rank):endpoint(rank))

    !Rebalance Ez
    CALL MPI_ALLGATHERV(Ez(1:nx),nx,mpireal,Temp_Field,nx_each_rank,displs,mpireal,comm,errcode)
    DEALLOCATE(Ez)
    ALLOCATE(Ez(-2:nx_each_rank_new(rank)+2))
    Ez(1:nx_each_rank_new(rank))=Temp_Field(startpoint(rank):endpoint(rank))

    !Rebalance Bx
    CALL MPI_ALLGATHERV(Bx(1:nx),nx,mpireal,Temp_Field,nx_each_rank,displs,mpireal,comm,errcode)
    DEALLOCATE(Bx)
    ALLOCATE(Bx(-2:nx_each_rank_new(rank)+2))
    Bx(1:nx_each_rank_new(rank))=Temp_Field(startpoint(rank):endpoint(rank))

    !Rebalance By
    CALL MPI_ALLGATHERV(By(1:nx),nx,mpireal,Temp_Field,nx_each_rank,displs,mpireal,comm,errcode)
    DEALLOCATE(By)
    ALLOCATE(By(-2:nx_each_rank_new(rank)+2))
    By(1:nx_each_rank_new(rank))=Temp_Field(startpoint(rank):endpoint(rank))

    !Rebalance Bz
    CALL MPI_ALLGATHERV(Bz(1:nx),nx,mpireal,Temp_Field,nx_each_rank,displs,mpireal,comm,errcode)
    DEALLOCATE(Bz)
    ALLOCATE(Bz(-2:nx_each_rank_new(rank)+2))
    Bz(1:nx_each_rank_new(rank))=Temp_Field(startpoint(rank):endpoint(rank))

    !Rebalance Bz
    CALL MPI_ALLGATHERV(Bz(1:nx),nx,mpireal,Temp_Field,nx_each_rank,displs,mpireal,comm,errcode)
    DEALLOCATE(Bz)
    ALLOCATE(Bz(-2:nx_each_rank_new(rank)+2))
    Bz(1:nx_each_rank_new(rank))=Temp_Field(startpoint(rank):endpoint(rank))

    !Rebalance X
    DEALLOCATE(x)
    ALLOCATE(x(-1:nx_each_rank_new(rank)+2))
    x(1:nx_each_rank_new(rank))=x_global(startpoint(rank):endpoint(rank))
    x(0)=x(1)-dx
    x(nx_each_rank_new(rank)+1)=x(nx_each_rank_new(rank))+dx
    !Recalculate the start and end points of each processor
    DO iproc=0,nproc-1
       x_starts(iproc)=x_global(startpoint(iproc))
       x_ends(iproc)=x_global(endpoint(iproc)+1)
    ENDDO


    DEALLOCATE(jx,jy,jz)
    ALLOCATE(jx(-2:nx_each_rank_new(rank)+3),jy(-2:nx_each_rank_new(rank)+3),jz(-2:nx_each_rank_new(rank)+3))
    jx=0.0_num
    jy=0.0_num
    jz=0.0_num

    DEALLOCATE(ek_bar,ekbar_sum,ct)
    ALLOCATE(ek_bar(-2:nx_each_rank_new(rank)+3,1:nspecies),ekbar_sum(-2:nx_each_rank_new(rank)+3,1:nspecies)&
         ,ct(-2:nx_each_rank_new(rank)+3,1:nspecies))

    x_start_local=x_starts(rank)
    x_end_local=x_ends(rank)

    !Now that the system knows about the new domain lengths, you can just distribute the particles using the normal routine
    CALL Distribute_Particles

    CALL MPI_BARRIER(comm,errcode)

    DEALLOCATE(displs)
    nx=nx_each_rank_new(rank)

  END SUBROUTINE Balance_Workload

  FUNCTION GetParticleProcessor(aParticle)

    TYPE(Particle),INTENT(IN) :: aParticle
    INTEGER :: GetParticleProcessor
    INTEGER :: CurrentLoc, iproc
    GetParticleProcessor = -1

    !This could be replaced by a bisection method, but for the moment I just don't care

    DO iproc=0,nproc-1
       IF (aParticle%Part_Pos .GE. x_starts(iproc) .AND. aParticle%Part_Pos .LE. x_ends(iproc)) THEN
          GetParticleProcessor=iproc
          RETURN
       ENDIF
    ENDDO


  END FUNCTION GetParticleProcessor


  !This subroutine is used to rearrange particles over processors
  SUBROUTINE Distribute_Particles

    TYPE(ParticlePointer),DIMENSION(:),ALLOCATABLE :: Pointers_Send, Pointers_Recv
    TYPE(Particle), POINTER :: Current, Next
    INTEGER :: part_proc, iproc_recv,iproc_send

    ALLOCATE(Pointers_Send(0:nprocx-1), Pointers_Recv(0:nprocx-1))    
    Current=>MainRoot%Head

    DO iproc_send = 0, nprocx-1
       CALL Create_Empty_PartList(Pointers_Send(iproc_send))
       CALL Create_Empty_PartList(Pointers_Recv(iproc_send))
    ENDDO


    DO WHILE(ASSOCIATED(Current))
       Next=>Current%Next
       part_proc=GetParticleProcessor(Current)
       IF (part_proc .LT. 0) THEN
          PRINT *,"Unlocatable particle on processor",rank,Current%Part_Pos
          CALL MPI_BARRIER(comm,errcode)
          STOP
          !PRINT *,x_starts,x_ends
       ENDIF
       IF (part_proc .NE. rank) THEN
          CALL Remove_Particle_From_PartList(MainRoot, Current)
          CALL Add_Particle_To_PartList(Pointers_Send(part_proc),Current)
       ENDIF
       Current=>Next
    ENDDO

    DO iproc_send = 0, nprocx-1
       DO iproc_recv=0, nprocx-1
          IF (iproc_send .NE. iproc_recv) THEN
             IF (rank .EQ. iproc_send) THEN 
                CALL PartList_Send(Pointers_Send(iproc_recv),iproc_recv)
                CALL Destroy_PartList(Pointers_Send(iproc_recv))
             ENDIF
             IF (rank .EQ. iproc_recv) CALL PartList_Recv(Pointers_Recv(iproc_send),iproc_send)
          ENDIF
       ENDDO
    ENDDO


    DO iproc_recv = 0, nprocx-1
       CALL Append_PartList(MainRoot,Pointers_Recv(iproc_recv))
    ENDDO

    npart=MainRoot%Count
    DEALLOCATE(Pointers_Send,Pointers_Recv)


  END SUBROUTINE Distribute_Particles

  SUBROUTINE CreateSubtypes

    INTEGER, DIMENSION(3) :: length,disp,type

    ! Create the subarray for the fields in this problem: subtype decribes where this
    ! process's data fits into the global picture.
    nx_each_rank=0
    CALL MPI_ALLGATHER(nx,1,MPI_INTEGER,nx_each_rank,1,MPI_INTEGER,comm,errcode)
    length(1)=1
    length(2)=nx
    length(3)=1
    disp(1)=0
    disp(2)=0
    DO ix=1,rank
       disp(2)=disp(2)+nx_each_rank(ix)*num
    ENDDO
    disp(3)= nx_global * num
    type(1)=MPI_LB
    type(2)=mpireal
    type(3)=MPI_UB
    IF (subtype_field /= 0) CALL MPI_TYPE_FREE(subtype_field,errcode)
    CALL MPI_TYPE_STRUCT(3,length,disp,type,subtype_field,errcode)
    CALL MPI_TYPE_COMMIT(subtype_field,errcode)

    ! Create the subarray for the particles in this problem: subtype decribes where this
    ! process's data fits into the global picture.
    CALL MPI_ALLGATHER(npart,1,MPI_INTEGER8,npart_each_rank,1,MPI_INTEGER8,comm,errcode)
    length(1)=1
    length(2)=npart * 1 !1D
    length(3)=1
    disp(1)=0
    disp(2)=0
    DO ix=1,rank
       disp(2)=disp(2)+npart_each_rank(ix)*num * 1 !1D
    ENDDO
    disp(3)= npart_global * num * 1 !1D
    type(1)=MPI_LB
    type(2)=mpireal
    type(3)=MPI_UB
    IF (subtype_particle /= 0) CALL MPI_TYPE_FREE(subtype_particle,errcode)
    CALL MPI_TYPE_STRUCT(3,length,disp,type,subtype_particle,errcode)
    CALL MPI_TYPE_COMMIT(subtype_particle,errcode)


    ! Create the subarray for the particles in this problem: subtype decribes where this
    ! process's data fits into the global picture.
    length(1)=1
    length(2)=npart
    length(3)=1
    disp(1)=0
    disp(2)=0
    DO ix=1,rank
       disp(2)=disp(2)+npart_each_rank(ix)*num
    ENDDO
    disp(3)= npart_global * num
    type(1)=MPI_LB
    type(2)=mpireal
    type(3)=MPI_UB
    IF (subtype_particle_var /= 0) CALL MPI_TYPE_FREE(subtype_particle_var,errcode)
    CALL MPI_TYPE_STRUCT(3,length,disp,type,subtype_particle_var,errcode)
    CALL MPI_TYPE_COMMIT(subtype_particle_var,errcode)


    ! Create the subarray for the integer properties of particles in this problem: subtype decribes where this
    ! process's data fits into the global picture.
    length(1)=1
    length(2)=npart
    length(3)=1
    disp(1)=0
    disp(2)=0
    DO ix=1,rank
       disp(2)=disp(2)+npart_each_rank(ix)*4
    ENDDO
    disp(3)= npart_global * num
    type(1)=MPI_LB
    type(2)=MPI_INTEGER
    type(3)=MPI_UB
    IF (subtype_particle_int /= 0) CALL MPI_TYPE_FREE(subtype_particle_int,errcode)
    CALL MPI_TYPE_STRUCT(3,length,disp,type,subtype_particle_int,errcode)
    CALL MPI_TYPE_COMMIT(subtype_particle_int,errcode)

  END SUBROUTINE CreateSubtypes


END MODULE balance
