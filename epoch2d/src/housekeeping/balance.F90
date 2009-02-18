MODULE balance

  USE shared_data
  USE partlist
  USE boundary

  IMPLICIT NONE

CONTAINS

  FUNCTION Get_Total_Local_Particles()

    !This subroutine describes the total number of particles on the current processor
    !It simply sums over every particle species

    INTEGER(KIND=8) :: Get_Total_Local_Particles
    INTEGER :: iSpecies
    Get_Total_Local_Particles=0
    DO iSpecies = 1,nspecies
       Get_Total_Local_Particles=Get_Total_Local_Particles+ParticleSpecies(iSpecies)%AttachedList%Count
    ENDDO

  END FUNCTION Get_Total_Local_Particles

  FUNCTION Get_Total_Local_Dumped_Particles(Force_Restart)

    !This subroutine describes the total number of particles on the current processor
    !which are members of species with the Dump=T attribute in the input deck
    !If FORCE_RESTART=.TRUE. then the subroutine simply counts all the particles

    LOGICAL,INTENT(IN) :: Force_Restart
    INTEGER(KIND=8) :: Get_Total_Local_Dumped_Particles
    INTEGER :: iSpecies

    Get_Total_Local_Dumped_Particles=0
    DO iSpecies = 1,nspecies
       IF (ParticleSpecies(iSpecies)%Dump .OR. Force_Restart) THEN
          Get_Total_Local_Dumped_Particles=Get_Total_Local_Dumped_Particles+ParticleSpecies(iSpecies)%AttachedList%Count
       ENDIF
    ENDDO

  END FUNCTION Get_Total_Local_Dumped_Particles

  SUBROUTINE Balance_Workload(OverRide)

    !This subroutine determines whether or not the code needs rebalancing, calculates where to split the domain
    !and calls other subroutines to actually rearrange the fields and particles onto the new processors

    !This is really, really hard to do properly
    !So cheat

    LOGICAL,INTENT(IN) :: OverRide
    INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE :: npart_each_rank
    INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE :: Density_x,Density_y
    INTEGER,DIMENSION(:),ALLOCATABLE,TARGET :: Starts_x,ends_x,Starts_y,ends_y
    INTEGER :: new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end
    REAL(num) :: balance_frac,balance_frac_x,balance_frac_y
    INTEGER(KIND=8) :: Max_x,Max_y,wk,Min_x,Min_y,npart_local
    INTEGER :: iProc,iSpecies,iLaser
#ifdef PART_DEBUG
    TYPE(Particle),POINTER :: Current
#endif

    !On one processor do nothing to save time
    IF (nproc .EQ. 1) RETURN

    !This parameter allows selecting the mode of the autobalancing
    !Between leftsweep, rightsweep, auto(best of leftsweep and rightsweep) or both
    Balance_Mode=LB_BOTH

    !Count particles
    npart_local=Get_Total_Local_Particles()

    !The OverRide flag allows the code to force a load balancing sweep at t=0
    IF (.NOT. OverRide) THEN
       ALLOCATE(npart_each_rank(1:nproc))
       !Get npart for each rank
       CALL MPI_ALLGATHER(npart_local,1,MPI_INTEGER8,npart_each_rank,1,MPI_INTEGER8,comm,errcode)
       !Determine ratio of npart on between most loaded and least loaded processor
       !Maybe this can be replaced by and MPI_ALLREDUCE to find min/max?
       balance_frac=REAL(MINVAL(npart_each_rank),num)/REAL(MAXVAL(npart_each_rank),num)
       IF (balance_frac .GT. dlb_threshold) RETURN
       IF (rank .EQ. 0) PRINT *,"Load balancing with fraction",balance_frac
       DEALLOCATE(npart_each_rank)
    ENDIF

    ALLOCATE(starts_x(1:nprocx),ends_x(1:nprocx))
    ALLOCATE(starts_y(1:nprocy),ends_y(1:nprocy))

    !Sweep in X
    IF (IAND(Balance_Mode,LB_X) .NE. 0 .OR. IAND(Balance_Mode,LB_AUTO) .NE. 0) THEN
       !Rebalancing in X
       ALLOCATE(Density_x(0:nx_global+1))
       CALL GetDensityInX(Density_x)
       CALL CalculateBreaks(Density_x,nprocx,starts_x,ends_x)
    ELSE
       !Just keep the original lengths
       starts_x=cell_x_start
       ends_x=cell_x_end
    ENDIF

    !Sweep in Y
    IF (IAND(Balance_Mode,LB_Y) .NE. 0 .OR. IAND(Balance_Mode,LB_AUTO) .NE. 0) THEN
       !Rebalancing in Y
       ALLOCATE(Density_y(0:ny_global+1))
       CALL GetDensityInY(Density_y)
       CALL CalculateBreaks(Density_y,nprocy,starts_y,ends_y)
    ELSE
       !Just keep the original lengths
       starts_y=cell_y_start
       ends_y=cell_y_end
    ENDIF


    !In the autobalancer then determine whether to balance in X or Y
    !Is this worth keeping?
    IF (IAND(Balance_Mode,LB_AUTO) .NE. 0 ) THEN

       !Code is auto load balancing
       Max_x=0
       Min_x=npart_global
       DO iproc=1,nprocx
          wk=SUM(Density_x(starts_x(iproc):ends_x(iproc)))
          IF (wk .GT. Max_x) Max_x=wk
          IF (wk .LT. Min_x) Min_x=wk
       ENDDO
       Max_y=0
       Min_y=npart_global
       DO iproc=1,nprocy
          wk=SUM(Density_y(starts_y(iproc):ends_y(iproc)))
          IF (wk .GT. Max_y) Max_y=wk
          IF (wk .LT. Min_y) Min_y=wk
       ENDDO

       balance_frac_x=REAL(Min_x,num)/REAL(Max_x,num)
       balance_frac_y=REAL(Min_y,num)/REAL(Max_y,num)

       IF (balance_frac_y .LT. balance_frac_x) THEN
          starts_y=cell_y_start
          ends_y=cell_y_end
       ELSE
          starts_x=cell_x_start
          ends_x=cell_x_end
       ENDIF

    ENDIF

    !Now need to calculate the start and end points for the new domain on the current processor
    new_cell_x_start=starts_x(coordinates(2)+1)
    new_cell_x_end=ends_x(coordinates(2)+1)

    new_cell_y_start=starts_y(coordinates(1)+1)
    new_cell_y_end=ends_y(coordinates(1)+1)

    !Redeistribute the field variables
    CALL Redistribute_Fields(new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end)

    !Copy the new lengths into the permanent variables
    cell_x_start=starts_x
    cell_y_start=starts_y
    cell_x_end=ends_x
    cell_y_end=ends_y

    !Set the new nx and ny
    nx=new_cell_x_end-new_cell_x_start+1
    ny=new_cell_y_end-new_cell_y_start+1

    !Do X and Y arrays separatly because we already have global copies of X and Y
    DEALLOCATE(x,y)
    ALLOCATE(x(-2:nx+3),y(-2:ny+3))
    x(0:nx+1)=x_global(new_cell_x_start-1:new_cell_x_end+1)
    y(0:ny+1)=y_global(new_cell_y_start-1:new_cell_y_end+1)

    !Reallocate Currents (don't need to keep current timestep values for current so just reallocate)
    DEALLOCATE(Jx,Jy,Jz)
    ALLOCATE(Jx(-2:nx+3,-2:ny+3),Jy(-2:nx+3,-2:ny+3),Jz(-2:nx+3,-2:ny+3))
    Jx=0.0_num
    Jy=0.0_num
    Jz=0.0_num

    !Reallocate the kinetic energy calculation
    DEALLOCATE(ekbar,ekbar_sum,ct)
    ALLOCATE(ekbar(1:nx,1:ny,1:nspecies),ekbar_sum(-2:nx+3,-2:ny+3,1:nspecies))
    ALLOCATE(ct(-2:nx+3,-2:ny+3,1:nspecies))

    !Recalculate x_starts and y_starts so that rebalancing works next time
    DO iproc=0,nprocx-1
       x_starts(iproc)=x_global(cell_x_start(iproc+1)) 
       x_ends(iproc)=x_global(cell_x_end(iproc+1))
    ENDDO
    !Same for y
    DO iproc=0,nprocy-1
       y_starts(iproc)=y_global(cell_y_start(iproc+1))
       y_ends(iproc)=y_global(cell_y_end(iproc+1))
    ENDDO

    !Set the lengths of the current domain so that the particle balancer works properly
    x_start_local=x_starts(coordinates(2))
    x_end_local=x_ends(coordinates(2))
    y_start_local=y_starts(coordinates(1))
    y_end_local=y_ends(coordinates(1))

    !Redistribute the particles onto their new processors
    CALL Distribute_Particles


    !If running with particle debugging then set the t=0 processor if
    !Override =true
#ifdef PART_DEBUG
    IF (OverRide) THEN
       DO iSpecies=1,nSpecies
          Current=>ParticleSpecies(iSpecies)%AttachedList%Head
          DO WHILE(ASSOCIATED(Current))
             Current%Processor_At_T0=rank
             Current=>Current%Next
          ENDDO
       ENDDO
    ENDIF
#endif

  END SUBROUTINE Balance_Workload

  SUBROUTINE Redistribute_Fields(new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end)

    !This subroutine redistributes the 2D field variables over the new processor layout
    !If using a 2D field of your own then se the Redistribute_Field subroutine to implement it
    !1D fields, you're on your own (have global copies and use those to repopulate?)
    INTEGER, INTENT(IN) :: new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end
    INTEGER :: nx_new,ny_new,iLaser
    REAL(num),DIMENSION(:,:), ALLOCATABLE :: temp
    REAL(num),DIMENSION(:), ALLOCATABLE :: temp1d
    TYPE(Laser_Block),POINTER :: Current

    nx_new=new_cell_x_end-new_cell_x_start+1
    ny_new=new_cell_y_end-new_cell_y_start+1

    ALLOCATE(temp(-2:nx_new+3,-2:ny_new+3))

    temp=0.0_num
    CALL Redistribute_Field(new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end,Ex,temp)
    DEALLOCATE(Ex)
    ALLOCATE(Ex(-2:nx_new+3,-2:ny_new+3))
    Ex=temp

    temp=0.0_num
    CALL Redistribute_Field(new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end,Ey,temp)
    DEALLOCATE(Ey)
    ALLOCATE(Ey(-2:nx_new+3,-2:ny_new+3))
    Ey=temp

    temp=0.0_num
    CALL Redistribute_Field(new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end,Ez,temp)
    DEALLOCATE(Ez)
    ALLOCATE(Ez(-2:nx_new+3,-2:ny_new+3))
    Ez=temp

    temp=0.0_num
    CALL Redistribute_Field(new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end,Bx,temp)
    DEALLOCATE(Bx)
    ALLOCATE(Bx(-2:nx_new+3,-2:ny_new+3))
    Bx=temp

    temp=0.0_num
    CALL Redistribute_Field(new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end,By,temp)
    DEALLOCATE(By)
    ALLOCATE(By(-2:nx_new+3,-2:ny_new+3))
    By=temp

    temp=0.0_num
    CALL Redistribute_Field(new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end,Bz,temp)
    DEALLOCATE(Bz)
    ALLOCATE(Bz(-2:nx_new+3,-2:ny_new+3))
    Bz=temp

    DEALLOCATE(temp)

#ifndef FULL_LASER_SETTINGS
    ALLOCATE(temp1d(-2:ny_new+3))
    Current=>Laser_Left
    DO WHILE(ASSOCIATED(Current))
       temp1d=0.0_num
       CALL Redistribute_Field_1D(new_cell_y_start,new_cell_y_end,cell_y_start(coordinates(1)+1)&
            ,cell_y_end(coordinates(1)+1),ny_global,Current%Profile,temp1d,DIR_X)
       DEALLOCATE(Current%Profile)
       ALLOCATE(Current%Profile(-2:ny_new+3))
       Current%Profile=temp1d
       temp1d=0.0_num
       CALL Redistribute_Field_1D(new_cell_y_start,new_cell_y_end,cell_y_start(coordinates(1)+1)&
            ,cell_y_end(coordinates(1)+1),ny_global,Current%Phase,temp1d,DIR_X)
       DEALLOCATE(Current%Phase)
       ALLOCATE(Current%Phase(-2:ny_new+3))
       Current%Phase=temp1d

       Current=>Current%Next
    ENDDO
    Current=>Laser_Right
    DO WHILE(ASSOCIATED(Current))
       temp1d=0.0_num
       CALL Redistribute_Field_1D(new_cell_y_start,new_cell_y_end,cell_y_start(coordinates(1)+1)&
            ,cell_y_end(coordinates(1)+1),ny_global,Current%Profile,temp1d,DIR_X)
       DEALLOCATE(Current%Profile)
       ALLOCATE(Current%Profile(-2:ny_new+3))
       Current%Profile=temp1d
       temp1d=0.0_num
       CALL Redistribute_Field_1D(new_cell_y_start,new_cell_y_end,cell_y_start(coordinates(1)+1)&
            ,cell_y_end(coordinates(1)+1),ny_global,Current%Phase,temp1d,DIR_X)
       DEALLOCATE(Current%Profile)
       ALLOCATE(Current%Profile(-2:ny_new+3))
       Current%Profile=temp1d

       Current=>Current%Next
    ENDDO
#endif


  END SUBROUTINE Redistribute_Fields

  SUBROUTINE Redistribute_Field(new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end,Field_In,Field_Out)

    !This subroutine redistributes the fields over the new processor layout
    !The current version works by producing a global copy on each processor
    !And then extracting the required part for the local processor.
    !This is not in general a good idea
    INTEGER, INTENT(IN) :: new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end
    REAL(num),DIMENSION(-2:,-2:),INTENT(IN) :: Field_In
    REAL(num),DIMENSION(-2:,-2:),INTENT(OUT) :: Field_Out
    REAL(num),DIMENSION(:,:),ALLOCATABLE :: Field_New,Field_Temp
    INTEGER :: nx_new,ny_new
    INTEGER :: comm_new,iproc,color

    nx_new=new_cell_x_end-new_cell_x_start+1
    ny_new=new_cell_y_end-new_cell_y_start+1

!    PRINT *,rank,new_cell_x_start,new_cell_y_start," "

    !This is a horrible, horrible way of doing this, I MUST think of a better way

    !Create a global copy of the whole array
    ALLOCATE(Field_New(1:nx_global,1:ny_global),Field_Temp(1:nx_global,1:ny_global))
    Field_New=0.0_num
    Field_New(cell_x_start(coordinates(2)+1):cell_x_end(coordinates(2)+1),&
         cell_y_start(coordinates(1)+1):cell_y_end(coordinates(1)+1))=Field_In(1:nx,1:ny)
    CALL MPI_ALLREDUCE(Field_New,Field_Temp,nx_global*ny_global,mpireal,MPI_SUM,comm,errcode)

    Field_Out(1:nx_new,1:ny_new)=Field_Temp(new_cell_x_start:new_cell_x_end,new_cell_y_start:new_cell_y_end)
    DEALLOCATE(Field_Temp)

    !Call boundary conditions (this does not include any special BCS)
    CALL Do_Field_MPI_With_Lengths(Field_out,nx_new,ny_new)

  END SUBROUTINE Redistribute_Field

  SUBROUTINE Redistribute_Field_1D(new_start,new_end,old_start,old_end,npts_global,Field_In,Field_Out,direction)
    !This subroutine redistributes a 1D field over the new processor layout
    !The current version works by producing a global copy on each processor
    !And then extracting the required part for the local processor.
    !in 1D, this is probably OK
    INTEGER, INTENT(IN) :: new_start,new_end,npts_global,old_start,old_end,direction
    REAL(num),DIMENSION(-2:),INTENT(IN) :: Field_In
    REAL(num),DIMENSION(-2:),INTENT(OUT) :: Field_Out
    REAL(num),DIMENSION(:),ALLOCATABLE :: Field_New,Field_Temp
    INTEGER :: new_pts,old_pts
    INTEGER :: New_Comm, color


    IF (IAND(Direction,DIR_X) .EQ. 0) color=coordinates(1)
    IF (IAND(Direction,DIR_Y) .EQ. 0) color=coordinates(2)

    CALL MPI_COMM_SPLIT(comm,color,rank,new_comm,errcode)

    new_pts=new_end-new_start+1
    old_pts=old_end-old_start+1

    !Create a global copy of the whole array
    ALLOCATE(Field_New(1:npts_global),Field_Temp(0:npts_global+1))
    Field_New=0.0_num
    Field_New(old_start:old_end)=Field_In(1:old_pts)
    CALL MPI_ALLREDUCE(Field_New,Field_Temp(1:npts_global),npts_global,mpireal,MPI_SUM,New_comm,errcode)
    DEALLOCATE(Field_New)

    Field_Out(0:new_pts+1)=Field_Temp(new_start-1:new_end+1)
    DEALLOCATE(Field_Temp)
    CALL MPI_COMM_FREE(new_comm,errcode)

  END SUBROUTINE Redistribute_Field_1D

  SUBROUTINE GetDensityInX(Density)

    !Calculate total particle density across the X direction
    !Summed in the Y direction

    INTEGER(KIND=8),DIMENSION(:),INTENT(INOUT) :: Density
    INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE :: Temp
    TYPE(Particle),POINTER :: Current
    REAL(num) :: part_x
    INTEGER :: cell_x1,iSpecies
    Density=0.0_num

    DO iSpecies=1,nspecies
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO WHILE(ASSOCIATED(Current))
          !Want global position, so not x_start, NOT x_start_local
          part_x=Current%Part_Pos(1)-x_start
          cell_x1=NINT(part_x/dx)+1
          Density(cell_x1)=Density(cell_x1)+1
          Current=>Current%Next
       ENDDO
    ENDDO
    !Now have local densities, so add using MPI
    ALLOCATE(Temp(0:nx_global+1))
    CALL MPI_ALLREDUCE(Density,Temp,nx_global+2,MPI_INTEGER8,MPI_SUM,comm,errcode)
    Density=Temp
    DEALLOCATE(Temp)

  END SUBROUTINE GetDensityInX

  SUBROUTINE GetDensityInY(Density)

    !Calculate total particle density across the Y direction
    !Summed in the X direction

    INTEGER(KIND=8),DIMENSION(:),INTENT(INOUT) :: Density
    INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE :: Temp
    TYPE(Particle),POINTER :: Current
    REAL(num) :: part_y
    INTEGER :: cell_y1,iSpecies
    Density=0.0_num

    DO iSpecies=1,nspecies
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO WHILE(ASSOCIATED(Current))
          !Want global position, so not x_start, NOT x_start_local
          part_y=Current%Part_Pos(2)-y_start
          cell_y1=NINT(part_y/dy)+1
          Density(cell_y1)=Density(cell_y1)+1
          Current=>Current%Next
       ENDDO
    ENDDO
    !Now have local densities, so add using MPI
    ALLOCATE(Temp(0:ny_global+1))
    CALL MPI_ALLREDUCE(Density,Temp,ny_global+2,MPI_INTEGER8,MPI_SUM,comm,errcode)
    Density=Temp

    DEALLOCATE(TEMP)
  END SUBROUTINE GetDensityInY

  SUBROUTINE CalculateBreaks(Density,nproc,starts,ends)

    !This subroutine calculates the places in a given density profile to split
    !The domain to give the most even subdivision possible

    INTEGER(KIND=8),INTENT(IN),DIMENSION(:) :: Density
    INTEGER,INTENT(IN) :: nproc
    INTEGER,DIMENSION(:),INTENT(OUT) :: starts,ends
    INTEGER :: Sz, iDim, partition
    INTEGER(KIND=8) :: Total
    INTEGER :: npart_per_proc_ideal

    !-2 because of ghost cells at each end
    Sz=SIZE(Density)-2
    IF (nproc .EQ. 1) THEN
       starts=1
       ends=Sz
    ENDIF

    npart_per_proc_ideal=SUM(Density)/nproc
    partition=2
    starts(1)=1
    total=0
    DO iDim=1,Sz
       IF (partition .GT. nproc) EXIT
       IF (total .GE. npart_per_proc_ideal .OR. ABS(total + Density(iDim) -npart_per_proc_ideal) .GT. ABS(total-npart_per_proc_ideal)  .OR. iDim .EQ. Sz) THEN
          total=Density(iDim)
          starts(partition)=iDim-1
          partition=partition+1
          !If you've reached the last processor, have already done the best you can, so just leave
       ELSE
          total=total+Density(iDim)
       ENDIF
    ENDDO
    ends(nproc)=Sz
    DO iDim=1,nproc-1
       ends(iDim)=starts(iDim+1)-1
    ENDDO

  END SUBROUTINE CalculateBreaks

  FUNCTION GetParticleProcessor(aParticle)

    !This subroutine calculates which processor a given particles resides on

    TYPE(Particle),INTENT(IN) :: aParticle
    INTEGER :: GetParticleProcessor
    INTEGER :: CurrentLoc, iproc,coords(2)
    GetParticleProcessor = -1
    coords=-1

    !This could be replaced by a bisection method, but for the moment I just don't care


    DO iproc=0,nprocx-1
       IF (aParticle%Part_Pos(1) .GE. x_starts(iproc) - dx/2.0_num .AND. aParticle%Part_Pos(1) .LE. x_ends(iproc) + dx/2.0_num) THEN
          coords(2)=iproc
          EXIT
       ENDIF
    ENDDO

    DO iproc=0,nprocy-1
       IF (aParticle%Part_Pos(2) .GE. y_starts(iproc) -dy/2.0_num .AND. aParticle%Part_Pos(2) .LE. y_ends(iproc) + dy/2.0_num) THEN
          coords(1)=iproc
          EXIT
       ENDIF
    ENDDO


    IF (MINVAL(coords) .LT. 0) PRINT *,"UNLOCATABLE PARTICLE"
    IF (MINVAL(coords) .LT. 0) RETURN
    CALL MPI_CART_RANK(comm,coords,GetParticleProcessor,errcode)
    !    IF (GetParticleProcessor .NE. rank) PRINT *,

  END FUNCTION GetParticleProcessor


  !This subroutine is used to rearrange particles over processors
  SUBROUTINE Distribute_Particles

    !This subroutine actually moves particles which are on the wrong processor
    !And moves then to the correct processor.

    TYPE(ParticleList),DIMENSION(:),ALLOCATABLE :: Pointers_Send, Pointers_Recv
    TYPE(Particle), POINTER :: Current, Next
    INTEGER :: part_proc, iproc_recv,iproc_send,iSpecies


    ALLOCATE(Pointers_Send(0:nproc-1), Pointers_Recv(0:nproc-1))  
    DO iSpecies=1,nspecies
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO iproc_send = 0, nproc-1
          CALL Create_Empty_PartList(Pointers_Send(iproc_send))
          CALL Create_Empty_PartList(Pointers_Recv(iproc_send))
       ENDDO


       DO WHILE(ASSOCIATED(Current))
          Next=>Current%Next
          part_proc=GetParticleProcessor(Current)
          IF (part_proc .LT. 0) THEN
             PRINT *,"Unlocatable particle on processor",rank,Current%Part_Pos,dx,dy
             CALL MPI_BARRIER(comm,errcode)
             STOP
          ENDIF
#ifdef PART_DEBUG
          Current%Processor=part_proc
#endif
          IF (part_proc .NE. rank) THEN
!!$          PRINT *,Current%Part_Pos,part_proc
             CALL Remove_Particle_From_PartList(ParticleSpecies(iSpecies)%AttachedList, Current)
             CALL Add_Particle_To_PartList(Pointers_Send(part_proc),Current)
          ENDIF
          Current=>Next
       ENDDO

       DO iproc_send = 0, nproc-1
          DO iproc_recv=0, nproc-1
             IF (iproc_send .NE. iproc_recv) THEN
                IF (rank .EQ. iproc_send) THEN 
                   CALL PartList_Send(Pointers_Send(iproc_recv),iproc_recv)
                   CALL Destroy_PartList(Pointers_Send(iproc_recv))
                ENDIF
                IF (rank .EQ. iproc_recv) CALL PartList_Recv(Pointers_Recv(iproc_send),iproc_send)
             ENDIF
          ENDDO
       ENDDO


       DO iproc_recv = 0, nproc-1
          CALL Append_PartList(ParticleSpecies(iSpecies)%AttachedList,Pointers_Recv(iproc_recv))
       ENDDO
    ENDDO

    DEALLOCATE(Pointers_Send,Pointers_Recv)

  END SUBROUTINE Distribute_Particles

  SUBROUTINE CreateSubtypes(Force_Restart)

    !This subroutines creates the MPI types which represent the data for the field and
    !particles data. It is used when writing data
    LOGICAL,INTENT(IN) :: Force_Restart
    INTEGER(KIND=8), DIMENSION(3) :: length,disp,type
    INTEGER :: ndims=2
    INTEGER,DIMENSION(:),ALLOCATABLE :: lengths, starts

    INTEGER(KIND=8) :: npart_local

    npart_local = Get_Total_Local_Dumped_Particles(Force_Restart)

    ! Create the subarray for the fields in this problem: subtype decribes where this
    ! process's data fits into the global picture.

    ALLOCATE(lengths(1:ny),starts(1:ny))

    lengths=nx
    DO iy=0,ny-1
       starts(iy+1)=(cell_y_start(coordinates(1)+1)+iy-1) * nx_global + cell_x_start(coordinates(2)+1) -1
    ENDDO

    CALL MPI_TYPE_INDEXED(ny,lengths,starts,mpireal,subtype_field,errcode)
    CALL MPI_TYPE_COMMIT(subtype_field,errcode)
    DEALLOCATE(lengths,starts)

    subtype_particle_var=Create_Particle_Subtype(npart_local)

  END SUBROUTINE CreateSubtypes

  FUNCTION Create_Particle_Subtype(nPart_Local)

    INTEGER(KIND=8),INTENT(IN) :: nPart_Local
    INTEGER :: Create_Particle_Subtype

    INTEGER,DIMENSION(:),ALLOCATABLE :: lengths, starts


    ! Create the subarray for the particles in this problem: subtype decribes where this
    ! process's data fits into the global picture.
    CALL MPI_ALLGATHER(npart_local,1,MPI_INTEGER8,npart_each_rank,1,MPI_INTEGER8,comm,errcode)

    ALLOCATE(lengths(1),starts(1))
    lengths=npart_local
    starts=0
    DO ix=1,rank
       starts=starts+npart_each_rank(ix)
    ENDDO

    CALL MPI_TYPE_INDEXED(1,lengths,starts,mpireal,Create_Particle_Subtype,errcode)
    CALL MPI_TYPE_COMMIT(Create_Particle_Subtype,errcode)

    DEALLOCATE(lengths,starts)

  END FUNCTION Create_Particle_Subtype

  SUBROUTINE CreateSubtypesForLoad(npart_local)

    !This subroutines creates the MPI types which represent the data for the field and
    !particles data. It is used when reading data. To this end, it takes npart_local
    !rather than determining it from the data structures

    INTEGER, DIMENSION(:),ALLOCATABLE :: lengths,starts
    INTEGER(KIND=8),INTENT(IN) :: npart_local


    ALLOCATE(lengths(1:ny),starts(1:ny))

    lengths=nx
    DO iy=0,ny-1
       starts(iy+1)=(cell_y_start(coordinates(1)+1)+iy-1) * nx_global + cell_x_start(coordinates(2)+1) -1
    ENDDO

    CALL MPI_TYPE_INDEXED(ny,lengths,starts,mpireal,subtype_field,errcode)
    CALL MPI_TYPE_COMMIT(subtype_field,errcode)
    DEALLOCATE(lengths,starts)

    ! Create the subarray for the particles in this problem: subtype decribes where this
    ! process's data fits into the global picture.
    CALL MPI_ALLGATHER(npart_local,1,MPI_INTEGER8,npart_each_rank,1,MPI_INTEGER8,comm,errcode)


    ALLOCATE(lengths(1),starts(1))
    lengths=npart_local
    starts=0
    DO ix=1,rank
       starts=starts+npart_each_rank(ix)
    ENDDO

    CALL MPI_TYPE_INDEXED(1,lengths,starts,mpireal,subtype_particle_var,errcode)
    CALL MPI_TYPE_COMMIT(subtype_particle_var,errcode)

    CALL MPI_TYPE_INDEXED(1,lengths,starts,MPI_INTEGER,subtype_particle_int,errcode)
    CALL MPI_TYPE_COMMIT(subtype_particle_int,errcode)
    DEALLOCATE(lengths,starts)


  END SUBROUTINE CreateSubtypesForLoad


END MODULE balance
