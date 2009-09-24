MODULE balance

  USE shared_data
  USE partlist
  USE boundary
  USE mpi_subtype_control

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Balance_Workload(OverRide)

    !This subroutine determines whether or not the code needs rebalancing, calculates where to split the domain
    !and calls other subroutines to actually rearrange the fields and particles onto the new processors

    !This is really, really hard to do properly
    !So cheat

    LOGICAL,INTENT(IN) :: OverRide
    INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE :: npart_each_rank
    INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE :: Density_x
    INTEGER,DIMENSION(:),ALLOCATABLE,TARGET :: Starts_x,ends_x
    INTEGER :: new_cell_x_start,new_cell_x_end
    REAL(num) :: balance_frac,balance_frac_x
    INTEGER(KIND=8) :: Max_x,wk,Min_x,npart_local
    INTEGER :: iProc,iSpecies,iLaser
    INTEGER, DIMENSION(2) :: Domain
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



    !Now need to calculate the start and end points for the new domain on the current processor
    new_cell_x_start=starts_x(coordinates(1)+1)
    new_cell_x_end=ends_x(coordinates(1)+1)

    !Redeistribute the field variables
    domain=(/new_cell_x_start,new_cell_x_end/)
    CALL Redistribute_Fields(domain)

    !Copy the new lengths into the permanent variables
    cell_x_start=starts_x
    cell_x_end=ends_x

    !Set the new nx and ny
    nx=new_cell_x_end-new_cell_x_start+1

    !Do X and Y arrays separatly because we already have global copies of X and Y
    DEALLOCATE(x)
    ALLOCATE(x(-2:nx+3))
    x(0:nx+1)=x_global(new_cell_x_start-1:new_cell_x_end+1)

    !Reallocate the kinetic energy calculation
    DEALLOCATE(ekbar,ekbar_sum,ct)
    ALLOCATE(ekbar(1:nx,1:nspecies),ekbar_sum(-2:nx+3,1:nspecies))
    ALLOCATE(ct(-2:nx+3,1:nspecies))

    !Recalculate x_starts and y_starts so that rebalancing works next time
    DO iproc=0,nprocx-1
       x_starts(iproc)=x_global(cell_x_start(iproc+1)) 
       x_ends(iproc)=x_global(cell_x_end(iproc+1))
    ENDDO

    !Set the lengths of the current domain so that the particle balancer works properly
    x_start_local=x_starts(coordinates(1))
    x_end_local=x_ends(coordinates(1))

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

  SUBROUTINE Redistribute_Fields(new_domain)

    !This subroutine redistributes the 2D field variables over the new processor layout
    !If using a 2D field of your own then se the Redistribute_Field subroutine to implement it
    !1D fields, you're on your own (have global copies and use those to repopulate?)
    INTEGER :: nx_new,iLaser
    INTEGER,DIMENSION(2),INTENT(IN) :: new_domain
    REAL(num),DIMENSION(:), ALLOCATABLE :: temp
    REAL(num),DIMENSION(:), ALLOCATABLE :: temp1d
    TYPE(Laser_Block),POINTER :: Current

    nx_new=new_domain(2)-new_domain(1)+1

    ALLOCATE(temp(-2:nx_new+3))

    temp=0.0_num
    CALL Redistribute_Field(new_domain,Ex,temp)
    DEALLOCATE(Ex)
    ALLOCATE(Ex(-2:nx_new+3))
    Ex=temp


    temp=0.0_num
    CALL Redistribute_Field(new_domain,Ey,temp)
    DEALLOCATE(Ey)
    ALLOCATE(Ey(-2:nx_new+3))
    Ey=temp

    temp=0.0_num
    CALL Redistribute_Field(new_domain,Ez,temp)
    DEALLOCATE(Ez)
    ALLOCATE(Ez(-2:nx_new+3))
    Ez=temp

    temp=0.0_num
    CALL Redistribute_Field(new_domain,Bx,temp)
    DEALLOCATE(Bx)
    ALLOCATE(Bx(-2:nx_new+3))
    Bx=temp


    temp=0.0_num
    CALL Redistribute_Field(new_domain,By,temp)
    DEALLOCATE(By)
    ALLOCATE(By(-2:nx_new+3))
    By=temp

    temp=0.0_num
    CALL Redistribute_Field(new_domain,Bz,temp)
    DEALLOCATE(Bz)
    ALLOCATE(Bz(-2:nx_new+3))
    Bz=temp

	temp=0.0_num
	CALL Redistribute_Field(new_domain,Jx,temp)
	DEALLOCATE(Jx)
	ALLOCATE(Jx(-2:nx_new+3))
	Jx=temp

	temp=0.0_num
	CALL Redistribute_Field(new_domain,Jy,temp)
	DEALLOCATE(Jy)
	ALLOCATE(Jy(-2:nx_new+3))
	Jy=temp

	temp=0.0_num
	CALL Redistribute_Field(new_domain,Jz,temp)
	DEALLOCATE(Jz)
	ALLOCATE(Jz(-2:nx_new+3))
	Jz=temp

    DEALLOCATE(temp)

!No need to rebalance lasers in 1D, lasers are just a point!


  END SUBROUTINE Redistribute_Fields

  SUBROUTINE Redistribute_Field(domain,Field,NewField)

    !This subroutine redistributes the fields over the new processor layout
    !The current version works by writing the field to a file and then each processor
    !Loads back in it's own part. This is better than the previous version where
    !Each processor produced it's own copy of the global array and then took
    !It's own subsection
    INTEGER,DIMENSION(2),INTENT(IN) :: domain
    REAL(num),DIMENSION(-2:),INTENT(IN) :: Field
    REAL(num),DIMENSION(-2:),INTENT(OUT) :: NewField
    INTEGER :: nx_new,ixd,iyd,izd
    INTEGER :: subtype_write, subtype_read, fh
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset=0
    CHARACTER(LEN=9+Data_Dir_Max_Length+n_zeros) :: filename

    WRITE(filename, '(a,"/balance.dat")') TRIM(data_dir)

    nx_new=domain(2)-domain(1)+1

    CALL MPI_FILE_OPEN(comm,TRIM(Filename),MPI_MODE_RDWR+MPI_MODE_CREATE,MPI_INFO_NULL,fh,errcode)
    subtype_write = Create_Current_Field_Subtype()
    subtype_read  = Create_Field_Subtype(nx_new,domain(1))

    CALL MPI_FILE_SET_VIEW(fh,offset,mpireal,subtype_write,"native",MPI_INFO_NULL,errcode)
    CALL MPI_FILE_WRITE_ALL(fh,Field(1:nx),nx,mpireal,status,errcode)
    CALL MPI_BARRIER(comm,errcode)
    CALL MPI_FILE_SEEK(fh,offset,MPI_SEEK_SET,errcode)
    CALL MPI_FILE_SET_VIEW(fh,offset,mpireal,subtype_read,"native",MPI_INFO_NULL,errcode)
    CALL MPI_FILE_READ_ALL(fh,NewField(1:nx_new),nx_new,mpireal,status,errcode)
    CALL MPI_FILE_CLOSE(fh,errcode)
    CALL MPI_BARRIER(comm,errcode)

    CALL MPI_TYPE_FREE(subtype_write,errcode)
    CALL MPI_TYPE_FREE(subtype_read,errcode)

    CALL Do_Field_MPI_With_Lengths(NewField,nx_new)

  END SUBROUTINE Redistribute_Field

!!$  SUBROUTINE Redistribute_Field(new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end,Field_In,Field_Out)
!!$
!!$    !This subroutine redistributes the fields over the new processor layout
!!$    !The current version works by producing a global copy on each processor
!!$    !And then extracting the required part for the local processor.
!!$    !This is not in general a good idea
!!$    INTEGER, INTENT(IN) :: new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end
!!$    REAL(num),DIMENSION(-2:,-2:),INTENT(IN) :: Field_In
!!$    REAL(num),DIMENSION(-2:,-2:),INTENT(OUT) :: Field_Out
!!$    REAL(num),DIMENSION(:,:),ALLOCATABLE :: Field_New,Field_Temp
!!$    INTEGER :: nx_new,ny_new
!!$    INTEGER :: comm_new,iproc,color
!!$
!!$    nx_new=new_cell_x_end-new_cell_x_start+1
!!$    ny_new=new_cell_y_end-new_cell_y_start+1
!!$
!!$!    PRINT *,rank,new_cell_x_start,new_cell_y_start," "
!!$
!!$    !This is a horrible, horrible way of doing this, I MUST think of a better way
!!$
!!$    !Create a global copy of the whole array
!!$    ALLOCATE(Field_New(1:nx_global,1:ny_global),Field_Temp(1:nx_global,1:ny_global))
!!$    Field_New=0.0_num
!!$    Field_New(cell_x_start(coordinates(2)+1):cell_x_end(coordinates(2)+1),&
!!$         cell_y_start(coordinates(1)+1):cell_y_end(coordinates(1)+1))=Field_In(1:nx,1:ny)
!!$    CALL MPI_ALLREDUCE(Field_New,Field_Temp,nx_global*ny_global,mpireal,MPI_SUM,comm,errcode)
!!$
!!$    Field_Out(1:nx_new,1:ny_new)=Field_Temp(new_cell_x_start:new_cell_x_end,new_cell_y_start:new_cell_y_end)
!!$    DEALLOCATE(Field_Temp)
!!$
!!$    !Call boundary conditions (this does not include any special BCS)
!!$    CALL Do_Field_MPI_With_Lengths(Field_out,nx_new,ny_new)
!!$
!!$  END SUBROUTINE Redistribute_Field

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
          part_x=Current%Part_Pos-x_start
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
    INTEGER :: CurrentLoc, iproc,coords(1)
    GetParticleProcessor = -1
    coords=-1

    !This could be replaced by a bisection method, but for the moment I just don't care


    DO iproc=0,nprocx-1
       IF (aParticle%Part_Pos .GE. x_starts(iproc) - dx/2.0_num .AND. aParticle%Part_Pos .LE. x_ends(iproc) + dx/2.0_num) THEN
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
             PRINT *,"Unlocatable particle on processor",rank,Current%Part_Pos,dx
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

END MODULE balance
