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
!!$    INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE :: npart_each_rank
    INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE :: Density_x, Density_y, Density_z
    INTEGER,DIMENSION(:),ALLOCATABLE,TARGET :: Starts_x, ends_x, Starts_y, ends_y, starts_z, ends_z
    INTEGER :: new_cell_x_start, new_cell_x_end
    INTEGER :: new_cell_y_start, new_cell_y_end
    INTEGER :: new_cell_z_start, new_cell_z_end
    REAL(num) :: balance_frac,balance_frac_x,balance_frac_y,balance_frac_z
    INTEGER(KIND=8) :: Max_x,Max_y,wk,Min_x,Min_y,npart_local,Max_z,Min_z
    INTEGER(KIND=8) :: max_npart,min_npart
    INTEGER :: iProc
    INTEGER,DIMENSION(3,2) :: domain
#ifdef PART_DEBUG
    TYPE(Particle),POINTER :: Current
#endif
    RETURN
    !On one processor do nothing to save time
    IF (nproc .EQ. 1) RETURN

    !This parameter allows selecting the mode of the autobalancing
    !Between leftsweep, rightsweep, auto(best of leftsweep and rightsweep) or both
    Balance_Mode=LB_ALL

    !Count particles
    npart_local=Get_Total_Local_Particles()

    !The OverRide flag allows the code to force a load balancing sweep at t=0
    IF (.NOT. OverRide) THEN
       CALL MPI_ALLREDUCE(npart_local,max_npart,1,MPI_INTEGER8,MPI_MAX,comm,errcode)
       CALL MPI_ALLREDUCE(npart_local,min_npart,1,MPI_INTEGER8,MPI_MIN,comm,errcode)
       balance_frac=REAL(min_npart,num)/REAL(max_npart,num)
       IF (balance_frac .GT. dlb_threshold) RETURN
       IF (rank .EQ. 0) PRINT *,"Load balancing with fraction",balance_frac
    ENDIF

    ALLOCATE(starts_x(1:nprocx),ends_x(1:nprocx))
    ALLOCATE(starts_y(1:nprocy),ends_y(1:nprocy))
    ALLOCATE(starts_z(1:nprocz),ends_z(1:nprocz))

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

    !Sweep in Z
    IF (IAND(Balance_Mode,LB_Z) .NE. 0 .OR. IAND(Balance_Mode,LB_AUTO) .NE. 0) THEN
       !Rebalancing in Z
       ALLOCATE(Density_z(0:nz_global+1))
       CALL GetDensityInZ(Density_z)
       CALL CalculateBreaks(Density_z,nprocz,starts_z,ends_z)
    ELSE
       !Just keep the original lengths
       starts_z=cell_z_start
       ends_z=cell_z_end
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
       Max_z=0
       Min_z=npart_global
       DO iproc=1,nprocz
          wk=SUM(Density_z(starts_z(iproc):ends_z(iproc)))
          IF (wk .GT. Max_z) Max_z=wk
          IF (wk .LT. Min_z) Min_z=wk
       ENDDO

       balance_frac_x=REAL(Min_x,num)/REAL(Max_x,num)
       balance_frac_y=REAL(Min_y,num)/REAL(Max_y,num)
       balance_frac_z=REAL(Min_z,num)/REAL(Max_z,num)

       IF (balance_frac_y .LT. balance_frac_x) THEN
          starts_y=cell_y_start
          ends_y=cell_y_end
       ELSE
          starts_x=cell_x_start
          ends_x=cell_x_end
       ENDIF

    ENDIF

    !Now need to calculate the start and end points for the new domain on the current processor
    new_cell_x_start=starts_x(coordinates(3)+1)
    new_cell_x_end=ends_x(coordinates(3)+1)

    new_cell_y_start=starts_y(coordinates(2)+1)
    new_cell_y_end=ends_y(coordinates(2)+1)

    new_cell_z_start=starts_z(coordinates(1)+1)
    new_cell_z_end=ends_z(coordinates(1)+1)

    domain(1,:)=(/new_cell_x_start,new_cell_x_end/)
    domain(2,:)=(/new_cell_y_start,new_cell_y_end/)
    domain(3,:)=(/new_cell_z_start,new_cell_z_end/)

    !Redeistribute the field variables
    CALL Redistribute_Fields(domain)

    !Copy the new lengths into the permanent variables
    cell_x_start=starts_x
    cell_y_start=starts_y
    cell_z_start=starts_z
    cell_x_end=ends_x
    cell_y_end=ends_y
    cell_z_end=ends_z

    !Set the new nx and ny
    nx=new_cell_x_end-new_cell_x_start+1
    ny=new_cell_y_end-new_cell_y_start+1
    nz=new_cell_z_end-new_cell_z_start+1

    !Do X and Y arrays separatly because we already have global copies of X and Y
    DEALLOCATE(x,y,z)
    ALLOCATE(x(-2:nx+3),y(-2:ny+3),z(-2:nz+3))
    x(0:nx+1)=x_global(new_cell_x_start-1:new_cell_x_end+1)
    y(0:ny+1)=y_global(new_cell_y_start-1:new_cell_y_end+1)
    z(0:nz+1)=z_global(new_cell_z_start-1:new_cell_y_end+1)

    !Reallocate Currents (don't need to keep current timestep values for current so just reallocate)
    DEALLOCATE(Jx,Jy,Jz)
    ALLOCATE(Jx(-2:nx+3,-2:ny+3,-2:nz+3),Jy(-2:nx+3,-2:ny+3,-2:nz+3),Jz(-2:nx+3,-2:ny+3,-2:nz+3))
    Jx=0.0_num
    Jy=0.0_num
    Jz=0.0_num

    !Reallocate the kinetic energy calculation
    DEALLOCATE(ekbar,ekbar_sum,ct)
    ALLOCATE(ekbar(1:nx,1:ny,1:nz,1:nspecies),ekbar_sum(-2:nx+3,-2:ny+3,-2:nz+3,1:nspecies))
    ALLOCATE(ct(-2:nx+3,-2:ny+3,-2:nz+3,1:nspecies))

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
    !Same for z
    DO iproc=0,nprocz-1
       z_starts(iproc)=z_global(cell_z_start(iproc+1))
       z_ends(iproc)=z_global(cell_z_end(iproc+1))
    ENDDO

    !Set the lengths of the current domain so that the particle balancer works properly
    x_start_local=x_starts(coordinates(3))
    x_end_local=x_ends(coordinates(3))
    y_start_local=y_starts(coordinates(2))
    y_end_local=y_ends(coordinates(2))
    z_start_local=z_starts(coordinates(1))
    z_end_local=z_ends(coordinates(1))

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
!!$    INTEGER, INTENT(IN) :: new_cell_x_start,new_cell_x_end,new_cell_y_start,new_cell_y_end
!!$    INTEGER, INTENT(IN) :: new_cell_z_start, new_cell_z_end
    INTEGER,DIMENSION(3,2),INTENT(IN) :: new_domain
    INTEGER,DIMENSION(2,2) :: new_domain_2d
    INTEGER :: nx_new, ny_new, nz_new
    TYPE(Laser_Block),POINTER :: Current
    REAL(num),DIMENSION(:,:,:),ALLOCATABLE :: Temp
    REAL(num),DIMENSION(:,:),ALLOCATABLE :: Temp2D

    nx_new=new_domain(1,2)-new_domain(1,1)+1
    ny_new=new_domain(2,2)-new_domain(2,1)+1
    nz_new=new_domain(3,2)-new_domain(3,1)+1

    !Domain is of form (direction,start/end)

    !This is used this way because some compilers (notably PGI)
    !Don't support the deallocation/reallocation of arrays when they
    !Are passed as parameters. This may change in future!

    ALLOCATE(Temp(-2:nx_new+3,-2:ny_new+3,-2:nz_new+3))
    Temp=0.0_num
    CALL Redistribute_Field(new_domain,Ex,Temp)
    DEALLOCATE(Ex)
    ALLOCATE(Ex(-2:nx_new+3,-2:ny_new+3,-2:nz_new+3))
    Ex=Temp

    Temp=0.0_num
    CALL Redistribute_Field(new_domain,Ey,Temp)
    DEALLOCATE(Ey)
    ALLOCATE(Ey(-2:nx_new+3,-2:ny_new+3,-2:nz_new+3))
    Ey=Temp

    Temp=0.0_num
    CALL Redistribute_Field(new_domain,Ez,Temp)
    DEALLOCATE(Ez)
    ALLOCATE(Ez(-2:nx_new+3,-2:ny_new+3,-2:nz_new+3))
    Ez=Temp

    Temp=0.0_num
    CALL Redistribute_Field(new_domain,Bx,Temp)
    DEALLOCATE(Bx)
    ALLOCATE(Bx(-2:nx_new+3,-2:ny_new+3,-2:nz_new+3))
    Bx=Temp

    Temp=0.0_num
    CALL Redistribute_Field(new_domain,By,Temp)
    DEALLOCATE(By)
    ALLOCATE(By(-2:nx_new+3,-2:ny_new+3,-2:nz_new+3))
    By=Temp

    Temp=0.0_num
    CALL Redistribute_Field(new_domain,Bz,Temp)
    DEALLOCATE(Bz)
    ALLOCATE(Bz(-2:nx_new+3,-2:ny_new+3,-2:nz_new+3))
    Bz=Temp
    DEALLOCATE(Temp)

    ALLOCATE(temp2d(-2:ny_new+3,-2:nz_new+3))
    new_domain_2d(1,:)=new_domain(2,:)
    new_domain_2d(2,:)=new_domain(3,:)
    Current=>Laser_Left
    DO WHILE(ASSOCIATED(Current))
       temp2d=0.0_num
       CALL Redistribute_Field_2D(new_domain_2D,DIR_X,Current%Profile,temp2d)
       DEALLOCATE(Current%Profile)
       ALLOCATE(Current%Profile(-2:ny_new+3,-2:nz_new+3))
       Current%Profile=temp2d
       temp2d=0.0_num
       CALL Redistribute_Field_2D(new_domain_2D,DIR_X,Current%Phase,temp2d)
       DEALLOCATE(Current%Phase)
       ALLOCATE(Current%Phase(-2:ny_new+3,-2:nz_new+3))
       Current%Phase=temp2d
       Current=>Current%Next
    ENDDO
    Current=>Laser_Right
    DO WHILE(ASSOCIATED(Current))
       temp2d=0.0_num
       CALL Redistribute_Field_2D(new_domain_2D,DIR_X,Current%Profile,temp2d)
       DEALLOCATE(Current%Profile)
       ALLOCATE(Current%Profile(-2:ny_new+3,-2:nz_new+3))
       Current%Profile=temp2d
       temp2d=0.0_num
       CALL Redistribute_Field_2D(new_domain_2D,DIR_X,Current%Phase,temp2d)
       DEALLOCATE(Current%Phase)
       ALLOCATE(Current%Phase(-2:ny_new+3,-2:nz_new+3))
       Current%Phase=temp2d
       Current=>Current%Next
    ENDDO
    new_domain_2d(1,:)=new_domain(1,:)
    new_domain_2d(2,:)=new_domain(3,:)
    Current=>Laser_Up
    DO WHILE(ASSOCIATED(Current))
       temp2d=0.0_num
       CALL Redistribute_Field_2D(new_domain_2D,DIR_Y,Current%Profile,temp2d)
       DEALLOCATE(Current%Profile)
       ALLOCATE(Current%Profile(-2:nx_new+3,-2:nz_new+3))
       Current%Profile=temp2d
       temp2d=0.0_num
       CALL Redistribute_Field_2D(new_domain_2D,DIR_Y,Current%Phase,temp2d)
       DEALLOCATE(Current%Phase)
       ALLOCATE(Current%Phase(-2:nx_new+3,-2:nz_new+3))
       Current%Phase=temp2d
       Current=>Current%Next
    ENDDO
    Current=>Laser_Down
    DO WHILE(ASSOCIATED(Current))
       temp2d=0.0_num
       CALL Redistribute_Field_2D(new_domain_2D,DIR_Y,Current%Profile,temp2d)
       DEALLOCATE(Current%Profile)
       ALLOCATE(Current%Profile(-2:nx_new+3,-2:nz_new+3))
       Current%Profile=temp2d
       temp2d=0.0_num
       CALL Redistribute_Field_2D(new_domain_2D,DIR_Y,Current%Phase,temp2d)
       DEALLOCATE(Current%Phase)
       ALLOCATE(Current%Phase(-2:nx_new+3,-2:nz_new+3))
       Current%Phase=temp2d
       Current=>Current%Next
    ENDDO
    new_domain_2d(1,:)=new_domain(2,:)
    new_domain_2d(2,:)=new_domain(3,:)
    Current=>Laser_Up
    DO WHILE(ASSOCIATED(Current))
       temp2d=0.0_num
       CALL Redistribute_Field_2D(new_domain_2D,DIR_Z,Current%Profile,temp2d)
       DEALLOCATE(Current%Profile)
       ALLOCATE(Current%Profile(-2:nx_new+3,-2:ny_new+3))
       Current%Profile=temp2d
       temp2d=0.0_num
       CALL Redistribute_Field_2D(new_domain_2D,DIR_Z,Current%Phase,temp2d)
       DEALLOCATE(Current%Phase)
       ALLOCATE(Current%Phase(-2:nx_new+3,-2:ny_new+3))
       Current%Phase=temp2d
       Current=>Current%Next
    ENDDO
    Current=>Laser_Down
    DO WHILE(ASSOCIATED(Current))
       temp2d=0.0_num
       CALL Redistribute_Field_2D(new_domain_2D,DIR_Z,Current%Profile,temp2d)
       DEALLOCATE(Current%Profile)
       ALLOCATE(Current%Profile(-2:nx_new+3,-2:ny_new+3))
       Current%Profile=temp2d
       temp2d=0.0_num
       CALL Redistribute_Field_2D(new_domain_2D,DIR_Z,Current%Phase,temp2d)
       DEALLOCATE(Current%Phase)
       ALLOCATE(Current%Phase(-2:nx_new+3,-2:ny_new+3))
       Current%Phase=temp2d
       Current=>Current%Next
    ENDDO

    DEALLOCATE(temp2d)

  END SUBROUTINE Redistribute_Fields

  SUBROUTINE Redistribute_Field(domain,Field,NewField)

    !This subroutine redistributes the fields over the new processor layout
    !The current version works by writing the field to a file and then each processor
    !Loads back in it's own part. This is better than the previous version where
    !Each processor produced it's own copy of the global array and then took
    !It's own subsection
    INTEGER,DIMENSION(3,2),INTENT(IN) :: domain
    REAL(num),DIMENSION(-2:,-2:,-2:),INTENT(IN) :: Field
    REAL(num),DIMENSION(-2:,-2:,-2:),INTENT(OUT) :: NewField
    INTEGER :: nx_new,ny_new,nz_new
    INTEGER :: subtype_write, subtype_read, fh
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset=0
    CHARACTER(LEN=9+Data_Dir_Max_Length+n_zeros) :: filename

    WRITE(filename, '(a,"/temp.dat")') TRIM(data_dir)

    nx_new=domain(1,2)-domain(1,1)+1
    ny_new=domain(2,2)-domain(2,1)+1
    nz_new=domain(3,2)-domain(3,1)+1

    CALL MPI_FILE_OPEN(comm,TRIM(Filename),MPI_MODE_RDWR+MPI_MODE_CREATE,MPI_INFO_NULL,fh,errcode)
    subtype_write = Create_Current_Field_Subtype()
    subtype_read  = Create_Field_Subtype(nx_new,ny_new,nz_new,domain(1,1),domain(2,1),domain(3,1))


    CALL MPI_FILE_SET_VIEW(fh,offset,mpireal,subtype_write,"native",MPI_INFO_NULL,errcode)
    CALL MPI_FILE_WRITE_ALL(fh,Field(1:nx,1:ny,1:nz),nx*ny*nz,mpireal,status,errcode)
    CALL MPI_BARRIER(comm,errcode)
    CALL MPI_FILE_SEEK(fh,offset,MPI_SEEK_SET,errcode)
    CALL MPI_FILE_SET_VIEW(fh,offset,mpireal,subtype_read,"native",MPI_INFO_NULL,errcode)
    CALL MPI_FILE_READ_ALL(fh,NewField(1:nx_new,1:ny_new,1:nz_new),nx_new*ny_new*nz_new,mpireal,status,errcode)
    CALL MPI_FILE_CLOSE(fh,errcode)
    CALL MPI_BARRIER(comm,errcode)

    CALL MPI_TYPE_FREE(subtype_write,errcode)
    CALL MPI_TYPE_FREE(subtype_read,errcode)

    CALL Do_Field_MPI_With_Lengths(NewField,nx_new,ny_new,nz_new)

  END SUBROUTINE Redistribute_Field

  SUBROUTINE Redistribute_Field_2D(domain,Direction,Field,NewField)

    !This subroutine redistributes the fields over the new processor layout
    !The current version works by writing the field to a file and then each processor
    !Loads back in it's own part. This is better than the previous version where
    !Each processor produced it's own copy of the global array and then took
    !It's own subsection
    INTEGER,DIMENSION(2,2),INTENT(IN) :: domain
    INTEGER, INTENT(IN) :: Direction
    REAL(num),DIMENSION(-2:,-2:),INTENT(IN) :: Field
    REAL(num),DIMENSION(-2:,-2:),INTENT(OUT) :: NewField
    INTEGER :: n1, n2, n1_new, n2_new, n1_global, n2_global
    INTEGER :: n1_start, n2_start
    INTEGER :: subtype_write, subtype_read, fh
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset=0
    CHARACTER(LEN=9+Data_Dir_Max_Length+n_zeros) :: filename

    IF (Direction .EQ. DIR_X) THEN
       n1=ny
       n2=nz
       n1_global=ny_global
       n2_global=nz_global
       n1_start=cell_y_start(coordinates(2)+1)
       n2_start=cell_z_start(coordinates(1)+1)
    ENDIF

    IF (Direction .EQ. DIR_Y) THEN
       n1=nx
       n2=nz
       n1_global=nx_global
       n2_global=nz_global
       n1_start=cell_x_start(coordinates(3)+1)
       n2_start=cell_z_start(coordinates(1)+1)
    ENDIF

    IF (Direction .EQ. DIR_Z) THEN
       n1=nx
       n2=ny
       n1_global=nx_global
       n2_global=ny_global
       n1_start=cell_x_start(coordinates(3)+1)
       n2_start=cell_y_start(coordinates(2)+1)
    ENDIF


    WRITE(filename, '(a,"/temp.dat")') TRIM(data_dir)

    n1_new=domain(1,2)-domain(1,1)+1
    n2_new=domain(2,2)-domain(2,1)+1

    CALL MPI_FILE_OPEN(comm,TRIM(Filename),MPI_MODE_RDWR+MPI_MODE_CREATE,MPI_INFO_NULL,fh,errcode)
    subtype_write = Create_2D_Array_Subtype((/n1,n2/),(/n1_global,n2_global/),&
         (/domain(1,1),domain(2,1)/))
    subtype_read  = Create_2D_Array_Subtype((/n1_new,n2_new/),(/n1_global,n2_global/),&
         (/domain(1,1),domain(2,1)/))


    CALL MPI_FILE_SET_VIEW(fh,offset,mpireal,subtype_write,"native",MPI_INFO_NULL,errcode)
    CALL MPI_FILE_WRITE_ALL(fh,Field(1:n1,1:n2),n1*n2,mpireal,status,errcode)
    CALL MPI_BARRIER(comm,errcode)
    CALL MPI_FILE_SEEK(fh,offset,MPI_SEEK_SET,errcode)
    CALL MPI_FILE_SET_VIEW(fh,offset,mpireal,subtype_read,"native",MPI_INFO_NULL,errcode)
    CALL MPI_FILE_READ_ALL(fh,NewField(1:n1_new,1:n2_new),n1_new*n2_new,mpireal,status,errcode)
    CALL MPI_FILE_CLOSE(fh,errcode)
    CALL MPI_BARRIER(comm,errcode)

    CALL MPI_TYPE_FREE(subtype_write,errcode)
    CALL MPI_TYPE_FREE(subtype_read,errcode)

  END SUBROUTINE Redistribute_Field_2D



  SUBROUTINE GetDensityInX(Density)

    !Calculate total particle density across the X direction
    !Summed in the Y and Z directions

    INTEGER(KIND=8),DIMENSION(:),INTENT(INOUT) :: Density
    INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE :: Temp
    TYPE(Particle),POINTER :: Current
    REAL(num) :: part_x
    INTEGER :: cell_x1,iSpecies
    Density=0.0_num

    DO iSpecies=1,nspecies
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO WHILE(ASSOCIATED(Current))
          !Want global position, so x_start, NOT x_start_local
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
    !Summed in the X and Z directions

    INTEGER(KIND=8),DIMENSION(:),INTENT(INOUT) :: Density
    INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE :: Temp
    TYPE(Particle),POINTER :: Current
    REAL(num) :: part_y
    INTEGER :: cell_y1,iSpecies
    Density=0.0_num

    DO iSpecies=1,nspecies
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO WHILE(ASSOCIATED(Current))
          !Want global position, so y_start, NOT y_start_local
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

  SUBROUTINE GetDensityInZ(Density)

    !Calculate total particle density across the Z direction
    !Summed in the X and Y directions

    INTEGER(KIND=8),DIMENSION(:),INTENT(INOUT) :: Density
    INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE :: Temp
    TYPE(Particle),POINTER :: Current
    REAL(num) :: part_z
    INTEGER :: cell_z1,iSpecies
    Density=0.0_num

    DO iSpecies=1,nspecies
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO WHILE(ASSOCIATED(Current))
          !Want global position, so z_start, NOT z_start_local
          part_z=Current%Part_Pos(3)-z_start
          cell_z1=NINT(part_z/dz)+1
          Density(cell_z1)=Density(cell_z1)+1
          Current=>Current%Next
       ENDDO
    ENDDO
    !Now have local densities, so add using MPI
    ALLOCATE(Temp(0:nz_global+1))
    CALL MPI_ALLREDUCE(Density,Temp,nz_global+2,MPI_INTEGER8,MPI_SUM,comm,errcode)
    Density=Temp
    DEALLOCATE(TEMP)

  END SUBROUTINE GetDensityInZ

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
    INTEGER :: iproc,coords(3)
    GetParticleProcessor = -1
    coords=-1

    !This could be replaced by a bisection method, but for the moment I just don't care


    DO iproc=0,nprocx-1
       IF (aParticle%Part_Pos(1) .GT. x_starts(iproc) - dx/2.0_num .AND. &
            aParticle%Part_Pos(1) .LE. x_ends(iproc) + dx/2.0_num) THEN
          coords(3)=iproc
          EXIT
       ENDIF
    ENDDO

    DO iproc=0,nprocy-1
       IF (aParticle%Part_Pos(2) .GT. y_starts(iproc) - dy/2.0_num .AND. &
            aParticle%Part_Pos(2) .LE. y_ends(iproc) + dy/2.0_num) THEN
          coords(2)=iproc
          EXIT
       ENDIF
    ENDDO
    DO iproc=0,nprocz-1
       IF (aParticle%Part_Pos(3) .GT. z_starts(iproc) - dz/2.0_num .AND. &
            aParticle%Part_Pos(3) .LE. z_ends(iproc) + dz/2.0_num) THEN
          coords(1)=iproc
          EXIT
       ENDIF
    ENDDO


    IF (MINVAL(coords) .LT. 0) PRINT *,"UNLOCATABLE PARTICLE",coords,dz
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
             PRINT *,"Unlocatable particle on processor",rank,Current%Part_Pos
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
