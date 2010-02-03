MODULE dist_fn

  USE shared_data
  USE partlist
  USE output_cartesian 
  USE output
  USE iocontrol
  IMPLICIT NONE

CONTAINS

  SUBROUTINE Attach_Dist_Fn(Block)
    TYPE(Distribution_Function_Block),POINTER :: Block
    TYPE(Distribution_Function_Block),POINTER :: Current

    Current=>Dist_Fns
    IF (.NOT. ASSOCIATED(Current)) THEN
       !This is the first distribution function to add
       Dist_Fns=>Block
       RETURN
    ENDIF
    DO WHILE(ASSOCIATED(Current%Next))
       Current=>Current%Next
    ENDDO
    Current%Next=>Block

  END SUBROUTINE Attach_Dist_Fn

  SUBROUTINE Setup_Dist_Fn(Block)
    TYPE(Distribution_Function_Block),POINTER :: Block

    NULLIFY(Block%Next)
    ALLOCATE(Block%Use_Species(1:nSpecies))
    Block%Use_Species=.FALSE.
    Block%Ranges=1.0_num
    Block%Use_Restrictions=.FALSE.
    Block%nDims=-1
  END SUBROUTINE Setup_Dist_Fn

  SUBROUTINE Clean_Dist_Fns

    TYPE(Distribution_Function_Block),POINTER :: Current,Next

    Current=>Dist_Fns
    DO WHILE(ASSOCIATED(Current))
       Next=>Current%Next
       DEALLOCATE(Current)
       Current=>Next
    ENDDO

  END SUBROUTINE Clean_Dist_Fns

  SUBROUTINE Write_Dist_Fns(code)
    INTEGER,INTENT(IN) :: code

    INTEGER :: iSpecies
    REAL(num),DIMENSION(3,2) :: ranges
    LOGICAL,DIMENSION(5) :: use_restrictions
    REAL(num),DIMENSION(5,2) :: restrictions
    INTEGER,DIMENSION(3) :: resolution=(/100,100,100/)
    TYPE(Distribution_Function_Block),POINTER :: Current

    !Write the distribution functions
    Current=>Dist_Fns
    DO WHILE(ASSOCIATED(Current))
       IF (IAND(Current%Dumpmask,code) .NE. 0) THEN
          DO iSpecies=1,nSpecies
             IF (.NOT. Current%Use_Species(iSpecies)) CYCLE
             ranges=Current%Ranges
             resolution=Current%Resolution
             restrictions=Current%restrictions
             use_restrictions=Current%Use_Restrictions

             IF (Current%nDims .EQ. 2) THEN
                CALL general_2d_dist_fn(Current%Name,Current%Directions(1:2),ranges(1:2,:),resolution(1:2),iSpecies,restrictions,use_restrictions)
             ELSE
                CALL general_3d_dist_fn(Current%Name,Current%Directions,ranges,resolution,iSpecies,restrictions,use_restrictions)
             ENDIF
             IF (Current%Store_Ranges) Current%Ranges=ranges
          ENDDO
       ENDIF
       Current=>Current%Next
    ENDDO
  END SUBROUTINE Write_Dist_Fns


  SUBROUTINE general_3d_dist_fn(name,direction,range,resolution,species,restrictions,use_restrictions)

    CHARACTER(LEN=*),INTENT(IN) :: name
    INTEGER,DIMENSION(3),INTENT(IN) :: direction
    REAL(num),DIMENSION(3,2),INTENT(INOUT) :: range
    INTEGER,DIMENSION(3),INTENT(INOUT) :: resolution
    INTEGER,INTENT(IN) :: species
    REAL(num),DIMENSION(5,2),INTENT(IN) :: restrictions
    LOGICAL,DIMENSION(5),INTENT(IN) :: use_restrictions

    REAL(num),DIMENSION(:,:,:),ALLOCATABLE :: Data, Data2
    REAL(num),DIMENSION(:),ALLOCATABLE :: Grid1,Grid2,Grid3
    LOGICAL,DIMENSION(3) :: Parallel
    REAL(num),DIMENSION(3) :: dgrid
    REAL(num) :: CurrentData, temp_data
    INTEGER :: iDim,iDir
    LOGICAL,DIMENSION(3) :: Calc_Range
    LOGICAL :: Calc_Ranges

    LOGICAL :: Use_x,need_reduce
    INTEGER,DIMENSION(3) :: Start_Local,global_resolution
    INTEGER :: color
    INTEGER :: comm_new, type_new

    LOGICAL,DIMENSION(3,4) :: Use_Direction
    LOGICAL,DIMENSION(3) :: Calc_Mod
    INTEGER,DIMENSION(3) :: P_Count
    INTEGER,DIMENSION(3) :: l_direction
    REAL(num),DIMENSION(3) :: conv 
    INTEGER,DIMENSION(3) :: cell
    LOGICAL :: UseThis
    REAL(num) :: RealSpace_Area


    TYPE(Particle),POINTER :: Current
    CHARACTER(LEN=ENTRYLENGTH) :: Grid_Name, Norm_Grid_Name, Var_Name
    REAL(num),DIMENSION(3) :: stagger=0.0_num
    REAL(num),DIMENSION(4) :: Particle_Data

    REAL(num) :: Max_P_Conv

    INTEGER :: ind

    Use_x=.FALSE.
    Need_Reduce=.TRUE.
    color=0
    global_resolution=resolution
    Parallel=.FALSE.
    start_local=1
    Calc_Range=.FALSE.
    calc_ranges=.FALSE.
    P_Count=0
    Use_Direction=.FALSE.
    l_Direction=0

    RealSpace_Area=1.0_num


    DO iDim=1,3
       IF (direction(iDim) .EQ. DIR_X) THEN
          Use_x=.TRUE.
          resolution(iDim) = nx
          range(iDim,1)=x_start_local
          range(iDim,2)=x_end_local
          start_local(iDim)=cell_x_start(rank+1)
          global_resolution(iDim)=nx_global
          dgrid(iDim)=dx
          Parallel(iDim)=.TRUE.
          Use_Direction(iDim,1)=.TRUE.
          l_Direction(iDim)=1
          conv(iDim)=length_x
          CYCLE
       ENDIF
       conv(iDim)=c*m0
       !If we're here then this must be a momentum space direction
       !So determine which momentum space directions are needed
       IF (range(iDim,1) .EQ. range(iDim,2)) THEN
          calc_range(iDim) = .TRUE.
          calc_ranges=.TRUE.
       ENDIF
       IF (IAND(direction(iDim),DIR_PX) .NE. 0) THEN
          Use_Direction(iDim,2)=.TRUE.
          P_Count(iDim)=P_Count(iDim)+1
          l_Direction(iDim)=2
       ENDIF
       IF (IAND(direction(iDim),DIR_PY) .NE. 0) THEN
          Use_Direction(iDim,3)=.TRUE.
          P_Count(iDim)=P_Count(iDim)+1
          l_Direction(iDim)=3
       ENDIF
       IF (IAND(direction(iDim),DIR_PZ) .NE. 0) THEN
          Use_Direction(iDim,4)=.TRUE.
          P_Count(iDim)=P_Count(iDim)+1
          l_Direction(iDim)=4
       ENDIF
    ENDDO

    IF (.NOT. use_x) RealSpace_Area=RealSpace_Area*dx

    DO iDim=1,3
       IF (P_Count(iDim) .GT. 1) calc_mod(iDim)=.TRUE.
    ENDDO
    !Calculate ranges for directions where needed
    IF (Calc_Ranges) THEN
       DO iDim=1,3
          IF (calc_range(iDim)) THEN
             range(iDim,1)=1.0e6_num
             range(iDim,2)=-1.0e6_num
          ENDIF
       ENDDO
       Current=>ParticleSpecies(Species)%AttachedList%Head
       ind =0
       DO WHILE(ASSOCIATED(Current))
          Particle_Data(1)=Current%part_pos
          Particle_Data(2:4)=Current%Part_P
          DO iDim=1,3
             IF (calc_range(iDim)) THEN
                IF (calc_mod(iDim)) THEN
                   DO iDir=1,4
                      IF (Use_Direction(iDim,iDir)) CurrentData=CurrentData+Particle_Data(iDir)
                   ENDDO
                   CurrentData=SQRT(CurrentData)
                ELSE
                   CurrentData=Particle_Data(l_Direction(iDim))
                ENDIF
                UseThis=.TRUE.
                DO iDir=1,4
                   IF (Use_Restrictions(iDir) .AND. (Particle_Data(iDir) .LT. Restrictions(iDir,1) &
                        .OR. Particle_Data(iDir) .GT. Restrictions(iDir,2))) UseThis=.FALSE.
                ENDDO
                IF (UseThis) THEN
                   IF (CurrentData .LT. range(iDim,1)) range(iDim,1)=CurrentData
                   IF (CurrentData .GT. range(iDim,2)) range(iDim,2)=CurrentData
                ENDIF
             ENDIF
          ENDDO
          ind=ind+1
          Current=>Current%Next
       ENDDO
    ENDIF

    Max_P_Conv=-10.0_num
    DO iDim=1,3
       IF (.NOT. Parallel(iDim)) THEN
          !If not parallel then this is a momentum dimension
          CALL MPI_ALLREDUCE(range(iDim,1),temp_data,1,mpireal,MPI_MIN,comm,errcode)
          range(iDim,1)=temp_data
          CALL MPI_ALLREDUCE(range(iDim,2),temp_data,1,mpireal,MPI_MAX,comm,errcode)
          range(iDim,2)=temp_data
       ENDIF
       !Fix so that if distribution function is zero then it picks an arbitrary scale in that direction
       IF (range(iDim,1) .EQ. range(iDim,2)) THEN
          range(iDim,1)=-1.0_num
          range(iDim,2)=1.0_num
       ENDIF
       !Calculate the maxmium range of a momentum direction
       IF (range(iDim,2)-range(iDim,1) .GT. Max_P_Conv .AND. .NOT. Parallel(iDim)) Max_P_Conv=range(iDim,2)-range(iDim,1)
    ENDDO

    DO iDim=1,3
       IF (.NOT. Parallel(iDim)) conv(iDim)=Max_P_Conv
    ENDDO

    !Setup MPI
    IF (Use_x) need_reduce=.FALSE.
    !In 1D, either fully parallel or fully reduced, no need for split communicator
    comm_new=comm

    Type_new=Create_3D_Field_Subtype(resolution,global_resolution,start_local)
    !Create grids
    DO iDim=1,3
       IF (.NOT. Parallel(iDim)) dgrid(iDim)=(range(iDim,2)-range(iDim,1))/REAL(resolution(iDim)-1,num)
    ENDDO
    ALLOCATE(grid1(0:global_resolution(1)),grid2(0:global_resolution(2)))
    ALLOCATE(grid3(0:global_resolution(3)))

    grid1(0)=-dgrid(1) + range(1,1)
    DO iDir=1,global_resolution(1)
       grid1(iDir)=Grid1(iDir-1)+dgrid(1)
    ENDDO

    grid2(0)=-dgrid(2) + range(2,1)
    DO iDir=1,global_resolution(2)
       grid2(iDir)=Grid2(iDir-1)+dgrid(2)
    ENDDO

    grid3(0)=-dgrid(3) + range(3,1)
    DO iDir=1,global_resolution(3)
       grid3(iDir)=Grid3(iDir-1)+dgrid(3)
    ENDDO

    ALLOCATE(Data(1:resolution(1),1:resolution(2),1:resolution(3)))
    Data=0.0_num

    Current=>ParticleSpecies(Species)%AttachedList%Head
    DO WHILE(ASSOCIATED(Current))
       Particle_Data(1)=Current%part_pos
       Particle_Data(2:4)=Current%Part_P
       UseThis=.TRUE.
       DO iDir=1,4
          IF (Use_Restrictions(iDir) .AND. (Particle_Data(iDir) .LT. Restrictions(iDir,1) &
               .OR. Particle_Data(iDir) .GT. Restrictions(iDir,2))) UseThis=.FALSE.
       ENDDO
       IF (UseThis) THEN
          DO iDim=1,3
             IF (calc_mod(iDim)) THEN
                DO iDir=1,4
                   IF (Use_Direction(iDim,iDir)) CurrentData=CurrentData+Particle_Data(iDir)
                ENDDO
                CurrentData=SQRT(CurrentData)
             ELSE
                CurrentData=Particle_Data(l_Direction(iDim))
             ENDIF
             cell(iDim)=NINT((CurrentData-range(iDim,1))/dgrid(iDim))+1
             IF (cell(iDim) .LT. 1 .OR. cell(iDim) .GT. resolution(iDim)) UseThis=.FALSE.
          ENDDO
          IF (UseThis) Data(cell(1),cell(2),cell(3))=Data(cell(1),cell(2),cell(3))+Current%Weight * realspace_area
       ENDIF
       Current=>Current%Next
    ENDDO

    IF (need_reduce) THEN
       ALLOCATE(Data2(1:resolution(1),1:resolution(2),1:resolution(3)))
       Data2=0.0_num
       CALL MPI_ALLREDUCE(Data,Data2,resolution(1)*resolution(2)*resolution(3),mpireal,MPI_SUM,comm_new,errcode)
       Data=Data2
       DEALLOCATE(Data2)
    ENDIF

    Grid_Name="Grid_"//TRIM(Name)//"_"//TRIM(ParticleSpecies(Species)%Name)
    Norm_Grid_Name="Norm_Grid_"//TRIM(Name)//"_"//TRIM(ParticleSpecies(Species)%Name)
    Var_Name=TRIM(Name)//"_"//TRIM(ParticleSpecies(Species)%Name)

    CALL cfd_Write_3D_Cartesian_Grid(TRIM(Grid_Name),"Grid",grid1(1:global_resolution(1)),grid2(1:global_resolution(2))&
         ,Grid3(1:global_resolution(3)),0)
    IF (Use_Offset_Grid) THEN
       IF (Parallel(1)) grid1=grid1-range(1,1)
       IF (Parallel(2)) grid2=grid2-range(2,1)
       IF (Parallel(3)) grid3=grid3-range(3,1)
    ENDIF
    CALL cfd_Write_3D_Cartesian_Grid(TRIM(Norm_Grid_Name),"Grid",grid1(1:global_resolution(1))/conv(1),grid2(1:global_resolution(2))/conv(2)&
         ,Grid3(1:global_resolution(3))/conv(3),0)

    CALL cfd_Write_3D_Cartesian_Variable_Parallel(TRIM(Var_Name),"dist_fn",global_resolution,Stagger,TRIM(Norm_Grid_Name),"Grid"&
         ,Data,Type_new)
    CALL MPI_Type_Free(Type_new,errcode)

    DEALLOCATE(Data)
    DEALLOCATE(Grid1,Grid2,Grid3)

  END SUBROUTINE general_3d_dist_fn

  SUBROUTINE general_2d_dist_fn(name,direction,range,resolution,species,restrictions,use_restrictions)

    CHARACTER(LEN=*),INTENT(IN) :: name
    INTEGER,DIMENSION(2),INTENT(IN) :: direction
    REAL(num),DIMENSION(2,2),INTENT(INOUT) :: range
    INTEGER,DIMENSION(2),INTENT(INOUT) :: resolution
    INTEGER,INTENT(IN) :: species
    REAL(num),DIMENSION(5,2),INTENT(IN) :: restrictions
    LOGICAL,DIMENSION(5),INTENT(IN) :: use_restrictions

    REAL(num),DIMENSION(:,:),ALLOCATABLE :: Data, Data2
    REAL(num),DIMENSION(:),ALLOCATABLE :: Grid1,Grid2
    LOGICAL,DIMENSION(2) :: Parallel
    REAL(num),DIMENSION(2) :: dgrid
    REAL(num) :: CurrentData,temp_data
    INTEGER :: iDim,iDir
    LOGICAL,DIMENSION(2) :: Calc_Range
    LOGICAL :: Calc_Ranges

    LOGICAL :: Use_x,need_reduce
    INTEGER,DIMENSION(2) :: Start_Local,global_resolution
    INTEGER :: color
    INTEGER :: comm_new, type_new

    LOGICAL,DIMENSION(2,4) :: Use_Direction
    LOGICAL,DIMENSION(2) :: Calc_Mod
    INTEGER,DIMENSION(2) :: P_Count
    INTEGER,DIMENSION(2) :: l_direction
    REAL(num),DIMENSION(2) :: conv 
    INTEGER,DIMENSION(2) :: cell
    LOGICAL :: UseThis
    REAL(num) :: RealSpace_Area


    TYPE(Particle),POINTER :: Current
    CHARACTER(LEN=ENTRYLENGTH) :: Grid_Name, Norm_Grid_Name, Var_Name
    REAL(num),DIMENSION(2) :: stagger=0.0_num
    REAL(num),DIMENSION(4) :: Particle_Data
    REAL(num) :: Max_P_Conv

    INTEGER :: ind


    Use_x=.FALSE.
    Need_Reduce=.TRUE.
    color=0
    global_resolution=resolution
    Parallel=.FALSE.
    start_local=1
    Calc_Range=.FALSE.
    calc_ranges=.FALSE.
    P_Count=0
    Use_Direction=.FALSE.
    l_Direction=0
    RealSpace_Area=1.0_num


    DO iDim=1,2
       IF (direction(iDim) .EQ. DIR_X) THEN
          Use_x=.TRUE.
          resolution(iDim) = nx
          range(iDim,1)=x_start_local
          range(iDim,2)=x_end_local
          !IF (Use_Offset_Grid) range(iDim,:)=range(iDim,:)-x_start
          start_local(iDim)=cell_x_start(rank+1)
          global_resolution(iDim)=nx_global
          dgrid(iDim)=dx
          Parallel(iDim)=.TRUE.
          Use_Direction(iDim,1)=.TRUE.
          l_Direction(iDim)=1
          conv(iDim)=length_x
          CYCLE
       ENDIF
       conv(iDim)=c*m0
       !If we're here then this must be a momentum space direction
       !So determine which momentum space directions are needed
       IF (range(iDim,1) .EQ. range(iDim,2)) THEN
          calc_range(iDim) = .TRUE.
          calc_ranges=.TRUE.
       ENDIF
       IF (IAND(direction(iDim),DIR_PX) .NE. 0) THEN
          Use_Direction(iDim,2)=.TRUE.
          P_Count(iDim)=P_Count(iDim)+1
          l_Direction(iDim)=2
       ENDIF
       IF (IAND(direction(iDim),DIR_PY) .NE. 0) THEN
          Use_Direction(iDim,3)=.TRUE.
          P_Count(iDim)=P_Count(iDim)+1
          l_Direction(iDim)=3
       ENDIF
       IF (IAND(direction(iDim),DIR_PZ) .NE. 0) THEN
          Use_Direction(iDim,4)=.TRUE.
          P_Count(iDim)=P_Count(iDim)+1
          l_Direction(iDim)=4
       ENDIF
    ENDDO

    IF (.NOT. use_x) RealSpace_Area=RealSpace_Area*dx

    DO iDim=1,2
       IF (P_Count(iDim) .GT. 1) calc_mod(iDim)=.TRUE.
    ENDDO
    !Calculate ranges for directions where needed
    IF (Calc_Ranges) THEN
       DO iDim=1,2
          IF (calc_range(iDim)) THEN
             range(iDim,1)=1.0e6_num
             range(iDim,2)=-1.0e6_num
          ENDIF
       ENDDO
       Current=>ParticleSpecies(Species)%AttachedList%Head
       ind =0
       DO WHILE(ASSOCIATED(Current))
          Particle_Data(1)=Current%part_pos
          Particle_Data(2:4)=Current%Part_P
          DO iDim=1,2
             IF (calc_range(iDim)) THEN
                IF (calc_mod(iDim)) THEN
                   DO iDir=1,4
                      IF (Use_Direction(iDim,iDir)) CurrentData=CurrentData+Particle_Data(iDir)
                   ENDDO
                   CurrentData=SQRT(CurrentData)
                ELSE
                   CurrentData=Particle_Data(l_Direction(iDim))
                ENDIF
                UseThis=.TRUE.
                DO iDir=1,4
                   IF (Use_Restrictions(iDir) .AND. (Particle_Data(iDir) .LT. Restrictions(iDir,1) &
                        .OR. Particle_Data(iDir) .GT. Restrictions(iDir,2))) UseThis=.FALSE.
                ENDDO
                IF (UseThis) THEN
                   IF (CurrentData .LT. range(iDim,1)) range(iDim,1)=CurrentData
                   IF (CurrentData .GT. range(iDim,2)) range(iDim,2)=CurrentData
                ENDIF
             ENDIF
          ENDDO
          ind=ind+1
          Current=>Current%Next
       ENDDO
    ENDIF

    Max_P_Conv=-10.0_num
    DO iDim=1,2
       IF (.NOT. Parallel(iDim)) THEN
          !If not parallel then this is a momentum dimension
          CALL MPI_ALLREDUCE(range(iDim,1),temp_data,1,mpireal,MPI_MIN,comm,errcode)
          range(iDim,1)=temp_data
          CALL MPI_ALLREDUCE(range(iDim,2),temp_data,1,mpireal,MPI_MAX,comm,errcode)
          range(iDim,2)=temp_data

          !Calculate the maxmium range of a momentum direction
          IF (range(iDim,2)-range(iDim,1) .GT. Max_P_Conv .AND. .NOT. Parallel(iDim)) Max_P_Conv=range(iDim,2)-range(iDim,1)
       ENDIF
       !Fix so that if distribution function is zero then it picks an arbitrary scale in that direction
       IF (range(iDim,1) .EQ. range(iDim,2)) THEN
          range(iDim,1)=-1.0_num
          range(iDim,2)=1.0_num
       ENDIF
    ENDDO

    DO iDim=1,2
       IF (.NOT. Parallel(iDim)) conv(iDim)=Max_P_Conv
    ENDDO


    !Setup MPI
    IF (Use_x) need_reduce=.FALSE.
    !In 1D don't need comm split since either fully parallel or fully reduced
    comm_new=comm

    Type_new=Create_2D_Field_Subtype(resolution,global_resolution,start_local)
    !Create grids
    DO iDim=1,2
       IF (.NOT. Parallel(iDim)) dgrid(iDim)=(range(iDim,2)-range(iDim,1))/REAL(resolution(iDim)-1,num)
    ENDDO
    ALLOCATE(grid1(0:global_resolution(1)),grid2(0:global_resolution(2)))

    grid1(0)=-dgrid(1) + range(1,1)
    DO iDir=1,global_resolution(1)
       grid1(iDir)=Grid1(iDir-1)+dgrid(1)
    ENDDO

    grid2(0)=-dgrid(2) + range(2,1)
    DO iDir=1,global_resolution(2)
       grid2(iDir)=Grid2(iDir-1)+dgrid(2)
    ENDDO

    ALLOCATE(Data(1:resolution(1),1:resolution(2)))
    Data=0.0_num

    Current=>ParticleSpecies(Species)%AttachedList%Head
    DO WHILE(ASSOCIATED(Current))
       Particle_Data(1)=Current%part_pos
       Particle_Data(2:4)=Current%Part_P
       UseThis=.TRUE.
       DO iDir=1,4
          IF (Use_Restrictions(iDir) .AND. (Particle_Data(iDir) .LT. Restrictions(iDir,1) &
               .OR. Particle_Data(iDir) .GT. Restrictions(iDir,2))) UseThis=.FALSE.
       ENDDO
       IF (UseThis) THEN
          DO iDim=1,2
             IF (calc_mod(iDim)) THEN
                DO iDir=1,4
                   IF (Use_Direction(iDim,iDir)) CurrentData=CurrentData+Particle_Data(iDir)
                ENDDO
                CurrentData=SQRT(CurrentData)
             ELSE
                CurrentData=Particle_Data(l_Direction(iDim))
             ENDIF
             cell(iDim)=NINT((CurrentData-range(iDim,1))/dgrid(iDim))+1
             IF (cell(iDim) .LT. 1 .OR. cell(iDim) .GT. resolution(iDim)) UseThis=.FALSE.
          ENDDO

          IF (UseThis)Data(cell(1),cell(2))=Data(cell(1),cell(2))+Current%Weight*RealSpace_Area
       ENDIF
       Current=>Current%Next
    ENDDO


    IF (need_reduce) THEN
       ALLOCATE(Data2(1:resolution(1),1:resolution(2)))
       Data2=0.0_num
       CALL MPI_ALLREDUCE(Data,Data2,resolution(1)*resolution(2),mpireal,MPI_SUM,comm_new,errcode)
       Data=Data2
       DEALLOCATE(Data2)
    ENDIF

    Grid_Name="Grid_"//TRIM(Name)//"_"//TRIM(ParticleSpecies(Species)%Name)
    Norm_Grid_Name="Norm_Grid_"//TRIM(Name)//"_"//TRIM(ParticleSpecies(Species)%Name)
    Var_Name=TRIM(Name)//"_"//TRIM(ParticleSpecies(Species)%Name)

    CALL cfd_Write_2D_Cartesian_Grid(TRIM(Grid_Name),"Grid",grid1(1:global_resolution(1)),grid2(1:global_resolution(2))&
         ,0)
    IF (Use_Offset_Grid) THEN
       IF (Parallel(1)) grid1=grid1-range(1,1)
       IF (Parallel(2)) grid2=grid2-range(2,1)
    ENDIF
    CALL cfd_Write_2D_Cartesian_Grid(TRIM(Norm_Grid_Name),"Grid",grid1(1:global_resolution(1))/conv(1),grid2(1:global_resolution(2))/conv(2)&
         ,0)

    CALL cfd_Write_2D_Cartesian_Variable_Parallel(TRIM(Var_Name),"dist_fn",global_resolution,Stagger,TRIM(Norm_Grid_Name),"Grid"&
         ,Data,Type_new)
    CALL MPI_Type_Free(Type_new,errcode)

    DEALLOCATE(Data)
    DEALLOCATE(Grid1,Grid2)

  END SUBROUTINE general_2d_dist_fn

  FUNCTION Create_3D_Field_Subtype(n_local,n_global,start)
    INTEGER,DIMENSION(3),INTENT(IN) :: n_local
    INTEGER,DIMENSION(3),INTENT(IN) :: n_global
    INTEGER,DIMENSION(3),INTENT(IN) :: start
    INTEGER, DIMENSION(:),ALLOCATABLE :: lengths,starts
    INTEGER :: iPoint,iy,iz
    INTEGER :: Create_3D_Field_Subtype
    ALLOCATE(lengths(1:n_local(2) * n_local(3)),starts(1:n_local(2) * n_local(3)))
    lengths=n_local(1)
    ipoint=0
    DO iz=0,n_local(3)-1
       DO iy=0,n_local(2)-1
          ipoint=ipoint+1
          starts(ipoint)=(start(3)+iz-1) * n_global(1) * n_global(2) + (start(2)+iy-1) * n_global(1) + start(1) -1
       ENDDO
    ENDDO

    CALL MPI_TYPE_INDEXED(n_local(2)*n_local(3),lengths,starts,mpireal,Create_3D_Field_Subtype,errcode)
    CALL MPI_TYPE_COMMIT(Create_3D_Field_Subtype,errcode)
    DEALLOCATE(lengths,starts)
  END FUNCTION Create_3D_Field_Subtype

  FUNCTION Create_2D_Field_Subtype(n_local,n_global,start)
    INTEGER,DIMENSION(2),INTENT(IN) :: n_local
    INTEGER,DIMENSION(2),INTENT(IN) :: n_global
    INTEGER,DIMENSION(2),INTENT(IN) :: start
    INTEGER, DIMENSION(:),ALLOCATABLE :: lengths,starts
    INTEGER :: iPoint
    INTEGER :: Create_2D_Field_Subtype
    ALLOCATE(lengths(1:n_local(2)),starts(1:n_local(2)))
    lengths=n_local(1)
    ipoint=0
    DO iy=0,n_local(2)-1
       ipoint=ipoint+1
       starts(ipoint)=(start(2)+iy-1) * n_global(1) + start(1) -1
    ENDDO

    CALL MPI_TYPE_INDEXED(n_local(2),lengths,starts,mpireal,Create_2D_Field_Subtype,errcode)
    CALL MPI_TYPE_COMMIT(Create_2D_Field_Subtype,errcode)
    DEALLOCATE(lengths,starts)
  END FUNCTION Create_2D_Field_Subtype

END MODULE dist_fn
