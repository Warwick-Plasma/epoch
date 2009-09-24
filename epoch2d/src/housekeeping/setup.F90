MODULE setup 
  USE shared_data
  USE input_cartesian
  USE input_particle
  USE input_arb
  USE input
  USE iocontrol
  USE inputfunctions
  USE Strings
  USE balance
  USE partlist
  USE mpi_subtype_control
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: after_control,minimal_init,restart_data
  PUBLIC :: open_files, close_files
  PUBLIC :: setup_species

  SAVE
  TYPE(ParticleList) :: MainRoot
  INTEGER, DIMENSION(:),ALLOCATABLE :: Species_ID

CONTAINS

  SUBROUTINE minimal_init

    IF (num .EQ. 4) mpireal=MPI_REAL
    DumpMask=0
    Comm=MPI_COMM_NULL

    Window_Shift=0.0_num
    npart_global=-1
	 Smooth_Currents=.FALSE.
	
	 Ex=0.0_num
	 Ey=0.0_num
	 Ez=0.0_num
	 Bx=0.0_num
	 By=0.0_num
	 Bz=0.0_num
	

    NULLIFY(Laser_Left)
    NULLIFY(Laser_Right)
    NULLIFY(Laser_Up)
    NULLIFY(Laser_Down)

    NULLIFY(Dist_Fns)
  END SUBROUTINE minimal_init

  SUBROUTINE after_control

    INTEGER :: iproc

    length_x=x_end-x_start
    dx = length_x / REAL(nx_global-1, num)
    length_y=y_end-y_start
    dy = length_y / REAL(ny_global-1, num)


    !Setup global grid
    x_global(0)=x_start-dx
    DO ix=1,nx_global+1
       x_global(ix)=x_global(ix-1)+dx
       x_offset_global(ix)=x_global(ix)
    ENDDO
    y_global(0)=y_start-dy
    DO iy=1,ny_global+1
       y_global(iy)=y_global(iy-1)+dy
       y_offset_global(iy)=y_global(iy)
    ENDDO

    DO iproc=0,nprocx-1
       x_starts(iproc)=x_global(iproc*nx+1) 
       x_ends(iproc)=x_global((iproc+1)*nx) 
    ENDDO
    DO iproc=0,nprocy-1
       y_starts(iproc)=y_global(iproc*ny+1) 
       y_ends(iproc)=y_global((iproc+1)*ny) 
    ENDDO
    x_start_local=x_starts(coordinates(2))
    x_end_local=x_ends(coordinates(2))
    y_start_local=y_starts(coordinates(1))
    y_end_local=y_ends(coordinates(1))      

    !Setup local grid
    x(-1)=x_start_local-dx*2.0_num
    DO ix=0,nx+2
       x(ix)=x(ix-1)+dx
    ENDDO
    y(-1)=y_start_local-dy*2.0_num
    DO iy=0,ny+2
       y(iy)=y(iy-1)+dy
    ENDDO

    CALL set_initial_values
	 CALL setup_data_averaging

  END SUBROUTINE after_control

  SUBROUTINE setup_data_averaging()

	INTEGER :: ioutput
   DO iOutput = 1, num_vars_to_dump
		IF(IAND(dumpmask(iOutput),IO_AVERAGED) .NE. 0) THEN
			ALLOCATE(averaged_data(iOutput)%data(-2:nx+3,-2:ny+3))
		ENDIF
	ENDDO

  END SUBROUTINE setup_data_averaging

  SUBROUTINE Setup_Species
    INTEGER :: iSpecies

    ALLOCATE(ParticleSpecies(1:nSpecies))

    DO iSpecies=1,nSpecies
       ParticleSpecies(iSpecies)%name=blank
       ParticleSpecies(iSpecies)%mass=-1.0_num
       ParticleSpecies(iSpecies)%Dump=.TRUE.
       ParticleSpecies(iSpecies)%Count=-1
#ifdef SPLIT_PARTICLES_AFTER_PUSH
       ParticleSpecies(iSpecies)%Split=.FALSE.
       ParticleSpecies(iSpecies)%npart_max=0.0_num     
#endif
#ifdef PART_IONISE
       ParticleSpecies(iSpecies)%ionise=.FALSE.
       ParticleSpecies(iSpecies)%ionise_to_species=.FALSE.
       ParticleSpecies(iSpecies)%release_species=-1
       ParticleSpecies(iSpecies)%critical_field=0.0_num
#endif
#ifdef TRACER_PARTICLES
       ParticleSpecies(iSpecies)%Tracer=.FALSE.
#endif
#ifdef PARTICLE_PROBES
       NULLIFY(ParticleSpecies(iSpecies)%AttachedProbes)
#endif
    ENDDO

  END SUBROUTINE Setup_Species

  SUBROUTINE open_files

    CHARACTER(LEN=11+Data_Dir_Max_Length) :: file2
    CHARACTER(LEN=7+Data_Dir_Max_Length) :: file3
    INTEGER :: Errcode

    IF (rank == 0) THEN
       WRITE(file2, '(a,"/epoch2d.dat")') TRIM(data_dir)
       OPEN(unit=20, STATUS = 'REPLACE',FILE = file2, IOSTAT=errcode)
       IF (errcode .NE. 0) THEN
             PRINT *,"***ERROR*** Cannot create epoch2d.dat output file. The most common cause of this problem is that the ouput directory does not exist"
             CALL MPI_ABORT(comm,errcode)
       ENDIF
    END IF

  END SUBROUTINE open_files


  SUBROUTINE close_files

    CLOSE(unit=20)
    CLOSE(unit=30)

  END SUBROUTINE close_files


  SUBROUTINE set_initial_values



  END SUBROUTINE set_initial_values

  SUBROUTINE restart_data

    CHARACTER(LEN=20+Data_Dir_Max_Length) :: filename
    CHARACTER(Len=MaxStringLen) :: Name,Class,MeshName,MeshClass
    INTEGER :: Type,nd
    INTEGER :: filehandle,sof
    INTEGER(KIND=8) :: npart_l,npart_per_it=10000,npart_local,npart_test
    REAL(num), DIMENSION(2) :: extents,stagger
    INTEGER,DIMENSION(1) :: dims
    REAL(KIND=8) :: time_d
    INTEGER :: snap,coord_type
    LOGICAL,DIMENSION(3) :: PisV
    REAL(num) :: v2
    TYPE(Particle), POINTER :: Current,Next
    LOGICAL :: Const_Weight
    INTEGER(KIND=8) :: npart

    npart=npart_global/nproc
    Const_Weight=.FALSE.
    CALL Create_SubTypes_For_Load(npart)

    ! Create the filename for the last snapshot
    WRITE(filename, '("nfs:",a,"/",i4.4,".cfd")') TRIM(data_dir), restart_snapshot
    CALL cfd_Open(filename,rank,comm,MPI_MODE_RDONLY)
    ! Open the file
    nBlocks=cfd_Get_nBlocks()

    Ex=0.0_num
    Ey=0.0_num
    Ez=0.0_num

    Bx=0.0_num
    By=0.0_num
    Bz=0.0_num

    IF (rank .EQ. 0) PRINT *,"Input file contains",nBlocks,"blocks"
    DO ix=1,nBlocks
       CALL cfd_Get_Next_Block_Info_All(Name,Class,Type)
       !IF (rank .EQ. 0) PRINT *,"Loading block",ix,Name,Type
       IF (Type .EQ. TYPE_SNAPSHOT) THEN
          CALL cfd_Get_Snapshot(time_d,snap)
          time=time_d
          IF (rank .EQ. 0) PRINT *,"Loading snapshot for time",time
       ENDIF
       SELECT CASE(Type)
       CASE(TYPE_MESH_VARIABLE)
          CALL cfd_Get_Common_MeshType_MetaData_All(Type,nd,sof)
          IF (sof .NE. num) THEN
             IF (rank .EQ. 0) PRINT *,"Precision does not match, recompile code so that sizeof(real) = ",sof
             CALL MPI_ABORT(comm,errcode)
          ENDIF

          IF (nd .NE. DIMENSION_2D .AND. nd .NE. DIMENSION_IRRELEVANT ) THEN
             IF (rank .EQ. 0) PRINT *,"Dimensionality does not match, file is ",nd,"D"
             CALL MPI_ABORT(comm,errcode)
          ENDIF
          SELECT CASE(Type)
          CASE(VAR_CARTESIAN)
             !Grid variables
             CALL cfd_Get_nD_Cartesian_Variable_MetaData_All(nd,dims,extents,stagger,meshname,meshclass)
             IF (dims(1) .NE. nx_global) THEN
                IF (rank .EQ. 0) PRINT *,"Number of gridpoints does not match, gridpoints in file is",dims(1)
                CALL MPI_ABORT(comm,errcode)
             ENDIF
             IF (StrCmp(Name(1:2),"Ex")) CALL cfd_Get_2D_Cartesian_Variable_Parallel(Ex(1:nx,1:ny),Subtype_Field)
             IF (StrCmp(Name(1:2),"Ey")) CALL cfd_Get_2D_Cartesian_Variable_Parallel(Ey(1:nx,1:ny),Subtype_Field)
             IF (StrCmp(Name(1:2),"Ez")) CALL cfd_Get_2D_Cartesian_Variable_Parallel(Ez(1:nx,1:ny),SubType_Field)

             IF (StrCmp(Name(1:2),"Bx")) CALL cfd_Get_2D_Cartesian_Variable_Parallel(Bx(1:nx,1:ny),Subtype_Field)
             IF (StrCmp(Name(1:2),"By")) CALL cfd_Get_2D_Cartesian_Variable_Parallel(By(1:nx,1:ny),Subtype_Field)
             IF (StrCmp(Name(1:2),"Bz")) CALL cfd_Get_2D_Cartesian_Variable_Parallel(Bz(1:nx,1:ny),Subtype_Field)

          CASE(VAR_PARTICLE)
             CALL cfd_Get_nD_Particle_Variable_MetaData_All(npart_l,extents,meshname,meshclass)
             IF (npart_l .NE. npart_global) THEN
                IF (rank .EQ. 0) PRINT *,"Number of particles does not match, changing npart to match",npart_l
                npart=npart_l/nproc
                CALL Create_SubTypes_For_Load(npart)
                CALL Create_Allocated_PartList(MainRoot,npart)
                ALLOCATE(Species_ID(1:npart))
                Current=>MainRoot%Head
                DO ipart=1,MainRoot%Count
                   Current%Part_P=0.0_num
                   Current%Part_Pos=0.0_num
                   Current=>Current%Next
                ENDDO
                Current=>MainRoot%Head
                npart_global=npart_l
             ENDIF
             !Particle variables
             IF (StrCmp(Name(1:2),"Px")) CALL cfd_Get_nD_Particle_Variable_Parallel_With_Iterator(npart,npart_per_it,subtype_particle_var,It_Px)
             IF (StrCmp(Name(1:2),"Py")) CALL cfd_Get_nD_Particle_Variable_Parallel_With_Iterator(npart,npart_per_it,subtype_particle_var,It_Py)
             IF (StrCmp(Name(1:2),"Pz")) CALL cfd_Get_nD_Particle_Variable_Parallel_With_Iterator(npart,npart_per_it,subtype_particle_var,It_Pz)
#ifdef PER_PARTICLE_WEIGHT
             IF (StrCmp(Name(1:6),"Weight")) CALL cfd_Get_nd_Particle_Variable_Parallel_With_Iterator(npart,npart_per_it,subtype_particle_var,It_Weight)
#else
             IF (StrCmp(Name(1:6),"Weight")) THEN
                IF (rank .EQ. 0) PRINT *,"Cannot load dump file with per particle weight if the code is compiled without per particle weights. Code terminates"
                CALL MPI_ABORT(comm,errcode)
             ENDIF
#endif
             IF (StrCmp(Name(1:7),"Species")) CALL cfd_Get_nD_Particle_Variable_Parallel_With_Iterator(npart,npart_per_it,subtype_particle_var,It_Species)
          END SELECT
       CASE(TYPE_MESH)
          CALL cfd_Get_Common_MeshType_MetaData_All(Type,nd,sof)
          IF (Type .EQ. MESH_PARTICLE) THEN
             CALL cfd_Get_nD_Particle_Grid_MetaData_All(nd,coord_type,npart_l,extents)
             IF (npart_l .NE. npart_global) THEN
                npart=npart_l/nproc
                IF (npart * nproc .NE. npart_l) THEN
                   IF (rank .EQ. 0) PRINT *,"Cannot evenly subdivide particles over",nproc,"processors. Trying to fix"
                   IF (rank .LT. npart_l-(npart*nproc)) npart=npart+1
                ENDIF
                CALL Create_SubTypes_For_Load(npart)
                CALL Create_Allocated_PartList(MainRoot,npart)
                ALLOCATE(Species_ID(1:npart))
                Current=>MainRoot%Head
                npart_global=npart_l
             ENDIF
             CALL cfd_Get_nD_Particle_Grid_Parallel_With_Iterator(nd,MainRoot%Count,npart_l,npart_per_it,sof,subtype_particle_var,It_Part)
          ENDIF
       CASE(TYPE_CONSTANT)
          CALL cfd_Get_Real_Constant(Weight)
          Const_Weight=.TRUE.
       END SELECT
       CALL cfd_Skip_Block()
    ENDDO
    CALL cfd_Close()


    Current=>MainRoot%Head
    ipart=1
    DO WHILE(ASSOCIATED(Current))
       Next=>Current%Next
       CALL Remove_Particle_From_PartList(MainRoot,Current)
       CALL Add_Particle_To_PartList(ParticleSpecies(Species_ID(ipart))%AttachedList,Current)
       Current=>Next
       ipart=ipart+1
    ENDDO

    DEALLOCATE(Species_ID)

#ifdef PER_PARTICLE_WEIGHT
    IF (Const_Weight) THEN
       Current=>MainRoot%Head
       DO WHILE(ASSOCIATED(Current))
          Current%Weight=Weight
          Current=>Current%Next
       ENDDO
    ENDIF
#endif

  END SUBROUTINE restart_data

  SUBROUTINE It_Part(Data,npart_this_it,Start,direction)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    INTEGER(KIND=8),INTENT(INOUT) :: npart_this_it
    LOGICAL,INTENT(IN) :: Start
    INTEGER,INTENT(IN) :: direction

    INTEGER(KIND=8) :: ipart
    TYPE(Particle),POINTER,Save :: Cur

    IF (Start) THEN 
       Cur=>MainRoot%Head
    ENDIF
    DO ipart=1,npart_this_it
       Cur%part_pos(direction) = Data(ipart)
       Cur=>Cur%Next
    ENDDO

  END SUBROUTINE It_Part

  SUBROUTINE It_Px(Data,npart_this_it,Start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    INTEGER(KIND=8),INTENT(INOUT) :: npart_this_it
    LOGICAL,INTENT(IN) :: Start
    INTEGER(KIND=8) :: ipart
    TYPE(Particle),POINTER,Save :: Cur

    IF (Start) Cur=>MainRoot%Head
    DO ipart=1,npart_this_it
       Cur%part_p(1) = Data(ipart)
       Cur=>Cur%Next
    ENDDO

  END SUBROUTINE It_Px

  SUBROUTINE It_Py(Data,npart_this_it,Start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    INTEGER(KIND=8),INTENT(INOUT) :: npart_this_it
    LOGICAL,INTENT(IN) :: Start
    INTEGER(KIND=8) :: ipart
    TYPE(Particle),POINTER,Save :: Cur

    IF (Start) Cur=>MainRoot%Head
    DO ipart=1,npart_this_it
       Cur%part_p(2) = Data(ipart)
       Cur=>Cur%Next
    ENDDO

  END SUBROUTINE It_Py

  SUBROUTINE It_Pz(Data,npart_this_it,Start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    INTEGER(KIND=8),INTENT(INOUT) :: npart_this_it
    LOGICAL,INTENT(IN) :: Start
    INTEGER(KIND=8) :: ipart
    TYPE(Particle),POINTER,Save :: Cur

    IF (Start) Cur=>MainRoot%Head
    DO ipart=1,npart_this_it
       Cur%part_p(3) = Data(ipart)
       Cur=>Cur%Next
    ENDDO

  END SUBROUTINE It_Pz

#ifdef PER_PARTICLE_WEIGHT
  SUBROUTINE It_Weight(Data,npart_this_it,Start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    INTEGER(KIND=8),INTENT(INOUT) :: npart_this_it
    LOGICAL,INTENT(IN) :: Start
    INTEGER(KIND=8) :: ipart
    TYPE(Particle),POINTER,Save :: Cur

    IF (Start) Cur=>MainRoot%Head
    DO ipart=1,npart_this_it
       Cur%weight = Data(ipart)
       Cur=>Cur%Next
    ENDDO

  END SUBROUTINE It_weight
#endif

  SUBROUTINE It_species(Data,npart_this_it,Start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    INTEGER(KIND=8),INTENT(INOUT) :: npart_this_it
    LOGICAL,INTENT(IN) :: Start
    INTEGER(KIND=8) :: ipart
    INTEGER(KIND=8), SAVE :: ipart_total

    IF (Start) ipart_total=1
    DO ipart=1,npart_this_it
       Species_ID(ipart_total) = NINT(Data(ipart))
       ipart_total=ipart_total+1
    ENDDO
  END SUBROUTINE It_species

END MODULE setup




