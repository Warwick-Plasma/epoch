MODULE setup 
  USE shared_data
  USE balance
  USE input_cartesian
  USE input_particle
  USE input_arb
  USE input
  USE iocontrol
  USE inputfunctions
  USE Strings
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: after_control,minimal_init,restart_data
  PUBLIC :: open_files, close_files

CONTAINS

  SUBROUTINE minimal_init

    IF (num .EQ. 4) mpireal=MPI_REAL

  END SUBROUTINE minimal_init

  SUBROUTINE after_control

    INTEGER :: iproc

    length_x=x_end-x_start
    dx = length_x / REAL(nx_global, num)


    x_global=0.0_num
    x_global(0)=x_start-dx
    DO ix=1,nx_global+1
       x_global(ix)=x_global(ix-1)+dx
    ENDDO

    IF (domain == DO_DECOMPOSED) THEN
       DO iproc=0,nproc-1
          x_starts(iproc)=x_global(iproc*nx+1) !x_start+length_x*REAL(iproc,num)/REAL(nprocx,num)-dx
          x_ends(iproc)=x_global((iproc+1)*nx+1) !x_start+length_x*REAL(iproc+1,num)/REAL(nprocx,num)
       ENDDO
       x_start_local=x_starts(rank)
       x_end_local=x_ends(rank)
    ELSE
       x_starts=x_start
       x_ends=x_end
       x_start_local=x_start
       x_end_local=x_end
    ENDIF


    !Setup grid
    x(0)=x_start_local-dx
    DO ix=1,nx
       x(ix)=x(ix-1)+dx
    ENDDO

    CALL set_initial_values

  END SUBROUTINE after_control

  SUBROUTINE open_files


    CHARACTER(LEN=11+Data_Dir_Max_Length) :: file2
    CHARACTER(LEN=7+Data_Dir_Max_Length) :: file3

    IF (rank == 0) THEN
       WRITE(file2, '(a,"/pic1d.dat")') TRIM(data_dir)
       OPEN(unit=20, STATUS = 'REPLACE',FILE = file2)
       WRITE(file3, '(a,"/en.dat")') TRIM(data_dir)
       OPEN(unit=30, STATUS = 'REPLACE',FILE = file3,FORM="binary")
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
    INTEGER(KIND=8) :: npart_l,npart_per_it=10000
    REAL(num), DIMENSION(2) :: extents,stagger
    INTEGER,DIMENSION(1) :: dims
    REAL(KIND=8) :: time_d
    INTEGER :: snap,coord_type
    LOGICAL,DIMENSION(3) :: PisV
    REAL(num) :: v2
    TYPE(Particle), POINTER :: Current
    LOGICAL :: Const_Weight

    Const_Weight=.FALSE.

    CALL CreateSubTypes

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
    PisV=.FALSE.

    !    PRINT *,"Input file contains",nBlocks,"blocks"
    DO ix=1,nBlocks
       CALL cfd_Get_Next_Block_Info_All(Name,Class,Type)
       !       IF (rank .EQ. 0) PRINT *,ix,Name,Type
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

          IF (nd .NE. DIMENSION_1D .AND. nd .NE. DIMENSION_IRRELEVANT ) THEN
             IF (rank .EQ. 0) PRINT *,"Dimensionality does not match, file is ",nd,"D"
             CALL MPI_ABORT(comm,errcode)
          ENDIF
          SELECT CASE(Type)
          CASE(VAR_CARTESIAN)
             !Grid variables
             CALL cfd_Get_nD_Cartesian_Variable_MetaData_All(nd,dims,extents,stagger,meshname,meshclass)
             IF (dims(1) .NE. nx) THEN
                IF (rank .EQ. 0) PRINT *,"Number of gridpoints does not match, gridpoints in file is",dims(1)
                CALL MPI_ABORT(comm,errcode)
             ENDIF
             IF (StrCmp(Name(1:2),"Ex")) CALL cfd_Get_1D_Cartesian_Variable_Parallel(Ex(1:nx),Subtype_Field)
             IF (StrCmp(Name(1:2),"Ey")) CALL cfd_Get_1D_Cartesian_Variable_Parallel(Ey(1:nx),Subtype_Field)
             IF (StrCmp(Name(1:2),"Ez")) CALL cfd_Get_1D_Cartesian_Variable_Parallel(Ez(1:nx),SubType_Field)

             IF (StrCmp(Name(1:2),"Bx")) CALL cfd_Get_1D_Cartesian_Variable_Parallel(Bx(1:nx),Subtype_Field)
             IF (StrCmp(Name(1:2),"By")) CALL cfd_Get_1D_Cartesian_Variable_Parallel(By(1:nx),Subtype_Field)
             IF (StrCmp(Name(1:2),"Bz")) CALL cfd_Get_1D_Cartesian_Variable_Parallel(Bz(1:nx),Subtype_Field)

          CASE(VAR_PARTICLE)
             CALL cfd_Get_nD_Particle_Variable_MetaData_All(npart_l,extents,meshname,meshclass)
             IF (npart_l .NE. npart_global) THEN
                IF (rank .EQ. 0) PRINT *,"Number of particles does not match, particles in file is",npart
                CALL MPI_ABORT(comm,errcode)
             ENDIF
             !Particle variables
             IF (StrCmp(Name(1:2),"Px")) CALL cfd_Get_nD_Particle_Variable_Parallel_With_Iterator(npart_l/nproc,npart_per_it,subtype_particle,It_Px)
             IF (StrCmp(Name(1:2),"Py")) CALL cfd_Get_nD_Particle_Variable_Parallel_With_Iterator(npart_l/nproc,npart_per_it,subtype_particle,It_Py)
             IF (StrCmp(Name(1:2),"Pz")) CALL cfd_Get_nD_Particle_Variable_Parallel_With_Iterator(npart_l/nproc,npart_per_it,subtype_particle,It_Pz)
#ifdef PER_PARTICLE_WEIGHT
             IF (StrCmp(Name(1:6),"weight")) CALL cfd_Get_nd_Particle_Variable_Parallel_With_Iterator(npart_l/nproc,npart_per_it,subtype_particle,It_Weight)
#else
             IF (StrCmp(Name(1:6),"weight")) THEN
                IF (rank .EQ. 0) PRINT *,"Cannot load dump file with per particle weight if the code is compiled without per particle weights. Code terminates"
                CALL MPI_ABORT(comm,errcode)
             ENDIF
#endif

          END SELECT
       CASE(TYPE_MESH)
          CALL cfd_Get_Common_MeshType_MetaData_All(Type,nd,sof)
          IF (Type .EQ. MESH_PARTICLE) THEN
             CALL cfd_Get_nD_Particle_Grid_MetaData_All(nd,coord_type,npart_l,extents)
             CALL cfd_Get_nD_Particle_Grid_Parallel_With_Iterator(nd,npart_l/nproc,npart_per_it,subtype_particle,It_Part)
          ENDIF
       CASE(TYPE_ARB_DB)
          CALL cfd_Get_Arb_Block(ReadSpecies)
       CASE(TYPE_CONSTANT)
          CALL cfd_Get_Real_Constant(Weight)
          Const_Weight=.TRUE.
       END SELECT
       CALL cfd_Skip_Block()
    ENDDO
    CALL cfd_Close()


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

  SUBROUTINE It_Part(Data,npart_this_it,Start)

    REAL(num),DIMENSION(:,:),INTENT(INOUT) :: Data
    INTEGER(KIND=8),INTENT(INOUT) :: npart_this_it
    LOGICAL,INTENT(IN) :: Start
    INTEGER(KIND=8) :: ipart
    TYPE(Particle),POINTER,Save :: Cur

    IF (Start) Cur=>MainRoot%Head
    !Only 1D, so simple. In higher dimensions, data is packed in x,y,z triplets
    DO ipart=1,npart_this_it
       Cur%part_pos = Data(ipart,1)
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
       Cur%part_p(2) = Data(ipart)
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

  SUBROUTINE ReadSpecies(filehandle,current_displacement,Generator_Name)

    INTEGER,INTENT(IN) :: filehandle
    INTEGER(KIND=MPI_OFFSET_KIND),INTENT(IN) :: current_displacement
    CHARACTER(LEN=*),INTENT(IN) :: Generator_Name
    INTEGER, PARAMETER :: npart_per_it = 100000
    INTEGER,DIMENSION(npart_per_it) :: Data
    INTEGER(KIND=8):: npart_this_it, npart_left
    TYPE(Particle),POINTER,Save :: Cur
    INTEGER(KIND=8) :: ipart

    CALL MPI_FILE_SET_VIEW(filehandle, current_displacement, MPI_INTEGER, subtype_particle_int,&
         "native", MPI_INFO_NULL,cfd_errcode)

    npart_left=npart
    npart_this_it=MIN(npart_per_it,npart_left)
    DO WHILE(npart_this_it .NE. 0)
       npart_this_it=MIN(npart_per_it,npart_left)
       CALL MPI_FILE_READ_ALL(filehandle, Data, npart_this_it, MPI_INTEGER, status, errcode)
       npart_left=npart_left-npart_this_it
       DO ipart=1,npart_this_it
          Cur%part_species=Data(ipart)
          Cur=>Cur%Next
       ENDDO
    ENDDO

  END SUBROUTINE ReadSpecies

END MODULE setup




