MODULE setup 
  USE shared_data
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
    DumpMask=.FALSE.

  END SUBROUTINE minimal_init

  SUBROUTINE after_control

    length_x=x_end-x_start
    dx = length_x / REAL(nx, num)

    !Setup grid
    x(0)=x_start-dx
    DO ix=1,nx
       x(ix)=x(ix-1)+dx
    ENDDO


    length_y=y_end-y_start
    dy = length_y / REAL(ny, num)

    !Setup grid
    y(0)=y_start-dy
    DO iy=1,ny
       y(iy)=y(iy-1)+dy
    ENDDO


    CALL set_initial_values

  END SUBROUTINE after_control

  SUBROUTINE open_files


    CHARACTER(LEN=11+Data_Dir_Max_Length) :: file2
    CHARACTER(LEN=7+Data_Dir_Max_Length) :: file3

    IF (rank == 0) THEN
       WRITE(file2, '(a,"/pic2d.dat")') TRIM(data_dir)
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

    Part_Pos=0.0_num
    Part_P=0.0_num
    Part_Species=0
    PisV=.FALSE.

    IF (nBlocks .EQ. 0) THEN
       IF (rank .EQ. 0) PRINT *,"This file is an empty or corrupted file"
       CALL MPI_ABORT(comm,errcode)
    ENDIF

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

          IF (nd .NE. DIMENSION_2D .AND. nd .NE. DIMENSION_IRRELEVANT ) THEN
             IF (rank .EQ. 0) PRINT *,"Dimensionality does not match, file is ",nd,"D"
             CALL MPI_ABORT(comm,errcode)
          ENDIF
          SELECT CASE(Type)
          CASE(VAR_CARTESIAN)
             !Grid variables
             CALL cfd_Get_nD_Cartesian_Variable_MetaData_All(nd,dims,extents,stagger,meshname,meshclass)
             IF (dims(1) .NE. nx) THEN
                IF (rank .EQ. 0) PRINT *,"Number of gridpoints does not match, gridpoints in file is",dims(1),"x",dims(2)
                CALL MPI_ABORT(comm,errcode)
             ENDIF
             IF (StrCmp(Name(1:2),"Ex")) CALL cfd_Get_2D_Cartesian_Variable_All(Ex(1:nx,1:ny))
             IF (StrCmp(Name(1:2),"Ey")) CALL cfd_Get_2D_Cartesian_Variable_All(Ey(1:nx,1:ny))
             IF (StrCmp(Name(1:2),"Ez")) CALL cfd_Get_2D_Cartesian_Variable_All(Ez(1:nx,1:ny))

             IF (StrCmp(Name(1:2),"Bx")) CALL cfd_Get_2D_Cartesian_Variable_All(Bx(1:nx,1:ny))
             IF (StrCmp(Name(1:2),"By")) CALL cfd_Get_2D_Cartesian_Variable_All(By(1:nx,1:ny))
             IF (StrCmp(Name(1:2),"Bz")) CALL cfd_Get_2D_Cartesian_Variable_All(Bz(1:nx,1:ny))

          CASE(VAR_PARTICLE)
             CALL cfd_Get_nD_Particle_Variable_MetaData_All(npart_l,extents,meshname,meshclass)
             IF (npart_l .NE. npart_global) THEN
                IF (rank .EQ. 0) PRINT *,"Number of particles does not match, particles in file is",npart
                CALL MPI_ABORT(comm,errcode)
             ENDIF
             !Particle variables
             IF (StrCmp(Name(1:3),"P_x")) CALL cfd_Get_nD_Particle_Variable_Parallel_With_Iterator(npart,npart_per_it,subtype,It_Px)
             IF (StrCmp(Name(1:3),"P_y")) CALL cfd_Get_nD_Particle_Variable_Parallel_With_Iterator(npart,npart_per_it,subtype,It_Py)
             IF (StrCmp(Name(1:3),"P_z")) CALL cfd_Get_nD_Particle_Variable_Parallel_With_Iterator(npart,npart_per_it,subtype,It_Pz)
             IF (StrCmp(Name(1:3),"V_x")) THEN
                CALL cfd_Get_nD_Particle_Variable_Parallel_With_Iterator(npart,npart_per_it,subtype,It_Px)
                PisV(1)=.TRUE.
             ENDIF
             IF (StrCmp(Name(1:3),"V_y")) THEN
                CALL cfd_Get_nD_Particle_Variable_Parallel_With_Iterator(npart,npart_per_it,subtype,It_Py)
                PisV(2)=.TRUE.
             ENDIF
             IF (StrCmp(Name(1:3),"V_z")) THEN
                CALL cfd_Get_nD_Particle_Variable_Parallel_With_Iterator(npart,npart_per_it,subtype,It_Pz)
                PisV(3)=.TRUE.
             ENDIF
          END SELECT
       CASE(TYPE_MESH)
          CALL cfd_Get_Common_MeshType_MetaData_All(Type,nd,sof)
          IF (Type .EQ. MESH_PARTICLE) THEN
             CALL cfd_Get_nD_Particle_Grid_MetaData_All(nd,coord_type,npart_l,extents)
             CALL cfd_Get_nD_Particle_Grid_Parallel(nd,npart,Part_Pos,subtype_particle_mesh)
          ENDIF
       CASE(TYPE_ARB_DB)
          CALL cfd_Get_Arb_Block(ReadSpecies)
       CASE(TYPE_CONSTANT)
          CALL cfd_Get_Real_Constant(Weight)
       END SELECT
       CALL cfd_Skip_Block()
    ENDDO
    CALL cfd_Close()

    DO ipart=1,npart
       v2=Part_p(ipart,1)**2+Part_p(ipart,2)**2+Part_p(ipart,3)**2
       IF (PisV(1)) THEN
          Part_p(ipart,1)=Part_p(ipart,1) * Species(part_species(ipart),2)/SQRT(1.0_num - v2/c**2)
       ENDIF
       IF (PisV(2)) THEN
          Part_p(ipart,2)=Part_p(ipart,2) * Species(part_species(ipart),2)/SQRT(1.0_num - v2/c**2)
       ENDIF
       IF (PisV(3)) THEN
          Part_p(ipart,3)=Part_p(ipart,3) * Species(part_species(ipart),2)/SQRT(1.0_num - v2/c**2)
       ENDIF
    ENDDO

  END SUBROUTINE restart_data

  SUBROUTINE It_Px(Data,npart_this_it,Start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    INTEGER(KIND=8),INTENT(INOUT) :: npart_this_it
    LOGICAL,INTENT(IN) :: Start
    INTEGER(KIND=8),SAVE :: ipart_start

    IF (Start) ipart_start=1
    Part_P(ipart_start:ipart_start+npart_this_it-1,1) = Data(1:npart_this_it)
    ipart_start=ipart_start+npart_this_it

  END SUBROUTINE It_Px

  SUBROUTINE It_Py(Data,npart_this_it,Start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    INTEGER(KIND=8),INTENT(INOUT) :: npart_this_it
    LOGICAL,INTENT(IN) :: Start
    INTEGER(KIND=8),SAVE :: ipart_start

    IF (Start) ipart_start=1
    Part_P(ipart_start:ipart_start+npart_this_it-1,2) = Data(1:npart_this_it)
    ipart_start=ipart_start+npart_this_it

  END SUBROUTINE It_Py

  SUBROUTINE It_Pz(Data,npart_this_it,Start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    INTEGER(KIND=8),INTENT(INOUT) :: npart_this_it
    LOGICAL,INTENT(IN) :: Start
    INTEGER(KIND=8),SAVE :: ipart_start

    IF (Start) ipart_start=1
    Part_P(ipart_start:ipart_start+npart_this_it-1,3) = Data(1:npart_this_it)
    ipart_start=ipart_start+npart_this_it

  END SUBROUTINE It_Pz

  SUBROUTINE ReadSpecies(filehandle,current_displacement,Generator_Name)

    INTEGER,INTENT(IN) :: filehandle
    INTEGER(KIND=MPI_OFFSET_KIND),INTENT(IN) :: current_displacement
    CHARACTER(LEN=*),INTENT(IN) :: Generator_Name

    CALL MPI_FILE_SET_VIEW(filehandle, current_displacement, MPI_INTEGER, subtype_int,&
         "native", MPI_INFO_NULL,cfd_errcode)
    CALL MPI_FILE_READ_ALL(filehandle, Part_Species,npart , MPI_INTEGER, status, errcode)

  END SUBROUTINE ReadSpecies

END MODULE setup




