MODULE input_particle

  USE shared_data
  USE iocommon
  USE inputfunctions

  IMPLICIT NONE

  SAVE

CONTAINS


  !Grid loading functions
  SUBROUTINE cfd_Get_nD_Particle_Grid_MetaData_All(ndims,coord_type,npart,extents)

    INTEGER,INTENT(OUT) :: coord_type
    INTEGER(8),INTENT(OUT) :: npart
    REAL(num),DIMENSION(:),INTENT(OUT) :: extents
    INTEGER,INTENT(IN) :: ndims
    !this subroutine MUST be called after the call to Get_Common_Mesh_MetaData_All or it will break everything
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, coord_type, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    current_displacement = current_displacement +  SoI 

    CALL MPI_FILE_READ_ALL(cfd_filehandle,npart,1,MPI_INTEGER8,cfd_status,cfd_errcode)
    current_displacement = current_displacement + SoI8

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, extents, ndims*2, mpireal, cfd_status, cfd_errcode)

    !After this subroutine, all the metadata should be read in, so to make sure, just jump to known
    !Start of Data

    CALL cfd_Skip_Block_MetaData()

  END SUBROUTINE cfd_Get_nD_Particle_Grid_MetaData_All

  SUBROUTINE cfd_Get_nD_Particle_Grid_All(ndims,npart,Data)

    INTEGER, INTENT(IN) :: ndims
    INTEGER(KIND=8), INTENT(IN) :: npart
    REAL(num),DIMENSION(:,:),INTENT(INOUT) :: Data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, Data,ndims*npart,mpireal,cfd_status,cfd_errcode)

    CALL cfd_Skip_Block

  END SUBROUTINE cfd_Get_nD_Particle_Grid_All

  SUBROUTINE cfd_Get_nD_Particle_Grid_Parallel(ndims,npart,Data,SubType)

    INTEGER, INTENT(IN) :: ndims
    INTEGER, INTENT(IN) :: SubType
    INTEGER(KIND=8), INTENT(IN) :: npart
    REAL(num),DIMENSION(:,:),INTENT(INOUT) :: Data


    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, SubType,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, Data,ndims*npart,mpireal,cfd_status,cfd_errcode)
    CALL cfd_Skip_Block

  END SUBROUTINE cfd_Get_nD_Particle_Grid_Parallel

  SUBROUTINE cfd_Get_nD_Particle_Grid_Parallel_With_Iterator(ndims,npart_local,npart_lglobal,npart_per_it,sof,SubType,Iterator)


    INTEGER, INTENT(IN) :: SubType
    INTEGER, INTENT(IN) :: ndims
    INTEGER, INTENT(IN) :: sof
    INTEGER(KIND=8), INTENT(IN) :: npart_local,npart_per_it,npart_lglobal
    INTEGER(KIND=8) :: npart_this_it,npart_remain,size_of_stencil
    INTEGER :: Direction
    LOGICAL :: Start
    REAL(num),DIMENSION(:),ALLOCATABLE :: Data

    INTERFACE
       SUBROUTINE iterator(data,npart_it,start,direction)
         USE shared_data
         REAL(num),DIMENSION(:),INTENT(INOUT) :: data
         INTEGER(8),INTENT(INOUT) :: npart_it
         LOGICAL,INTENT(IN) :: start
         INTEGER,INTENT(IN) :: Direction
       END SUBROUTINE iterator
    END INTERFACE

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, SubType,&
         "native", MPI_INFO_NULL, cfd_errcode)

    ALLOCATE(Data(1:npart_per_it))
    DO Direction=1,ndims
       Start=.TRUE.
       npart_remain=npart_local
       npart_this_it=MIN(npart_remain,npart_per_it)
       CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, SubType,&
            "native", MPI_INFO_NULL, cfd_errcode)
       DO WHILE (npart_this_it .GT. 0)
          CALL MPI_FILE_READ(cfd_filehandle, Data,npart_this_it,mpireal,cfd_status,cfd_errcode)
          npart_remain=npart_remain-npart_this_it
          CALL Iterator(Data,npart_this_it,Start,direction)
          Start=.FALSE.
          npart_this_it=MIN(npart_remain,npart_per_it)
       ENDDO
       current_displacement=current_displacement+npart_lglobal*sof
    ENDDO
    DEALLOCATE(Data)
    CALL MPI_BARRIER(cfd_comm,cfd_errcode)
    CALL cfd_Skip_Block

  END SUBROUTINE cfd_Get_nD_Particle_Grid_Parallel_With_Iterator

  !Grid loading functions
  SUBROUTINE cfd_Get_nD_Particle_Variable_MetaData_All(npart,range,mesh,mesh_class)

    INTEGER(8),INTENT(OUT) :: npart
    REAL(num),DIMENSION(2),INTENT(OUT) :: range
    CHARACTER(Len=MaxStringLen), INTENT(OUT) :: mesh,mesh_class
    !this subroutine MUST be called after the call to Get_Common_Mesh_MetaData_All or it will break everything

    CALL MPI_FILE_READ_ALL(cfd_filehandle,npart,1,MPI_INTEGER8,cfd_status,cfd_errcode)
    current_displacement = current_displacement + SoI8

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, range, 2, mpireal, cfd_status, cfd_errcode)
    current_displacement=current_displacement+2*num

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, mesh, MaxStringLen, MPI_CHARACTER, cfd_status, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, mesh_class, MaxStringLen, MPI_CHARACTER, cfd_status, cfd_errcode)

    !After this subroutine, all the metadata should be read in, so to make sure, just jump to known
    !Start of Data

    CALL cfd_Skip_Block_MetaData()

  END SUBROUTINE cfd_Get_nD_Particle_Variable_MetaData_All

  SUBROUTINE cfd_Get_nD_Particle_Variable_All(npart,Data)

    INTEGER(KIND=8), INTENT(IN) :: npart
    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, Data,npart,mpireal,cfd_status,cfd_errcode)

    CALL cfd_Skip_Block

  END SUBROUTINE cfd_Get_nD_Particle_Variable_All

  SUBROUTINE cfd_Get_nD_Particle_Variable_Parallel(npart_local,Data,SubType)

    INTEGER, INTENT(IN) :: SubType
    INTEGER(KIND=8), INTENT(IN) :: npart_local
    REAL(num),DIMENSION(:,:),INTENT(INOUT) :: Data

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, SubType,&
         "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, Data,npart_local,mpireal,cfd_status,cfd_errcode)

    CALL cfd_Skip_Block

  END SUBROUTINE cfd_Get_nD_Particle_Variable_Parallel

  SUBROUTINE cfd_Get_nD_Particle_Variable_Parallel_With_Iterator(npart_local,npart_per_it,SubType,Iterator)

    INTEGER, INTENT(IN) :: SubType
    INTEGER(KIND=8), INTENT(IN) :: npart_local,npart_per_it
    INTEGER(KIND=8) :: npart_this_it,npart_remain
    LOGICAL :: Start
    REAL(num),DIMENSION(:),ALLOCATABLE :: Data

    INTERFACE
       SUBROUTINE iterator(data,npart_it,start)
         USE shared_data
         REAL(num),DIMENSION(:),INTENT(INOUT) :: data
         INTEGER(8),INTENT(INOUT) :: npart_it
         LOGICAL,INTENT(IN) :: start
       END SUBROUTINE iterator
    END INTERFACE


    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, SubType,&
         "native", MPI_INFO_NULL, cfd_errcode)

    Start=.TRUE.
    ALLOCATE(Data(1:npart_per_it))
    npart_remain=npart_local
    npart_this_it=MIN(npart_remain,npart_per_it)
    DO WHILE (npart_this_it .GT. 0)
       npart_this_it=MIN(npart_remain,npart_per_it)
       CALL MPI_FILE_READ(cfd_filehandle, Data,npart_this_it,mpireal,cfd_status,cfd_errcode)
       npart_remain=npart_remain-npart_this_it
       CALL Iterator(Data,npart_this_it,Start)
       Start=.FALSE.
    ENDDO
    CALL MPI_BARRIER(cfd_comm,cfd_errcode)
    DEALLOCATE(Data)
    CALL cfd_Skip_Block

  END SUBROUTINE cfd_Get_nD_Particle_Variable_Parallel_With_Iterator

END MODULE input_particle
