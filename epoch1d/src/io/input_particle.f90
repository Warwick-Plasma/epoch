MODULE input_particle

  USE shared_data
  USE iocommon
  USE input_functions

  IMPLICIT NONE

  SAVE

CONTAINS

  ! Grid loading functions
  SUBROUTINE cfd_get_nd_particle_grid_metadata_all(ndims, coord_type, npart, extents)

    INTEGER, INTENT(OUT) :: coord_type
    INTEGER(8), INTENT(OUT) :: npart
    REAL(num), DIMENSION(:), INTENT(OUT) :: extents
    INTEGER, INTENT(IN) :: ndims

    ! this subroutine MUST be called after the call to Get_Common_Mesh_MetaData_All or it will break everything
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, coord_type, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    current_displacement = current_displacement +  soi

    CALL MPI_FILE_READ_ALL(cfd_filehandle, npart, 1, MPI_INTEGER8, cfd_status, cfd_errcode)
    current_displacement = current_displacement + soi8

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, extents, ndims*2, mpireal, cfd_status, cfd_errcode)

    ! After this subroutine, all the metadata should be read in, so to make sure, just jump to known
    ! start of data

    CALL cfd_skip_block_metadata()

  END SUBROUTINE cfd_get_nd_particle_grid_metadata_all



  SUBROUTINE cfd_get_nd_particle_grid_all(ndims, npart, DATA)

    INTEGER, INTENT(IN) :: ndims
    INTEGER(KIND=8), INTENT(IN) :: npart
    REAL(num), DIMENSION(:, :), INTENT(INOUT) :: DATA

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, DATA, ndims*npart, mpireal, cfd_status, cfd_errcode)

    CALL cfd_skip_block

  END SUBROUTINE cfd_get_nd_particle_grid_all



  SUBROUTINE cfd_get_nd_particle_grid_parallel(ndims, npart, DATA, subtype)

    INTEGER, INTENT(IN) :: ndims
    INTEGER, INTENT(IN) :: subtype
    INTEGER(KIND=8), INTENT(IN) :: npart
    REAL(num), DIMENSION(:, :), INTENT(INOUT) :: DATA

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, subtype, "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, DATA, ndims*npart, mpireal, cfd_status, cfd_errcode)
    CALL cfd_skip_block

  END SUBROUTINE cfd_get_nd_particle_grid_parallel



  SUBROUTINE cfd_get_nd_particle_grid_parallel_with_iterator(ndims, npart_local, npart_lglobal, npart_per_it, sof, subtype, iterator)

    INTEGER, INTENT(IN) :: subtype
    INTEGER, INTENT(IN) :: ndims
    INTEGER, INTENT(IN) :: sof
    INTEGER(KIND=8), INTENT(IN) :: npart_local, npart_per_it, npart_lglobal
    INTEGER(KIND=8) :: npart_this_it, npart_remain
    INTEGER :: direction
    LOGICAL :: start
    REAL(num), DIMENSION(:), ALLOCATABLE :: DATA

    INTERFACE
      SUBROUTINE iterator(DATA, npart_it, start, direction)

        USE shared_data
        REAL(num), DIMENSION(:), INTENT(INOUT) :: DATA
        INTEGER(8), INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN) :: direction

      END SUBROUTINE iterator
    END INTERFACE

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, subtype, "native", MPI_INFO_NULL, cfd_errcode)

    ALLOCATE(DATA(1:npart_per_it))
    DO direction = 1, ndims
      start = .TRUE.
      npart_remain = npart_local
      npart_this_it = MIN(npart_remain, npart_per_it)
      CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, subtype, "native", MPI_INFO_NULL, cfd_errcode)
      DO WHILE (npart_this_it .GT. 0)
        CALL MPI_FILE_READ(cfd_filehandle, DATA, npart_this_it, mpireal, cfd_status, cfd_errcode)
        npart_remain = npart_remain-npart_this_it
        CALL iterator(DATA, npart_this_it, start, direction)
        start = .FALSE.
        npart_this_it = MIN(npart_remain, npart_per_it)
      ENDDO
      current_displacement = current_displacement+npart_lglobal*sof
    ENDDO
    DEALLOCATE(DATA)
    CALL MPI_BARRIER(cfd_comm, cfd_errcode)
    CALL cfd_skip_block

  END SUBROUTINE cfd_get_nd_particle_grid_parallel_with_iterator



  ! Grid loading functions
  SUBROUTINE cfd_get_nd_particle_variable_metadata_all(npart, range, mesh, mesh_class)

    INTEGER(8), INTENT(OUT) :: npart
    REAL(num), DIMENSION(2), INTENT(OUT) :: range
    CHARACTER(LEN=max_string_len), INTENT(OUT) :: mesh, mesh_class
    ! this subroutine MUST be called after the call to Get_Common_Mesh_MetaData_All or it will break everything

    CALL MPI_FILE_READ_ALL(cfd_filehandle, npart, 1, MPI_INTEGER8, cfd_status, cfd_errcode)
    current_displacement = current_displacement + soi8

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, range, 2, mpireal, cfd_status, cfd_errcode)
    current_displacement = current_displacement+2*num

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, mesh, max_string_len, MPI_CHARACTER, cfd_status, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, mesh_class, max_string_len, MPI_CHARACTER, cfd_status, cfd_errcode)

    ! After this subroutine, all the metadata should be read in, so to make sure, just jump to known
    ! start of data

    CALL cfd_skip_block_metadata()

  END SUBROUTINE cfd_get_nd_particle_variable_metadata_all



  SUBROUTINE cfd_get_nd_particle_variable_all(npart, DATA)

    INTEGER(KIND=8), INTENT(IN) :: npart
    REAL(num), DIMENSION(:), INTENT(INOUT) :: DATA

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, DATA, npart, mpireal, cfd_status, cfd_errcode)

    CALL cfd_skip_block

  END SUBROUTINE cfd_get_nd_particle_variable_all



  SUBROUTINE cfd_get_nd_particle_variable_parallel(npart_local, DATA, subtype)

    INTEGER, INTENT(IN) :: subtype
    INTEGER(KIND=8), INTENT(IN) :: npart_local
    REAL(num), DIMENSION(:, :), INTENT(INOUT) :: DATA

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, subtype, "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, DATA, npart_local, mpireal, cfd_status, cfd_errcode)

    CALL cfd_skip_block

  END SUBROUTINE cfd_get_nd_particle_variable_parallel



  SUBROUTINE cfd_get_nd_particle_variable_parallel_with_iterator(npart_local, npart_per_it, subtype, iterator)

    INTEGER, INTENT(IN) :: subtype
    INTEGER(KIND=8), INTENT(IN) :: npart_local, npart_per_it
    INTEGER(KIND=8) :: npart_this_it, npart_remain
    LOGICAL :: start
    REAL(num), DIMENSION(:), ALLOCATABLE :: DATA

    INTERFACE
      SUBROUTINE iterator(DATA, npart_it, start)
        USE shared_data
        REAL(num), DIMENSION(:), INTENT(INOUT) :: DATA
        INTEGER(8), INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
      END SUBROUTINE iterator
    END INTERFACE

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, subtype, "native", MPI_INFO_NULL, cfd_errcode)

    start = .TRUE.
    ALLOCATE(DATA(1:npart_per_it))
    npart_remain = npart_local
    npart_this_it = MIN(npart_remain, npart_per_it)
    DO WHILE (npart_this_it .GT. 0)
      npart_this_it = MIN(npart_remain, npart_per_it)
      CALL MPI_FILE_READ(cfd_filehandle, DATA, npart_this_it, mpireal, cfd_status, cfd_errcode)
      npart_remain = npart_remain-npart_this_it
      CALL iterator(DATA, npart_this_it, start)
      start = .FALSE.
    ENDDO
    CALL MPI_BARRIER(cfd_comm, cfd_errcode)
    DEALLOCATE(DATA)
    CALL cfd_skip_block

  END SUBROUTINE cfd_get_nd_particle_variable_parallel_with_iterator

END MODULE input_particle
