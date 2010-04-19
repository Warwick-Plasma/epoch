MODULE input_particle

  USE input_functions

  IMPLICIT NONE

  SAVE

CONTAINS

  ! Grid loading functions
  SUBROUTINE cfd_get_nd_particle_grid_metadata_all(ndims, coord_type, npart, &
      extents)

    INTEGER, INTENT(IN) :: ndims
    INTEGER, INTENT(OUT) :: coord_type
    INTEGER(8), INTENT(OUT) :: npart
    REAL(num), DIMENSION(:), INTENT(OUT) :: extents
    INTEGER(4) :: coord_type4
    INTEGER :: errcode

    ! This subroutine MUST be called after the call to
    ! get_common_mesh_metadata_all or it will break everything
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, coord_type4, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    coord_type = coord_type4

    current_displacement = current_displacement +  soi

    CALL MPI_FILE_READ_ALL(cfd_filehandle, npart, 1, MPI_INTEGER8, &
        MPI_STATUS_IGNORE, errcode)

    current_displacement = current_displacement + soi8

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, extents, ndims * 2, mpireal, &
        MPI_STATUS_IGNORE, errcode)

    ! After this subroutine, all the metadata should be read in, so to make
    ! sure, just jump to known start of data

    CALL cfd_skip_block_metadata()

  END SUBROUTINE cfd_get_nd_particle_grid_metadata_all



  SUBROUTINE cfd_get_nd_particle_grid_all(ndims, npart, data)

    INTEGER, INTENT(IN) :: ndims
    INTEGER(8), INTENT(IN) :: npart
    REAL(num), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER :: errcode

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, data, ndims * npart, mpireal, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block

  END SUBROUTINE cfd_get_nd_particle_grid_all



  SUBROUTINE cfd_get_nd_particle_grid_parallel(ndims, npart, data, subtype)

    INTEGER, INTENT(IN) :: ndims
    INTEGER(8), INTENT(IN) :: npart
    REAL(num), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER, INTENT(IN) :: subtype
    INTEGER :: errcode

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        subtype, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, data, ndims * npart, mpireal, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block

  END SUBROUTINE cfd_get_nd_particle_grid_parallel



  SUBROUTINE cfd_get_nd_particle_grid_parallel_with_iterator(ndims, &
      npart_local, npart_lglobal, npart_per_it, sof, subtype, iterator)

    INTEGER, INTENT(IN) :: ndims
    INTEGER(8), INTENT(IN) :: npart_local, npart_lglobal, npart_per_it
    INTEGER, INTENT(IN) :: sof
    INTEGER, INTENT(IN) :: subtype
    INTEGER(8) :: npart_this_it, npart_remain
    INTEGER :: direction, errcode
    LOGICAL :: start
    REAL(num), DIMENSION(:), ALLOCATABLE :: data

    INTERFACE
      SUBROUTINE iterator(data, npart_it, start, direction)
        USE shared_data
        REAL(num), DIMENSION(:), INTENT(INOUT) :: data
        INTEGER(8), INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN) :: direction
      END SUBROUTINE iterator
    END INTERFACE

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        subtype, "native", MPI_INFO_NULL, errcode)

    ALLOCATE(data(1:npart_per_it))

    DO direction = 1, ndims
      start = .TRUE.
      npart_remain = npart_local
      npart_this_it = MIN(npart_remain, npart_per_it)

      CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
          subtype, "native", MPI_INFO_NULL, errcode)

      DO WHILE (npart_this_it .GT. 0)
        CALL MPI_FILE_READ(cfd_filehandle, data, npart_this_it, mpireal, &
            MPI_STATUS_IGNORE, errcode)

        npart_remain = npart_remain - npart_this_it
        CALL iterator(data, npart_this_it, start, direction)
        start = .FALSE.
        npart_this_it = MIN(npart_remain, npart_per_it)
      ENDDO

      current_displacement = current_displacement + npart_lglobal * sof
    ENDDO

    DEALLOCATE(data)

    CALL MPI_BARRIER(cfd_comm, errcode)
    CALL cfd_skip_block

  END SUBROUTINE cfd_get_nd_particle_grid_parallel_with_iterator



  ! Grid loading functions
  SUBROUTINE cfd_get_nd_particle_variable_metadata_all(npart, range, mesh, &
      mesh_class)

    INTEGER(8), INTENT(OUT) :: npart
    REAL(num), DIMENSION(2), INTENT(OUT) :: range
    CHARACTER(LEN=max_string_len), INTENT(OUT) :: mesh, mesh_class
    INTEGER :: errcode

    ! This subroutine MUST be called after the call to
    ! get_common_mesh_metadata_all or it will break everything

    CALL MPI_FILE_READ_ALL(cfd_filehandle, npart, 1, MPI_INTEGER8, &
        MPI_STATUS_IGNORE, errcode)

    current_displacement = current_displacement + soi8

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, range, 2, mpireal, &
        MPI_STATUS_IGNORE, errcode)

    current_displacement = current_displacement + 2 * num

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, mesh, max_string_len, &
        MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, mesh_class, max_string_len, &
        MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)

    ! After this subroutine, all the metadata should be read in, so to make
    ! sure, just jump to known start of data

    CALL cfd_skip_block_metadata()

  END SUBROUTINE cfd_get_nd_particle_variable_metadata_all



  SUBROUTINE cfd_get_nd_particle_variable_all(npart, data)

    INTEGER(8), INTENT(IN) :: npart
    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER :: errcode

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, data, npart, mpireal, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block

  END SUBROUTINE cfd_get_nd_particle_variable_all



  SUBROUTINE cfd_get_nd_particle_variable_parallel(npart_local, data, subtype)

    INTEGER(8), INTENT(IN) :: npart_local
    REAL(num), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER, INTENT(IN) :: subtype
    INTEGER :: errcode

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        subtype, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, data, npart_local, mpireal, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block

  END SUBROUTINE cfd_get_nd_particle_variable_parallel



  SUBROUTINE cfd_get_nd_particle_variable_parallel_with_iterator(npart_local, &
      npart_per_it, subtype, iterator)

    INTEGER(8), INTENT(IN) :: npart_local, npart_per_it
    INTEGER, INTENT(IN) :: subtype
    INTEGER(8) :: npart_this_it, npart_remain
    INTEGER :: errcode
    LOGICAL :: start
    REAL(num), DIMENSION(:), ALLOCATABLE :: data

    INTERFACE
      SUBROUTINE iterator(data, npart_it, start)
        USE shared_data
        REAL(num), DIMENSION(:), INTENT(INOUT) :: data
        INTEGER(8), INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
      END SUBROUTINE iterator
    END INTERFACE

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        subtype, "native", MPI_INFO_NULL, errcode)

    start = .TRUE.
    ALLOCATE(data(1:npart_per_it))
    npart_remain = npart_local
    npart_this_it = MIN(npart_remain, npart_per_it)

    DO WHILE (npart_this_it .GT. 0)
      npart_this_it = MIN(npart_remain, npart_per_it)
      CALL MPI_FILE_READ(cfd_filehandle, data, npart_this_it, mpireal, &
          MPI_STATUS_IGNORE, errcode)

      npart_remain = npart_remain - npart_this_it
      CALL iterator(data, npart_this_it, start)
      start = .FALSE.
    ENDDO

    CALL MPI_BARRIER(cfd_comm, errcode)
    DEALLOCATE(data)
    CALL cfd_skip_block

  END SUBROUTINE cfd_get_nd_particle_variable_parallel_with_iterator

END MODULE input_particle
