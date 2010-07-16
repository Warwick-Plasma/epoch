MODULE cfd_input_particle

  USE cfd_common
  USE cfd_input_functions
  USE mpi

  IMPLICIT NONE

  SAVE

CONTAINS

  ! Grid loading functions
  SUBROUTINE cfd_get_nd_particle_grid_metadata_all(h, ndims, coord_type, &
      npart_global, extents)

    TYPE(cfd_file_handle) :: h
    INTEGER, INTENT(IN) :: ndims
    INTEGER, INTENT(OUT) :: coord_type
    INTEGER(8), INTENT(OUT) :: npart_global
    REAL(num), DIMENSION(:), INTENT(OUT) :: extents
    INTEGER(4) :: coord_type4
    INTEGER :: errcode

    ! This subroutine MUST be called after the call to
    ! get_common_mesh_metadata_all or it will break everything
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, coord_type4, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    coord_type = coord_type4

    h%current_displacement = h%current_displacement +  soi

    CALL MPI_FILE_READ_ALL(h%filehandle, npart_global, 1, MPI_INTEGER8, &
        MPI_STATUS_IGNORE, errcode)

    h%current_displacement = h%current_displacement + soi8

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        h%mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, extents, ndims * 2, h%mpireal, &
        MPI_STATUS_IGNORE, errcode)

    ! After this subroutine, all the metadata should be read in, so to make
    ! sure, just jump to known start of data

    CALL cfd_skip_block_metadata(h)

  END SUBROUTINE cfd_get_nd_particle_grid_metadata_all



  SUBROUTINE cfd_get_nd_particle_grid_parallel(h, ndims, npart, array, subtype)

    TYPE(cfd_file_handle) :: h
    INTEGER, INTENT(IN) :: ndims
    INTEGER(8), INTENT(IN) :: npart
    REAL(num), DIMENSION(:,:), INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: subtype
    INTEGER :: errcode

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        subtype, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, array, ndims * npart, h%mpireal, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block(h)

  END SUBROUTINE cfd_get_nd_particle_grid_parallel



  SUBROUTINE cfd_get_nd_particle_grid_parallel_with_iterator(h, ndims, &
      npart_local, npart_global, npart_per_it, sof, subtype, iterator)

    TYPE(cfd_file_handle) :: h
    INTEGER, INTENT(IN) :: ndims
    INTEGER(8), INTENT(IN) :: npart_local, npart_global, npart_per_it
    INTEGER, INTENT(IN) :: sof
    INTEGER, INTENT(IN) :: subtype
    INTEGER(8) :: npart_this_it, npart_remain
    INTEGER :: direction, errcode
    LOGICAL :: start
    REAL(num), DIMENSION(:), ALLOCATABLE :: array

    INTERFACE
      SUBROUTINE iterator(array, npart_it, start, direction)
        USE cfd_common
        REAL(num), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER(8), INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN) :: direction
      END SUBROUTINE iterator
    END INTERFACE

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        subtype, "native", MPI_INFO_NULL, errcode)

    ALLOCATE(array(1:npart_per_it))

    DO direction = 1, ndims
      start = .TRUE.
      npart_remain = npart_local
      npart_this_it = MIN(npart_remain, npart_per_it)

      CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
          h%mpireal, subtype, "native", MPI_INFO_NULL, errcode)

      DO WHILE (npart_this_it .GT. 0)
        CALL MPI_FILE_READ(h%filehandle, array, npart_this_it, h%mpireal, &
            MPI_STATUS_IGNORE, errcode)

        npart_remain = npart_remain - npart_this_it
        CALL iterator(array, npart_this_it, start, direction)
        start = .FALSE.
        npart_this_it = MIN(npart_remain, npart_per_it)
      ENDDO

      h%current_displacement = h%current_displacement + npart_global * sof
    ENDDO

    DEALLOCATE(array)

    CALL cfd_skip_block(h)

  END SUBROUTINE cfd_get_nd_particle_grid_parallel_with_iterator



  ! Grid loading functions
  SUBROUTINE cfd_get_nd_particle_variable_metadata_all(h, npart_global, &
      range, mesh, mesh_class)

    TYPE(cfd_file_handle) :: h
    INTEGER(8), INTENT(OUT) :: npart_global
    REAL(num), DIMENSION(2), INTENT(OUT) :: range
    CHARACTER(LEN=h%max_string_len), INTENT(OUT) :: mesh, mesh_class
    INTEGER :: errcode

    ! This subroutine MUST be called after the call to
    ! get_common_mesh_metadata_all or it will break everything

    CALL MPI_FILE_READ_ALL(h%filehandle, npart_global, 1, MPI_INTEGER8, &
        MPI_STATUS_IGNORE, errcode)

    h%current_displacement = h%current_displacement + soi8

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        h%mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, range, 2, h%mpireal, &
        MPI_STATUS_IGNORE, errcode)

    h%current_displacement = h%current_displacement + 2 * num

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, mesh, h%max_string_len, &
        MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, mesh_class, h%max_string_len, &
        MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)

    ! After this subroutine, all the metadata should be read in, so to make
    ! sure, just jump to known start of data

    CALL cfd_skip_block_metadata(h)

  END SUBROUTINE cfd_get_nd_particle_variable_metadata_all



  SUBROUTINE cfd_get_nd_particle_variable_parallel(h, npart_local, array, &
      subtype)

    TYPE(cfd_file_handle) :: h
    INTEGER(8), INTENT(IN) :: npart_local
    REAL(num), DIMENSION(:,:), INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: subtype
    INTEGER :: errcode

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        subtype, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, array, npart_local, h%mpireal, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block(h)

  END SUBROUTINE cfd_get_nd_particle_variable_parallel



  SUBROUTINE cfd_get_nd_particle_variable_parallel_with_iterator(h, &
      npart_local, npart_per_it, subtype, iterator)

    TYPE(cfd_file_handle) :: h
    INTEGER(8), INTENT(IN) :: npart_local, npart_per_it
    INTEGER, INTENT(IN) :: subtype
    INTEGER(8) :: npart_this_it, npart_remain
    INTEGER :: errcode
    LOGICAL :: start
    REAL(num), DIMENSION(:), ALLOCATABLE :: array

    INTERFACE
      SUBROUTINE iterator(array, npart_it, start)
        USE cfd_common
        REAL(num), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER(8), INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
      END SUBROUTINE iterator
    END INTERFACE

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        subtype, "native", MPI_INFO_NULL, errcode)

    start = .TRUE.
    ALLOCATE(array(1:npart_per_it))
    npart_remain = npart_local
    npart_this_it = MIN(npart_remain, npart_per_it)

    DO WHILE (npart_this_it .GT. 0)
      npart_this_it = MIN(npart_remain, npart_per_it)
      CALL MPI_FILE_READ(h%filehandle, array, npart_this_it, h%mpireal, &
          MPI_STATUS_IGNORE, errcode)

      npart_remain = npart_remain - npart_this_it
      CALL iterator(array, npart_this_it, start)
      start = .FALSE.
    ENDDO

    DEALLOCATE(array)
    CALL cfd_skip_block(h)

  END SUBROUTINE cfd_get_nd_particle_variable_parallel_with_iterator

END MODULE cfd_input_particle
