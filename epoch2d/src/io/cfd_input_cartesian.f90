MODULE cfd_input_cartesian

  USE cfd_common
  USE cfd_input_functions
  USE mpi

  IMPLICIT NONE

  SAVE

CONTAINS

  ! Grid loading functions

  SUBROUTINE cfd_get_nd_cartesian_grid_metadata_all(h, ndims, dims, extents)

    TYPE(cfd_file_handle) :: h
    INTEGER, INTENT(IN) :: ndims
    INTEGER, DIMENSION(:), INTENT(OUT) :: dims
    REAL(num), DIMENSION(:), INTENT(OUT) :: extents
    INTEGER(4), DIMENSION(:), ALLOCATABLE :: dims4
    INTEGER :: errcode

    ! This subroutine MUST be called after the call to
    ! get_common_mesh_metadata_all or it will break everything
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    IF (KIND(dims) .EQ. 4) THEN
      CALL MPI_FILE_READ_ALL(h%filehandle, dims, ndims, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      ALLOCATE(dims4(SIZE(dims)))

      CALL MPI_FILE_READ_ALL(h%filehandle, dims4, ndims, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      dims = dims4

      DEALLOCATE(dims4)
    ENDIF

    h%current_displacement = h%current_displacement + ndims * soi

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        h%mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, extents, ndims*2, h%mpireal, &
        MPI_STATUS_IGNORE, errcode)

    ! After this subroutine, all the metadata should be read in, so to make
    ! sure, just jump to known start of data

    CALL cfd_skip_block_metadata(h)

  END SUBROUTINE cfd_get_nd_cartesian_grid_metadata_all



  SUBROUTINE cfd_get_1d_cartesian_grid_all(h, x)

    TYPE(cfd_file_handle) :: h
    REAL(num), DIMENSION(:), INTENT(INOUT) :: x
    INTEGER :: errcode, nx

    nx = SIZE(x)

    CALL cfd_skip_block_metadata(h)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        h%mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, x, nx, h%mpireal, &
        MPI_STATUS_IGNORE, errcode)

    ! That should be it, so now skip to end of block
    CALL cfd_skip_block(h)

  END SUBROUTINE cfd_get_1d_cartesian_grid_all



  SUBROUTINE cfd_get_2d_cartesian_grid_all(h, x, y)

    TYPE(cfd_file_handle) :: h
    REAL(num), DIMENSION(:), INTENT(INOUT) :: x, y
    INTEGER :: errcode, nx, ny

    nx = SIZE(x)
    ny = SIZE(y)

    CALL cfd_skip_block_metadata(h)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        h%mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, x, nx, h%mpireal, &
        MPI_STATUS_IGNORE, errcode)
    CALL MPI_FILE_READ_ALL(h%filehandle, y, ny, h%mpireal, &
        MPI_STATUS_IGNORE, errcode)

    ! That should be it, so now skip to end of block
    CALL cfd_skip_block(h)

  END SUBROUTINE cfd_get_2d_cartesian_grid_all



  SUBROUTINE cfd_get_3d_cartesian_grid_all(h, x, y, z)

    TYPE(cfd_file_handle) :: h
    REAL(num), DIMENSION(:), INTENT(INOUT) :: x, y, z
    INTEGER :: errcode, nx, ny, nz

    nx = SIZE(x)
    ny = SIZE(y)
    nz = SIZE(z)

    CALL cfd_skip_block_metadata(h)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        h%mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, x, nx, h%mpireal, &
        MPI_STATUS_IGNORE, errcode)
    CALL MPI_FILE_READ_ALL(h%filehandle, y, ny, h%mpireal, &
        MPI_STATUS_IGNORE, errcode)
    CALL MPI_FILE_READ_ALL(h%filehandle, z, nz, h%mpireal, &
        MPI_STATUS_IGNORE, errcode)

    ! That should be it, so now skip to end of block
    CALL cfd_skip_block(h)

  END SUBROUTINE cfd_get_3d_cartesian_grid_all



  ! variable loading functions
  SUBROUTINE cfd_get_nd_cartesian_variable_metadata_all(h, ndims, dims, &
      extents, stagger, mesh_name, mesh_class)

    TYPE(cfd_file_handle) :: h
    INTEGER, INTENT(IN) :: ndims
    INTEGER, DIMENSION(:), INTENT(OUT) :: dims
    REAL(num), DIMENSION(:), INTENT(OUT) :: extents
    REAL(num), DIMENSION(:), INTENT(OUT) :: stagger
    CHARACTER(LEN=*), INTENT(INOUT) :: mesh_name, mesh_class
    CHARACTER(LEN=h%max_string_len) :: mesh_name_file, mesh_class_file
    INTEGER(4), DIMENSION(:), ALLOCATABLE :: dims4
    INTEGER :: len_name, len_class, errcode

    len_name = LEN(mesh_name)
    len_class = LEN(mesh_class)

    ! This subroutine MUST be called after the call to
    ! get_common_mesh_metadata_all or it will break everything
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    IF (KIND(dims) .EQ. 4) THEN
      CALL MPI_FILE_READ_ALL(h%filehandle, dims, ndims, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      ALLOCATE(dims4(SIZE(dims)))

      CALL MPI_FILE_READ_ALL(h%filehandle, dims4, ndims, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      dims = dims4

      DEALLOCATE(dims4)
    ENDIF

    h%current_displacement = h%current_displacement + ndims * soi

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        h%mpireal, "native", MPI_INFO_NULL, errcode)

    ! Read grid stagger
    CALL MPI_FILE_READ_ALL(h%filehandle, stagger, ndims, h%mpireal, &
        MPI_STATUS_IGNORE, errcode)

    ! Read data range
    CALL MPI_FILE_READ_ALL(h%filehandle, extents, 2, h%mpireal, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, mesh_name_file, h%max_string_len, &
        MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, mesh_class_file, h%max_string_len, &
        MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)

    mesh_name = mesh_name_file(1:MIN(h%max_string_len, len_name))
    mesh_class = mesh_class_file(1:MIN(h%max_string_len, len_class))

    ! After this subroutine, all the metadata should be read in, so to make
    ! sure, just jump to known start of data

    CALL cfd_skip_block_metadata(h)

  END SUBROUTINE cfd_get_nd_cartesian_variable_metadata_all



  SUBROUTINE cfd_get_1d_cartesian_variable_parallel(h, variable, subtype, &
      subarray)

    TYPE(cfd_file_handle) :: h
    REAL(num), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subtype, subarray
    INTEGER :: errcode

    CALL cfd_skip_block_metadata(h)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, subarray, &
        subtype, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block(h)

  END SUBROUTINE cfd_get_1d_cartesian_variable_parallel



  SUBROUTINE cfd_get_2d_cartesian_variable_parallel(h, variable, subtype, &
      subarray)

    TYPE(cfd_file_handle) :: h
    REAL(num), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subtype, subarray
    INTEGER :: errcode

    CALL cfd_skip_block_metadata(h)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, subarray, &
        subtype, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block(h)

  END SUBROUTINE cfd_get_2d_cartesian_variable_parallel



  SUBROUTINE cfd_get_3d_cartesian_variable_parallel(h, variable, subtype, &
      subarray)

    TYPE(cfd_file_handle) :: h
    REAL(num), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subtype, subarray
    INTEGER :: errcode

    CALL cfd_skip_block_metadata(h)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, subarray, &
        subtype, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block(h)

  END SUBROUTINE cfd_get_3d_cartesian_variable_parallel

END MODULE cfd_input_cartesian
