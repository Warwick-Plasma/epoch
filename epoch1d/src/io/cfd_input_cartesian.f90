MODULE input_cartesian

  USE input_functions

  IMPLICIT NONE

  SAVE

CONTAINS

  ! Grid loading functions

  SUBROUTINE cfd_get_nd_cartesian_grid_metadata_all(ndims, dims, extents)

    INTEGER, INTENT(IN) :: ndims
    INTEGER, DIMENSION(:), INTENT(OUT) :: dims
    REAL(num), DIMENSION(:), INTENT(OUT) :: extents
    INTEGER(4), DIMENSION(:), ALLOCATABLE :: dims4
    INTEGER :: errcode

    ! This subroutine MUST be called after the call to
    ! get_common_mesh_metadata_all or it will break everything
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    IF (KIND(dims) .EQ. 4) THEN
      CALL MPI_FILE_READ_ALL(cfd_filehandle, dims, ndims, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      ALLOCATE(dims4(SIZE(dims)))

      CALL MPI_FILE_READ_ALL(cfd_filehandle, dims4, ndims, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      dims = dims4

      DEALLOCATE(dims4)
    ENDIF

    current_displacement = current_displacement + ndims * soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, extents, ndims*2, mpireal, &
        MPI_STATUS_IGNORE, errcode)

    ! After this subroutine, all the metadata should be read in, so to make
    ! sure, just jump to known start of data

    CALL cfd_skip_block_metadata()

  END SUBROUTINE cfd_get_nd_cartesian_grid_metadata_all



  SUBROUTINE cfd_get_1d_cartesian_grid_all(x)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: x
    INTEGER :: errcode, nx

    nx = SIZE(x)

    CALL cfd_skip_block_metadata()

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, x, nx, mpireal, MPI_STATUS_IGNORE, &
        errcode)

    ! That should be it, so now skip to end of block
    CALL cfd_skip_block

  END SUBROUTINE cfd_get_1d_cartesian_grid_all



  SUBROUTINE cfd_get_2d_cartesian_grid_all(x, y)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: x, y
    INTEGER :: errcode, nx, ny

    nx = SIZE(x)
    ny = SIZE(y)

    CALL cfd_skip_block_metadata()

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, x, nx, mpireal, MPI_STATUS_IGNORE, &
        errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, y, ny, mpireal, MPI_STATUS_IGNORE, &
        errcode)

    ! That should be it, so now skip to end of block
    CALL cfd_skip_block

  END SUBROUTINE cfd_get_2d_cartesian_grid_all



  SUBROUTINE cfd_get_3d_cartesian_grid_all(x, y, z)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: x, y, z
    INTEGER :: errcode, nx, ny, nz

    nx = SIZE(x)
    ny = SIZE(y)
    nz = SIZE(z)

    CALL cfd_skip_block_metadata()

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, x, nx, mpireal, MPI_STATUS_IGNORE, &
        errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, y, ny, mpireal, MPI_STATUS_IGNORE, &
        errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, z, nz, mpireal, MPI_STATUS_IGNORE, &
        errcode)

    ! That should be it, so now skip to end of block
    CALL cfd_skip_block

  END SUBROUTINE cfd_get_3d_cartesian_grid_all



  ! variable loading functions
  SUBROUTINE cfd_get_nd_cartesian_variable_metadata_all(ndims, dims, &
      extents, stagger, mesh_name, mesh_class)

    INTEGER, INTENT(IN) :: ndims
    INTEGER, DIMENSION(:), INTENT(OUT) :: dims
    REAL(num), DIMENSION(:), INTENT(OUT) :: extents
    REAL(num), DIMENSION(:), INTENT(OUT) :: stagger
    CHARACTER(LEN=*), INTENT(INOUT) :: mesh_name, mesh_class
    CHARACTER(LEN=max_string_len) :: mesh_name_file, mesh_class_file
    INTEGER(4), DIMENSION(:), ALLOCATABLE :: dims4
    INTEGER :: len_name, len_class, errcode

    len_name = LEN(mesh_name)
    len_class = LEN(mesh_class)

    ! This subroutine MUST be called after the call to
    ! get_common_mesh_metadata_all or it will break everything
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    IF (KIND(dims) .EQ. 4) THEN
      CALL MPI_FILE_READ_ALL(cfd_filehandle, dims, ndims, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      ALLOCATE(dims4(SIZE(dims)))

      CALL MPI_FILE_READ_ALL(cfd_filehandle, dims4, ndims, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      dims = dims4

      DEALLOCATE(dims4)
    ENDIF

    current_displacement = current_displacement + ndims * soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, errcode)

    ! Read grid stagger
    CALL MPI_FILE_READ_ALL(cfd_filehandle, stagger, ndims, mpireal, &
        MPI_STATUS_IGNORE, errcode)

    ! Read data range
    CALL MPI_FILE_READ_ALL(cfd_filehandle, extents, 2, mpireal, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, mesh_name_file, max_string_len, &
        MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, mesh_class_file, max_string_len, &
        MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)

    mesh_name = mesh_name_file(1:MIN(max_string_len, len_name))
    mesh_class = mesh_class_file(1:MIN(max_string_len, len_class))

    ! After this subroutine, all the metadata should be read in, so to make
    ! sure, just jump to known start of data

    CALL cfd_skip_block_metadata()

  END SUBROUTINE cfd_get_nd_cartesian_variable_metadata_all



  SUBROUTINE cfd_get_1d_cartesian_variable_parallel(variable, subtype, subarray)

    REAL(num), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subtype, subarray
    INTEGER :: errcode

    CALL cfd_skip_block_metadata()

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, subarray, &
        subtype, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_1d_cartesian_variable_parallel



  SUBROUTINE cfd_get_1d_cartesian_variable_all(variable)

    REAL(num), DIMENSION(:), INTENT(IN) :: variable
    INTEGER :: errcode

    CALL cfd_skip_block_metadata()

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_1d_cartesian_variable_all



  SUBROUTINE cfd_get_2d_cartesian_variable_parallel(variable, subtype, subarray)

    REAL(num), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subtype, subarray
    INTEGER :: errcode

    CALL cfd_skip_block_metadata()

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, subarray, &
        subtype, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_2d_cartesian_variable_parallel



  SUBROUTINE cfd_get_2d_cartesian_variable_all(variable)

    REAL(num), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER :: errcode

    CALL cfd_skip_block_metadata()

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_2d_cartesian_variable_all



  SUBROUTINE cfd_get_3d_cartesian_variable_parallel(variable, subtype, subarray)

    REAL(num), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subtype, subarray
    INTEGER :: errcode

    CALL cfd_skip_block_metadata()

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, subarray, &
        subtype, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_3d_cartesian_variable_parallel



  SUBROUTINE cfd_get_3d_cartesian_variable_all(variable)

    REAL(num), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER :: errcode

    CALL cfd_skip_block_metadata()

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_3d_cartesian_variable_all

END MODULE input_cartesian
