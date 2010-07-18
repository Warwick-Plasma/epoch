MODULE output_cartesian

  USE iocommon
  USE output

  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian grid in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify nx
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_1d_cartesian_grid(name, class, x, rank_write)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), DIMENSION(:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: rank_write
    INTEGER, PARAMETER :: ndims = 1
    INTEGER(4), DIMENSION(ndims) :: dims4
    INTEGER(8) :: block_length, md_length
    INTEGER :: i, len_var

    dims4(1) = INT(SIZE(x),4)

    len_var = 0
    DO i = 1,ndims
      len_var = len_var + dims4(i)
    ENDDO

    ! Metadata is
    ! - meshtype (INTEGER(4)) All mesh blocks contain this
    ! - nd    INTEGER(4)
    ! - sof   INTEGER(4)
    ! Specific to cartesian grid
    ! - nx    INTEGER(4)
    ! - xmin  REAL(num)
    ! - xmax  REAL(num)

    md_length = meshtype_header_offset + ndims * soi + 2 * ndims * sof
    block_length = md_length + len_var * sof

    ! Write header
    CALL cfd_write_block_header(name, class, c_type_mesh, &
        block_length, md_length, rank_write)
    CALL cfd_write_meshtype_header(c_mesh_cartesian, c_dimension_1d, sof, &
        rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      DO i = 1,ndims
        CALL MPI_FILE_WRITE(cfd_filehandle, dims4(i), 1, MPI_INTEGER4, &
            cfd_status, cfd_errcode)
      ENDDO
    ENDIF

    current_displacement = current_displacement + ndims * soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, x(1), 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, x(dims4(1)), 1, mpireal, &
          cfd_status, cfd_errcode)

      ! Now write the real arrays
      CALL MPI_FILE_WRITE(cfd_filehandle, x, dims4(1), mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + (2 * ndims + len_var) * sof

  END SUBROUTINE cfd_write_1d_cartesian_grid



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian grid in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify nx, ny
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_2d_cartesian_grid(name, class, x, y, rank_write)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), DIMENSION(:), INTENT(IN) :: x, y
    INTEGER, INTENT(IN) :: rank_write
    INTEGER, PARAMETER :: ndims = 2
    INTEGER(4), DIMENSION(ndims) :: dims4
    INTEGER(8) :: block_length, md_length
    INTEGER :: i, len_var

    dims4(1) = INT(SIZE(x),4)
    dims4(2) = INT(SIZE(y),4)

    len_var = 0
    DO i = 1,ndims
      len_var = len_var + dims4(i)
    ENDDO

    ! Metadata is
    ! - meshtype (INTEGER(4)) All mesh blocks contain this
    ! - nd    INTEGER(4)
    ! - sof   INTEGER(4)
    ! Specific to cartesian grid
    ! - nx    INTEGER(4)
    ! - ny    INTEGER(4)
    ! - xmin  REAL(num)
    ! - xmax  REAL(num)
    ! - ymin  REAL(num)
    ! - ymax  REAL(num)

    md_length = meshtype_header_offset + ndims * soi + 2 * ndims * sof
    block_length = md_length + len_var * sof

    ! Write header
    CALL cfd_write_block_header(name, class, c_type_mesh, &
        block_length, md_length, rank_write)
    CALL cfd_write_meshtype_header(c_mesh_cartesian, c_dimension_2d, sof, &
        rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      DO i = 1,ndims
        CALL MPI_FILE_WRITE(cfd_filehandle, dims4(i), 1, MPI_INTEGER4, &
            cfd_status, cfd_errcode)
      ENDDO
    ENDIF

    current_displacement = current_displacement + ndims * soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, x(1), 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, x(dims4(1)), 1, mpireal, &
          cfd_status, cfd_errcode)

      CALL MPI_FILE_WRITE(cfd_filehandle, y(1), 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, y(dims4(2)), 1, mpireal, &
          cfd_status, cfd_errcode)

      ! Now write the real arrays
      CALL MPI_FILE_WRITE(cfd_filehandle, x, dims4(1), mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, y, dims4(2), mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + (2 * ndims + len_var) * sof

  END SUBROUTINE cfd_write_2d_cartesian_grid



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian grid in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify nx, ny, nz
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_3d_cartesian_grid(name, class, x, y, z, rank_write)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), DIMENSION(:), INTENT(IN) :: x, y, z
    INTEGER, INTENT(IN) :: rank_write
    INTEGER, PARAMETER :: ndims = 3
    INTEGER(4), DIMENSION(ndims) :: dims4
    INTEGER(8) :: block_length, md_length
    INTEGER :: i, len_var

    dims4(1) = INT(SIZE(x),4)
    dims4(2) = INT(SIZE(y),4)
    dims4(3) = INT(SIZE(z),4)

    len_var = 0
    DO i = 1,ndims
      len_var = len_var + dims4(i)
    ENDDO

    ! Metadata is
    ! - meshtype (INTEGER(4)) All mesh blocks contain this
    ! - nd    INTEGER(4)
    ! - sof   INTEGER(4)
    ! Specific to cartesian grid
    ! - nx    INTEGER(4)
    ! - ny    INTEGER(4)
    ! - nz    INTEGER(4)
    ! - xmin  REAL(num)
    ! - xmax  REAL(num)
    ! - ymin  REAL(num)
    ! - ymax  REAL(num)
    ! - zmin  REAL(num)
    ! - zmax  REAL(num)

    md_length = meshtype_header_offset + ndims * soi + 2 * ndims * sof
    block_length = md_length + len_var * sof

    ! Write header
    CALL cfd_write_block_header(name, class, c_type_mesh, &
        block_length, md_length, rank_write)
    CALL cfd_write_meshtype_header(c_mesh_cartesian, c_dimension_3d, sof, &
        rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      DO i = 1,ndims
        CALL MPI_FILE_WRITE(cfd_filehandle, dims4(i), 1, MPI_INTEGER4, &
            cfd_status, cfd_errcode)
      ENDDO
    ENDIF

    current_displacement = current_displacement + ndims * soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, x(1), 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, x(dims4(1)), 1, mpireal, &
          cfd_status, cfd_errcode)

      CALL MPI_FILE_WRITE(cfd_filehandle, y(1), 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, y(dims4(2)), 1, mpireal, &
          cfd_status, cfd_errcode)

      CALL MPI_FILE_WRITE(cfd_filehandle, z(1), 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, z(dims4(3)), 1, mpireal, &
          cfd_status, cfd_errcode)

      ! Now write the real arrays
      CALL MPI_FILE_WRITE(cfd_filehandle, x, dims4(1), mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, y, dims4(2), mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, z, dims4(3), mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + (2 * ndims + len_var) * sof

  END SUBROUTINE cfd_write_3d_cartesian_grid



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian variable in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global nx
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_1d_cartesian_variable_parallel(name, class, dims, &
      stagger, mesh_name, mesh_class, variable, distribution, subarray)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    REAL(num), DIMENSION(:), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_name, mesh_class
    REAL(num), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num) :: mn, mx, mn_global, mx_global
    INTEGER(8) :: block_length, md_length, len_var
    INTEGER, PARAMETER :: ndims = 1
    INTEGER(4), DIMENSION(ndims) :: dims4
    INTEGER :: i

    len_var = 1
    DO i = 1,ndims
      dims4(i) = INT(dims(i),4)
      len_var = len_var * dims(i)
    ENDDO

    ! Metadata is
    ! - meshtype (INTEGER(4)) All mesh blocks contain this
    ! - nd    INTEGER(4)
    ! - sof   INTEGER(4)
    ! Specific to cartesian variable
    ! - nx    INTEGER(4)
    ! - stx   REAL(num)
    ! - dmin  REAL(num)
    ! - dmax  REAL(num)
    ! - mesh  CHARACTER(max_string_len)
    ! - class CHARACTER(max_string_len)

    md_length = meshtype_header_offset + ndims * soi + (ndims + 2) * sof &
        + 2 * max_string_len
    block_length = md_length + len_var * sof

    ! Write header
    CALL cfd_write_block_header(name, class, c_type_mesh_variable, &
        block_length, md_length, default_rank)
    CALL cfd_write_meshtype_header(c_var_cartesian, c_dimension_1d, sof, &
        default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      DO i = 1,ndims
        CALL MPI_FILE_WRITE(cfd_filehandle, dims4(i), 1, MPI_INTEGER4, &
            cfd_status, cfd_errcode)
      ENDDO
    ENDIF

    current_displacement = current_displacement + ndims * soi

    ! Determine data ranges and write out
    mn = MINVAL(variable)
    mx = MAXVAL(variable)
    CALL MPI_ALLREDUCE(mn, mn_global, 1, mpireal, MPI_MIN, cfd_comm, &
        cfd_errcode)
    CALL MPI_ALLREDUCE(mx, mx_global, 1, mpireal, MPI_MAX, cfd_comm, &
        cfd_errcode)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      ! Write out grid stagger
      CALL MPI_FILE_WRITE(cfd_filehandle, stagger, ndims, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mn_global, 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx_global, 1, mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + (ndims + 2) * sof

    ! Write the mesh name and class
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF

    current_displacement = current_displacement + 2 * max_string_len

    ! Write the actual data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, subarray, &
        distribution, "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_WRITE_ALL(cfd_filehandle, variable, 1, subarray, &
        cfd_status, cfd_errcode)

    current_displacement = current_displacement + len_var * sof

  END SUBROUTINE cfd_write_1d_cartesian_variable_parallel



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian variable in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global nx, ny
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_2d_cartesian_variable_parallel(name, class, dims, &
      stagger, mesh_name, mesh_class, variable, distribution, subarray)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    REAL(num), DIMENSION(:), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_name, mesh_class
    REAL(num), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num) :: mn, mx, mn_global, mx_global
    INTEGER(8) :: block_length, md_length, len_var
    INTEGER, PARAMETER :: ndims = 2
    INTEGER(4), DIMENSION(ndims) :: dims4
    INTEGER :: i

    len_var = 1
    DO i = 1,ndims
      dims4(i) = INT(dims(i),4)
      len_var = len_var * dims(i)
    ENDDO

    ! Metadata is
    ! - meshtype (INTEGER(4)) All mesh blocks contain this
    ! - nd    INTEGER(4)
    ! - sof   INTEGER(4)
    ! Specific to cartesian variable
    ! - nx    INTEGER(4)
    ! - ny    INTEGER(4)
    ! - stx   REAL(num)
    ! - sty   REAL(num)
    ! - dmin  REAL(num)
    ! - dmax  REAL(num)
    ! - mesh  CHARACTER(max_string_len)
    ! - class CHARACTER(max_string_len)

    md_length = meshtype_header_offset + ndims * soi + (ndims + 2) * sof &
        + 2 * max_string_len
    block_length = md_length + len_var * sof

    ! Write header
    CALL cfd_write_block_header(name, class, c_type_mesh_variable, &
        block_length, md_length, default_rank)
    CALL cfd_write_meshtype_header(c_var_cartesian, c_dimension_2d, sof, &
        default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      DO i = 1,ndims
        CALL MPI_FILE_WRITE(cfd_filehandle, dims4(i), 1, MPI_INTEGER4, &
            cfd_status, cfd_errcode)
      ENDDO
    ENDIF

    current_displacement = current_displacement + ndims * soi

    ! Determine data ranges and write out
    mn = MINVAL(variable)
    mx = MAXVAL(variable)
    CALL MPI_ALLREDUCE(mn, mn_global, 1, mpireal, MPI_MIN, cfd_comm, &
        cfd_errcode)
    CALL MPI_ALLREDUCE(mx, mx_global, 1, mpireal, MPI_MAX, cfd_comm, &
        cfd_errcode)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      ! Write out grid stagger
      CALL MPI_FILE_WRITE(cfd_filehandle, stagger, ndims, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mn_global, 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx_global, 1, mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + (ndims + 2) * sof

    ! Write the mesh name and class
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF

    current_displacement = current_displacement + 2 * max_string_len

    ! Write the actual data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, subarray, &
        distribution, "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_WRITE_ALL(cfd_filehandle, variable, 1, subarray, &
        cfd_status, cfd_errcode)

    current_displacement = current_displacement + len_var * sof

  END SUBROUTINE cfd_write_2d_cartesian_variable_parallel



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian variable in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global nx, ny, nz
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_3d_cartesian_variable_parallel(name, class, dims, &
      stagger, mesh_name, mesh_class, variable, distribution, subarray)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    REAL(num), DIMENSION(:), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_name, mesh_class
    REAL(num), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num) :: mn, mx, mn_global, mx_global
    INTEGER(8) :: block_length, md_length, len_var
    INTEGER, PARAMETER :: ndims = 3
    INTEGER(4), DIMENSION(ndims) :: dims4
    INTEGER :: i

    len_var = 1
    DO i = 1,ndims
      dims4(i) = INT(dims(i),4)
      len_var = len_var * dims(i)
    ENDDO

    ! Metadata is
    ! - meshtype (INTEGER(4)) All mesh blocks contain this
    ! - nd    INTEGER(4)
    ! - sof   INTEGER(4)
    ! Specific to cartesian variable
    ! - nx    INTEGER(4)
    ! - ny    INTEGER(4)
    ! - nz    INTEGER(4)
    ! - stx   REAL(num)
    ! - sty   REAL(num)
    ! - stz   REAL(num)
    ! - dmin  REAL(num)
    ! - dmax  REAL(num)
    ! - mesh  CHARACTER(max_string_len)
    ! - class CHARACTER(max_string_len)

    md_length = meshtype_header_offset + ndims * soi + (ndims + 2) * sof &
        + 2 * max_string_len
    block_length = md_length + len_var * sof

    ! Write header
    CALL cfd_write_block_header(name, class, c_type_mesh_variable, &
        block_length, md_length, default_rank)
    CALL cfd_write_meshtype_header(c_var_cartesian, c_dimension_3d, sof, &
        default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      DO i = 1,ndims
        CALL MPI_FILE_WRITE(cfd_filehandle, dims4(i), 1, MPI_INTEGER4, &
            cfd_status, cfd_errcode)
      ENDDO
    ENDIF

    current_displacement = current_displacement + ndims * soi

    ! Determine data ranges and write out
    mn = MINVAL(variable)
    mx = MAXVAL(variable)
    CALL MPI_ALLREDUCE(mn, mn_global, 1, mpireal, MPI_MIN, cfd_comm, &
        cfd_errcode)
    CALL MPI_ALLREDUCE(mx, mx_global, 1, mpireal, MPI_MAX, cfd_comm, &
        cfd_errcode)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      ! Write out grid stagger
      CALL MPI_FILE_WRITE(cfd_filehandle, stagger, ndims, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mn_global, 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx_global, 1, mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + (ndims + 2) * sof

    ! Write the mesh name and class
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF

    current_displacement = current_displacement + 2 * max_string_len

    ! Write the actual data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, subarray, &
        distribution, "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_WRITE_ALL(cfd_filehandle, variable, 1, subarray, &
        cfd_status, cfd_errcode)

    current_displacement = current_displacement + len_var * sof

  END SUBROUTINE cfd_write_3d_cartesian_variable_parallel



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian variable in serial using node with
  ! rank {rank_write} for writing
  ! Serial operation, so no need for nx
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_1d_cartesian_variable(name, class, stagger, mesh_name, &
      mesh_class, variable, rank_write)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), DIMENSION(:), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_name, mesh_class
    REAL(num), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: rank_write
    REAL(num) :: mn, mx, mn_global, mx_global
    INTEGER(8) :: block_length, md_length, len_var
    INTEGER, PARAMETER :: ndims = 1
    INTEGER(4), DIMENSION(ndims) :: dims4
    INTEGER :: i

    dims4 = SHAPE(variable)
    len_var = 1
    DO i = 1,ndims
      len_var = len_var * dims4(i)
    ENDDO

    ! Metadata is
    ! - meshtype (INTEGER(4)) All mesh blocks contain this
    ! - nd    INTEGER(4)
    ! - sof   INTEGER(4)
    ! Specific to cartesian variable
    ! - nx    INTEGER(4)
    ! - stx   REAL(num)
    ! - dmin  REAL(num)
    ! - dmax  REAL(num)
    ! - mesh  CHARACTER(max_string_len)
    ! - class CHARACTER(max_string_len)

    md_length = meshtype_header_offset + ndims * soi + (ndims + 2) * sof &
        + 2 * max_string_len
    block_length = md_length + len_var * sof

    ! Write header
    CALL cfd_write_block_header(name, class, c_type_mesh_variable, &
        block_length, md_length, rank_write)
    CALL cfd_write_meshtype_header(c_var_cartesian, c_dimension_1d, sof, &
        rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      DO i = 1,ndims
        CALL MPI_FILE_WRITE(cfd_filehandle, dims4(i), 1, MPI_INTEGER4, &
            cfd_status, cfd_errcode)
      ENDDO
    ENDIF

    current_displacement = current_displacement + ndims * soi

    ! Determine data ranges and write out
    mn = MINVAL(variable)
    mx = MAXVAL(variable)
    CALL MPI_ALLREDUCE(mn, mn_global, 1, mpireal, MPI_MIN, cfd_comm, &
        cfd_errcode)
    CALL MPI_ALLREDUCE(mx, mx_global, 1, mpireal, MPI_MAX, cfd_comm, &
        cfd_errcode)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      ! Write out grid stagger
      CALL MPI_FILE_WRITE(cfd_filehandle, stagger, ndims, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mn_global, 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx_global, 1, mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + (ndims + 2) * sof

    ! Write the mesh name and class
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF

    current_displacement = current_displacement + 2 * max_string_len

    ! Write the actual data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, variable, len_var, mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + len_var * sof

  END SUBROUTINE cfd_write_1d_cartesian_variable



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian variable in serial using node with
  ! rank {rank_write} for writing
  ! Serial operation, so no need for nx, ny
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_2d_cartesian_variable(name, class, stagger, mesh_name, &
      mesh_class, variable, rank_write)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), DIMENSION(:), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_name, mesh_class
    REAL(num), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: rank_write
    REAL(num) :: mn, mx, mn_global, mx_global
    INTEGER(8) :: block_length, md_length, len_var
    INTEGER, PARAMETER :: ndims = 2
    INTEGER(4), DIMENSION(ndims) :: dims4
    INTEGER :: i

    dims4 = SHAPE(variable)
    len_var = 1
    DO i = 1,ndims
      len_var = len_var * dims4(i)
    ENDDO

    ! Metadata is
    ! - meshtype (INTEGER(4)) All mesh blocks contain this
    ! - nd    INTEGER(4)
    ! - sof   INTEGER(4)
    ! Specific to cartesian variable
    ! - nx    INTEGER(4)
    ! - ny    INTEGER(4)
    ! - stx   REAL(num)
    ! - sty   REAL(num)
    ! - dmin  REAL(num)
    ! - dmax  REAL(num)
    ! - mesh  CHARACTER(max_string_len)
    ! - class CHARACTER(max_string_len)

    md_length = meshtype_header_offset + ndims * soi + (ndims + 2) * sof &
        + 2 * max_string_len
    block_length = md_length + len_var * sof

    ! Write header
    CALL cfd_write_block_header(name, class, c_type_mesh_variable, &
        block_length, md_length, rank_write)
    CALL cfd_write_meshtype_header(c_var_cartesian, c_dimension_2d, sof, &
        rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      DO i = 1,ndims
        CALL MPI_FILE_WRITE(cfd_filehandle, dims4(i), 1, MPI_INTEGER4, &
            cfd_status, cfd_errcode)
      ENDDO
    ENDIF

    current_displacement = current_displacement + ndims * soi

    ! Determine data ranges and write out
    mn = MINVAL(variable)
    mx = MAXVAL(variable)
    CALL MPI_ALLREDUCE(mn, mn_global, 1, mpireal, MPI_MIN, cfd_comm, &
        cfd_errcode)
    CALL MPI_ALLREDUCE(mx, mx_global, 1, mpireal, MPI_MAX, cfd_comm, &
        cfd_errcode)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      ! Write out grid stagger
      CALL MPI_FILE_WRITE(cfd_filehandle, stagger, ndims, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mn_global, 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx_global, 1, mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + (ndims + 2) * sof

    ! Write the mesh name and class
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF

    current_displacement = current_displacement + 2 * max_string_len

    ! Write the actual data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, variable, len_var, mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + len_var * sof

  END SUBROUTINE cfd_write_2d_cartesian_variable



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian variable in serial using node with
  ! rank {rank_write} for writing
  ! Serial operation, so no need for nx, ny, nz
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_3d_cartesian_variable(name, class, stagger, mesh_name, &
      mesh_class, variable, rank_write)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), DIMENSION(:), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_name, mesh_class
    REAL(num), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: rank_write
    REAL(num) :: mn, mx, mn_global, mx_global
    INTEGER(8) :: block_length, md_length, len_var
    INTEGER, PARAMETER :: ndims = 3
    INTEGER(4), DIMENSION(ndims) :: dims4
    INTEGER :: i

    dims4 = SHAPE(variable)
    len_var = 1
    DO i = 1,ndims
      len_var = len_var * dims4(i)
    ENDDO

    ! Metadata is
    ! - meshtype (INTEGER(4)) All mesh blocks contain this
    ! - nd    INTEGER(4)
    ! - sof   INTEGER(4)
    ! Specific to cartesian variable
    ! - nx    INTEGER(4)
    ! - ny    INTEGER(4)
    ! - nz    INTEGER(4)
    ! - stx   REAL(num)
    ! - sty   REAL(num)
    ! - stz   REAL(num)
    ! - dmin  REAL(num)
    ! - dmax  REAL(num)
    ! - mesh  CHARACTER(max_string_len)
    ! - class CHARACTER(max_string_len)

    md_length = meshtype_header_offset + ndims * soi + (ndims + 2) * sof &
        + 2 * max_string_len
    block_length = md_length + len_var * sof

    ! Write header
    CALL cfd_write_block_header(name, class, c_type_mesh_variable, &
        block_length, md_length, rank_write)
    CALL cfd_write_meshtype_header(c_var_cartesian, c_dimension_3d, sof, &
        rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      DO i = 1,ndims
        CALL MPI_FILE_WRITE(cfd_filehandle, dims4(i), 1, MPI_INTEGER4, &
            cfd_status, cfd_errcode)
      ENDDO
    ENDIF

    current_displacement = current_displacement + ndims * soi

    ! Determine data ranges and write out
    mn = MINVAL(variable)
    mx = MAXVAL(variable)
    CALL MPI_ALLREDUCE(mn, mn_global, 1, mpireal, MPI_MIN, cfd_comm, &
        cfd_errcode)
    CALL MPI_ALLREDUCE(mx, mx_global, 1, mpireal, MPI_MAX, cfd_comm, &
        cfd_errcode)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      ! Write out grid stagger
      CALL MPI_FILE_WRITE(cfd_filehandle, stagger, ndims, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mn_global, 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx_global, 1, mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + (ndims + 2) * sof

    ! Write the mesh name and class
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF

    current_displacement = current_displacement + 2 * max_string_len

    ! Write the actual data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, variable, len_var, mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + len_var * sof

  END SUBROUTINE cfd_write_3d_cartesian_variable

END MODULE output_cartesian
