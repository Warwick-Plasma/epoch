MODULE output_cartesian

  USE shared_data
  USE iocommon
  USE output

  IMPLICIT NONE

CONTAINS

  ! -------------------------------------------------------------------------------------------------------------------
  ! Code to write a 1D Cartesian grid in serial from the node with rank {rank_write}
  ! Serial operation, so no need to specify nx, ny
  ! -------------------------------------------------------------------------------------------------------------------
  SUBROUTINE cfd_write_1d_cartesian_grid(name, class, x, rank_write)

    REAL(num), DIMENSION(:), INTENT(IN) :: x
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(4) :: nx
    INTEGER(8) :: block_length, md_length

    nx = SIZE(x)

    ! Metadata is
    ! * ) MeshType (INTEGER(4)) All mesh blocks contain this
    ! * ) nd   INTEGER(4)
    ! * ) sof  INTEGER(4)
    ! Specific to Cartesian Grid
    ! 1 ) nx   INTEGER(4)
    ! 2 ) xmin REAL(num)
    ! 3 ) xmax REAL(num)

    md_length = meshtype_header_offset + 1 * soi + 2 * num ! 1 ints, 2 reals + MeshType Header
    block_length = md_length + nx*num

    ! Now written header, write metadata
    CALL cfd_write_block_header(name, class, c_type_mesh, block_length, md_length, rank_write)

    CALL cfd_write_meshtype_header(c_mesh_cartesian, c_dimension_1d, num, rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, nx, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 1 * soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, MINVAL(x), 1, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, MAXVAL(x), 1, mpireal, cfd_status, cfd_errcode)

      ! Now write the real arrays
      CALL MPI_FILE_WRITE(cfd_filehandle, x, nx, mpireal, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 2*num + nx * num

  END SUBROUTINE cfd_write_1d_cartesian_grid



  ! -------------------------------------------------------------------------------------------------------------------
  ! Code to write a 2D Cartesian grid in serial from the node with rank {rank_write}
  ! Serial operation, so no need to specify nx, ny
  ! -------------------------------------------------------------------------------------------------------------------
  SUBROUTINE cfd_write_2d_cartesian_grid(name, class, x, y, rank_write)

    REAL(num), DIMENSION(:), INTENT(IN) :: x, y
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(4) :: nx, ny
    INTEGER(8) :: block_length, md_length

    nx = SIZE(x)
    ny = SIZE(y)

    ! Metadata is
    ! * ) MeshType (INTEGER(4)) All mesh blocks contain this
    ! * ) nd   INTEGER(4)
    ! * ) sof  INTEGER(4)
    ! Specific to Cartesian Grid
    ! 1 ) nx   INTEGER(4)
    ! 2 ) ny   INTEGER(4)
    ! 3 ) xmin REAL(num)
    ! 4 ) xmax REAL(num)
    ! 5 ) ymin REAL(num)
    ! 6 ) ymax REAL(num)

    md_length = meshtype_header_offset + 2 * soi + 4 * num ! 2 ints, 6 reals + MeshType Header
    block_length = md_length + (nx+ny)*num

    ! Now written header, write metadata
    CALL cfd_write_block_header(name, class, c_type_mesh, block_length, md_length, rank_write)

    CALL cfd_write_meshtype_header(c_mesh_cartesian, c_dimension_2d, num, rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, nx, 1, MPI_INTEGER, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, ny, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 2 * soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, MINVAL(x), 1, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, MAXVAL(x), 1, mpireal, cfd_status, cfd_errcode)

      CALL MPI_FILE_WRITE(cfd_filehandle, MINVAL(y), 1, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, MAXVAL(y), 1, mpireal, cfd_status, cfd_errcode)

      ! Now write the real arrays
      CALL MPI_FILE_WRITE(cfd_filehandle, x, nx, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, y, ny, mpireal, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 4*num + (nx + ny)*num

  END SUBROUTINE cfd_write_2d_cartesian_grid



  ! -------------------------------------------------------------------------------------------------------------------
  ! Code to write a 3D Cartesian grid in serial from the node with rank {rank_write}
  ! Serial operation, so no need to specify nx, ny, nz
  ! -------------------------------------------------------------------------------------------------------------------
  SUBROUTINE cfd_write_3d_cartesian_grid(name, class, x, y, z, rank_write)

    REAL(num), DIMENSION(:), INTENT(IN) :: x, y, z
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(4) :: nx, ny, nz
    INTEGER(8) :: block_length, md_length

    nx = SIZE(x)
    ny = SIZE(y)
    nz = SIZE(z)

    ! Metadata is
    ! * ) MeshType (INTEGER(4)) All mesh blocks contain this
    ! * ) nd   INTEGER(4)
    ! * ) sof  INTEGER(4)
    ! Specific to Cartesian Grid
    ! 2 ) nx   INTEGER(4)
    ! 3 ) ny   INTEGER(4)
    ! 4 ) nz   INTEGER(4)
    ! 5 ) xmin REAL(num)
    ! 6 ) xmax REAL(num)
    ! 7 ) ymin REAL(num)
    ! 8 ) ymax REAL(num)
    ! 9 ) zmin REAL(num)
    ! 10) zmax REAL(num)

    md_length = meshtype_header_offset + 3 * soi + 6 * num ! 3 ints, 6 reals + MeshType Header
    block_length = md_length + (nx+ny+nz)*num

    ! Now written header, write metadata
    CALL cfd_write_block_header(name, class, c_type_mesh, block_length, md_length, rank_write)

    CALL cfd_write_meshtype_header(c_mesh_cartesian, c_dimension_3d, num, rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, nx, 1, MPI_INTEGER, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, ny, 1, MPI_INTEGER, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, nz, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 3 * soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, MINVAL(x), 1, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, MAXVAL(x), 1, mpireal, cfd_status, cfd_errcode)

      CALL MPI_FILE_WRITE(cfd_filehandle, MINVAL(y), 1, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, MAXVAL(y), 1, mpireal, cfd_status, cfd_errcode)

      CALL MPI_FILE_WRITE(cfd_filehandle, MINVAL(z), 1, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, MAXVAL(z), 1, mpireal, cfd_status, cfd_errcode)

      ! Now write the real arrays
      CALL MPI_FILE_WRITE(cfd_filehandle, x, nx, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, y, ny, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, z, nz, mpireal, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 6*num + (nx + ny + nz)*num

  END SUBROUTINE cfd_write_3d_cartesian_grid



  ! -------------------------------------------------------------------------------------------------------------------
  ! Code to write a 3D Cartesian variable in parallel using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution
  ! parallel operation, so need global nx, ny, nz
  ! -------------------------------------------------------------------------------------------------------------------

  SUBROUTINE cfd_write_3d_cartesian_variable_parallel(name, class, dims, stagger, mesh_name, mesh_class, variable, distribution)

    CHARACTER(LEN=*), INTENT(IN) :: name, class, mesh_name, mesh_class
    REAL(num), DIMENSION(:, :, :), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution
    INTEGER, INTENT(IN), DIMENSION(3) :: dims
    REAL(num), INTENT(IN), DIMENSION(3) :: stagger
    REAL(num) :: mn, mx, mn_global, mx_global
    INTEGER(8) :: block_length, md_length, len_var
    INTEGER :: nx, ny, nz

    ! *) VariableType (INTEGER(4)) All variable blocks contain this
    ! These are specific to a cartesian variable
    ! * ) nd   INTEGER(4)
    ! * ) sof  INTEGER(4)
    ! 1 ) nx   INTEGER(4)
    ! 2 ) ny   INTEGER(4)
    ! 3 ) nz   INTEGER(4)
    ! 4 ) stx  REAL(num)
    ! 5 ) sty  REAL(num)
    ! 6 ) stz  REAL(num)
    ! 7 ) dmin REAL(num)
    ! 8 ) dmax REAL(num)
    ! 9 ) Mesh CHARACTER(max_string_len)
    ! 10) class CHARACTER(max_string_len)

    len_var = SIZE(variable)
    nx = dims(1)
    ny = dims(2)
    nz = dims(3)

    md_length = meshtype_header_offset + 3 * soi + 5 * num + 2 * max_string_len ! 3 INTs 5 REALs 2STRINGs
    block_length = md_length + num * nx * ny * nz

    ! Write the common stuff
    CALL cfd_write_block_header(name, class, c_type_mesh_variable, block_length, md_length, default_rank)
    CALL cfd_write_meshtype_header(c_var_cartesian, c_dimension_3d, num, default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, nx, 1, MPI_INTEGER, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, ny, 1, MPI_INTEGER, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, nz, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 3 * soi

    ! Set the file view
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    ! Determine data ranges and write out
    mn = MINVAL(variable)
    mx = MAXVAL(variable)
    CALL MPI_ALLREDUCE(mn, mn_global, 1, mpireal, MPI_MIN, cfd_comm, cfd_errcode)
    CALL MPI_ALLREDUCE(mx, mx_global, 1, mpireal, MPI_MAX, cfd_comm, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
      ! Write out grid stagger
      CALL MPI_FILE_WRITE(cfd_filehandle, stagger, 3, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mn_global, 1, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx_global, 1, mpireal, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 5 * num

    ! Write the mesh name and class
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF
    current_displacement = current_displacement + 2 * max_string_len

    ! Write the actual data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, distribution, "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_WRITE_ALL(cfd_filehandle, variable, len_var, mpireal, cfd_status, cfd_errcode)
    current_displacement = current_displacement + num*nx*ny*nz

  END SUBROUTINE cfd_write_3d_cartesian_variable_parallel



  ! -------------------------------------------------------------------------------------------------------------------
  ! Code to write a 2D Cartesian variable in parallel using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution
  ! parallel operation, so need global nx, ny
  ! -------------------------------------------------------------------------------------------------------------------

  SUBROUTINE cfd_write_2d_cartesian_variable_parallel(name, class, dims, stagger, mesh_name, mesh_class, variable, distribution)

    CHARACTER(LEN=*), INTENT(IN) :: name, class, mesh_name, mesh_class
    REAL(num), DIMENSION(:, :), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution
    INTEGER, INTENT(IN), DIMENSION(2) :: dims
    REAL(num), INTENT(IN), DIMENSION(2) :: stagger
    INTEGER :: nx, ny
    REAL(num) :: mn, mx, mn_global, mx_global
    INTEGER(8) :: block_length, md_length, len_var

    ! *) VariableType (INTEGER(4)) All variable blocks contain this
    ! These are specific to a cartesian variable
    ! * ) nd   INTEGER(4)
    ! * ) sof  INTEGER(4)
    ! 1 ) nx   INTEGER(4)
    ! 2 ) ny   INTEGER(4)
    ! 3 ) stx  REAL(num)
    ! 4 ) sty  REAL(num)
    ! 5 ) dmin REAL(num)
    ! 6 ) dmax REAL(num)
    ! 7 ) Mesh CHARACTER(max_string_len)
    ! 8 ) class CHARACTER(max_string_len)

    len_var = SIZE(variable)
    nx = dims(1)
    ny = dims(2)

    md_length = meshtype_header_offset + 2 * soi + 4 * num + 2 * max_string_len ! 3 INTs 2 REALs
    block_length = md_length + num * nx * ny

    CALL cfd_write_block_header(name, class, c_type_mesh_variable, block_length, md_length, default_rank)
    CALL cfd_write_meshtype_header(c_var_cartesian, c_dimension_2d, num, default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, nx, 1, MPI_INTEGER, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, ny, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 2 * soi

    ! Determine data ranges and write out
    mn = MINVAL(variable)
    mx = MAXVAL(variable)
    CALL MPI_ALLREDUCE(mn, mn_global, 1, mpireal, MPI_MIN, cfd_comm, cfd_errcode)
    CALL MPI_ALLREDUCE(mx, mx_global, 1, mpireal, MPI_MAX, cfd_comm, cfd_errcode)
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, stagger, 2, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mn_global, 1, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx_global, 1, mpireal, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 4 * num

    ! Write the mesh name and class
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF
    current_displacement = current_displacement + 2 * max_string_len

    ! Write the actual data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, distribution, "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_WRITE_ALL(cfd_filehandle, variable, len_var, mpireal, cfd_status, cfd_errcode)

    current_displacement = current_displacement + num*nx*ny

  END SUBROUTINE cfd_write_2d_cartesian_variable_parallel



  ! -------------------------------------------------------------------------------------------------------------------
  ! Code to write a 1D Cartesian variable in parallel using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution
  ! parallel operation, so need global nx
  ! -------------------------------------------------------------------------------------------------------------------

  SUBROUTINE cfd_write_1d_cartesian_variable_parallel(name, class, dims, stagger, mesh_name, mesh_class, variable, distribution)

    CHARACTER(LEN=*), INTENT(IN) :: name, class, mesh_name, mesh_class
    REAL(num), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution
    INTEGER, INTENT(IN) :: dims
    REAL(num), INTENT(IN) :: stagger
    INTEGER :: nx
    REAL(num) :: mn, mx, mn_global, mx_global
    INTEGER(8) :: block_length, md_length, len_var

    ! *) VariableType (INTEGER(4)) All variable blocks contain this
    ! These are specific to a cartesian variable
    ! * ) nd   INTEGER(4)
    ! * ) sof  INTEGER(4)
    ! 1 ) nx   INTEGER(4)
    ! 2 ) stx  REAL(num)
    ! 3 ) dmin REAL(num)
    ! 4 ) dmax REAL(num)
    ! 5 ) Mesh CHARACTER(max_string_len)
    ! 6 ) class CHARACTER(max_string_len)

    len_var = SIZE(variable)
    nx = dims
    md_length = meshtype_header_offset + 1 * soi + 3 * num + 2 * max_string_len ! 1 INTs 3 REALs 2 strings
    block_length = md_length + num * nx

    CALL cfd_write_block_header(name, class, c_type_mesh_variable, block_length, md_length, default_rank)
    CALL cfd_write_meshtype_header(c_var_cartesian, c_dimension_1d, num, default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, nx, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 1 * soi

    ! Determine data ranges and write out
    mn = MINVAL(variable)
    mx = MAXVAL(variable)
    CALL MPI_ALLREDUCE(mn, mn_global, 1, mpireal, MPI_MIN, cfd_comm, cfd_errcode)
    CALL MPI_ALLREDUCE(mx, mx_global, 1, mpireal, MPI_MAX, cfd_comm, cfd_errcode)
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, stagger, 1, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mn_global, 1, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx_global, 1, mpireal, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 3 * num

    ! Write the mesh name and class
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF
    current_displacement = current_displacement + 2 * max_string_len

    ! Write the actual data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, distribution, "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_WRITE_ALL(cfd_filehandle, variable, len_var, mpireal, cfd_status, cfd_errcode)

    current_displacement = current_displacement + num*nx

  END SUBROUTINE cfd_write_1d_cartesian_variable_parallel



  ! -------------------------------------------------------------------------------------------------------------------
  ! Code to write a 1D Cartesian variable in serial using node with rank {rank_write} for writing
  ! Serial operation, so no need for nx, ny
  ! -------------------------------------------------------------------------------------------------------------------

  SUBROUTINE cfd_write_1d_cartesian_variable(name, class, stagger, mesh_name, mesh_class, variable, rank_write)

    CHARACTER(LEN=*), INTENT(IN) :: name, class, mesh_name, mesh_class
    REAL(num), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: rank_write
    REAL(num), INTENT(IN), DIMENSION(1) :: stagger
    INTEGER :: nx
    REAL(num) :: mn, mx, mn_global, mx_global
    INTEGER, DIMENSION(1) :: dims
    INTEGER(8) :: block_length, md_length, len_var

    ! *) VariableType (INTEGER(4)) All variable blocks contain this
    ! These are specific to a cartesian variable
    ! * ) nd   INTEGER(4)
    ! * ) sof  INTEGER(4)
    ! 1 ) nx   INTEGER(4)
    ! 2 ) stx  REAL(num)
    ! 3 ) dmin REAL(num)
    ! 4 ) dmax REAL(num)
    ! 5 ) Mesh CHARACTER(max_string_len)
    ! 6 ) class CHARACTER(max_string_len)
    len_var = SIZE(variable)
    dims = SHAPE(variable)

    nx = dims(1)

    md_length = meshtype_header_offset + 1 * soi + 3 * num + 2 * max_string_len ! 1 INTs 3 REALs
    block_length = md_length + num * nx

    CALL cfd_write_block_header(name, class, c_type_mesh_variable, block_length, md_length, rank_write)
    CALL cfd_write_meshtype_header(c_var_cartesian, c_dimension_1d, num, rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)
    ! This is the serial version remember

    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, nx, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 1 * soi

    ! Determine data ranges and write out
    mn = MINVAL(variable)
    mx = MAXVAL(variable)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, stagger, 2, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mn_global, 1, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx_global, 1, mpireal, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 3 * num

    ! Write the mesh name and class
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF

    current_displacement = current_displacement + 2 * max_string_len

    ! Write the actual data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write) CALL MPI_FILE_WRITE(cfd_filehandle, variable, len_var, mpireal, cfd_status, cfd_errcode)

    current_displacement = current_displacement + num*nx

  END SUBROUTINE cfd_write_1d_cartesian_variable



  ! -------------------------------------------------------------------------------------------------------------------
  ! Code to write a 2D Cartesian variable in serial using node with rank {rank_write} for writing
  ! Serial operation, so no need for nx, ny
  ! -------------------------------------------------------------------------------------------------------------------

  SUBROUTINE cfd_write_2d_cartesian_variable(name, class, stagger, mesh_name, mesh_class, variable, rank_write)

    CHARACTER(LEN=*), INTENT(IN) :: name, class, mesh_name, mesh_class
    REAL(num), DIMENSION(:, :), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: rank_write
    REAL(num), INTENT(IN), DIMENSION(2) :: stagger
    INTEGER :: nx, ny
    REAL(num) :: mn, mx, mn_global, mx_global
    INTEGER, DIMENSION(2) :: dims
    INTEGER(8) :: block_length, md_length, len_var

    ! *) VariableType (INTEGER(4)) All variable blocks contain this
    ! These are specific to a cartesian variable
    ! * ) nd   INTEGER(4)
    ! * ) sof  INTEGER(4)
    ! 1 ) nx   INTEGER(4)
    ! 2 ) ny   INTEGER(4)
    ! 3 ) stx  REAL(num)
    ! 4 ) sty  REAL(num)
    ! 5 ) dmin REAL(num)
    ! 6 ) dmax REAL(num)
    ! 7 ) Mesh CHARACTER(max_string_len)
    ! 8 ) class CHARACTER(max_string_len)

    len_var = SIZE(variable)
    dims = SHAPE(variable)

    nx = dims(1)
    ny = dims(2)

    md_length = meshtype_header_offset + 2 * soi + 4 * num + 2 * max_string_len ! 2 INTs 4 REALs
    block_length = md_length + num * nx * ny

    CALL cfd_write_block_header(name, class, c_type_mesh_variable, block_length, md_length, rank_write)
    CALL cfd_write_meshtype_header(c_var_cartesian, c_dimension_2d, num, rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)
    ! This is the serial version remember

    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, nx, 1, MPI_INTEGER, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, ny, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 2 * soi

    ! Determine data ranges and write out
    mn = MINVAL(variable)
    mx = MAXVAL(variable)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, stagger, 2, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mn_global, 1, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx_global, 1, mpireal, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 4 * num

    ! Write the mesh name and class
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF

    current_displacement = current_displacement + 2 * max_string_len

    ! Write the actual data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write) CALL MPI_FILE_WRITE(cfd_filehandle, variable, len_var, mpireal, cfd_status, cfd_errcode)

    current_displacement = current_displacement + num*nx*ny

  END SUBROUTINE cfd_write_2d_cartesian_variable



  ! -------------------------------------------------------------------------------------------------------------------
  ! Code to write a 3D Cartesian variable in serial using node with rank {rank_write} for writing
  ! Serial operation, so no need for nx, ny, nz
  ! -------------------------------------------------------------------------------------------------------------------

  SUBROUTINE cfd_write_3d_cartesian_variable(name, class, stagger, mesh_name, mesh_class, variable, rank_write)

    CHARACTER(LEN=*), INTENT(IN) :: name, class, mesh_name, mesh_class
    REAL(num), DIMENSION(:, :, :), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: rank_write
    REAL(num), INTENT(IN), DIMENSION(3) :: stagger
    INTEGER :: nx, ny, nz
    INTEGER, DIMENSION(3) :: dims
    REAL(num) :: mn, mx, mn_global, mx_global
    INTEGER(8) :: block_length, md_length, len_var

    ! *) VariableType (INTEGER(4)) All variable blocks contain this
    ! These are specific to a cartesian variable
    ! * ) nd   INTEGER(4)
    ! * ) sof  INTEGER(4)
    ! 1 ) nx   INTEGER(4)
    ! 2 ) ny   INTEGER(4)
    ! 3 ) nz   INTEGER(4)
    ! 4 ) stx  REAL(num)
    ! 5 ) sty  REAL(num)
    ! 6 ) stz  REAL(num)
    ! 7 ) dmin REAL(num)
    ! 8 ) dmax REAL(num)
    ! 9 ) Mesh CHARACTER(max_string_len)
    ! 10) class CHARACTER(max_string_len)

    len_var = SIZE(variable)
    dims = SHAPE(variable)

    nx = dims(1)
    ny = dims(2)
    nz = dims(3)

    md_length = meshtype_header_offset + 3 * soi + 5 * num + 2 * max_string_len ! 3 INTs 5 REALs
    block_length = md_length + num * nx * ny * nz

    CALL cfd_write_block_header(name, class, c_type_mesh_variable, block_length, md_length, rank_write)
    CALL cfd_write_meshtype_header(c_var_cartesian, c_dimension_3d, num, rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)
    ! This is the serial version remember

    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, nx, 1, MPI_INTEGER, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, ny, 1, MPI_INTEGER, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, nz, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 3 * soi

    ! Determine data ranges and write out
    mn = MINVAL(variable)
    mx = MAXVAL(variable)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, stagger, 3, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mn_global, 1, mpireal, cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx_global, 1, mpireal, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 5 * num

    ! Write the mesh name and class
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF
    current_displacement = current_displacement + 2 * max_string_len

    ! Write the actual data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write) CALL MPI_FILE_WRITE(cfd_filehandle, variable, len_var, mpireal, cfd_status, cfd_errcode)
    current_displacement = current_displacement + num*nx*ny*nz

  END SUBROUTINE cfd_write_3d_cartesian_variable

END MODULE output_cartesian
