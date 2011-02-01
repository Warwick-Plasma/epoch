MODULE sdf_output_cartesian

  USE sdf_common
  USE sdf_output
  USE mpi

  IMPLICIT NONE

  INTERFACE sdf_write_srl_plain_mesh
    MODULE PROCEDURE &
        sdf_write_srl_1d_mesh, &
        sdf_write_srl_2d_mesh, &
        sdf_write_srl_3d_mesh
  END INTERFACE sdf_write_srl_plain_mesh

  INTERFACE sdf_write_plain_mesh
    MODULE PROCEDURE &
        sdf_write_1d_mesh, &
        sdf_write_2d_mesh, &
        sdf_write_3d_mesh
  END INTERFACE sdf_write_plain_mesh

  INTERFACE sdf_write_plain_variable
    MODULE PROCEDURE &
        sdf_write_1d_float, &
        sdf_write_2d_float, &
        sdf_write_3d_float, &
        sdf_write_1d_character, &
        sdf_write_2d_character, &
        sdf_write_3d_character, &
        sdf_write_1d_integer, &
        sdf_write_2d_integer, &
        sdf_write_3d_integer
  END INTERFACE sdf_write_plain_variable

  INTERFACE sdf_write_material
    MODULE PROCEDURE &
        sdf_write_1d_material, &
        sdf_write_2d_material, &
        sdf_write_3d_material
  END INTERFACE sdf_write_material

  INTERFACE sdf_write_matvar
    MODULE PROCEDURE &
        sdf_write_1d_matvar, &
        sdf_write_2d_matvar, &
        sdf_write_3d_matvar
  END INTERFACE sdf_write_matvar

  INTERFACE sdf_write_species
    MODULE PROCEDURE &
        sdf_write_1d_species, &
        sdf_write_2d_species, &
        sdf_write_3d_species
  END INTERFACE sdf_write_species

CONTAINS

  SUBROUTINE write_mesh_meta(h, id, name, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER :: ndims
    TYPE(sdf_block_type), POINTER :: b
    INTEGER :: i, errcode

    b => h%current_block

    b%blocktype = c_blocktype_plain_mesh
    ndims = b%ndims

    b%nelements = 0
    DO i = 1,ndims
      b%nelements = b%nelements + b%dims(i)
    ENDDO

    ! Metadata is
    ! - mults     REAL(r8), DIMENSION(ndims)
    ! - labels    CHARACTER(id_length), DIMENSION(ndims)
    ! - units     CHARACTER(id_length), DIMENSION(ndims)
    ! - geometry  INTEGER(i4)
    ! - minval    REAL(r8), DIMENSION(ndims)
    ! - maxval    REAL(r8), DIMENSION(ndims)
    ! - dims      INTEGER(i4), DIMENSION(ndims)

    b%info_length = h%block_header_length + (ndims + 1) * soi4 &
        + (3 * ndims) * sof8 + 2 * ndims * c_id_length
    b%data_length = b%nelements * b%type_size

    ! Write header

    IF (PRESENT(id)) THEN
      ALLOCATE(b%dim_labels(ndims), b%dim_units(ndims), b%dim_mults(ndims))

      IF (PRESENT(dim_labels)) THEN
        DO i = 1,ndims
          CALL safe_copy_string(dim_labels(i), b%dim_labels(i))
        ENDDO
      ELSE
        IF (ndims .GE. 1) CALL safe_copy_string('X', b%dim_labels(1))
        IF (ndims .GE. 2) CALL safe_copy_string('Y', b%dim_labels(2))
        IF (ndims .GE. 3) CALL safe_copy_string('Z', b%dim_labels(3))
      ENDIF

      IF (PRESENT(dim_units)) THEN
        DO i = 1,ndims
          CALL safe_copy_string(dim_units(i), b%dim_units(i))
        ENDDO
      ELSE
        DO i = 1,ndims
          CALL safe_copy_string('m', b%dim_units(i))
        ENDDO
      ENDIF

      IF (PRESENT(dim_mults)) THEN
        DO i = 1,ndims
          b%dim_mults(i) = REAL(dim_mults(i),r8)
        ENDDO
      ELSE
        DO i = 1,ndims
          b%dim_mults(i) = 1.d0
        ENDDO
      ENDIF

      CALL sdf_write_block_header(h, id, name)
    ELSE
      CALL write_block_header(h)
    ENDIF

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_WRITE(h%filehandle, b%dim_mults, ndims, MPI_REAL8, &
          MPI_STATUS_IGNORE, errcode)

      DO i = 1,ndims
        CALL sdf_safe_write_id(h, b%dim_labels(i))
      ENDDO

      DO i = 1,ndims
        CALL sdf_safe_write_id(h, b%dim_units(i))
      ENDDO

      CALL MPI_FILE_WRITE(h%filehandle, b%geometry, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%extents, 2 * ndims, MPI_REAL8, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%dims, ndims, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_location = b%block_start + b%info_length
    b%done_info = .TRUE.

  END SUBROUTINE write_mesh_meta



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian mesh in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_srl_1d_mesh(h, id, name, x, dim_labels, &
      dim_units, dim_mults, rank_write)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(num), DIMENSION(:), INTENT(IN) :: x
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: errcode, intn
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%geometry = c_geometry_cartesian
    b%ndims = ndims

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    b%dims(1) = INT(SIZE(x),i4)

    IF (h%rank .EQ. h%rank_master) THEN
      b%extents(1) = REAL(x(1),r8)
      b%extents(ndims+1) = REAL(x(b%dims(1)),r8)
    ENDIF

    ! Write header

    CALL write_mesh_meta(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    IF (h%rank .EQ. h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      intn = b%dims(1)
      CALL MPI_FILE_WRITE(h%filehandle, x, intn, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_srl_1d_mesh



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian mesh in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_srl_2d_mesh(h, id, name, x, y, dim_labels, &
      dim_units, dim_mults, rank_write)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(num), DIMENSION(:), INTENT(IN) :: x, y
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: errcode, intn
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%geometry = c_geometry_cartesian
    b%ndims = ndims

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    b%dims(1) = INT(SIZE(x),i4)
    b%dims(2) = INT(SIZE(y),i4)

    IF (h%rank .EQ. h%rank_master) THEN
      b%extents(1) = REAL(x(1),r8)
      b%extents(2) = REAL(y(1),r8)
      b%extents(ndims+1) = REAL(x(b%dims(1)),r8)
      b%extents(ndims+2) = REAL(y(b%dims(2)),r8)
    ENDIF

    ! Write header

    CALL write_mesh_meta(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    IF (h%rank .EQ. h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      intn = b%dims(1)
      CALL MPI_FILE_WRITE(h%filehandle, x, intn, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
      intn = b%dims(2)
      CALL MPI_FILE_WRITE(h%filehandle, y, intn, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_srl_2d_mesh



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian mesh in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_srl_3d_mesh(h, id, name, x, y, z, dim_labels, &
      dim_units, dim_mults, rank_write)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(num), DIMENSION(:), INTENT(IN) :: x, y, z
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: errcode, intn
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%geometry = c_geometry_cartesian
    b%ndims = ndims

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    b%dims(1) = INT(SIZE(x),i4)
    b%dims(2) = INT(SIZE(y),i4)
    b%dims(3) = INT(SIZE(z),i4)

    IF (h%rank .EQ. h%rank_master) THEN
      b%extents(1) = REAL(x(1),r8)
      b%extents(2) = REAL(y(1),r8)
      b%extents(3) = REAL(z(1),r8)
      b%extents(ndims+1) = REAL(x(b%dims(1)),r8)
      b%extents(ndims+2) = REAL(y(b%dims(2)),r8)
      b%extents(ndims+3) = REAL(z(b%dims(3)),r8)
    ENDIF

    ! Write header

    CALL write_mesh_meta(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    IF (h%rank .EQ. h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      intn = b%dims(1)
      CALL MPI_FILE_WRITE(h%filehandle, x, intn, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
      intn = b%dims(2)
      CALL MPI_FILE_WRITE(h%filehandle, y, intn, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
      intn = b%dims(3)
      CALL MPI_FILE_WRITE(h%filehandle, z, intn, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_srl_3d_mesh



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_1d_mesh(h, id, name, x, dims, xmin, xmax, &
      distribution, subarray, dim_labels, dim_units, &
      dim_mults)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(num), DIMENSION(:), INTENT(IN) :: x
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    REAL(num), INTENT(IN) :: xmin, xmax
    INTEGER, DIMENSION(:), INTENT(IN) :: distribution, subarray
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%geometry = c_geometry_cartesian
    b%ndims = ndims

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    ENDDO

    gmn(1) = xmin
    gmx(1) = xmax

    CALL MPI_ALLREDUCE(gmn, b%extents, ndims, MPI_REAL8, MPI_MIN, &
        h%comm, errcode)
    CALL MPI_ALLREDUCE(gmx, b%extents(ndims+1), ndims, MPI_REAL8, MPI_MAX, &
        h%comm, errcode)

    ! Write header

    CALL write_mesh_meta(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray(1), &
        distribution(1), "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, x, 1, subarray(1), &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_1d_mesh



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_2d_mesh(h, id, name, x, y, dims, xmin, xmax, &
      ymin, ymax, distribution, subarray, dim_labels, dim_units, &
      dim_mults)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(num), DIMENSION(:), INTENT(IN) :: x, y
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    REAL(num), INTENT(IN) :: xmin, xmax, ymin, ymax
    INTEGER, DIMENSION(:), INTENT(IN) :: distribution, subarray
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%geometry = c_geometry_cartesian
    b%ndims = ndims

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    ENDDO

    gmn(1) = xmin
    gmx(1) = xmax
    gmn(2) = ymin
    gmx(2) = ymax

    CALL MPI_ALLREDUCE(gmn, b%extents, ndims, MPI_REAL8, MPI_MIN, &
        h%comm, errcode)
    CALL MPI_ALLREDUCE(gmx, b%extents(ndims+1), ndims, MPI_REAL8, MPI_MAX, &
        h%comm, errcode)

    ! Write header

    CALL write_mesh_meta(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray(1), &
        distribution(1), "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, x, 1, subarray(1), &
        MPI_STATUS_IGNORE, errcode)

    h%current_location = h%current_location + b%dims(1) * b%type_size

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray(2), &
        distribution(2), "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, y, 1, subarray(2), &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_2d_mesh



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_3d_mesh(h, id, name, x, y, z, dims, xmin, xmax, &
      ymin, ymax, zmin, zmax, distribution, subarray, dim_labels, dim_units, &
      dim_mults)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(num), DIMENSION(:), INTENT(IN) :: x, y, z
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    REAL(num), INTENT(IN) :: xmin, xmax, ymin, ymax, zmin, zmax
    INTEGER, DIMENSION(:), INTENT(IN) :: distribution, subarray
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%geometry = c_geometry_cartesian
    b%ndims = ndims

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    ENDDO

    gmn(1) = xmin
    gmx(1) = xmax
    gmn(2) = ymin
    gmx(2) = ymax
    gmn(3) = zmin
    gmx(3) = zmax

    CALL MPI_ALLREDUCE(gmn, b%extents, ndims, MPI_REAL8, MPI_MIN, &
        h%comm, errcode)
    CALL MPI_ALLREDUCE(gmx, b%extents(ndims+1), ndims, MPI_REAL8, MPI_MAX, &
        h%comm, errcode)

    ! Write header

    CALL write_mesh_meta(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray(1), &
        distribution(1), "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, x, 1, subarray(1), &
        MPI_STATUS_IGNORE, errcode)

    h%current_location = h%current_location + b%dims(1) * b%type_size

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray(2), &
        distribution(2), "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, y, 1, subarray(2), &
        MPI_STATUS_IGNORE, errcode)

    h%current_location = h%current_location + b%dims(2) * b%type_size

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray(3), &
        distribution(3), "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, z, 1, subarray(3), &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_3d_mesh



  SUBROUTINE write_mesh_variable_meta(h, id, name, units, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name, units, mesh_id
    REAL(num), INTENT(IN), OPTIONAL :: mult
    INTEGER :: ndims
    TYPE(sdf_block_type), POINTER :: b
    INTEGER :: i, errcode

    b => h%current_block

    b%blocktype = c_blocktype_plain_variable
    ndims = b%ndims

    b%nelements = 1
    DO i = 1,ndims
      b%nelements = b%nelements * b%dims(i)
    ENDDO

    ! Metadata is
    ! - mult      REAL(r8)
    ! - units     CHARACTER(id_length)
    ! - meshid    CHARACTER(id_length)
    ! - dims      INTEGER(i4), DIMENSION(ndims)
    ! - stagger   INTEGER(i4)

    b%info_length = h%block_header_length + (ndims + 1) * soi4 + sof8 &
        + 2 * c_id_length
    b%data_length = b%nelements * b%type_size

    ! Write header

    IF (PRESENT(id)) THEN
      CALL safe_copy_string(units, b%units)
      CALL safe_copy_string(mesh_id, b%mesh_id)

      IF (PRESENT(mult)) THEN
        b%mult = REAL(mult,r8)
      ELSE
        b%mult = 1.d0
      ENDIF

      CALL sdf_write_block_header(h, id, name)
    ELSE
      CALL write_block_header(h)
    ENDIF

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_WRITE(h%filehandle, b%mult, 1, MPI_REAL8, &
          MPI_STATUS_IGNORE, errcode)

      CALL sdf_safe_write_id(h, b%units)

      CALL sdf_safe_write_id(h, b%mesh_id)

      CALL MPI_FILE_WRITE(h%filehandle, b%dims, ndims, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%stagger, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_location = b%block_start + b%info_length
    b%done_info = .TRUE.

  END SUBROUTINE write_mesh_variable_meta



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_1d_float_gen(h, id, name, units, ndims, dims, stagger, &
      mesh_id, variable, distribution, subarray, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, INTENT(IN) :: ndims
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), DIMENSION(1), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%ndims = ndims
    b%stagger = stagger

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    ENDDO

    ! Write header

    CALL write_mesh_variable_meta(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_1d_float_gen



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_2d_float_gen(h, id, name, units, ndims, dims, stagger, &
      mesh_id, variable, distribution, subarray, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, INTENT(IN) :: ndims
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), DIMENSION(1,1), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%ndims = ndims
    b%stagger = stagger

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    ENDDO

    ! Write header

    CALL write_mesh_variable_meta(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_2d_float_gen



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_3d_float_gen(h, id, name, units, ndims, dims, stagger, &
      mesh_id, variable, distribution, subarray, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, INTENT(IN) :: ndims
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), DIMENSION(1,1,1), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%ndims = ndims
    b%stagger = stagger

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    ENDDO

    ! Write header

    CALL write_mesh_variable_meta(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_3d_float_gen



  !----------------------------------------------------------------------------
  ! Code to write a 4D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_4d_float_gen(h, id, name, units, ndims, dims, stagger, &
      mesh_id, variable, distribution, subarray, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, INTENT(IN) :: ndims
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), DIMENSION(1,1,1,1), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%ndims = ndims
    b%stagger = stagger

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    ENDDO

    ! Write header

    CALL write_mesh_variable_meta(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_4d_float_gen



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_1d_float(h, id, name, units, dims, stagger, mesh_id, &
      variable, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult

    CALL sdf_write_1d_float_gen(h, id, name, units, ndims, dims, stagger, &
        mesh_id, variable, distribution, subarray, mult)

  END SUBROUTINE sdf_write_1d_float



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_2d_float(h, id, name, units, dims, stagger, mesh_id, &
      variable, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult

    CALL sdf_write_2d_float_gen(h, id, name, units, ndims, dims, stagger, &
        mesh_id, variable, distribution, subarray, mult)

  END SUBROUTINE sdf_write_2d_float



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_3d_float(h, id, name, units, dims, stagger, mesh_id, &
      variable, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult

    CALL sdf_write_3d_float_gen(h, id, name, units, ndims, dims, stagger, &
        mesh_id, variable, distribution, subarray, mult)

  END SUBROUTINE sdf_write_3d_float



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian integer variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_1d_integer(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER, DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%soi
    b%datatype = h%datatype_integer
    b%mpitype = h%mpitype_integer
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    ENDIF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    ENDDO

    ! Write header

    CALL write_mesh_variable_meta(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_1d_integer



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian integer variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_2d_integer(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER, DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%soi
    b%datatype = h%datatype_integer
    b%mpitype = h%mpitype_integer
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    ENDIF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    ENDDO

    ! Write header

    CALL write_mesh_variable_meta(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_2d_integer



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian integer variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_3d_integer(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER, DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%soi
    b%datatype = h%datatype_integer
    b%mpitype = h%mpitype_integer
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    ENDIF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    ENDDO

    ! Write header

    CALL write_mesh_variable_meta(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_3d_integer



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian character variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_1d_character(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 1
    b%datatype = c_datatype_character
    b%mpitype = MPI_CHARACTER
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    ENDIF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    ENDDO

    ! Write header

    CALL write_mesh_variable_meta(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_1d_character



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian character variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_2d_character(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 1
    b%datatype = c_datatype_character
    b%mpitype = MPI_CHARACTER
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    ENDIF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    ENDDO

    ! Write header

    CALL write_mesh_variable_meta(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_2d_character



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian character variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_3d_character(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 1
    b%datatype = c_datatype_character
    b%mpitype = MPI_CHARACTER
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    ENDIF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    ENDDO

    ! Write header

    CALL write_mesh_variable_meta(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, "native", MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_3d_character



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian variable component (with material index first)
  ! in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_1d_var_first(h, id, name, units, dims, stagger, &
      mesh_id, variable, idx, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), INTENT(IN) :: variable(dims(ndims+1),dims(1))
    INTEGER, INTENT(IN) :: idx, distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult

    CALL sdf_write_2d_float_gen(h, id, name, units, ndims, dims, stagger, &
        mesh_id, variable(idx,1), distribution, subarray, mult)

  END SUBROUTINE sdf_write_1d_var_first



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian variable component (with material index first)
  ! in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_2d_var_first(h, id, name, units, dims, stagger, &
      mesh_id, variable, idx, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), INTENT(IN) :: variable(dims(ndims+1),dims(1),dims(2))
    INTEGER, INTENT(IN) :: idx, distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult

    CALL sdf_write_3d_float_gen(h, id, name, units, ndims, dims, stagger, &
        mesh_id, variable(idx,1,1), distribution, subarray, mult)

  END SUBROUTINE sdf_write_2d_var_first



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian variable component (with material index first)
  ! in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_3d_var_first(h, id, name, units, dims, stagger, &
      mesh_id, variable, idx, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), INTENT(IN) :: variable(dims(ndims+1),dims(1),dims(2),dims(3))
    INTEGER, INTENT(IN) :: idx, distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult

    CALL sdf_write_4d_float_gen(h, id, name, units, ndims, dims, stagger, &
        mesh_id, variable(idx,1,1,1), distribution, subarray, mult)

  END SUBROUTINE sdf_write_3d_var_first



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian variable component (with material index last)
  ! in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_1d_var_last(h, id, name, units, dims, stagger, &
      mesh_id, variable, idx, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), INTENT(IN) :: variable(dims(1),dims(ndims+1))
    INTEGER, INTENT(IN) :: idx, distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult

    CALL sdf_write_2d_float_gen(h, id, name, units, ndims, dims, stagger, &
        mesh_id, variable(1,idx), distribution, subarray, mult)

  END SUBROUTINE sdf_write_1d_var_last



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian variable component (with material index last)
  ! in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_2d_var_last(h, id, name, units, dims, stagger, &
      mesh_id, variable, idx, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), INTENT(IN) :: variable(dims(1),dims(2),dims(ndims+1))
    INTEGER, INTENT(IN) :: idx, distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult

    CALL sdf_write_3d_float_gen(h, id, name, units, ndims, dims, stagger, &
        mesh_id, variable(1,1,idx), distribution, subarray, mult)

  END SUBROUTINE sdf_write_2d_var_last



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian variable component (with material index last)
  ! in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_3d_var_last(h, id, name, units, dims, stagger, &
      mesh_id, variable, idx, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), INTENT(IN) :: variable(dims(1),dims(2),dims(3),dims(ndims+1))
    INTEGER, INTENT(IN) :: idx, distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult

    CALL sdf_write_4d_float_gen(h, id, name, units, ndims, dims, stagger, &
        mesh_id, variable(1,1,1,idx), distribution, subarray, mult)

  END SUBROUTINE sdf_write_3d_var_last



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian material mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_1d_material(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_names, variable, distribution, subarray, mult, last_in)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(num), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    ENDIF

    ALLOCATE(variable_ids(nmat))

    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) .EQ. 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(material_names(i)), variable_ids(i))
      ENDIF
    ENDDO

    CALL sdf_safe_string_composite(h, name, 'All', temp_name)
    CALL sdf_write_stitched_material(h, id, temp_name, mesh_id, stagger, &
        material_names, variable_ids)

    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) .EQ. 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, material_names(i), temp_name)
      IF (last) THEN
        CALL sdf_write_1d_var_last(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ELSE
        CALL sdf_write_1d_var_first(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ENDIF
    ENDDO

    DEALLOCATE(variable_ids)

  END SUBROUTINE sdf_write_1d_material



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian material mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_2d_material(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_names, variable, distribution, subarray, mult, last_in)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(num), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    ENDIF

    ALLOCATE(variable_ids(nmat))

    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) .EQ. 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(material_names(i)), variable_ids(i))
      ENDIF
    ENDDO

    CALL sdf_safe_string_composite(h, name, 'All', temp_name)
    CALL sdf_write_stitched_material(h, id, temp_name, mesh_id, stagger, &
        material_names, variable_ids)

    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) .EQ. 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, material_names(i), temp_name)
      IF (last) THEN
        CALL sdf_write_2d_var_last(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ELSE
        CALL sdf_write_2d_var_first(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ENDIF
    ENDDO

    DEALLOCATE(variable_ids)

  END SUBROUTINE sdf_write_2d_material



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian material mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_3d_material(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_names, variable, distribution, subarray, mult, last_in)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(num), DIMENSION(:,:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    ENDIF

    ALLOCATE(variable_ids(nmat))

    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) .EQ. 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(material_names(i)), variable_ids(i))
      ENDIF
    ENDDO

    CALL sdf_safe_string_composite(h, name, 'All', temp_name)
    CALL sdf_write_stitched_material(h, id, temp_name, mesh_id, stagger, &
        material_names, variable_ids)

    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) .EQ. 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, material_names(i), temp_name)
      IF (last) THEN
        CALL sdf_write_3d_var_last(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ELSE
        CALL sdf_write_3d_var_first(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ENDIF
    ENDDO

    DEALLOCATE(variable_ids)

  END SUBROUTINE sdf_write_3d_material



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian material variable in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_1d_matvar(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_id, material_names, variable, distribution, &
      subarray, mult, last_in)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id, material_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(num), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    ENDIF

    ALLOCATE(variable_ids(nmat))

    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) .EQ. 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(material_names(i)), variable_ids(i))
      ENDIF
    ENDDO

    CALL sdf_safe_string_composite(h, name, 'All', temp_name)
    CALL sdf_write_stitched_matvar(h, id, temp_name, mesh_id, stagger, &
        material_id, variable_ids)

    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) .EQ. 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, material_names(i), temp_name)
      IF (last) THEN
        CALL sdf_write_1d_var_last(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ELSE
        CALL sdf_write_1d_var_first(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ENDIF
    ENDDO

    DEALLOCATE(variable_ids)

  END SUBROUTINE sdf_write_1d_matvar



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian material variable in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_2d_matvar(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_id, material_names, variable, distribution, &
      subarray, mult, last_in)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id, material_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(num), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    ENDIF

    ALLOCATE(variable_ids(nmat))

    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) .EQ. 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(material_names(i)), variable_ids(i))
      ENDIF
    ENDDO

    CALL sdf_safe_string_composite(h, name, 'All', temp_name)
    CALL sdf_write_stitched_matvar(h, id, temp_name, mesh_id, stagger, &
        material_id, variable_ids)

    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) .EQ. 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, material_names(i), temp_name)
      IF (last) THEN
        CALL sdf_write_2d_var_last(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ELSE
        CALL sdf_write_2d_var_first(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ENDIF
    ENDDO

    DEALLOCATE(variable_ids)

  END SUBROUTINE sdf_write_2d_matvar



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian material variable in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_3d_matvar(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_id, material_names, variable, distribution, &
      subarray, mult, last_in)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id, material_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(num), DIMENSION(:,:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    ENDIF

    ALLOCATE(variable_ids(nmat))

    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) .EQ. 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(material_names(i)), variable_ids(i))
      ENDIF
    ENDDO

    CALL sdf_safe_string_composite(h, name, 'All', temp_name)
    CALL sdf_write_stitched_matvar(h, id, temp_name, mesh_id, stagger, &
        material_id, variable_ids)

    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) .EQ. 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, material_names(i), temp_name)
      IF (last) THEN
        CALL sdf_write_3d_var_last(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ELSE
        CALL sdf_write_3d_var_first(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ENDIF
    ENDDO

    DEALLOCATE(variable_ids)

  END SUBROUTINE sdf_write_3d_matvar



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian material species in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_1d_species(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_id, material_name, specnames, variable, distribution, &
      subarray, mult, last_in)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id, material_id
    CHARACTER(LEN=*), INTENT(IN) :: material_name
    CHARACTER(LEN=*), INTENT(IN) :: specnames(:)
    REAL(num), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    ENDIF

    ALLOCATE(variable_ids(nmat))

    DO i = 1,nmat
      IF (LEN_TRIM(specnames(i)) .EQ. 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(specnames(i)), variable_ids(i))
      ENDIF
    ENDDO

    CALL sdf_safe_string_composite(h, name, 'All', temp_name)
    CALL sdf_write_stitched_species(h, id, temp_name, mesh_id, stagger, &
        material_id, material_name, specnames, variable_ids)

    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(specnames(i)) .EQ. 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, specnames(i), temp_name)
      IF (last) THEN
        CALL sdf_write_1d_var_last(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ELSE
        CALL sdf_write_1d_var_first(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ENDIF
    ENDDO

    DEALLOCATE(variable_ids)

  END SUBROUTINE sdf_write_1d_species



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian material species in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_2d_species(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_id, material_name, specnames, variable, distribution, &
      subarray, mult, last_in)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id, material_id
    CHARACTER(LEN=*), INTENT(IN) :: material_name
    CHARACTER(LEN=*), INTENT(IN) :: specnames(:)
    REAL(num), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    ENDIF

    ALLOCATE(variable_ids(nmat))

    DO i = 1,nmat
      IF (LEN_TRIM(specnames(i)) .EQ. 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(specnames(i)), variable_ids(i))
      ENDIF
    ENDDO

    CALL sdf_safe_string_composite(h, name, 'All', temp_name)
    CALL sdf_write_stitched_species(h, id, temp_name, mesh_id, stagger, &
        material_id, material_name, specnames, variable_ids)

    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(specnames(i)) .EQ. 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, specnames(i), temp_name)
      IF (last) THEN
        CALL sdf_write_2d_var_last(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ELSE
        CALL sdf_write_2d_var_first(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ENDIF
    ENDDO

    DEALLOCATE(variable_ids)

  END SUBROUTINE sdf_write_2d_species



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian material species in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_3d_species(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_id, material_name, specnames, variable, distribution, &
      subarray, mult, last_in)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id, material_id
    CHARACTER(LEN=*), INTENT(IN) :: material_name
    CHARACTER(LEN=*), INTENT(IN) :: specnames(:)
    REAL(num), DIMENSION(:,:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(num), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    ENDIF

    ALLOCATE(variable_ids(nmat))

    DO i = 1,nmat
      IF (LEN_TRIM(specnames(i)) .EQ. 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(specnames(i)), variable_ids(i))
      ENDIF
    ENDDO

    CALL sdf_safe_string_composite(h, name, 'All', temp_name)
    CALL sdf_write_stitched_species(h, id, temp_name, mesh_id, stagger, &
        material_id, material_name, specnames, variable_ids)

    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(specnames(i)) .EQ. 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, specnames(i), temp_name)
      IF (last) THEN
        CALL sdf_write_3d_var_last(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ELSE
        CALL sdf_write_3d_var_first(h, variable_ids(i), temp_name, units, &
            dims, stagger, mesh_id, variable, idx, distribution, subarray, mult)
      ENDIF
    ENDDO

    DEALLOCATE(variable_ids)

  END SUBROUTINE sdf_write_3d_species

END MODULE sdf_output_cartesian
