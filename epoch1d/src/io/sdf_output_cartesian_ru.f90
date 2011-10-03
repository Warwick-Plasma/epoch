MODULE sdf_output_cartesian_ru

  USE mpi
  USE sdf_common
  USE sdf_output

  IMPLICIT NONE

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
        distribution, 'native', MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
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
        distribution, 'native', MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
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
        distribution, 'native', MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
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
        distribution, 'native', MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
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
        distribution, 'native', MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
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
        distribution, 'native', MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_3d_character

END MODULE sdf_output_cartesian_ru
