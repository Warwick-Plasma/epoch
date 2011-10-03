MODULE sdf_input_cartesian_ru

  USE mpi
  USE sdf_common
  USE sdf_input

  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! Code to read common mesh info
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_info_init(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT. ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    h%current_location = b%block_start + h%block_header_length

    IF (.NOT. ASSOCIATED(h%buffer)) THEN
      CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
          MPI_BYTE, 'native', MPI_INFO_NULL, errcode)
    ENDIF

  END SUBROUTINE sdf_info_init



  ! Mesh loading functions

  SUBROUTINE sdf_read_plain_mesh_info(h, geometry, dims, extents, dim_labels, &
      dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(OUT), OPTIONAL :: geometry
    INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: dims
    REAL(num), DIMENSION(:), INTENT(OUT), OPTIONAL :: extents
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(OUT), OPTIONAL :: dim_mults
    INTEGER :: i, clen
    TYPE(sdf_block_type), POINTER :: b

    ! Metadata is
    ! - mults     REAL(r8), DIMENSION(ndims)
    ! - labels    CHARACTER(id_length), DIMENSION(ndims)
    ! - units     CHARACTER(id_length), DIMENSION(ndims)
    ! - geometry  INTEGER(i4)
    ! - minval    REAL(r8), DIMENSION(ndims)
    ! - maxval    REAL(r8), DIMENSION(ndims)
    ! - dims      INTEGER(i4), DIMENSION(ndims)

    CALL sdf_info_init(h)

    b => h%current_block

    IF (.NOT. b%done_info) THEN
      ALLOCATE(b%dim_labels(b%ndims))
      ALLOCATE(b%dim_units(b%ndims))
      ALLOCATE(b%dim_mults(b%ndims))

      CALL read_entry_array_real8(h, b%dim_mults, b%ndims)

      DO i = 1,b%ndims
        CALL read_entry_id(h, b%dim_labels(i))
      ENDDO

      DO i = 1,b%ndims
        CALL read_entry_id(h, b%dim_units(i))
      ENDDO

      CALL read_entry_int4(h, b%geometry)

      CALL read_entry_array_real8(h, b%extents, 2*b%ndims)

      CALL read_entry_array_int4(h, b%dims, b%ndims)
    ENDIF

    IF (PRESENT(geometry)) geometry = b%geometry
    IF (PRESENT(dims)) dims(1:b%ndims) = b%dims(1:b%ndims)
    IF (PRESENT(extents)) extents(1:2*b%ndims) = b%extents(1:2*b%ndims)
    IF (PRESENT(dim_mults)) dim_mults = b%dim_mults
    IF (PRESENT(dim_labels)) THEN
      DO i = 1,b%ndims
        clen = MIN(LEN(dim_labels(i)),c_id_length)
        dim_labels(i)(1:clen) = b%dim_labels(i)(1:clen)
      ENDDO
    ENDIF
    IF (PRESENT(dim_units)) THEN
      DO i = 1,b%ndims
        clen = MIN(LEN(dim_units(i)),c_id_length)
        dim_units(i)(1:clen) = b%dim_units(i)(1:clen)
      ENDDO
    ENDIF

    h%current_location = b%block_start + h%block_header_length
    b%done_info = .TRUE.

  END SUBROUTINE sdf_read_plain_mesh_info



  ! Variable loading functions

  SUBROUTINE sdf_read_plain_variable_info(h, dims, units, mesh_id, stagger, &
      mult)

    TYPE(sdf_file_handle) :: h
    INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: dims
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: units, mesh_id
    INTEGER, INTENT(OUT), OPTIONAL :: stagger
    REAL(num), INTENT(OUT), OPTIONAL :: mult
    INTEGER :: clen
    TYPE(sdf_block_type), POINTER :: b

    ! Metadata is
    ! - mult      REAL(r8)
    ! - units     CHARACTER(id_length)
    ! - meshid    CHARACTER(id_length)
    ! - dims      INTEGER(i4), DIMENSION(ndims)
    ! - stagger   INTEGER(i4)

    CALL sdf_info_init(h)

    b => h%current_block
    IF (.NOT. b%done_info) THEN
      CALL read_entry_real8(h, b%mult)

      CALL read_entry_id(h, b%units)

      CALL read_entry_id(h, b%mesh_id)

      CALL read_entry_array_int4(h, b%dims, b%ndims)

      CALL read_entry_int4(h, b%stagger)
    ENDIF

    IF (PRESENT(dims)) dims(1:b%ndims) = b%dims(1:b%ndims)
    IF (PRESENT(stagger)) stagger = b%stagger
    IF (PRESENT(mult)) mult = b%mult
    IF (PRESENT(units)) THEN
      clen = MIN(LEN(units),c_id_length)
      units(1:clen) = b%units(1:clen)
    ENDIF
    IF (PRESENT(mesh_id)) THEN
      clen = MIN(LEN(mesh_id),c_id_length)
      mesh_id(1:clen) = b%mesh_id(1:clen)
    ENDIF

    h%current_location = b%block_start + h%block_header_length
    b%done_info = .TRUE.

  END SUBROUTINE sdf_read_plain_variable_info



  !----------------------------------------------------------------------------
  ! Code to read a 1D cartesian integer variable in parallel
  ! using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_1d_integer(h, variable, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    INTEGER, DIMENSION(:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_plain_variable_info(h)

    ! Read the actual data

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, 'native', MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_1d_integer



  !----------------------------------------------------------------------------
  ! Code to read a 2D cartesian integer variable in parallel
  ! using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_2d_integer(h, variable, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_plain_variable_info(h)

    ! Read the actual data

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, 'native', MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_2d_integer



  !----------------------------------------------------------------------------
  ! Code to read a 3D cartesian integer variable in parallel
  ! using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_3d_integer(h, variable, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_plain_variable_info(h)

    ! Read the actual data

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, 'native', MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_3d_integer



  !----------------------------------------------------------------------------
  ! Code to read a 1D cartesian character variable in parallel
  ! using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_1d_character(h, variable, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=1), DIMENSION(:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_plain_variable_info(h)

    ! Read the actual data

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, 'native', MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_1d_character



  !----------------------------------------------------------------------------
  ! Code to read a 2D cartesian character variable in parallel
  ! using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_2d_character(h, variable, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=1), DIMENSION(:,:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_plain_variable_info(h)

    ! Read the actual data

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, 'native', MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_2d_character



  !----------------------------------------------------------------------------
  ! Code to read a 3D cartesian character variable in parallel
  ! using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_3d_character(h, variable, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=1), DIMENSION(:,:,:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_plain_variable_info(h)

    ! Read the actual data

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, subarray, &
        distribution, 'native', MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_3d_character



  ! Material mesh loading functions
  SUBROUTINE sdf_read_material_info(h, material_names)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), DIMENSION(:), INTENT(OUT) :: material_names
    INTEGER :: iloop
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    CALL sdf_read_stitched_material(h)
    b => h%current_block

    DO iloop = 1, b%ndims
      CALL safe_copy_string(b%material_names(iloop), material_names(iloop))
    ENDDO

  END SUBROUTINE sdf_read_material_info

END MODULE sdf_input_cartesian_ru
