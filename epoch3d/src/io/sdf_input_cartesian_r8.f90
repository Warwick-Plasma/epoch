MODULE sdf_input_cartesian_r8

  USE mpi
  USE sdf_common
  USE sdf_input
  USE sdf_input_cartesian_ru

  IMPLICIT NONE

CONTAINS

  SUBROUTINE sdf_read_srl_1d_mesh(h, x)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:), INTENT(OUT) :: x
    INTEGER :: errcode, intn
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_plain_mesh_info(h)

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        MPI_BYTE, 'native', MPI_INFO_NULL, errcode)

    intn = b%dims(1)
    CALL MPI_FILE_READ_ALL(h%filehandle, x, intn, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    ! That should be it, so now skip to end of block
    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_srl_1d_mesh



  SUBROUTINE sdf_read_srl_2d_mesh(h, x, y)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:), INTENT(OUT) :: x, y
    INTEGER :: errcode, intn
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_plain_mesh_info(h)

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        MPI_BYTE, 'native', MPI_INFO_NULL, errcode)

    intn = b%dims(1)
    CALL MPI_FILE_READ_ALL(h%filehandle, x, intn, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)
    intn = b%dims(2)
    CALL MPI_FILE_READ_ALL(h%filehandle, y, intn, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    ! That should be it, so now skip to end of block
    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_srl_2d_mesh



  SUBROUTINE sdf_read_srl_3d_mesh(h, x, y, z)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:), INTENT(OUT) :: x, y, z
    INTEGER :: errcode, intn
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_plain_mesh_info(h)

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        MPI_BYTE, 'native', MPI_INFO_NULL, errcode)

    intn = b%dims(1)
    CALL MPI_FILE_READ_ALL(h%filehandle, x, intn, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)
    intn = b%dims(2)
    CALL MPI_FILE_READ_ALL(h%filehandle, y, intn, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)
    intn = b%dims(3)
    CALL MPI_FILE_READ_ALL(h%filehandle, z, intn, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    ! That should be it, so now skip to end of block
    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_srl_3d_mesh



  !----------------------------------------------------------------------------
  ! Code to read a 1D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_1d_float(h, variable, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:), INTENT(OUT) :: variable
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

  END SUBROUTINE sdf_read_1d_float



  !----------------------------------------------------------------------------
  ! Code to read a 2D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_2d_float(h, variable, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(1,1), INTENT(OUT) :: variable
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

  END SUBROUTINE sdf_read_2d_float



  !----------------------------------------------------------------------------
  ! Code to read a 3D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_3d_float(h, variable, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(1,1,1), INTENT(OUT) :: variable
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

  END SUBROUTINE sdf_read_3d_float



  !----------------------------------------------------------------------------
  ! Code to read a 4D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_4d_float(h, variable, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(1,1,1,1), INTENT(OUT) :: variable
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

  END SUBROUTINE sdf_read_4d_float



  !----------------------------------------------------------------------------
  ! Code to read a 1D cartesian variable component (with material index first)
  ! in parallel using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_1d_var_first(h, nm, dims, variable, idx, distribution, &
      subarray)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: dims(:)
    REAL(num), INTENT(OUT) :: variable(nm, dims(1))
    INTEGER, INTENT(IN) :: idx, distribution, subarray

    CALL sdf_read_2d_float(h, variable(idx,1), distribution, subarray)

  END SUBROUTINE sdf_read_1d_var_first



  !----------------------------------------------------------------------------
  ! Code to read a 2D cartesian variable component (with material index first)
  ! in parallel using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_2d_var_first(h, nm, dims, variable, idx, distribution, &
      subarray)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: dims(:)
    REAL(num), INTENT(OUT) :: variable(nm, dims(1), dims(2))
    INTEGER, INTENT(IN) :: idx, distribution, subarray

    CALL sdf_read_3d_float(h, variable(idx,1,1), distribution, subarray)

  END SUBROUTINE sdf_read_2d_var_first



  !----------------------------------------------------------------------------
  ! Code to read a 3D cartesian variable component (with material index first)
  ! in parallel using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_3d_var_first(h, nm, dims, variable, idx, distribution, &
      subarray)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: dims(:)
    REAL(num), INTENT(OUT) :: variable(nm, dims(1), dims(2), dims(3))
    INTEGER, INTENT(IN) :: idx, distribution, subarray

    CALL sdf_read_4d_float(h, variable(idx,1,1,1), distribution, subarray)

  END SUBROUTINE sdf_read_3d_var_first



  !----------------------------------------------------------------------------
  ! Code to read a 3D cartesian variable component (with material index last)
  ! in parallel using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_1d_var_last(h, nm, dims, variable, idx, distribution, &
      subarray)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: dims(:)
    REAL(num), INTENT(OUT) :: variable(dims(1), nm)
    INTEGER, INTENT(IN) :: idx, distribution, subarray

    CALL sdf_read_2d_float(h, variable(1,idx), distribution, subarray)

  END SUBROUTINE sdf_read_1d_var_last



  !----------------------------------------------------------------------------
  ! Code to read a 2D cartesian variable component (with material index last)
  ! in parallel using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_2d_var_last(h, nm, dims, variable, idx, distribution, &
      subarray)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: dims(:)
    REAL(num), INTENT(OUT) :: variable(dims(1), dims(2), nm)
    INTEGER, INTENT(IN) :: idx, distribution, subarray

    CALL sdf_read_3d_float(h, variable(1,1,idx), distribution, subarray)

  END SUBROUTINE sdf_read_2d_var_last



  !----------------------------------------------------------------------------
  ! Code to read a 3D cartesian variable component (with material index last)
  ! in parallel using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_3d_var_last(h, nm, dims, variable, idx, distribution, &
      subarray)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: dims(:)
    REAL(num), INTENT(OUT) :: variable(dims(1), dims(2), dims(3), nm)
    INTEGER, INTENT(IN) :: idx, distribution, subarray

    CALL sdf_read_4d_float(h, variable(1,1,1,idx), distribution, subarray)

  END SUBROUTINE sdf_read_3d_var_last



  !----------------------------------------------------------------------------
  ! Code to read a 1D cartesian material mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_1d_material(h, variable, distribution, subarray, last_in)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:,:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: last_in
    INTEGER :: i, nm
    LOGICAL :: last, found
    TYPE(sdf_block_type), POINTER :: cur, b

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    ENDIF

    cur => h%current_block
    nm = cur%ndims

    IF (last) THEN
      DO i = 1,nm
        found = sdf_find_block(h, b, cur%variable_ids(i))
        IF (.NOT. found) RETURN
        h%current_block => b
        CALL sdf_read_1d_var_last(h, nm, b%dims, variable, i, distribution, &
            subarray)
      ENDDO
    ELSE
      DO i = 1,nm
        found = sdf_find_block(h, b, cur%variable_ids(i))
        IF (.NOT. found) RETURN
        h%current_block => b
        CALL sdf_read_1d_var_first(h, nm, b%dims, variable, i, distribution, &
            subarray)
      ENDDO
    ENDIF

    h%current_block => cur

  END SUBROUTINE sdf_read_1d_material



  !----------------------------------------------------------------------------
  ! Code to read a 2D cartesian material mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_2d_material(h, variable, distribution, subarray, last_in)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:,:,:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: last_in
    INTEGER :: i, nm
    LOGICAL :: last, found
    TYPE(sdf_block_type), POINTER :: cur, b

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    ENDIF

    cur => h%current_block
    nm = cur%ndims

    IF (last) THEN
      DO i = 1,nm
        found = sdf_find_block(h, b, cur%variable_ids(i))
        IF (.NOT. found) RETURN
        h%current_block => b
        CALL sdf_read_2d_var_last(h, nm, b%dims, variable, i, distribution, &
            subarray)
      ENDDO
    ELSE
      DO i = 1,nm
        found = sdf_find_block(h, b, cur%variable_ids(i))
        IF (.NOT. found) RETURN
        h%current_block => b
        CALL sdf_read_2d_var_first(h, nm, b%dims, variable, i, distribution, &
            subarray)
      ENDDO
    ENDIF

    h%current_block => cur

  END SUBROUTINE sdf_read_2d_material



  !----------------------------------------------------------------------------
  ! Code to read a 3D cartesian material mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_3d_material(h, variable, distribution, subarray, last_in)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:,:,:,:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: last_in
    INTEGER :: i, nm
    LOGICAL :: last, found
    TYPE(sdf_block_type), POINTER :: cur, b

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    ENDIF

    cur => h%current_block
    nm = cur%ndims

    IF (last) THEN
      DO i = 1,nm
        found = sdf_find_block(h, b, cur%variable_ids(i))
        IF (.NOT. found) RETURN
        h%current_block => b
        CALL sdf_read_3d_var_last(h, nm, b%dims, variable, i, distribution, &
            subarray)
      ENDDO
    ELSE
      DO i = 1,nm
        found = sdf_find_block(h, b, cur%variable_ids(i))
        IF (.NOT. found) RETURN
        h%current_block => b
        CALL sdf_read_3d_var_first(h, nm, b%dims, variable, i, distribution, &
            subarray)
      ENDDO
    ENDIF

    h%current_block => cur

  END SUBROUTINE sdf_read_3d_material



  !----------------------------------------------------------------------------
  ! Code to read either 2D cartesian variable or 1D cartesian multi-material
  ! in parallel
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_2d_variable(h, variable, distribution, subarray, last_in)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:,:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: last_in
    TYPE(sdf_block_type), POINTER :: b

    b => h%current_block
    IF (b%blocktype .EQ. c_blocktype_plain_variable) THEN
      CALL sdf_read_2d_float(h, variable, distribution, subarray)
    ELSE
      CALL sdf_read_1d_material(h, variable, distribution, subarray, last_in)
    ENDIF

  END SUBROUTINE sdf_read_2d_variable



  !----------------------------------------------------------------------------
  ! Code to read either 3D cartesian variable or 2D cartesian multi-material
  ! in parallel
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_read_3d_variable(h, variable, distribution, subarray, last_in)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:,:,:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: last_in
    TYPE(sdf_block_type), POINTER :: b

    b => h%current_block
    IF (b%blocktype .EQ. c_blocktype_plain_variable) THEN
      CALL sdf_read_3d_float(h, variable, distribution, subarray)
    ELSE
      CALL sdf_read_2d_material(h, variable, distribution, subarray, last_in)
    ENDIF

  END SUBROUTINE sdf_read_3d_variable

END MODULE sdf_input_cartesian_r8
