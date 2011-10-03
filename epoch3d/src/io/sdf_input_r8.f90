MODULE sdf_input_r8

  USE mpi
  USE sdf_common
  USE sdf_input_ru

  IMPLICIT NONE

CONTAINS

  SUBROUTINE sdf_read_constant_real(h, value)

    TYPE(sdf_file_handle) :: h
    REAL(num), INTENT(OUT) :: value
    REAL(r4) :: real4
    REAL(r8) :: real8
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_read_constant(h)

    b => h%current_block

    IF (b%datatype .EQ. c_datatype_real4) THEN
      real4 = TRANSFER(b%const_value, real4)
      value = REAL(real4,num)
    ELSE IF (b%datatype .EQ. c_datatype_real8) THEN
      real8 = TRANSFER(b%const_value, real8)
      value = REAL(real8,num)
    ENDIF

  END SUBROUTINE sdf_read_constant_real



  SUBROUTINE sdf_read_1d_array_real(h, values)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:), INTENT(OUT) :: values
    INTEGER, DIMENSION(c_maxdims) :: dims
    INTEGER :: errcode, n1
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_array_info(h, dims)

    h%current_location = b%data_location

    n1 = b%dims(1)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        MPI_BYTE, 'native', MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, values, n1, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_1d_array_real



  SUBROUTINE sdf_read_2d_array_real(h, values)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:,:), INTENT(OUT) :: values
    INTEGER, DIMENSION(c_maxdims) :: dims
    INTEGER :: errcode, i, n1, n2
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_array_info(h, dims)

    h%current_location = b%data_location

    n1 = b%dims(1)
    n2 = b%dims(2)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        MPI_BYTE, 'native', MPI_INFO_NULL, errcode)

    DO i = 1,n2
      CALL MPI_FILE_READ_ALL(h%filehandle, values(1,i), n1, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    ENDDO

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_2d_array_real

END MODULE sdf_input_r8
