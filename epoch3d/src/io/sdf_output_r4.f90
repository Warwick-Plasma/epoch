MODULE sdf_output_r4

  USE sdf_output_ru

  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: sof = 4
  INTEGER, PARAMETER, PRIVATE :: datatype_real = c_datatype_real4
  INTEGER, PARAMETER, PRIVATE :: mpitype_real = MPI_REAL4

CONTAINS

  SUBROUTINE write_constant_real_r4(h, id, name, value, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), INTENT(IN) :: value
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = sof
    b%datatype = datatype_real
    b%mpitype = mpitype_real

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    b%const_value(1:sof) = TRANSFER(value, b%const_value(1:sof))

    CALL write_constant_meta(h, id, name)

    h%rank_master = h%default_rank

  END SUBROUTINE write_constant_real_r4



  SUBROUTINE write_1d_array_real_spec_r4(h, id, name, n1, array, rank_write)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER, INTENT(IN) :: n1
    REAL(r4), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = sof
    b%datatype = datatype_real
    b%mpitype = mpitype_real
    b%ndims = ndims

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    b%dims(1) = n1

    ! Write header

    CALL write_array_meta(h, id, name)

    h%current_location = b%data_location

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      ! Actual array
      CALL MPI_FILE_WRITE(h%filehandle, array, n1, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_1d_array_real_spec_r4



  SUBROUTINE write_1d_array_real_r4(h, id, name, array, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: n1

    n1 = SIZE(array,1)
    CALL write_1d_array_real_spec_r4(h, id, name, n1, array, rank_write)

  END SUBROUTINE write_1d_array_real_r4



  SUBROUTINE write_2d_array_real_spec_r4(h, id, name, n1, n2, array, rank_write)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER, INTENT(IN) :: n1, n2
    REAL(r4), DIMENSION(:,:), INTENT(IN) :: array
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: errcode, var_len, i
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = sof
    b%datatype = datatype_real
    b%mpitype = mpitype_real
    b%ndims = ndims

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    b%dims(1) = n1
    b%dims(2) = n2

    ! Write header

    CALL write_array_meta(h, id, name)

    h%current_location = b%data_location

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      ! Actual array
      IF (n1 .EQ. SIZE(array,1)) THEN
        var_len = INT(b%nelements)
        CALL MPI_FILE_WRITE(h%filehandle, array, var_len, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
      ELSE
        DO i = 1,n2
          CALL MPI_FILE_WRITE(h%filehandle, array(1,i), n1, b%mpitype, &
              MPI_STATUS_IGNORE, errcode)
        ENDDO
      ENDIF
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_2d_array_real_spec_r4



  SUBROUTINE write_2d_array_real_r4(h, id, name, array, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:,:), INTENT(IN) :: array
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: n1, n2

    n1 = SIZE(array,1)
    n2 = SIZE(array,2)
    CALL write_2d_array_real_spec_r4(h, id, name, n1, n2, array, rank_write)

  END SUBROUTINE write_2d_array_real_r4

END MODULE sdf_output_r4
