MODULE sdf_output_station_r4

  USE sdf_output_station_ru

  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: sof = 4
  INTEGER, PARAMETER, PRIVATE :: datatype_real = c_datatype_real4
  INTEGER, PARAMETER, PRIVATE :: mpitype_real = MPI_REAL4

CONTAINS

  SUBROUTINE write_station_array_r8_r4(h, time, step, array)

    TYPE(sdf_file_handle) :: h
    REAL(r8), INTENT(IN) :: time
    INTEGER, INTENT(IN) :: step
    REAL(r4), DIMENSION(:), INTENT(IN) :: array
    TYPE(sdf_block_type), POINTER :: b
    REAL(r4) :: real4
    INTEGER :: errcode

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%print_warnings .AND. h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'SDF block cannot be found. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (b%blocktype .NE. c_blocktype_station) THEN
      IF (h%print_warnings .AND. h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'SDF unable to write station data. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    h%time = time
    h%step = step

    IF (h%rank .EQ. h%rank_master) THEN
      h%current_location = b%data_location + b%data_length
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      real4 = REAL(time,r4)
      CALL MPI_FILE_WRITE(h%filehandle, real4, 1, mpitype_real, &
            MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, array, b%nvariables-1, mpitype_real, &
            MPI_STATUS_IGNORE, errcode)
    ENDIF
    b%nelements = b%nelements + 1
    b%data_length = b%data_length + b%type_size

    CALL write_station_update(h)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_station_array_r8_r4



  SUBROUTINE write_station_array_r4_r4(h, time, step, array)

    TYPE(sdf_file_handle) :: h
    REAL(r4), INTENT(IN) :: time
    INTEGER, INTENT(IN) :: step
    REAL(r4), DIMENSION(:), INTENT(IN) :: array

    CALL write_station_array_r8_r4(h, REAL(time,r8), step, array)

  END SUBROUTINE write_station_array_r4_r4

END MODULE sdf_output_station_r4
