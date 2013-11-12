MODULE sdf_input_station_r8

  USE sdf_input_station_ru

  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: sof = 8
  INTEGER, PARAMETER, PRIVATE :: datatype_real = c_datatype_real8
  INTEGER, PARAMETER, PRIVATE :: mpitype_real = MPI_REAL8

CONTAINS

  SUBROUTINE read_station_array_r8_r8(h, time, step, array)

    TYPE(sdf_file_handle) :: h
    REAL(r8), INTENT(OUT) :: time
    INTEGER, INTENT(OUT) :: step
    REAL(r8), DIMENSION(:), INTENT(OUT) :: array
    TYPE(sdf_block_type), POINTER :: b
    INTEGER :: errcode
    INTEGER(i8) :: bstart, bstop

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_station_info(h)

    bstart = b%data_location
    bstop = bstart + b%data_length
    IF (.NOT. (h%current_location .GE. bstart &
        .AND. h%current_location .LT. bstop)) THEN
      h%current_location = b%data_location
      IF (h%rank .EQ. h%rank_master) THEN
        CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
            errcode)
      ENDIF
    ENDIF

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_READ(h%filehandle, time, 1, mpitype_real, &
            MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_READ(h%filehandle, array, b%nvariables-1, mpitype_real, &
            MPI_STATUS_IGNORE, errcode)
      h%current_location = h%current_location + b%type_size
    ENDIF

    step = h%step
    b%done_data = .TRUE.

  END SUBROUTINE read_station_array_r8_r8



  SUBROUTINE read_station_array_r4_r8(h, time, step, array)

    TYPE(sdf_file_handle) :: h
    REAL(r4), INTENT(OUT) :: time
    INTEGER, INTENT(OUT) :: step
    REAL(r8), DIMENSION(:), INTENT(OUT) :: array
    REAL(r8) :: real8

    IF (sdf_check_block_header(h)) RETURN

    CALL read_station_array_r8_r8(h, real8, step, array)
    time = REAL(real8,r4)

  END SUBROUTINE read_station_array_r4_r8

END MODULE sdf_input_station_r8
