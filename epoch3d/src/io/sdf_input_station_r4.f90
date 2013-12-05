MODULE sdf_input_station_r4

  USE sdf_input_station_ru

  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: sof = 4
  INTEGER, PARAMETER, PRIVATE :: datatype_real = c_datatype_real4
  INTEGER, PARAMETER, PRIVATE :: mpitype_real = MPI_REAL4

CONTAINS

  SUBROUTINE read_station_array_r4_r4(h, time, step, array)

    TYPE(sdf_file_handle) :: h
    REAL(r4), INTENT(OUT) :: time
    INTEGER, INTENT(OUT) :: step
    REAL(r4), DIMENSION(:), INTENT(OUT) :: array
    TYPE(sdf_block_type), POINTER :: b
    INTEGER :: errcode
    INTEGER(i8) :: bstart, bstop
    REAL(r8), DIMENSION(:), ALLOCATABLE :: r8array

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_station_info(h)

    bstart = b%data_location
    bstop = bstart + b%data_length
    IF (.NOT. (h%current_location .GT. bstart &
        .AND. h%current_location .LT. bstop)) THEN
      h%current_location = b%data_location
      IF (h%rank .EQ. h%rank_master) THEN
        CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
            errcode)
      ENDIF
    ENDIF

    IF (h%rank .EQ. h%rank_master) THEN
      IF (b%variable_types(1) .EQ. c_datatype_real4) THEN
         CALL MPI_FILE_READ(h%filehandle, time, 1, mpitype_real, &
               MPI_STATUS_IGNORE, errcode)
         CALL MPI_FILE_READ(h%filehandle, array, b%nvariables-1, mpitype_real, &
               MPI_STATUS_IGNORE, errcode)
      ELSE
         ALLOCATE(r8array(b%nvariables))
         CALL MPI_FILE_READ(h%filehandle, r8array, b%nvariables, MPI_REAL8, &
               MPI_STATUS_IGNORE, errcode)
         time = REAL(r8array(1), r4)
         array = REAL(r8array(2:), r4)
         DEALLOCATE(r8array)
      ENDIF
      h%current_location = h%current_location + b%type_size
    ENDIF

    step = h%step
    b%done_data = .TRUE.

  END SUBROUTINE read_station_array_r4_r4



  SUBROUTINE read_station_array_r8_r4(h, time, step, array)

    TYPE(sdf_file_handle) :: h
    REAL(r8), INTENT(OUT) :: time
    INTEGER, INTENT(OUT) :: step
    REAL(r4), DIMENSION(:), INTENT(OUT) :: array
    REAL(r4) :: real4

    IF (sdf_check_block_header(h)) RETURN

    CALL read_station_array_r4_r4(h, real4, step, array)
    time = REAL(real4,r8)

  END SUBROUTINE read_station_array_r8_r4

END MODULE sdf_input_station_r4
