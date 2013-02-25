MODULE sdf_control

  USE sdf_output_util

  IMPLICIT NONE

CONTAINS

  SUBROUTINE sdf_open(h, filename, sdf_comm_in, mode)

    TYPE(sdf_file_handle), TARGET :: h
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: sdf_comm_in, mode
    INTEGER :: errcode, ierr, i

    CALL initialise_file_handle(h)
    CALL sdf_set_default_rank(h, 0)

    h%comm = sdf_comm_in
    CALL MPI_COMM_RANK(sdf_comm_in, h%rank, errcode)

    ierr = KIND(errcode)
    IF (ierr .EQ. i4) THEN
      ! Should use MPI_SIZEOF() but this breaks on scalimpi
      h%soi = 4
      h%datatype_integer = c_datatype_integer4
      h%mpitype_integer = MPI_INTEGER4
    ELSE IF (ierr .EQ. i8) THEN
      h%soi = 8
      h%datatype_integer = c_datatype_integer8
      h%mpitype_integer = MPI_INTEGER8
    ELSE
      h%error_code = c_err_unsupported_datarep + 64 * h%nblocks
      h%handled_error = .TRUE.
      RETURN
    ENDIF

    IF (mode .EQ. c_sdf_write) THEN
      h%writing = .TRUE.
      h%mode = MPI_MODE_CREATE + MPI_MODE_WRONLY

      ! Delete file
      IF (h%rank .EQ. h%rank_master) &
          CALL MPI_FILE_DELETE(TRIM(filename), MPI_INFO_NULL, errcode)
    ELSE
      ! We're opening a file which already exists, so don't damage it
      h%writing = .FALSE.
      h%mode = MPI_MODE_RDONLY
    ENDIF

    CALL MPI_FILE_OPEN(h%comm, TRIM(filename), h%mode, MPI_INFO_NULL, &
        h%filehandle, errcode)
    IF (errcode .NE. 0) h%error_code = map_error_code(errcode)

    IF (h%rank .EQ. h%rank_master .AND. h%filehandle .NE. 0) THEN
      CALL MPI_FILE_CREATE_ERRHANDLER(error_handler, h%errhandler, errcode)
      CALL MPI_FILE_SET_ERRHANDLER(h%filehandle, h%errhandler, errcode)
      DO i = 1, max_handles
        IF (sdf_handles(i)%filehandle .EQ. 0) THEN
          sdf_handles(i)%filehandle = h%filehandle
          sdf_handles(i)%handle => h
          EXIT
        ENDIF
      ENDDO
    ENDIF

  END SUBROUTINE sdf_open



  SUBROUTINE sdf_close(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: errcode

    ! No open file
    IF (h%filehandle .EQ. -1) RETURN

    ! If writing
    IF (h%writing) THEN
      CALL sdf_write_summary(h)

      CALL sdf_flush(h)
    ENDIF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, &
        MPI_BYTE, 'native', MPI_INFO_NULL, errcode)

    CALL MPI_BARRIER(h%comm, errcode)

    CALL MPI_FILE_CLOSE(h%filehandle, errcode)

    CALL sdf_destroy_blocklist(h)
    CALL deallocate_file_handle(h)

  END SUBROUTINE sdf_close



  SUBROUTINE sdf_destroy_block(b)

    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT. ASSOCIATED(b)) RETURN

    CALL deallocate_block_type(b)
    DEALLOCATE(b)

  END SUBROUTINE sdf_destroy_block



  SUBROUTINE sdf_destroy_blocklist(h)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: b, next
    INTEGER :: i

    IF (.NOT.ASSOCIATED(h%blocklist)) RETURN

    b => h%blocklist
    DO i = 1,h%nblocks
      next => b%next_block
      CALL sdf_destroy_block(b)
      b => next
    ENDDO

    NULLIFY(h%blocklist)
    NULLIFY(h%current_block)

  END SUBROUTINE sdf_destroy_blocklist



  SUBROUTINE sdf_set_string_length(h, maxlen)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(IN) :: maxlen

    h%string_length = INT(maxlen,i4)

  END SUBROUTINE sdf_set_string_length



  SUBROUTINE sdf_set_default_rank(h, rank_in)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(IN) :: rank_in

    h%default_rank = INT(rank_in,i4)
    h%rank_master = h%default_rank

  END SUBROUTINE sdf_set_default_rank



  FUNCTION sdf_read_nblocks(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: sdf_read_nblocks

    sdf_read_nblocks = h%nblocks

  END FUNCTION sdf_read_nblocks



  FUNCTION sdf_read_jobid(h)

    TYPE(sdf_file_handle) :: h
    TYPE(jobid_type) :: sdf_read_jobid

    sdf_read_jobid = h%jobid

  END FUNCTION sdf_read_jobid



  FUNCTION sdf_errorcode(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: sdf_errorcode

    IF (h%handled_error) THEN
      sdf_errorcode = h%error_code
    ELSE
      sdf_errorcode = c_err_success
    ENDIF

  END FUNCTION sdf_errorcode

END MODULE sdf_control
