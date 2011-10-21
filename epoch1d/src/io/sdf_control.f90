MODULE sdf_control

  USE mpi
  USE sdf_common
  USE sdf_input
  USE sdf_output
  USE sdf_output_util

  IMPLICIT NONE

CONTAINS

  SUBROUTINE sdf_open(h, filename, sdf_rank_in, sdf_comm_in, mode)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: sdf_rank_in, sdf_comm_in, mode
    INTEGER :: errcode, ierr

    CALL initialise_file_handle(h)
    CALL sdf_set_default_rank(h, 0)

    h%comm = sdf_comm_in
    h%rank = sdf_rank_in

    ! We currently only support files written at the same precision
    ! as the 'num' kind
    IF (num .EQ. r4) THEN
      ! Should use MPI_SIZEOF() but this breaks on scalimpi
      h%sof = 4
      h%datatype_real = c_datatype_real4
      h%mpitype_real = MPI_REAL4
    ELSE IF (num .EQ. r8) THEN
      h%sof = 8
      h%datatype_real = c_datatype_real8
      h%mpitype_real = MPI_REAL8
    ELSE
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Error writing SDF output file - unknown datatype'
      ENDIF
      CALL MPI_ABORT(h%comm, errcode, ierr)
    ENDIF

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
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Error writing SDF output file - unknown datatype'
      ENDIF
      CALL MPI_ABORT(h%comm, errcode, ierr)
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

  END SUBROUTINE sdf_open



  SUBROUTINE sdf_close(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: errcode
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset

    ! No open file
    IF (h%filehandle .EQ. -1) RETURN

    ! If writing
    IF (h%writing) THEN
      CALL sdf_write_summary(h)

      ! Update summary and nblocks info
      IF (h%rank .EQ. h%rank_master) THEN
        offset = c_summary_offset
        CALL MPI_FILE_SEEK(h%filehandle, offset, MPI_SEEK_SET, &
            errcode)
        CALL MPI_FILE_WRITE(h%filehandle, h%summary_location, 1, MPI_INTEGER8, &
            MPI_STATUS_IGNORE, errcode)
        CALL MPI_FILE_WRITE(h%filehandle, h%summary_size, 1, MPI_INTEGER4, &
            MPI_STATUS_IGNORE, errcode)
        CALL MPI_FILE_WRITE(h%filehandle, h%nblocks, 1, MPI_INTEGER4, &
            MPI_STATUS_IGNORE, errcode)
      ENDIF
    ENDIF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, &
        MPI_BYTE, "native", MPI_INFO_NULL, errcode)

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

END MODULE sdf_control
