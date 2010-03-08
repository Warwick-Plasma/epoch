MODULE iocontrol

  USE input
  USE output

  IMPLICIT NONE

CONTAINS

  SUBROUTINE cfd_open(filename, cfd_rank_in, cfd_comm_in, mode)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: cfd_comm_in, cfd_rank_in, mode

    cfd_comm = cfd_comm_in
    cfd_rank = cfd_rank_in
    cfd_mode = mode

    cfd_writing = IOR(IAND(mode, MPI_MODE_RDWR), &
        IAND(mode, MPI_MODE_WRONLY)) .NE. 0
    cfd_reading = IOR(IAND(mode, MPI_MODE_RDWR), &
        IAND(mode, MPI_MODE_RDONLY)) .NE. 0

    IF (IAND(mode, MPI_MODE_CREATE) .NE. 0) THEN
      ! Creating a new file of the current version, so set the header offset
      ! to reflect current version
      header_offset = header_offset_this_version

      ! We are opening a file to be created, so use the destructive file
      ! opening command
      CALL cfd_open_clobber(filename)
    ELSE
      ! We're opening a file which already exists, so don't damage it
      CALL cfd_open_read(filename)
    ENDIF

  END SUBROUTINE cfd_open



  SUBROUTINE cfd_close

    ! No open file
    IF (cfd_filehandle .EQ. -1) RETURN

    ! If writing
    IF (cfd_writing) THEN
      ! Go to place where the empty value for nblocks is
      current_displacement = header_offset - 4
      CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
          MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)

      IF (cfd_rank .EQ. default_rank) &
          CALL MPI_FILE_WRITE(cfd_filehandle, nblocks, 1, MPI_INTEGER, &
              cfd_status, cfd_errcode)
    ENDIF

    CALL MPI_BARRIER(comm, cfd_errcode)

    CALL MPI_FILE_CLOSE(cfd_filehandle, cfd_errcode)

    ! Set cfd_filehandle to -1 to show that the file is closed
    cfd_filehandle = -1

  END SUBROUTINE cfd_close



  SUBROUTINE cfd_set_max_string_length(maxlen)

    INTEGER, INTENT(IN) :: maxlen

    max_string_len = maxlen

  END SUBROUTINE cfd_set_max_string_length



  SUBROUTINE cfd_set_default_rank(rank_in)

    INTEGER, INTENT(IN) :: rank_in

    default_rank = rank_in

  END SUBROUTINE cfd_set_default_rank



  FUNCTION cfd_get_nblocks()

    INTEGER :: cfd_get_nblocks

    cfd_get_nblocks = nblocks

  END FUNCTION cfd_get_nblocks

END MODULE iocontrol
