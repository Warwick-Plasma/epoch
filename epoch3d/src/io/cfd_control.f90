MODULE cfd_control

  USE cfd_common
  USE cfd_input
  USE cfd_output
  USE mpi

  IMPLICIT NONE

CONTAINS

  SUBROUTINE cfd_open(filename, cfd_rank_in, cfd_comm_in, mode, step, time, &
      jobid)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: cfd_rank_in, cfd_comm_in, mode
    INTEGER, OPTIONAL, INTENT(INOUT) :: step
    REAL(num), OPTIONAL, INTENT(INOUT) :: time
    TYPE(jobid_type), OPTIONAL, INTENT(IN) :: jobid
    REAL(num) :: dummy = 1.0_num
    INTEGER :: errcode, sof4, ostep = 0
    DOUBLE PRECISION :: otime = 0

    cfd_comm = cfd_comm_in
    cfd_rank = cfd_rank_in
    CALL MPI_SIZEOF(dummy, sof4, errcode)
    sof = sof4

    IF (mode .EQ. c_cfd_write) THEN
      cfd_mode = MPI_MODE_CREATE + MPI_MODE_WRONLY
      cfd_writing = .TRUE.

      ! Creating a new file of the current version, so set the header offset
      ! to reflect current version
      header_offset = header_offset_this_version

      IF (PRESENT(step)) ostep = step

      IF (PRESENT(time)) otime = DBLE(time)

      IF (PRESENT(jobid)) THEN
        cfd_jobid = jobid
      ELSE
        cfd_jobid%start_seconds = 0
        cfd_jobid%start_milliseconds = 0
      ENDIF

      ! We are opening a file to be created, so use the destructive file
      ! opening command
      CALL cfd_open_clobber(filename, ostep, otime)
    ELSE
      cfd_mode = MPI_MODE_RDONLY
      cfd_writing = .FALSE.

      ! We're opening a file which already exists, so don't damage it
      CALL cfd_open_read(filename, ostep, otime)

      IF (PRESENT(step)) step = ostep
      IF (PRESENT(time)) time = otime

    ENDIF

  END SUBROUTINE cfd_open



  SUBROUTINE cfd_close

    INTEGER :: errcode

    ! No open file
    IF (cfd_filehandle .EQ. -1) RETURN

    ! If writing
    IF (cfd_writing) THEN
      ! Go to place where the empty value for nblocks is
      current_displacement = nblocks_offset_this_version
      CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
          MPI_INTEGER4, MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

      IF (cfd_rank .EQ. default_rank) &
          CALL MPI_FILE_WRITE(cfd_filehandle, cfd_nblocks, 1, MPI_INTEGER4, &
              MPI_STATUS_IGNORE, errcode)
    ENDIF

    CALL MPI_BARRIER(comm, errcode)

    CALL MPI_FILE_CLOSE(cfd_filehandle, errcode)

    ! Set cfd_filehandle to -1 to show that the file is closed
    cfd_filehandle = -1

  END SUBROUTINE cfd_close



  SUBROUTINE cfd_set_max_string_length(maxlen)

    INTEGER, INTENT(IN) :: maxlen

    max_string_len = INT(maxlen,4)

  END SUBROUTINE cfd_set_max_string_length



  SUBROUTINE cfd_set_default_rank(rank_in)

    INTEGER, INTENT(IN) :: rank_in

    default_rank = INT(rank_in,4)

  END SUBROUTINE cfd_set_default_rank



  FUNCTION cfd_get_nblocks()

    INTEGER :: cfd_get_nblocks

    cfd_get_nblocks = cfd_nblocks

  END FUNCTION cfd_get_nblocks



  FUNCTION cfd_get_jobid()

    TYPE(jobid_type) :: cfd_get_jobid

    cfd_get_jobid = cfd_jobid

  END FUNCTION cfd_get_jobid

END MODULE cfd_control
