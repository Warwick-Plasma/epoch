MODULE cfd_control

  USE cfd_common
  USE cfd_input
  USE cfd_output
  USE mpi

  IMPLICIT NONE

CONTAINS

  SUBROUTINE cfd_open(h, filename, cfd_rank_in, cfd_comm_in, mode, step, &
      time, jobid)

    TYPE(cfd_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: cfd_rank_in, cfd_comm_in, mode
    INTEGER, OPTIONAL, INTENT(INOUT) :: step
    REAL(num), OPTIONAL, INTENT(INOUT) :: time
    TYPE(jobid_type), OPTIONAL, INTENT(IN) :: jobid
    INTEGER :: errcode, ierr, ostep = 0
    DOUBLE PRECISION :: otime = 0

    h%default_rank = 0
    h%filehandle = -1
    h%max_string_len = 64

    h%comm = cfd_comm_in
    h%rank = cfd_rank_in

    ! We currently only support files written at the same precision
    ! as the 'num' kind
    IF (num .EQ. KIND(1.0)) THEN
      ! Should use MPI_SIZEOF() but this breaks on scalimpi
      sof = 4
      h%mpireal = MPI_REAL
    ELSE IF (num .EQ. KIND(1.D0)) THEN
      sof = 8
      h%mpireal = MPI_DOUBLE_PRECISION
    ELSE
      IF (h%rank .EQ. h%default_rank) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Error writing CFD output file - unknown datatype'
      ENDIF
      CALL MPI_ABORT(h%comm, errcode, ierr)
    ENDIF

    IF (mode .EQ. c_cfd_write) THEN
      h%mode = MPI_MODE_CREATE + MPI_MODE_WRONLY
      h%writing = .TRUE.

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
      CALL cfd_open_clobber(h, filename, ostep, otime)
    ELSE
      h%mode = MPI_MODE_RDONLY
      h%writing = .FALSE.

      ! We're opening a file which already exists, so don't damage it
      CALL cfd_open_read(h, filename, ostep, otime)

      IF (PRESENT(step)) step = ostep
      IF (PRESENT(time)) time = otime
    ENDIF

  END SUBROUTINE cfd_open



  SUBROUTINE cfd_close(h)

    TYPE(cfd_file_handle) :: h
    INTEGER :: errcode

    ! No open file
    IF (h%filehandle .EQ. -1) RETURN

    ! If writing
    IF (h%writing) THEN
      ! Go to place where the empty value for nblocks is
      h%current_displacement = c_nblocks_offset
      CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
          MPI_INTEGER4, MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

      IF (h%rank .EQ. h%default_rank) &
          CALL MPI_FILE_WRITE(h%filehandle, h%nblocks, 1, MPI_INTEGER4, &
              MPI_STATUS_IGNORE, errcode)
    ENDIF

    CALL MPI_BARRIER(h%comm, errcode)

    CALL MPI_FILE_CLOSE(h%filehandle, errcode)

    ! Set cfd_filehandle to -1 to show that the file is closed
    h%filehandle = -1

  END SUBROUTINE cfd_close



  SUBROUTINE cfd_set_max_string_length(h, maxlen)

    TYPE(cfd_file_handle) :: h
    INTEGER, INTENT(IN) :: maxlen

    h%max_string_len = INT(maxlen,4)

  END SUBROUTINE cfd_set_max_string_length



  SUBROUTINE cfd_set_default_rank(h, rank_in)

    TYPE(cfd_file_handle) :: h
    INTEGER, INTENT(IN) :: rank_in

    h%default_rank = INT(rank_in,4)

  END SUBROUTINE cfd_set_default_rank



  FUNCTION cfd_get_nblocks(h)

    TYPE(cfd_file_handle) :: h
    INTEGER :: cfd_get_nblocks

    cfd_get_nblocks = h%nblocks

  END FUNCTION cfd_get_nblocks



  FUNCTION cfd_get_jobid(h)

    TYPE(cfd_file_handle) :: h
    TYPE(jobid_type) :: cfd_get_jobid

    cfd_get_jobid = cfd_jobid

  END FUNCTION cfd_get_jobid

END MODULE cfd_control
