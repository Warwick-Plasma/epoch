MODULE iocontrol

  USE shared_data
  USE iocommon
  USE input
  USE output

  IMPLICIT NONE

CONTAINS

  SUBROUTINE cfd_Open(filename,cfd_rank_in,cfd_comm_in,mode)

    CHARACTER(len=*),INTENT(IN) :: filename
    INTEGER,INTENT(IN) ::cfd_comm_in,cfd_rank_in,mode

    cfd_comm=cfd_comm_in
    cfd_rank=cfd_rank_in
    cfd_mode=mode

    cfd_writing=IOR(IAND(mode,MPI_MODE_RDWR) ,IAND(mode,MPI_MODE_WRONLY)) .NE. 0
    cfd_reading=IOR(IAND(mode,MPI_MODE_RDWR) ,IAND(mode,MPI_MODE_RDONLY)) .NE. 0

    IF (IAND(mode,MPI_MODE_CREATE) .NE. 0) THEN
      !Creating a new file of the current version, so set the header offset to reflect current version
      header_offset = header_offset_this_version
      !We are opening a file to be created, so use the destructive file opening command
      CALL cfd_Open_Clobber(filename)
    ELSE
      !We're opening a file which already exists, so don't damage it
      CALL cfd_Open_Read(filename)
    ENDIF

  END SUBROUTINE

  SUBROUTINE cfd_Close

    INTEGER :: a

    !No open file
    IF (cfd_filehandle == -1) RETURN

    !If writing
    IF (cfd_writing) THEN
    !Go to place where the empty value for nblocks is
       current_displacement=header_offset-4
       CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER,MPI_INTEGER,&
            "native", MPI_INFO_NULL, cfd_errcode)

       IF (cfd_rank == default_rank) &
            CALL MPI_FILE_WRITE(cfd_filehandle,nBlocks,1,MPI_INTEGER,cfd_status,cfd_errcode)
    ENDIF
    CALL MPI_BARRIER(comm,cfd_errcode)

    CALL MPI_FILE_CLOSE(cfd_filehandle,cfd_errcode)

    !Set cfd_filehandle to -1 to show that the file is closed
    cfd_filehandle=-1

  END SUBROUTINE cfd_Close

  SUBROUTINE cfd_Set_Max_String_Length(maxlen)

    INTEGER, INTENT(IN) :: maxlen

    MaxStringLen=maxlen

  END SUBROUTINE cfd_Set_Max_String_Length

  SUBROUTINE cfd_Set_Default_Rank(rank_in)

    INTEGER, INTENT(IN) :: rank_in

    default_rank=rank_in

  END SUBROUTINE cfd_Set_Default_Rank

  FUNCTION cfd_Get_nBlocks()

    INTEGER :: cfd_Get_nBlocks

    cfd_Get_nBlocks=nBlocks

  END FUNCTION

END MODULE iocontrol
