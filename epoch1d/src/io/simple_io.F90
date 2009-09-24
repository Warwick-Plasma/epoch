MODULE simple_io

  USE shared_data
  USE mpi_subtype_control
  USE boundary

  IMPLICIT NONE

CONTAINS

  !-----------------------------------------------------------------------
  !This subroutine opens a file containing an array the size of the entire
  !domain (1:nx_global,1:ny_global) and splits it up onto each processor
  !(-2:nx+3,-2:nx+3). If there are multiple variables in the file use
  !offset to specify where to start loading the requested variable from.
  !Returns errors in an input deck like fashion.
  !-----------------------------------------------------------------------
  SUBROUTINE Load_Single_Array_From_Data_File(FileName,Array,offset,ERR)
    CHARACTER(LEN=*),INTENT(IN) :: FileName
    REAL(num),DIMENSION(-2:),INTENT(INOUT) :: Array
    INTEGER(KIND=MPI_OFFSET_KIND),INTENT(IN) :: offset
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER :: subtype, fh

    CALL MPI_FILE_OPEN(comm,TRIM(Filename),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,errcode)
    IF (errcode .NE. 0) THEN
       IF (rank .EQ. 0) PRINT *,"File ",TRIM(FileName), " does not exist."
       ERR=IOR(ERR,ERR_BAD_VALUE)
       RETURN
    ENDIF
    subtype = Create_Current_Field_Subtype()
    CALL MPI_FILE_SET_VIEW(fh,offset,mpireal,subtype,"native",MPI_INFO_NULL,errcode)
    CALL MPI_FILE_READ_ALL(fh,Array(1:nx),nx,mpireal,status,errcode)
    CALL MPI_FILE_CLOSE(fh,errcode)
    CALL MPI_TYPE_FREE(subtype,errcode)

    CALL Field_BC(Array)
    CALL Field_Zero_Gradient(Array,.TRUE.)

  END SUBROUTINE Load_Single_Array_From_Data_File

END MODULE simple_io
