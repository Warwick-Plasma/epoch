MODULE simple_io

  USE mpi
  USE mpi_subtype_control
  USE boundary

  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! This subroutine opens a file containing an array the size of the entire
  ! domain (1:nx_global, 1:ny_global) and splits it up onto each processor
  ! (-2:nx+3, -2:nx+3). If there are multiple variables in the file use
  ! offset to specify where to start loading the requested variable from.
  ! Returns errors in an input deck like fashion.
  !----------------------------------------------------------------------------

  SUBROUTINE load_single_array_from_file(filename, array, offset, err)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: offset
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: subtype, subarray, fh, i

    CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_RDONLY, &
        MPI_INFO_NULL, fh, errcode)

    IF (errcode .NE. 0) THEN
      IF (rank .EQ. 0) PRINT *, 'file ', TRIM(filename), ' does not exist.'
      err = IOR(err, c_err_bad_value)
      RETURN
    ENDIF

    subtype = create_current_field_subtype()
    subarray = create_current_field_subarray(ng)
    CALL MPI_FILE_SET_VIEW(fh, offset, MPI_BYTE, subtype, 'native', &
        MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(fh, array, 1, subarray, MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_CLOSE(fh, errcode)
    CALL MPI_TYPE_FREE(subtype, errcode)

    CALL field_bc(array)
    DO i = 1, 2*c_ndims
      CALL field_zero_gradient(array, c_stagger_centre, i)
    ENDDO

  END SUBROUTINE load_single_array_from_file

END MODULE simple_io
