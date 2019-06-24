! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE simple_io

  USE boundary
  USE mpi_subtype_control

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
    REAL(num), DIMENSION(:,:), INTENT(INOUT) :: array
    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: offset
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: subtype, subarray, fh, i, oerr
#ifdef NO_MPI3
    INTEGER :: iasz
#else
    INTEGER(KIND=MPI_COUNT_KIND) :: asz
#endif
    INTEGER(KIND=MPI_OFFSET_KIND) :: sz
    INTEGER(KIND=i8) :: asz8, tsz

    CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_RDONLY, &
        MPI_INFO_NULL, fh, errcode)

    IF (errcode /= 0) THEN
      IF (rank == 0) PRINT *, 'file ', TRIM(filename), ' does not exist.'
      err = IOR(err, c_err_bad_value)
      RETURN
    END IF

    subtype = create_current_field_subtype()
    subarray = create_current_field_subarray(ng)

    ! Check that file size is divisible by the total size of a field array
#ifdef NO_MPI3
    CALL MPI_TYPE_SIZE(subtype, iasz, errcode)
    asz8 = INT(iasz, i8)
#else
    CALL MPI_TYPE_SIZE_X(subarray, asz, errcode)
    asz8 = INT(asz, i8)
#endif
    oerr = errcode
    CALL MPI_REDUCE(asz8, tsz, 1, MPI_INTEGER8, MPI_SUM, 0, comm, errcode)

    IF (rank == 0) THEN
      CALL MPI_FILE_GET_SIZE(fh, sz, errcode)
      IF (oerr == MPI_UNDEFINED) THEN
        PRINT*, '*** WARNING ***'
        PRINT*, 'Cannot automatically test size of file"' // TRIM(filename) &
            // '". Ensure that file is of correct size'
      ELSE
        IF (MOD(sz-offset, tsz) /= 0) THEN
          PRINT*, '*** WARNING ***'
          PRINT*, 'Binary input file "' // TRIM(filename) // '"', &
              ' does not appear to match the domain dimensions'
        END IF
      END IF
    END IF

    CALL MPI_FILE_SET_VIEW(fh, offset, MPI_BYTE, subtype, 'native', &
        MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(fh, array, 1, subarray, MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_CLOSE(fh, errcode)
    CALL MPI_TYPE_FREE(subtype, errcode)

    CALL field_bc(array, ng)
    DO i = 1, 2*c_ndims
      CALL field_zero_gradient(array, c_stagger_centre, i)
    END DO

  END SUBROUTINE load_single_array_from_file



  !----------------------------------------------------------------------------
  ! These subroutines allow loading of a simple 1D array of data from a raw
  ! binary file. The file is assumed to contain complete valid data after a
  ! specified header of length offset.
  !----------------------------------------------------------------------------

  FUNCTION load_1d_real_array(filename, array, offset, err) RESULT(records)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: offset
    REAL(num), DIMENSION(:), POINTER, INTENT(INOUT) :: array
    INTEGER, INTENT(INOUT) :: err
    INTEGER(KIND=MPI_OFFSET_KIND) :: filesize, disp
    INTEGER(KIND=MPI_OFFSET_KIND) :: total_records, remainder, tail
    INTEGER :: fh, typesize, records
    INTEGER :: stat(MPI_STATUS_SIZE)
    INTEGER :: datatype

    datatype = mpireal

    records = 0

    CALL MPI_TYPE_SIZE(datatype, typesize, errcode)

    IF (errcode /= 0) THEN
      IF (rank == 0) THEN
        PRINT *, '*** DEVELOPER WARNING ***'
        PRINT *, 'Unknown MPI_DATATYPE passed'
      END IF
      err = IOR(err, c_err_io_error)
      RETURN
    END IF

    CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_RDONLY, &
        MPI_INFO_NULL, fh, errcode)

    IF (errcode /= 0) THEN
      IF (rank == 0) PRINT *, 'file ', TRIM(filename), ' does not exist.'
      err = IOR(err, c_err_io_error)
      RETURN
    END IF

    CALL MPI_FILE_GET_SIZE(fh, filesize, errcode)

    tail = MOD(INT(filesize, i8) - offset, INT(typesize, i8))
    IF (rank == 0 .AND. tail /= 0) THEN
      PRINT *, '***WARNING***'
      PRINT *, 'Length (less offset) of ', TRIM(filename), &
          ' not an integer multiple of datasize'
      PRINT *, 'Corrupt data?'
    END IF

    total_records = (filesize - offset - tail) / typesize
    records = INT(total_records / nproc, i4)
    remainder = MOD(total_records, INT(nproc, i8))
    IF (rank < remainder) THEN
      records = records + 1
      disp = offset + rank * records * typesize
    ELSE
      disp = offset + (records + 1) * remainder * typesize &
          + records * (rank - remainder) * typesize
    END IF

    ALLOCATE(array(records))

    CALL MPI_FILE_SET_VIEW(fh, disp, datatype, datatype, 'native', &
        MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(fh, array, records, datatype, stat, errcode)

    CALL MPI_GET_COUNT(stat, datatype, records, errcode)

    CALL MPI_FILE_CLOSE(fh, errcode)

  END FUNCTION load_1d_real_array



  FUNCTION load_1d_integer4_array(filename, array, offset, err) RESULT(records)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: offset
    INTEGER(KIND=i4), DIMENSION(:), POINTER, INTENT(INOUT) :: array
    INTEGER, INTENT(INOUT) :: err
    INTEGER(KIND=MPI_OFFSET_KIND) :: filesize, disp
    INTEGER(KIND=MPI_OFFSET_KIND) :: total_records, remainder, tail
    INTEGER :: fh, typesize, records
    INTEGER :: stat(MPI_STATUS_SIZE)
    INTEGER :: datatype = MPI_INTEGER4

    records = 0

    CALL MPI_TYPE_SIZE(datatype, typesize, errcode)

    IF (errcode /= 0) THEN
      IF (rank == 0) THEN
        PRINT *, '*** DEVELOPER WARNING ***'
        PRINT *, 'Unknown MPI_DATATYPE passed'
      END IF
      err = IOR(err, c_err_io_error)
      RETURN
    END IF

    CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_RDONLY, &
        MPI_INFO_NULL, fh, errcode)

    IF (errcode /= 0) THEN
      IF (rank == 0) PRINT *, 'file ', TRIM(filename), ' does not exist.'
      err = IOR(err, c_err_io_error)
      RETURN
    END IF

    CALL MPI_FILE_GET_SIZE(fh, filesize, errcode)

    tail = MOD(INT(filesize, i8) - offset, INT(typesize, i8))
    IF (rank == 0 .AND. tail /= 0) THEN
      PRINT *, '***WARNING***'
      PRINT *, 'Length (less offset) of ', TRIM(filename), &
          ' not an integer multiple of datasize'
      PRINT *, 'Corrupt data?'
    END IF

    total_records = (filesize - offset - tail) / typesize
    records = INT(total_records / nproc, i4)
    remainder = MOD(total_records, INT(nproc, i8))
    IF (rank < remainder) THEN
      records = records + 1
      disp = offset + rank * records * typesize
    ELSE
      disp = offset + (records + 1) * remainder * typesize &
          + records * (rank - remainder) * typesize
    END IF

    ALLOCATE(array(records))

    CALL MPI_FILE_SET_VIEW(fh, disp, datatype, datatype, 'native', &
        MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(fh, array, records, datatype, stat, errcode)

    CALL MPI_GET_COUNT(stat, datatype, records, errcode)

    CALL MPI_FILE_CLOSE(fh, errcode)

  END FUNCTION load_1d_integer4_array



  FUNCTION load_1d_integer8_array(filename, array, offset, err) RESULT(records)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN) :: offset
    INTEGER(KIND=i8), DIMENSION(:), POINTER, INTENT(INOUT) :: array
    INTEGER, INTENT(INOUT) :: err
    INTEGER(KIND=MPI_OFFSET_KIND) :: filesize, disp
    INTEGER(KIND=MPI_OFFSET_KIND) :: total_records, remainder, tail
    INTEGER :: fh, typesize, records
    INTEGER :: stat(MPI_STATUS_SIZE)
    INTEGER :: datatype = MPI_INTEGER8

    records = 0

    CALL MPI_TYPE_SIZE(datatype, typesize, errcode)

    IF (errcode /= 0) THEN
      IF (rank == 0) THEN
        PRINT *, '*** DEVELOPER WARNING ***'
        PRINT *, 'Unknown MPI_DATATYPE passed'
      END IF
      err = IOR(err, c_err_io_error)
      RETURN
    END IF

    CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_RDONLY, &
        MPI_INFO_NULL, fh, errcode)

    IF (errcode /= 0) THEN
      IF (rank == 0) PRINT *, 'file ', TRIM(filename), ' does not exist.'
      err = IOR(err, c_err_io_error)
      RETURN
    END IF

    CALL MPI_FILE_GET_SIZE(fh, filesize, errcode)

    tail = MOD(INT(filesize, i8) - offset, INT(typesize, i8))
    IF (rank == 0 .AND. tail /= 0) THEN
      PRINT *, '***WARNING***'
      PRINT *, 'Length (less offset) of ', TRIM(filename), &
          ' not an integer multiple of datasize'
      PRINT *, 'Corrupt data?'
    END IF

    total_records = (filesize - offset - tail) / typesize
    records = INT(total_records / nproc, i4)
    remainder = MOD(total_records, INT(nproc, i8))
    IF (rank < remainder) THEN
      records = records + 1
      disp = offset + rank * records * typesize
    ELSE
      disp = offset + (records + 1) * remainder * typesize &
          + records * (rank - remainder) * typesize
    END IF

    ALLOCATE(array(records))

    CALL MPI_FILE_SET_VIEW(fh, disp, datatype, datatype, 'native', &
        MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(fh, array, records, datatype, stat, errcode)

    CALL MPI_GET_COUNT(stat, datatype, records, errcode)

    CALL MPI_FILE_CLOSE(fh, errcode)

  END FUNCTION load_1d_integer8_array

END MODULE simple_io
