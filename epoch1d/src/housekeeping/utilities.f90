! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
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

MODULE utilities

  USE constants

  IMPLICIT NONE

  INTERFACE grow_array
    MODULE PROCEDURE grow_real_array, grow_integer_array, grow_logical_array, &
                     grow_string_array, grow_real_array2d, grow_integer_array2d
  END INTERFACE grow_array

  PRIVATE
  PUBLIC :: abort_code, grow_array

CONTAINS

  SUBROUTINE grow_real_array(array, idx)

    REAL(num), DIMENSION(:), POINTER :: array
    INTEGER, INTENT(IN) :: idx
    REAL(num), DIMENSION(:), ALLOCATABLE :: tmp_array
    INTEGER :: old_size, new_size, i

    old_size = SIZE(array)
    IF (idx <= old_size) RETURN

    ALLOCATE(tmp_array(old_size))
    DO i = 1, old_size
      tmp_array(i) = array(i)
    END DO

    new_size = 2 * old_size
    IF (new_size < idx) THEN
      new_size = idx + 1
    END IF
    DEALLOCATE(array)
    ALLOCATE(array(new_size))

    DO i = 1, old_size
      array(i) = tmp_array(i)
    END DO

    DEALLOCATE(tmp_array)

  END SUBROUTINE grow_real_array



  SUBROUTINE grow_integer_array(array, idx)

    INTEGER, DIMENSION(:), POINTER :: array
    INTEGER, INTENT(IN) :: idx
    INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_array
    INTEGER :: old_size, new_size, i

    old_size = SIZE(array)
    IF (idx <= old_size) RETURN

    ALLOCATE(tmp_array(old_size))
    DO i = 1, old_size
      tmp_array(i) = array(i)
    END DO

    new_size = 2 * old_size
    IF (new_size < idx) THEN
      new_size = idx + 1
    END IF
    DEALLOCATE(array)
    ALLOCATE(array(new_size))

    DO i = 1, old_size
      array(i) = tmp_array(i)
    END DO

    DEALLOCATE(tmp_array)

  END SUBROUTINE grow_integer_array



  SUBROUTINE grow_logical_array(array, idx)

    LOGICAL, DIMENSION(:), POINTER :: array
    INTEGER, INTENT(IN) :: idx
    LOGICAL, DIMENSION(:), ALLOCATABLE :: tmp_array
    INTEGER :: old_size, new_size, i

    old_size = SIZE(array)
    IF (idx <= old_size) RETURN

    ALLOCATE(tmp_array(old_size))
    DO i = 1, old_size
      tmp_array(i) = array(i)
    ENDDO

    new_size = 2 * old_size
    DEALLOCATE(array)
    ALLOCATE(array(new_size))

    DO i = 1, old_size
      array(i) = tmp_array(i)
    ENDDO

    DEALLOCATE(tmp_array)

  END SUBROUTINE grow_logical_array



  SUBROUTINE grow_string_array(array, idx)

    CHARACTER(LEN=string_length), DIMENSION(:), POINTER :: array
    INTEGER, INTENT(IN) :: idx
    CHARACTER(LEN=string_length), DIMENSION(:), ALLOCATABLE :: tmp_array
    INTEGER :: old_size, new_size, i

    old_size = SIZE(array)
    IF (idx <= old_size) RETURN

    ALLOCATE(tmp_array(old_size))
    DO i = 1, old_size
      tmp_array(i) = array(i)
    END DO

    new_size = 2 * old_size
    IF (new_size < idx) THEN
      new_size = idx + 1
    END IF
    DEALLOCATE(array)
    ALLOCATE(array(new_size))

    DO i = 1, old_size
      array(i) = tmp_array(i)
    END DO

    DEALLOCATE(tmp_array)

  END SUBROUTINE grow_string_array



  SUBROUTINE grow_real_array2d(array, idx, idy)

    REAL(num), DIMENSION(:,:), POINTER :: array
    INTEGER, INTENT(IN) :: idx, idy
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: tmp_array
    INTEGER :: old_size(2), new_size(2), i, j, idxy, idir

    old_size = SHAPE(array)
    IF (idx /= old_size(1)) THEN
      IF (idy /= old_size(2)) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'grow_real_array2d can only grow an array in one dimension'
        RETURN
      END IF
      idir = 1
      idxy = idx
    ELSE
      idir = 2
      idxy = idy
    END IF

    IF (idxy <= old_size(idir)) RETURN

    new_size = old_size
    ALLOCATE(tmp_array(old_size(1),old_size(2)))
    DO j = 1, old_size(2)
    DO i = 1, old_size(1)
      tmp_array(i,j) = array(i,j)
    END DO
    END DO

    new_size(idir) = 2 * old_size(idir)
    IF (new_size(idir) < idxy) THEN
      new_size(idir) = idxy + 1
    END IF
    DEALLOCATE(array)
    ALLOCATE(array(new_size(1),new_size(2)))

    DO j = 1, old_size(2)
    DO i = 1, old_size(1)
      array(i,j) = tmp_array(i,j)
    END DO
    END DO

    DEALLOCATE(tmp_array)

  END SUBROUTINE grow_real_array2d



  SUBROUTINE grow_integer_array2d(array, idx, idy)

    INTEGER, DIMENSION(:,:), POINTER :: array
    INTEGER, INTENT(IN) :: idx, idy
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: tmp_array
    INTEGER :: old_size(2), new_size(2), i, j, idxy, idir

    old_size = SHAPE(array)
    IF (idx /= old_size(1)) THEN
      IF (idy /= old_size(2)) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'grow_integer_array2d can only grow an array in one dimension'
        RETURN
      END IF
      idir = 1
      idxy = idx
    ELSE
      idir = 2
      idxy = idy
    END IF

    IF (idxy <= old_size(idir)) RETURN

    new_size = old_size
    ALLOCATE(tmp_array(old_size(1),old_size(2)))
    DO j = 1, old_size(2)
    DO i = 1, old_size(1)
      tmp_array(i,j) = array(i,j)
    END DO
    END DO

    new_size(idir) = 2 * old_size(idir)
    IF (new_size(idir) < idxy) THEN
      new_size(idir) = idxy + 1
    END IF
    DEALLOCATE(array)
    ALLOCATE(array(new_size(1),new_size(2)))

    DO j = 1, old_size(2)
    DO i = 1, old_size(1)
      array(i,j) = tmp_array(i,j)
    END DO
    END DO

    DEALLOCATE(tmp_array)

  END SUBROUTINE grow_integer_array2d



  SUBROUTINE abort_code(errcode)

    USE mpi

    INTEGER, INTENT(IN) :: errcode
    INTEGER :: i, newcode, ierr

    ! Translate error code bitmask into 0-255 range.
    ! Only the lowest bit is kept
    newcode = errcode
    DO i = 0, 255
      IF (newcode == 0) THEN
        newcode = i
        EXIT
      END IF
      newcode = newcode / 2
    END DO

    CALL MPI_ABORT(MPI_COMM_WORLD, newcode, ierr)

  END SUBROUTINE abort_code

END MODULE utilities
