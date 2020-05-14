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

MODULE utilities

  USE constants
  USE shared_data

  IMPLICIT NONE

  INTERFACE grow_array
    MODULE PROCEDURE grow_real_array, grow_integer_array, grow_logical_array, &
                     grow_string_array, grow_real_array2d, grow_integer_array2d
  END INTERFACE grow_array

  PRIVATE

  PUBLIC :: erf_func
  PUBLIC :: abort_code, grow_array, get_free_lun
  PUBLIC :: setup_grid_x, setup_grid_y

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
    END DO

    new_size = 2 * old_size
    DEALLOCATE(array)
    ALLOCATE(array(new_size))

    DO i = 1, old_size
      array(i) = tmp_array(i)
    END DO

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



  !----------------------------------------------------------------------------
  ! Approximation to ERF from Abramowitz and Stegun
  ! "Handbook Of Mathematical Functions" - 7.1.28
  ! Should be accurate to a few parts in 10^7
  !----------------------------------------------------------------------------

  FUNCTION erf_func(val)

   REAL(num), INTENT(IN) :: val
   REAL(num), PARAMETER :: a1 = 0.0705230784_num
   REAL(num), PARAMETER :: a2 = 0.0422820123_num
   REAL(num), PARAMETER :: a3 = 0.0092705272_num
   REAL(num), PARAMETER :: a4 = 0.0001520143_num
   REAL(num), PARAMETER :: a5 = 0.0002765672_num
   REAL(num), PARAMETER :: a6 = 0.0000430638_num
   REAL(num), PARAMETER :: unity = 1.0_num
   REAL(num) :: y, denom
   REAL(num) :: erf_func

   y = ABS(val)
   denom = unity + y * (a1 + y * (a2 + y * (a3 + y * (a4 + y * (a5 + y * a6)))))

   ! Use ERF(-ABS(X)) = -ERF(ABS(X))
   erf_func = SIGN(unity - unity / denom**16, val)

  END FUNCTION erf_func



  FUNCTION get_free_lun()

    ! This subroutine simply cycles round until it finds a free lun between
    ! min_lun and max_lun
    INTEGER :: get_free_lun
    INTEGER :: lun
    INTEGER, PARAMETER :: min_lun = 10, max_lun = 20
    LOGICAL :: is_open

    is_open = .TRUE.

    lun = min_lun
    DO
      INQUIRE(unit=lun, opened=is_open)
      IF (.NOT. is_open) EXIT
      lun = lun+1
      IF (lun > max_lun) THEN
        WRITE(*,*) '*** ERROR ***'
        WRITE(*,*) 'Unable to open lun for input deck read'
        CALL abort_code(c_err_io_error)
      END IF
    END DO

    get_free_lun = lun

  END FUNCTION get_free_lun



  SUBROUTINE setup_grid_x

    INTEGER :: iproc
    REAL(num) :: boundary_shift

    ! Setup local grid
    x(1-ng:nx+ng) = x_global(nx_global_min-ng:nx_global_max+ng)
    xb(1-ng:nx+ng) = xb_global(nx_global_min-ng:nx_global_max+ng)

    DO iproc = 0, nprocx-1
      x_grid_mins(iproc) = x_global(cell_x_min(iproc+1))
      x_grid_maxs(iproc) = x_global(cell_x_max(iproc+1))
    END DO

    x_grid_min_local = x_grid_mins(x_coords)
    x_grid_max_local = x_grid_maxs(x_coords)

    x_min_local = x_grid_min_local + (cpml_x_min_offset - 0.5_num) * dx
    x_max_local = x_grid_max_local - (cpml_x_max_offset - 0.5_num) * dx

    boundary_shift = REAL((1 + png + cpml_thickness) / 2, num)

    dir_d(1) = dx
    dir_min(1) = x_min
    dir_max(1) = x_max
    dir_grid_min(1) = x_grid_min
    dir_grid_max(1) = x_grid_max
    dir_min_local(1) = x_min_local
    dir_max_local(1) = x_max_local

  END SUBROUTINE setup_grid_x



  SUBROUTINE setup_grid_y

    INTEGER :: iproc
    REAL(num) :: boundary_shift

    ! Setup local grid
    y(1-ng:ny+ng) = y_global(ny_global_min-ng:ny_global_max+ng)
    yb(1-ng:ny+ng) = yb_global(ny_global_min-ng:ny_global_max+ng)

    DO iproc = 0, nprocy-1
      y_grid_mins(iproc) = y_global(cell_y_min(iproc+1))
      y_grid_maxs(iproc) = y_global(cell_y_max(iproc+1))
    END DO

    y_grid_min_local = y_grid_mins(y_coords)
    y_grid_max_local = y_grid_maxs(y_coords)

    y_min_local = y_grid_min_local + (cpml_y_min_offset - 0.5_num) * dy
    y_max_local = y_grid_max_local - (cpml_y_max_offset - 0.5_num) * dy

    boundary_shift = REAL((1 + png + cpml_thickness) / 2, num)

    dir_d(2) = dy
    dir_min(2) = y_min
    dir_max(2) = y_max
    dir_grid_min(2) = y_grid_min
    dir_grid_max(2) = y_grid_max
    dir_min_local(2) = y_min_local
    dir_max_local(2) = y_max_local

  END SUBROUTINE setup_grid_y

END MODULE utilities
