! Copyright (C) 2010-2013 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009      Chris Brady <C.S.Brady@warwick.ac.uk>
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

MODULE current_smooth

#ifdef HIGH_ORDER_SMOOTHING
  USE shape_functions
#else
  USE shared_data
#endif
  USE boundary

  IMPLICIT NONE

CONTAINS

  SUBROUTINE smooth_current()

    ! A very simple current smoothing routine

    ! First copy in values to ghost cells
    CALL field_bc(jx, jng)
    CALL field_bc(jy, jng)
    CALL field_bc(jz, jng)

    CALL smooth_array(jx)
    CALL smooth_array(jy)
    CALL smooth_array(jz)

  END SUBROUTINE smooth_current



  SUBROUTINE smooth_array(array)

    REAL(num), DIMENSION(1-jng:,1-jng:), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: wk_array
    INTEGER :: ix, iy
#ifdef HIGH_ORDER_SMOOTHING
    INTEGER :: isubx, isuby
    REAL(num), DIMENSION(sf_min:sf_max) :: weight_fn
    REAL(num) :: val, w1, w2
#endif

    ALLOCATE(wk_array(nx, ny))

#ifdef HIGH_ORDER_SMOOTHING
    CALL particle_to_grid(0.0_num, weight_fn)

    DO iy = 1, ny
    DO ix = 1, nx
      val = 0.0_num
      DO isuby = sf_min, sf_max
        w2 = weight_fn(isuby)
        DO isubx = sf_min, sf_max
          w1 = w2 * weight_fn(isubx)
          val = val + array(ix+isubx, iy+isuby) * w1
        ENDDO
      ENDDO
      wk_array(ix, iy) = val
    ENDDO
    ENDDO
#else
    DO iy = 1, ny
    DO ix = 1, nx
      wk_array(ix, iy) = 0.25_num * array(ix, iy) &
          + (array(ix-1, iy) + array(ix+1, iy) &
          + array(ix, iy-1) + array(ix, iy+1)) * 0.125_num &
          + (array(ix-1, iy-1) + array(ix+1, iy-1) &
          + array(ix-1, iy+1) + array(ix+1, iy+1)) * 0.0625_num
    ENDDO
    ENDDO
#endif
    array(1:nx, 1:ny) = wk_array

    DEALLOCATE(wk_array)

  END SUBROUTINE smooth_array

END MODULE current_smooth
