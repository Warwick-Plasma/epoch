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

    REAL(num), DIMENSION(1-jng:), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:), ALLOCATABLE :: wk_array
    INTEGER :: ix
#ifdef HIGH_ORDER_SMOOTHING
    INTEGER :: isubx
    REAL(num), DIMENSION(sf_min:sf_max) :: weight_fn
    REAL(num) :: val, w1
#endif

    ALLOCATE(wk_array(nx))

#ifdef HIGH_ORDER_SMOOTHING
    CALL particle_to_grid(0.0_num, weight_fn)

    DO ix = 1, nx
      val = 0.0_num
      DO isubx = sf_min, sf_max
        w1 = weight_fn(isubx)
        val = val + array(ix+isubx) * w1
      ENDDO
      wk_array(ix) = val
    ENDDO
#else
    DO ix = 1, nx
      wk_array(ix) = 0.5_num * array(ix) &
          + (array(ix-1) + array(ix+1)) * 0.25_num
    ENDDO
#endif
    array(1:nx) = wk_array

    DEALLOCATE(wk_array)

  END SUBROUTINE smooth_array

END MODULE current_smooth
