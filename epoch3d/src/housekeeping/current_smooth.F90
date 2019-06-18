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

MODULE current_smooth

#ifdef HIGH_ORDER_SMOOTHING
  USE shape_functions
#else
  USE constants
#endif
  USE boundary

  IMPLICIT NONE

CONTAINS

  SUBROUTINE current_finish

    CALL current_bcs

    CALL field_bc(jx, jng)
    CALL field_bc(jy, jng)
    CALL field_bc(jz, jng)

    IF (smooth_currents) CALL smooth_current

    IF (use_current_correction) THEN
      jx = jx - initial_jx
      jy = jy - initial_jy
      jz = jz - initial_jz
    END IF

  END SUBROUTINE current_finish



  SUBROUTINE smooth_current

    ! Implements strided compensated binomial filtering

    CALL smooth_array(jx, smooth_its, smooth_comp_its, smooth_strides)
    CALL smooth_array(jy, smooth_its, smooth_comp_its, smooth_strides)
    CALL smooth_array(jz, smooth_its, smooth_comp_its, smooth_strides)

  END SUBROUTINE smooth_current



  SUBROUTINE smooth_array(array, its, comp_its, stride)

    REAL(num), DIMENSION(1-jng:,1-jng:,1-jng:), INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: its
    INTEGER, INTENT(IN) :: comp_its
    INTEGER, INTENT(IN), DIMENSION(:), ALLOCATABLE :: stride
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: wk_array
    INTEGER :: ix, iy, iz
#ifdef HIGH_ORDER_SMOOTHING
    INTEGER :: isubx, isuby, isubz
    REAL(num), DIMENSION(sf_min:sf_max) :: weight_fn
    REAL(num) :: val, w1, w2, w3
#else
    INTEGER, DIMENSION(:), ALLOCATABLE :: stride_inner
    INTEGER :: ng_l, it, istride, cstride
    REAL(num) :: alpha, beta
#endif

#ifdef HIGH_ORDER_SMOOTHING
    CALL particle_to_grid(0.0_num, weight_fn)

    ALLOCATE(wk_array(nx,ny,nz))

    DO iz = 1, nz
    DO iy = 1, ny
    DO ix = 1, nx
      val = 0.0_num
      DO isubz = sf_min, sf_max
        w3 = weight_fn(isubz)
        DO isuby = sf_min, sf_max
          w2 = w3 * weight_fn(isuby)
          DO isubx = sf_min, sf_max
            w1 = w2 * weight_fn(isubx)
            val = val + array(ix+isubx,iy+isuby,iz+isubz) * w1
          END DO
        END DO
      END DO
      wk_array(ix,iy,iz) = val
    END DO
    END DO
    END DO

    array(1:nx,1:ny,1:nz) = wk_array(:,:,:)

    DEALLOCATE(wk_array)
#else
    IF (ALLOCATED(stride)) THEN
      ALLOCATE(stride_inner(SIZE(stride)), SOURCE=stride)
    ELSE
      ALLOCATE(stride_inner(1), SOURCE=[1])
    END IF

    ng_l = MAX(sng, jng)
    alpha = 0.5_num
    beta = (1.0_num - alpha) / 6.0_num

    ALLOCATE(wk_array(1-ng_l:nx+ng_l,1-ng_l:ny+ng_l,1-ng_l:nz+ng_l))

    wk_array = 0.0_num
    wk_array(1-jng:nx+jng,1-jng:ny+jng,1-jng:nz+jng) = &
        array(1-jng:nx+jng,1-jng:ny+jng,1-jng:nz+jng)

    DO it = 1, its + comp_its
      DO istride = 1, SIZE(stride_inner)
        CALL field_bc(wk_array, ng_l)
        cstride = stride_inner(istride)
        DO iz = 1, nz
        DO iy = 1, ny
        DO ix = 1, nx
          array(ix,iy,iz) = alpha * wk_array(ix,iy,iz) &
              + (wk_array(ix-cstride,iy,iz) + wk_array(ix+cstride,iy,iz) &
              +  wk_array(ix,iy-cstride,iz) + wk_array(ix,iy+cstride,iz) &
              +  wk_array(ix,iy,iz-cstride) + wk_array(ix,iy,iz+cstride)) &
              * beta
        END DO
        END DO
        END DO
        wk_array(1:nx,1:ny,1:nz) = array(1:nx,1:ny,1:nz)
      END DO
      IF (it > its) THEN
        alpha = REAL(its, num) * 0.5_num + 1.0_num
      END IF
    END DO

    array(1:nx,1:ny,1:nz) = wk_array(1:nx,1:ny,1:nz)

    DEALLOCATE(wk_array)
    DEALLOCATE(stride_inner)
#endif

  END SUBROUTINE smooth_array

END MODULE current_smooth
