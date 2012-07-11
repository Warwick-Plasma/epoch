MODULE current_smooth

  USE shape_functions

  IMPLICIT NONE

CONTAINS

  SUBROUTINE smooth_current()

    ! A very simple current smoothing routine

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
      wk_array(ix, iy) = (4.0_num * array(ix, iy) &
          + array(ix-1, iy) + array(ix+1, iy) &
          + array(ix, iy-1) + array(ix, iy+1)) * 0.125_num
    ENDDO
    ENDDO
#endif
    array(1:nx, 1:ny) = wk_array

    DEALLOCATE(wk_array)

  END SUBROUTINE smooth_array

END MODULE current_smooth
