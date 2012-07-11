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
      wk_array(ix) = (2.0_num * array(ix) &
          + array(ix-1) + array(ix+1)) * 0.25_num
    ENDDO
#endif
    array(1:nx) = wk_array

    DEALLOCATE(wk_array)

  END SUBROUTINE smooth_array

END MODULE current_smooth
