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

    REAL(num), DIMENSION(-2:), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:), ALLOCATABLE :: wk_array
    INTEGER :: ix
#ifdef HIGH_ORDER_SMOOTHING
    INTEGER :: isubx
    REAL(num), DIMENSION(sf_min:sf_max) :: weight_fn
    REAL(num) :: val
#endif

    ALLOCATE(wk_array(1:nx))

#ifdef HIGH_ORDER_SMOOTHING
    CALL particle_to_grid(0.0_num, weight_fn)

    DO ix = 1, nx
      val = 0.0_num
      DO isubx = sf_min, sf_max
        val = val + array(ix+isubx) * weight_fn(isubx)
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
