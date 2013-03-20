MODULE current_smooth

#ifdef HIGH_ORDER_SMOOTHING
  USE shape_functions
#else
  USE shared_data
#endif

  IMPLICIT NONE

CONTAINS

  SUBROUTINE smooth_current()

    ! A very simple current smoothing routine

    CALL smooth_array(jx)
    CALL smooth_array(jy)
    CALL smooth_array(jz)

  END SUBROUTINE smooth_current



  SUBROUTINE smooth_array(array)

    REAL(num), DIMENSION(1-jng:,1-jng:,1-jng:), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: wk_array
    INTEGER :: ix, iy, iz
#ifdef HIGH_ORDER_SMOOTHING
    INTEGER :: isubx, isuby, isubz
    REAL(num), DIMENSION(sf_min:sf_max) :: weight_fn
    REAL(num) :: val, w1, w2, w3
#endif

    ALLOCATE(wk_array(nx, ny, nz))

#ifdef HIGH_ORDER_SMOOTHING
    CALL particle_to_grid(0.0_num, weight_fn)

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
            val = val + array(ix+isubx, iy+isuby, iz+isubz) * w1
          ENDDO
        ENDDO
      ENDDO
      wk_array(ix, iy, iz) = val
    ENDDO
    ENDDO
    ENDDO
#else
    DO iz = 1, nz
    DO iy = 1, ny
    DO ix = 1, nx
      wk_array(ix, iy, iz) = (6.0_num * array(ix, iy, iz) &
          + array(ix-1, iy, iz) + array(ix+1, iy, iz) &
          + array(ix, iy-1, iz) + array(ix, iy+1, iz) &
          + array(ix, iy, iz-1) + array(ix, iy, iz+1)) / 12.0_num
    ENDDO
    ENDDO
    ENDDO
#endif
    array(1:nx, 1:ny, 1:nz) = wk_array

    DEALLOCATE(wk_array)

  END SUBROUTINE smooth_array

END MODULE current_smooth
