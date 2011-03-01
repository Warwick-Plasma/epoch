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

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: wk_array
    INTEGER :: ix, iy, iz
#ifdef HIGH_ORDER_SMOOTHING
    INTEGER :: isubx, isuby, isubz
    REAL(num), DIMENSION(sf_min:sf_max) :: weight_fn
    REAL(num) :: val, wf, wz
#endif

    ALLOCATE(wk_array(1:nx, 1:ny, 1:nz))

#ifdef HIGH_ORDER_SMOOTHING
    CALL particle_to_grid(0.0_num, weight_fn)

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          val = 0.0_num
          DO isubz = sf_min, sf_max
            wz = weight_fn(isubz)
            DO isuby = sf_min, sf_max
              wf = wz * weight_fn(isuby)
              DO isubx = sf_min, sf_max
                val = val + array(ix+isubx, iy+isubx, iz+isubz) &
                    * weight_fn(isubx) * wf
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
