MODULE shape_functions

  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE particle_to_grid(cell_frac, gf)

    REAL(num), INTENT(IN) :: cell_frac
    REAL(num), DIMENSION(sf_min:sf_max), INTENT(OUT) :: gf
#ifndef PARTICLE_SHAPE_TOPHAT
    REAL(num) :: cf2
#endif
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num), PARAMETER :: third = 1.0_num / 3.0_num
    REAL(num), PARAMETER :: fac1 = 0.125_num * third
    REAL(num), PARAMETER :: fac2 = 0.5_num * third
    REAL(num), PARAMETER :: fac3 = 7.1875_num * third

    cf2 = cell_frac**2
    gf(-2) = fac1 * (0.5_num + cell_frac)**4
    gf(-1) = fac2 * (1.1875_num + 2.75_num * cell_frac &
        + cf2 * (1.5_num - cell_frac - cf2))
    gf( 0) = 0.25_num * (fac3 + cf2 * (cf2 - 2.5_num))
    gf( 1) = fac2 * (1.1875_num - 2.75_num * cell_frac &
        + cf2 * (1.5_num + cell_frac - cf2))
    gf( 2) = fac1 * (0.5_num - cell_frac)**4
#elif  PARTICLE_SHAPE_TOPHAT
    gf( 0) = 0.5_num + cell_frac
    gf( 1) = 0.5_num - cell_frac
#else
    cf2 = cell_frac**2
    gf(-1) = 0.5_num * (0.25_num + cf2 + cell_frac)
    gf( 0) = 0.75_num - cf2
    gf( 1) = 0.5_num * (0.25_num + cf2 - cell_frac)
#endif

  END SUBROUTINE particle_to_grid



  SUBROUTINE grid_to_particle(cell_frac, gf)

    REAL(num), INTENT(IN) :: cell_frac
    REAL(num), DIMENSION(sf_min:sf_max), INTENT(OUT) :: gf

    CALL particle_to_grid(cell_frac, gf)

  END SUBROUTINE grid_to_particle

END MODULE shape_functions
