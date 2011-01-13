MODULE shape_functions

  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE particle_to_grid(cell_frac, output)

    REAL(num), INTENT(IN) :: cell_frac
    REAL(num), DIMENSION(-2:2), INTENT(INOUT) :: output
    REAL(num), PARAMETER :: fac1 = 1.0_num / 24.0_num
    REAL(num), PARAMETER :: fac2 = 1.0_num / 6.0_num
    REAL(num) :: cf2

#ifdef PARTICLE_SHAPE_BSPLINE3
    output(-2) = fac1 * (0.5_num + cell_frac)**4
    output(-1) = fac2 * (1.1875_num + cell_frac * (2.75_num &
        + cell_frac * (1.5_num - cell_frac - cell_frac**2)))
    output( 0) = 0.25_num * (115.0_num / 48.0_num + cell_frac**2 &
        * (cell_frac**2 - 2.5_num))
    output( 1) = fac2 * (1.1875_num - cell_frac * (2.75_num &
        - cell_frac * (1.5_num + cell_frac - cell_frac**2)))
    output( 2) = fac1 * (0.5_num - cell_frac)**4
#elif  PARTICLE_SHAPE_TOPHAT
    output( 0) = 0.5_num + cell_frac
    output( 1) = 0.5_num - cell_frac
#else
    cf2 = cell_frac**2
    output(-1) = 0.5_num * (0.25_num + cf2 + cell_frac)
    output( 0) = 0.75_num - cf2
    output( 1) = 0.5_num * (0.25_num + cf2 - cell_frac)
#endif

  END SUBROUTINE particle_to_grid



  SUBROUTINE grid_to_particle(cell_frac, output)

    REAL(num), INTENT(IN) :: cell_frac
    REAL(num), DIMENSION(-2:2), INTENT(INOUT) :: output
    REAL(num), PARAMETER :: fac1 = 1.0_num / 24.0_num
    REAL(num), PARAMETER :: fac2 = 1.0_num / 6.0_num
    REAL(num) :: cf2

#ifdef PARTICLE_SHAPE_BSPLINE3
    output(-2) = fac1 * (0.5_num + cell_frac)**4
    output(-1) = fac2 * (1.1875_num + cell_frac * (2.75_num &
        + cell_frac * (1.5_num - cell_frac - cell_frac**2)))
    output( 0) = 0.25_num * (115.0_num / 48.0_num + cell_frac**2 &
        * (cell_frac**2 - 2.5_num))
    output( 1) = fac2 * (1.1875_num - cell_frac * (2.75_num &
        - cell_frac * (1.5_num + cell_frac - cell_frac**2)))
    output( 2) = fac1 * (0.5_num - cell_frac)**4
#elif  PARTICLE_SHAPE_TOPHAT
    output( 0) = 0.5_num + cell_frac
    output( 1) = 0.5_num - cell_frac
#else
    cf2 = cell_frac**2
    output(-1) = 0.5_num * (0.25_num + cf2 + cell_frac)
    output( 0) = 0.75_num - cf2
    output( 1) = 0.5_num * (0.25_num + cf2 - cell_frac)
#endif

  END SUBROUTINE grid_to_particle

END MODULE shape_functions
