    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy, gz
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    REAL(num) :: cell_z_r, cell_frac_z
    INTEGER :: cell_x, cell_y, cell_z
#ifndef PARTICLE_SHAPE_TOPHAT
    REAL(num) :: cx2, cy2, cz2
#endif
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num), PARAMETER :: third = 1.0_num / 3.0_num
    REAL(num), PARAMETER :: fac1 = 0.125_num * third
    REAL(num), PARAMETER :: fac2 = 0.5_num * third
    REAL(num), PARAMETER :: fac3 = 7.1875_num * third
#endif
