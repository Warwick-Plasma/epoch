    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: gx, gy

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y

    REAL(num) :: cell_x_r, cell_y_r
    INTEGER :: cell_x, cell_y
#ifndef PARTICLE_SHAPE_TOPHAT
    REAL(num) :: cx2, cy2
#endif
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num), PARAMETER :: third = 1.0_num / 3.0_num
    REAL(num), PARAMETER :: fac1 = 0.125_num * third
    REAL(num), PARAMETER :: fac2 = 0.5_num * third
    REAL(num), PARAMETER :: fac3 = 7.1875_num * third
#endif
