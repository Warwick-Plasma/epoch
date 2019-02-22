    ! Temporary variables
    REAL(num) :: idx

    ! Contains the integer cell position of the particle in x, y, z
    ! cell_x3 is unused here, but is always used later where this is included
    INTEGER :: cell_x1, cell_x2, cell_x3

! in particle_head.inc
!    ! Contains the floating point version of the cell number (never actually
!    ! used)
!    REAL(num) :: cell_x_r

! in particle_head.inc
!    ! The fraction of a cell between the particle position and cell boundary
!    REAL(num) :: cell_frac_x

! in particle_head.inc
!    ! Weighting factors as Eqn 4.77 page 25 of manual
!    ! Eqn 4.77 would be written as
!    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
!    ! Defined at the particle position
!    REAL(num), DIMENSION(sf_min:sf_max) :: gx

#include "particle_head.inc"

    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(sf_min:sf_max) :: hx
    ! Temporary variables
    INTEGER :: dcellx
    REAL(num) :: ex_part, ey_part, ez_part
    REAL(num) :: bx_part, by_part, bz_part

#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = 1.0_num
#else
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

    ! Temporary variables
    REAL(num) :: idtx

    ! Unvarying multiplication factors
    idx = 1.0_num / dx
