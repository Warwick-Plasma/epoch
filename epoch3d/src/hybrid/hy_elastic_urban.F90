! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------
!
! hy_elastic_urban.F90
!
! This module contains a PIC implementation of the elastic scatter functions in
! Geant4. Here we will be using the Urban MSC method, as it can be used to
! calculate the elastic scatter angles for arbitrary charged particles below
! 100 MeV.

MODULE hy_elastic_urban
#ifdef HYBRID

  USE particles
  USE hy_shared

  IMPLICIT NONE

  ! Variables which can only be seen by functions and subroutines in this module

  ! Pre-calculated particle variables to prevent recalculation in functions
  REAL(num), PRIVATE :: part_x, part_y, part_z
  REAL(num), PRIVATE :: part_p, part_m, part_q, part_beta2, part_rad, part_v
  REAL(num), PRIVATE :: part_mc2
  INTEGER, PRIVATE :: ispecies

  ! Table declarations for Urban MSC (only useful in this file)
  REAL(num), PRIVATE :: Z_tab(15), sig_tab(15), high_tab(15)
  REAL(num), PRIVATE :: KE_tab(22), beta2_tab(22)
  REAL(num), PRIVATE :: el_tab(22,15), pos_tab(22,15)

  ! Switch variables
  REAL(num), PARAMETER, PRIVATE :: KE_10MeV = 10.0e6*q0
  REAL(num), PARAMETER, PRIVATE :: beta2_10MeV = KE_10MeV &
      * (KE_10MeV + 2.0_num*mc2) / (KE_10MeV + mc2)**2
  REAL(num), PARAMETER, PRIVATE :: beta_gam2_10MeV = KE_10MeV &
      * (KE_10MeV + 2.0_num*mc2) /  mc2**2

  ! Constants in Geant4 units
  REAL(num), PARAMETER, PRIVATE :: mec2_MeV = mc2 / (1.0e6_num * q0) ! [MeV]
  REAL(num), PARAMETER, PRIVATE :: one_MeV = 1.0e6_num * q0 ! [J]
  REAL(num), PARAMETER, PRIVATE :: r_Bohr_mm = h_bar/(mc0 * alpha) &
      * 1.0e3_num ! [mm]
  REAL(num), PARAMETER, PRIVATE :: hbarc_MeVmm = h_bar*c/(1.0e6_num*q0) &
      * 1.0e3_num ! [MeV.mm]
  REAL(num), PARAMETER, PRIVATE :: epsfac = 2.0_num * mec2_MeV**2 &
      * r_Bohr_mm**2 / hbarc_MeVmm**2 ! [Dimensionless]
  REAL(num), PARAMETER, PRIVATE :: sig_const = 2 * pi * r_Bohr_mm**2 * alpha**4
  REAL(num), PARAMETER, PRIVATE :: c_highland_J = 13.6_num * 1.0e6 * q0 ! [J]

CONTAINS

  SUBROUTINE setup_urban_scatter

    ! Runs the initialisation routines for the elastic scatter process. These
    ! include reading the material tables needed for the Urban algorithm, and
    ! caching the interpolations for the solids present in the simulation

    INTEGER :: i_sol

    CALL setup_urban_tables

    DO i_sol = 1, solid_count
      CALL setup_urban_cache(i_sol)
    END DO

    CALL set_Z_cache

  END SUBROUTINE setup_urban_scatter



  SUBROUTINE setup_urban_tables

    ! Initialises the tables for the Urban MSC algorithm

    ! Table of target atomic number values, Z
    Z_tab = (/ 4.0_num, 6.0_num, 13.0_num, 20.0_num, 26.0_num, 29.0_num, &
        32.0_num, 38.0_num, 47.0_num, 50.0_num, 56.0_num, 64.0_num, 74.0_num, &
        79.0_num, 82.0_num /)

    ! Table of incident particle kinetic energies [MeV/c]
    KE_tab = (/ 10.0e-6_num, 200.0e-6_num, 400.0e-6_num, 700.0e-6_num, &
        1.0e-3_num, 2.0e-3_num, 4.0e-3_num, 7.0e-3_num, 10.0e-3_num, &
        20.0e-3_num, 40.0e-3_num, 70.0e-3_num, 100.0e-3_num, 200.0e-3_num, &
        400.0e-3_num, 700.0e-3_num, 1.0_num, 2.0_num, 4.0_num, 7.0_num, &
        10.0_num, 20.0_num /)

    ! Convert the KE table into a beta² table
    beta2_tab = KE_tab * (KE_tab + 2.0_num*mec2_MeV)/(KE_tab + mec2_MeV)**2

    ! Cross sections at each Z [mm²]
    sig_tab = (/0.2672_num, 0.5922_num, 2.653_num, 6.235_num, 11.69_num, &
        13.24_num, 16.12_num, 23.00_num, 35.13_num, 39.95_num, 50.85_num, &
        67.19_num, 91.15_num, 104.4_num, 113.1_num/)
    sig_tab = sig_tab * 1.0e-22_num

    ! High energy correction factors for each Z [mm²]
    high_tab = (/120.70_num, 117.50_num, 105.00_num, 92.92_num, 79.23_num, &
        74.510_num, 68.29_num, 57.39_num, 41.97_num, 36.14_num, 24.53_num, &
        10.21_num, -7.855_num, -16.84_num, -22.30_num/)

    ! Electron factors for these KE and Z vals (KE vals, Z vals) (dimensionless)
    el_tab = RESHAPE((/ 1.125_num,1.072_num,1.051_num,1.047_num,1.047_num, &
        1.050_num,1.052_num,1.054_num,1.054_num,1.057_num,1.062_num,1.069_num, &
        1.075_num,1.090_num,1.105_num,1.111_num,1.112_num,1.108_num,1.100_num, &
        1.093_num,1.089_num,1.087_num, &
        1.408_num,1.246_num,1.143_num,1.096_num,1.077_num,1.059_num,1.053_num, &
        1.051_num,1.052_num,1.053_num,1.058_num,1.065_num,1.072_num,1.087_num, &
        1.101_num,1.108_num,1.109_num,1.105_num,1.097_num,1.090_num,1.086_num, &
        1.082_num, &
        2.833_num,2.268_num,1.861_num,1.612_num,1.486_num,1.309_num,1.204_num, &
        1.156_num,1.136_num,1.114_num,1.106_num,1.106_num,1.109_num,1.119_num, &
        1.129_num,1.132_num,1.131_num,1.124_num,1.113_num,1.104_num,1.099_num, &
        1.098_num, &
        3.879_num,3.016_num,2.380_num,2.007_num,1.818_num,1.535_num,1.340_num, &
        1.236_num,1.190_num,1.133_num,1.107_num,1.099_num,1.098_num,1.103_num, &
        1.110_num,1.113_num,1.112_num,1.105_num,1.096_num,1.089_num,1.085_num, &
        1.098_num, &
        6.937_num,4.330_num,2.886_num,2.256_num,1.987_num,1.628_num,1.395_num, &
        1.265_num,1.203_num,1.122_num,1.080_num,1.065_num,1.061_num,1.063_num, &
        1.070_num,1.073_num,1.073_num,1.070_num,1.064_num,1.059_num,1.056_num, &
        1.056_num, &
        9.616_num,5.708_num,3.424_num,2.551_num,2.204_num,1.762_num,1.485_num, &
        1.330_num,1.256_num,1.155_num,1.099_num,1.077_num,1.070_num,1.068_num, &
        1.072_num,1.074_num,1.074_num,1.070_num,1.063_num,1.059_num,1.056_num, &
        1.052_num, &
        11.72_num,6.364_num,3.811_num,2.806_num,2.401_num,1.884_num,1.564_num, &
        1.386_num,1.300_num,1.180_num,1.112_num,1.082_num,1.073_num,1.066_num, &
        1.068_num,1.069_num,1.068_num,1.064_num,1.059_num,1.054_num,1.051_num, &
        1.050_num, &
        18.08_num,8.601_num,4.569_num,3.183_num,2.662_num,2.025_num,1.646_num, &
        1.439_num,1.339_num,1.195_num,1.108_num,1.068_num,1.053_num,1.040_num, &
        1.039_num,1.039_num,1.039_num,1.037_num,1.034_num,1.031_num,1.030_num, &
        1.036_num, &
        18.22_num,10.48_num,5.333_num,3.713_num,3.115_num,2.367_num,1.898_num, &
        1.631_num,1.498_num,1.301_num,1.171_num,1.105_num,1.077_num,1.048_num, &
        1.036_num,1.033_num,1.031_num,1.028_num,1.024_num,1.022_num,1.021_num, &
        1.024_num, &
        14.14_num,10.65_num,5.710_num,3.929_num,3.266_num,2.453_num,1.951_num, &
        1.669_num,1.528_num,1.319_num,1.178_num,1.106_num,1.075_num,1.040_num, &
        1.027_num,1.022_num,1.020_num,1.017_num,1.015_num,1.013_num,1.013_num, &
        1.020_num, &
        14.11_num,11.73_num,6.312_num,4.240_num,3.478_num,2.566_num,2.022_num, &
        1.720_num,1.569_num,1.342_num,1.186_num,1.102_num,1.065_num,1.022_num, &
        1.003_num,0.997_num,0.995_num,0.993_num,0.993_num,0.993_num,0.993_num, &
        1.011_num, &
        22.76_num,20.01_num,8.835_num,5.287_num,4.144_num,2.901_num,2.219_num, &
        1.855_num,1.677_num,1.410_num,1.224_num,1.121_num,1.073_num,1.014_num, &
        0.986_num,0.976_num,0.974_num,0.972_num,0.973_num,0.974_num,0.975_num, &
        0.987_num, &
        50.77_num,40.85_num,14.13_num,7.184_num,5.284_num,3.435_num,2.520_num, &
        2.059_num,1.837_num,1.512_num,1.283_num,1.153_num,1.091_num,1.010_num, &
        0.969_num,0.954_num,0.950_num,0.947_num,0.949_num,0.952_num,0.954_num, &
        0.963_num, &
        65.87_num,59.06_num,15.87_num,7.570_num,5.567_num,3.650_num,2.682_num, &
        2.182_num,1.939_num,1.579_num,1.325_num,1.178_num,1.108_num,1.014_num, &
        0.965_num,0.947_num,0.941_num,0.938_num,0.940_num,0.944_num,0.946_num, &
        0.954_num, &
        55.60_num,47.34_num,15.92_num,7.810_num,5.755_num,3.767_num,2.760_num, &
        2.239_num,1.985_num,1.609_num,1.343_num,1.188_num,1.113_num,1.013_num, &
        0.960_num,0.939_num,0.933_num,0.930_num,0.933_num,0.936_num,0.939_num, &
        0.949_num /), SHAPE(el_tab))

    ! Positron factors for these KE and Z values (KE_vals, Z_vals)
    ! (dimensionless)
    pos_tab = RESHAPE((/ 2.589_num,2.044_num,1.658_num,1.446_num,1.347_num, &
        1.217_num,1.144_num,1.110_num,1.097_num,1.083_num,1.080_num,1.086_num, &
        1.092_num,1.108_num,1.123_num,1.131_num,1.131_num,1.126_num,1.117_num, &
        1.108_num,1.103_num,1.100_num, &
        3.904_num,2.794_num,2.079_num,1.710_num,1.543_num,1.325_num,1.202_num, &
        1.145_num,1.122_num,1.096_num,1.089_num,1.092_num,1.098_num,1.114_num, &
        1.130_num,1.137_num,1.138_num,1.132_num,1.122_num,1.113_num,1.108_num, &
        1.102_num, &
        7.970_num,6.080_num,4.442_num,3.398_num,2.872_num,2.127_num,1.672_num, &
        1.451_num,1.357_num,1.246_num,1.194_num,1.179_num,1.178_num,1.188_num, &
        1.201_num,1.205_num,1.203_num,1.190_num,1.173_num,1.159_num,1.151_num, &
        1.145_num, &
        9.714_num,7.607_num,5.747_num,4.493_num,3.815_num,2.777_num,2.079_num, &
        1.715_num,1.553_num,1.353_num,1.253_num,1.219_num,1.211_num,1.214_num, &
        1.225_num,1.228_num,1.225_num,1.210_num,1.191_num,1.175_num,1.166_num, &
        1.174_num, &
        17.97_num,12.95_num,8.628_num,6.065_num,4.849_num,3.222_num,2.275_num, &
        1.820_num,1.624_num,1.382_num,1.259_num,1.214_num,1.202_num,1.202_num, &
        1.214_num,1.219_num,1.217_num,1.203_num,1.184_num,1.169_num,1.160_num, &
        1.151_num, &
        24.83_num,17.06_num,10.84_num,7.355_num,5.767_num,3.707_num,2.546_num, &
        1.996_num,1.759_num,1.465_num,1.311_num,1.252_num,1.234_num,1.228_num, &
        1.238_num,1.241_num,1.237_num,1.222_num,1.201_num,1.184_num,1.174_num, &
        1.159_num, &
        23.26_num,17.15_num,11.52_num,8.049_num,6.375_num,4.114_num,2.792_num, &
        2.155_num,1.880_num,1.535_num,1.353_num,1.281_num,1.258_num,1.247_num, &
        1.254_num,1.256_num,1.252_num,1.234_num,1.212_num,1.194_num,1.183_num, &
        1.170_num, &
        22.33_num,18.01_num,12.86_num,9.212_num,7.336_num,4.702_num,3.117_num, &
        2.348_num,2.015_num,1.602_num,1.385_num,1.297_num,1.268_num,1.251_num, &
        1.256_num,1.258_num,1.254_num,1.237_num,1.214_num,1.195_num,1.185_num, &
        1.179_num, &
        33.91_num,24.13_num,15.71_num,10.80_num,8.507_num,5.467_num,3.692_num, &
        2.808_num,2.407_num,1.873_num,1.564_num,1.425_num,1.374_num,1.330_num, &
        1.324_num,1.320_num,1.312_num,1.288_num,1.258_num,1.235_num,1.221_num, &
        1.205_num, &
        32.14_num,24.11_num,16.30_num,11.40_num,9.015_num,5.782_num,3.868_num, &
        2.917_num,2.490_num,1.925_num,1.596_num,1.447_num,1.391_num,1.342_num, &
        1.332_num,1.327_num,1.320_num,1.294_num,1.264_num,1.240_num,1.226_num, &
        1.214_num, &
        29.51_num,24.07_num,17.19_num,12.28_num,9.766_num,6.238_num,4.112_num, &
        3.066_num,2.602_num,1.995_num,1.641_num,1.477_num,1.414_num,1.356_num, &
        1.342_num,1.336_num,1.328_num,1.302_num,1.270_num,1.245_num,1.231_num, &
        1.233_num, &
        38.19_num,30.85_num,21.76_num,15.35_num,12.07_num,7.521_num,4.812_num, &
        3.498_num,2.926_num,2.188_num,1.763_num,1.563_num,1.484_num,1.405_num, &
        1.382_num,1.371_num,1.361_num,1.330_num,1.294_num,1.267_num,1.251_num, &
        1.239_num, &
        49.71_num,39.80_num,27.96_num,19.63_num,15.36_num,9.407_num,5.863_num, &
        4.155_num,3.417_num,2.478_num,1.944_num,1.692_num,1.589_num,1.480_num, &
        1.441_num,1.423_num,1.409_num,1.372_num,1.330_num,1.298_num,1.280_num, &
        1.258_num, &
        59.25_num,45.08_num,30.36_num,20.83_num,16.15_num,9.834_num,6.166_num, &
        4.407_num,3.641_num,2.648_num,2.064_num,1.779_num,1.661_num,1.531_num, &
        1.482_num,1.459_num,1.442_num,1.400_num,1.354_num,1.319_num,1.299_num, &
        1.272_num, &
        56.38_num,44.29_num,30.50_num,21.18_num,16.51_num,10.11_num,6.354_num, &
        4.542_num,3.752_num,2.724_num,2.116_num,1.817_num,1.692_num,1.554_num, &
        1.499_num,1.474_num,1.456_num,1.412_num,1.364_num,1.328_num,1.307_num, &
        1.282_num /), SHAPE(pos_tab))

  END SUBROUTINE setup_urban_tables



  SUBROUTINE setup_urban_cache(i_sol)

    ! The tables in setup_urban_tables give values of cross section sig_tab,
    ! high energy correction high_tab and charge-specific variables el_tab and
    ! pos_tab. We can save time by interpolating these tables to the atomic
    ! number of the present solids. This subroutine interpolates in Z² and gets
    ! these tables and values for the current solid, i_sol

    INTEGER, INTENT(IN) :: i_sol
    REAL(num) :: frac
    INTEGER :: iZ1, iZ2

    ! Pull out the Z indices above and below the current element Z
    IF (solid_array(i_sol)%Z < Z_tab(2)) THEN
      iZ2 = 2
    ELSE IF (solid_array(i_sol)%Z >= Z_tab(14)) THEN
      iZ2 = 15
    ELSE
      DO iZ2 = 3,14
        IF (Z_tab(iZ2) > solid_array(i_sol)%Z) THEN
          EXIT
        END IF
      END DO
    END IF
    iZ1 = iZ2 - 1

    ! Interpolate in Z²
    IF (solid_array(i_sol)%Z < Z_tab(1)) THEN
      ! Z < 4
      frac = (solid_array(i_sol)%Z / Z_tab(1))**2
      solid_array(i_sol)%urb_sig = sig_tab(1) * frac
      solid_array(i_sol)%urb_high = high_tab(1) * solid_array(i_sol)%urb_sig
      solid_array(i_sol)%urb_el = el_tab(:,1) *frac
      solid_array(i_sol)%urb_pos = pos_tab(:,1) *frac

    ELSE IF (solid_array(i_sol)%Z > Z_tab(15)) THEN
      ! Z > 82
      frac = (solid_array(i_sol)%Z / Z_tab(15))**2
      solid_array(i_sol)%urb_sig = sig_tab(15) * frac
      solid_array(i_sol)%urb_high = high_tab(15) * solid_array(i_sol)%urb_sig
      solid_array(i_sol)%urb_el = el_tab(:,15) * frac
      solid_array(i_sol)%urb_pos = pos_tab(:,15) * frac

    ELSE
      ! 4 <= Z <= 82
      frac = (solid_array(i_sol)%Z**2 - Z_tab(iZ1)**2) &
          / (Z_tab(iZ2)**2 - Z_tab(iZ1)**2)
      solid_array(i_sol)%urb_sig = sig_tab(iZ1) &
          + frac*(sig_tab(iZ2) - sig_tab(iZ1))
      solid_array(i_sol)%urb_high = high_tab(iZ1)*sig_tab(iZ1) &
          + frac*(high_tab(iZ2)*sig_tab(iZ2) - high_tab(iZ1)*sig_tab(iZ1))
      solid_array(i_sol)%urb_el = el_tab(:,iZ1) &
          + frac*(el_tab(:,iZ2) - el_tab(:,iZ1))
      solid_array(i_sol)%urb_pos = pos_tab(:,iZ1) &
          + frac*(pos_tab(:,iZ2) - pos_tab(:,iZ1))
    END IF

  END SUBROUTINE setup_urban_cache



  SUBROUTINE set_Z_cache

    ! Some variables will depend on the atomic numbers and number densities of
    ! materials, and so can be different in each cell. These can be quite
    ! computationally expensive to calculate, so pre-calculate these and
    ! ensure the load balancer handles them right

    REAL(num), ALLOCATABLE :: Z16(:,:,:)
    REAL(num) :: sum_ni, ni_max
    INTEGER :: ix, iy, iz, isol

    ! Allocate arrays
    ALLOCATE(urb_d1(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(urb_d2(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(urb_d3(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(urb_d4(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(urb_Zeff(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(urb_iZ23(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(urb_theta1(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(urb_theta2(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(urb_rad_len(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))

    ! Allocate temporary arrays
    ALLOCATE(Z16(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))

    ! Effective Z is the average atomic number of the solid, weighted by number
    ! density. Also calculate the radiation length in each cell. Hybrid routines
    ! only expect one radiation length per cell, so assume this belongs to the
    ! solid with the highest number density in this cell
    urb_Zeff = 0.0_num
    DO iz = 1-ng,nz+ng
      DO iy = 1-ng,ny+ng
        DO ix = 1-ng,nx+ng
          sum_ni = 0.0_num
          ni_max = 0.0_num
          DO isol = 1,solid_count
            ! Effective Z
            urb_Zeff(ix,iy,iz) = urb_Zeff(ix,iy,iz) &
            + solid_array(isol)%ion_density(ix,iy,iz) &
            * solid_array(isol)%Z
            sum_ni = sum_ni + solid_array(isol)%ion_density(ix,iy,iz)

            ! Radiation length
            IF (solid_array(isol)%ion_density(ix,iy,iz) > ni_max) THEN
              urb_rad_len(ix,iy,iz) = solid_array(isol)%rad_len
              ni_max = solid_array(isol)%ion_density(ix,iy,iz)
            END IF
          END DO
          urb_Zeff(ix,iy,iz) = urb_Zeff(ix,iy,iz)/sum_ni
        END DO
      END DO
    END DO

    ! Powers of effective Z
    Z16 = urb_Zeff**(1.0_num/6.0_num)
    urb_iZ23 = 1.0_num/Z16**4

    ! Theta constants
    urb_theta1 = (0.990395_num - 0.168386_num*Z16 + 0.093286_num*Z16**2) &
        * (1.0_num - 0.08778_num/urb_Zeff)
    urb_theta2 = (0.990395_num - 0.168386_num*Z16 + 0.093286_num*Z16**2) &
        * (4.078e-2_num + 1.7315e-4_num*urb_Zeff)

    ! Tail constants
    urb_d1 = 2.3785_num - 0.41981_num*Z16**2 + 0.0631_num*Z16**4
    urb_d2 = 0.47526_num + 1.7694_num*Z16**2 - 0.33885_num*Z16**4
    urb_d3 = 0.23683_num - 1.8111_num*Z16**2 + 0.32774_num*Z16**4
    urb_d4 = 0.017888_num + 0.019659_num*Z16**2 - 0.002664_num*Z16**4

    ! Deallocate temporary arrays
    DEALLOCATE(Z16)

  END SUBROUTINE set_Z_cache



  SUBROUTINE Urban_elastic_scatter

    ! Cycles through all charged particles and calculates the elastic scatter
    ! according to the Urban MSC algorithm in Geant4

    REAL(num) :: part_p2, part_KE, e_KE, lambda1, KE_MeV
    REAL(num) :: step_lim, step_size, small_step, theta0, costheta, phi
    TYPE(particle), POINTER :: current

    ! Update the optical depth for each particle species
    DO ispecies = 1, n_species

      ! No elastic scatter for photons
      IF (species_list(ispecies)%species_type == c_species_id_photon) CYCLE

      part_m = species_list(ispecies)%mass
      part_mc2 = part_m * c**2
      part_q = species_list(ispecies)%charge / q0

      current => species_list(ispecies)%attached_list%head
      DO WHILE (ASSOCIATED(current))

        ! Particle position on grid
        part_x = current%part_pos(1) - x_grid_min_local
        part_y = current%part_pos(2) - y_grid_min_local
        part_z = current%part_pos(3) - z_grid_min_local

        ! Calculate effective KE (electron equilvalent), skip immobile particles
        part_p2 = current%part_p(1)**2 + current%part_p(2)**2 &
            + current%part_p(3)**2
        IF (part_p2 < c_tiny) THEN
          current=>current%next
          CYCLE
        END IF
        part_p = SQRT(part_p2)
        part_KE = SQRT(part_p2*c**2 + part_mc2**2) - part_mc2
        part_v = part_p * c**2 / (part_KE + part_mc2)
        IF (species_list(ispecies)%species_type == c_species_id_electron &
            .OR. species_list(ispecies)%species_type == c_species_id_positron) &
            THEN
          e_KE = part_KE
        ELSE
          e_KE = compute_effective_KE(part_KE, part_m)
        END IF

        ! Calculate the first transport mean free path
        part_beta2 = part_KE*(part_KE + 2.0_num*part_mc2)&
            / (part_KE + part_mc2)**2
        lambda1 = compute_lambda1(e_KE)

        ! Calcualte the theta0 parameter to describe angular distribution
        ! We use the higher of current step size (vdt) or small_step
        ! In EPOCH, small_step is measured in [m], (different to Geant4 in [mm])
        KE_MeV = part_KE / one_MeV
        step_lim = MAX(1.0e-2_num*lambda1/(KE_MeV*(10.0_num + KE_MeV)), &
            1.0e-11_num)
        small_step = MIN(step_lim, 1.0e-3_num)
        step_size = SQRT(part_beta2) * c * dt
        IF (step_size > small_step) THEN
          theta0 = compute_theta0(part_KE, step_size)
        ELSE
          theta0 = compute_theta0(part_KE, small_step) &
              * SQRT(step_size/small_step)
        END IF

        ! Sample cos(theta)
        costheta = sample_cos_theta(step_size, small_step, lambda1, theta0)

        ! Apply rotation
        phi = 2.0_num * pi * random()
        CALL rotate_p(current, costheta, phi, part_p)

        current => current%next
      END DO
    END DO

  END SUBROUTINE Urban_elastic_scatter



  FUNCTION compute_effective_KE(KE, mass)

    ! The Urban algorithm depends on the product of momentum and beta only, and
    ! the tables correspond to particles with an electron/positron mass. For a
    ! given particle, the effective KE is the kinetic energy an electron would
    ! need such that p*beta for the electron matches p*beta for the particle
    !
    ! Inputs: KE [J], mass [kg], outputs: compute_effective_KE [J]

    REAL(num), INTENT(IN) :: KE, mass
    REAL(num) :: compute_effective_KE
    REAL(num) :: tau_1, tau_2

    tau_1 = KE/(mass*c**2)
    tau_2 = (mass/m0) * tau_1 * (tau_1 + 2.0_num) / (tau_1 + 1.0_num)
    compute_effective_KE = 0.5_num * mc2 &
        * (tau_2-2.0_num + sqrt((tau_2-2.0_num)**2 + 4.0_num*tau_2))

  END FUNCTION compute_effective_KE



  FUNCTION compute_lambda1(KE)

    ! Calculates the first transport mean free path lambda1

    REAL(num), INTENT(IN) :: KE
    REAL(num) :: compute_lambda1
    REAL(num) :: beta_gam2
    REAL(num) :: eps, sigma_1, sigma_2, sigma, corr, part_ni, part_iZ23, frac
    INTEGER :: iB21, iB22, isol

    ! Get beta²gamma² and beta²
    beta_gam2 = KE*(KE + 2.0_num*part_mc2)/part_mc2**2

    ! Find the energy indices surrounding the effective KE of the particle if
    ! its under 10 MeV
    IF (KE < KE_10MeV) THEN

      ! Find upper index
      IF (part_beta2 <= beta2_tab(2)) THEN
        iB22 = 2
      ELSE
        DO iB22 = 1,22
          IF (beta2_tab(iB22) > part_beta2) THEN
            EXIT
          END IF
        END DO
      END IF

      ! Lower index and interpolation fraction (interpolate in beta²)
      iB21 = iB22 - 1
      frac = (part_beta2 - beta2_tab(iB21))/(beta2_tab(iB22) - beta2_tab(iB21))

      ! Cross section per unit volume prefactor
      CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, part_iZ23, &
          urb_iZ23)
      eps = epsfac * beta_gam2 * part_iZ23
      IF (eps < 1.0e-4_num) THEN
        sigma_1 = 2.0_num*eps**2
      ELSE IF(eps > 1.0e10_num) THEN
        sigma_1 = LOG(1.0_num+2.0_num*eps) - 2.0_num*eps/(1.0_num+2.0_num*eps)
      ELSE
        sigma_1 = LOG(2.0_num*eps) - 1.0_num + 1.0_num/eps
      END IF

      ! Sum cross section contributions from each solid
      sigma = 0.0_num
      sigma_2 = sigma_1 * part_q**2 / (part_beta2*beta_gam2) * sig_const
      DO isol = 1, solid_count
        IF (part_q < 0.0_num) THEN
          ! Treat negative charges like electrons
          corr = solid_array(isol)%urb_el(iB21) &
              + frac*(solid_array(isol)%urb_el(iB22) &
              - solid_array(isol)%urb_el(iB21))
        ELSE
          ! Treat positive charges like positrons
          corr = solid_array(isol)%urb_pos(iB21) &
              + frac*(solid_array(isol)%urb_pos(iB22) &
              - solid_array(isol)%urb_pos(iB21))
        END IF

        CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, part_ni, &
            solid_array(isol)%ion_density)
        part_ni = part_ni * 1.0e-9_num  ! Convert to 1/mm³

        sigma = sigma + part_ni * solid_array(isol)%Z**2 * sigma_2 / corr
      END DO

    ! If particle is over 10 MeV effective KE
    ELSE
        sigma = 0
        DO isol = 1,solid_count
          CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, &
              part_ni, solid_array(isol)%ion_density)
          part_ni = part_ni * 1.0e-9_num  ! Convert to 1/mm³

          sigma = sigma + part_ni * (beta_gam2_10MeV / beta_gam2) &
              * (solid_array(isol)%urb_sig + (part_beta2 - beta2_10MeV) &
              * solid_array(isol)%urb_high)
        END DO
    END IF

    ! Lambda1 = 1/sigma if sigma is positive, 1.0e-3 changes [mm] to [m]
    IF (sigma > 0) THEN
        compute_lambda1 = 1.0e-3_num/sigma
    ELSE
        compute_lambda1 = c_largest_number
    END IF

  END FUNCTION compute_lambda1



  FUNCTION compute_theta0(KE, step_size)

    ! Calculates theta_0, a parameter from the Highland-Lynch-Dahl formula which
    ! describes the width of the angular distribution

    REAL(num), INTENT(IN) :: KE, step_size
    REAL(num) :: compute_theta0
    REAL(num) :: part_Zeff, part_th1, part_th2, inv_vp
    REAL(num) :: y_val, tau, corr, y_low, y_high, y0, y1
    REAL(num) :: x_low, x_high, x_val, x1, x2, x3, x4

    ! Ratio of step size to radiation length
    CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, part_rad, &
        urb_rad_len)
    y_val = step_size / part_rad

    inv_vp = 1.0_num / (part_v * part_p)

    ! Correction for positrons
    IF (species_list(ispecies)%species_type == c_species_id_positron) THEN
      x_low = 0.6_num
      x_high = 0.9_num

      CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, part_Zeff, &
          urb_Zeff)

      tau = KE/(part_m*c**2)
      x_val = SQRT(tau * (tau + 2.0_num)/(tau + 1.0_num)**2)
      x1 = 0.994_num - 4.08e-3_num*part_Zeff
      x2 = 7.16_num + (52.6_num + 365.0_num/part_Zeff)/part_Zeff
      x3 = 1.000_num - 4.47e-3_num*part_Zeff
      x4 = 1.21e-3_num*part_Zeff

      IF (x_val < x_low) THEN
        corr = x1 * (1.0_num - EXP(-x2*x_val))
      ELSE IF (x_val > x_high) THEN
        corr = x3 + x4*EXP(113.0_num*(x_val - 1.0_num))
      ELSE
        y_low = x1*(1.0_num - EXP(-x2*x_low))
        y_high = x3 + x4*EXP(113.0_num*(x_high - 1.0_num))
        y0 = (y_high - y_low)/(x_high - x_low)
        y1 = y_low - y0*x_low
        corr = y0*x_val + y1
      END IF

      y_val = y_val * corr*(1.41125_num &
          + part_Zeff*(1.84035e-4_num*part_Zeff - 1.86427e-2_num))
    END IF

    compute_theta0 = c_highland_J * ABS(part_q) * SQRT(y_val) * inv_vp

    ! Correction factor from e- scattering data
    CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, part_th1, &
        urb_theta1)
    CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, part_th2, &
        urb_theta2)
    compute_theta0 = compute_theta0 * (part_th1 + part_th2*LOG(y_val))

  END FUNCTION compute_theta0



  FUNCTION sample_cos_theta(step_size, small_step, lambda1, theta0)

    ! When sampling the scattering cos(theta), we require two variables to act
    ! as switches between the different forms of cos(theta).
    !
    ! This function is a recreation of Geant4's:
    ! G4UrbanMscModel::SampleCosineTheta()
    ! We duplicate the comments found in the Geant4 source code, but make no
    ! attempt to explain how it works

    REAL(num), INTENT(IN) :: step_size, small_step, lambda1, theta0
    REAL(num) :: sample_cos_theta
    REAL(num) :: tau, u_val, x_val, d1, d2, d3, d4, xsi, c_fac, d_val
    REAL(num) :: iexp_xsi, cmin2, one_min_d
    REAL(num) :: xmean1, xmean2, theta_mean, theta_mean2
    REAL(num) :: p_fac, prob_p, prob_q, xr, xr1, xr2, xr3, xr4, v1, v2

    ! Protection for very small angles
    IF (theta0 < 1e-16_num) THEN
      sample_cos_theta = 1.0_num
      RETURN
    END IF

    tau = step_size / lambda1

    IF (tau >= 8.0_num) THEN
      sample_cos_theta = 2.0_num*random() - 1.0_num
      RETURN
    END IF

    IF (step_size > small_step) THEN
      u_val = tau**(1.0_num/6.0_num)
    ELSE
      u_val = (small_step/lambda1)**(1.0_num/6.0_num)
    END IF

    IF (theta0 < 0.1) THEN
      x_val = theta0**2 * (1.0_num - theta0**2/12.0_num)
    ELSE
      x_val = 4.0_num * SIN(0.5_num*theta0)**2;
    END IF

    ! Parameter for tail, but should not be too big (1.9 minimum)
    CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, d1, urb_d1)
    CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, d2, urb_d2)
    CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, d3, urb_d3)
    CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, d4, urb_d4)
    xsi = MAX(d1 + d2*u_val + d3*u_val**2 + d4*LOG(lambda1/part_rad), 1.9_num)

    IF (ABS(3.0_num - xsi) < 0.001) THEN
        c_fac = 3.001_num
    ELSE IF (ABS(2.0_num - xsi) < 0.001) THEN
        c_fac = 2.001_num
    ELSE
        c_fac = xsi
    END IF

    ! From continuity of derivatives, b = 1 + (c-xsi)*x
    d_val = (c_fac*x_val / ((c_fac - xsi)*x_val + 2.0_num))**(c_fac-1.0_num)
    iexp_xsi = EXP(-xsi)
    cmin2 = c_fac - 2.0_num
    one_min_d = 1.0_num - d_val
    xmean1 = 1.0_num - x_val*(1.0_num - (1.0_num + xsi)*iexp_xsi) &
        / (1.0_num - iexp_xsi)
    xmean2 = ((1.0_num - xsi*x_val + d_val)*cmin2 - c_fac*x_val &
        + (2.0_num + (c_fac - xsi)*x_val)*d_val) / (one_min_d*cmin2)

    IF (tau < 0.01_num) THEN
        theta_mean = 1.0_num - tau*(1.0_num - 0.5_num*tau)
        theta_mean2 = EXP(-tau)
    ELSE
        theta_mean = 1.0_num - 1.0_num/3.0_num * tau * (5.0_num - 6.25_num*tau)
        theta_mean2 = 1.0_num/3.0_num * (1.0_num + 2.0_num*EXP(-2.5_num*tau))
    END IF

    IF (theta0 > pi/6.0_num .OR. xmean1 <= 0.999_num*theta_mean) THEN
      sample_cos_theta = simple_scatter(theta_mean, theta_mean2)
      RETURN
    END IF

    p_fac = (c_fac - 1.0_num)*(1.0_num - iexp_xsi)
    prob_p = p_fac / (p_fac + iexp_xsi*c_fac*one_min_d)
    prob_q = theta_mean/(prob_p*xmean1 + (1.0_num - prob_p)*xmean2)

    ! Sample cos(theta)
    xr = random()
    v1 = one_min_d * xr
    v2 = v1 / (d_val*(c_fac - 1.0_num))

    xr1 = random()
    xr2 = random()
    xr3 = random()
    IF (xr1 < prob_q) THEN
      IF (xr2 < prob_p) THEN
        sample_cos_theta = 1.0_num &
            + LOG(iexp_xsi + xr3*(1.0_num - iexp_xsi))*x_val
      ELSE
        IF (v1 < 0.01_num*d_val) THEN
          sample_cos_theta = -1.0_num + v2 &
              * (1.0_num - 0.5_num*v2*c_fac)*(2.0_num + (c_fac - xsi)*x_val)
        ELSE
          sample_cos_theta = 1.0_num + x_val &
              * (c_fac - xsi - c_fac*(v1 + d_val)**(-1.0_num/(c_fac - 1.0_num)))
        END IF
      END IF
    ELSE
      sample_cos_theta = -1.0_num + 2.0_num*xr2;
    END IF

  END FUNCTION sample_cos_theta



  FUNCTION simple_scatter(theta_mean, theta_mean2)

    ! Simple sampling formula for cos(theta), using <cos(theta)>, <cos²(theta)>

    REAL(num), INTENT(IN) :: theta_mean, theta_mean2
    REAL(num) :: simple_scatter
    REAL(num) :: a_temp, prob_temp, xr1, xr2

    a_temp = (2.0_num*theta_mean + 9.0_num*theta_mean2 - 3.0_num) &
        / (2.0_num*theta_mean - 3.0_num*theta_mean2 + 1.0_num)
    prob_temp = (a_temp - 2.0_num)*theta_mean/a_temp
    xr1 = random()
    xr2 = random()
    IF (xr1 < prob_temp) THEN
        simple_scatter = 2.0_num*xr2**(a_temp+1.0_num) - 1.0_num
    ELSE
        simple_scatter = 2.0_num*xr2 - 1.0_num
    END IF

  END FUNCTION simple_scatter



  SUBROUTINE clear_Z_cache

    ! Deallocate pre-calculated variable arrays

    DEALLOCATE(urb_d1, urb_d2, urb_d3, urb_d4)
    DEALLOCATE(urb_iZ23, urb_theta1, urb_theta2, urb_Zeff)

  END SUBROUTINE clear_Z_cache

#endif
END MODULE hy_elastic_urban
