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
! hy_resistivity.F90
!
! This module holds all the subroutines related to resistivity, and also the
! subroutines for updating the ionisation state and Coulomb logarithm (as these
! are both used in the reduced Lee-More resistivity model)

MODULE hy_resistivity
#ifdef HYBRID

  USE hy_shared

  IMPLICIT NONE

  REAL(num), PRIVATE, PARAMETER :: rlm_c_hot = &
      128.0_num*epsilon0**2/q0**4 * SQRT(0.5_num*pi*m0)
  REAL(num), PRIVATE, PARAMETER :: rlm_c_ve = 3.0_num*kb/m0
  REAL(num), PRIVATE, PARAMETER :: rlm_c_sphere = 4.0_num*pi/3.0_num

CONTAINS

  SUBROUTINE update_ionisation

    ! Loops over all cells and updates the ionisation charge state of the
    ! background solid. In the case of compound solids, the ion properties are
    ! averaged over species

    INTEGER :: ix

    DO ix = 1-ng, nx+ng
      ion_charge(ix) = thomas_fermi_ionisation(ix)
    END DO

  END SUBROUTINE update_ionisation



  SUBROUTINE update_coulomb_logarithm

    ! Loops over all cells and updates the Coulomb logarithm of the background
    ! solid. In the case of compound solids, the ion properties are averaged
    ! over species

    INTEGER :: ix

    DO ix = 1-ng, nx+ng
      ion_cou_log(ix) = lee_more_coulomb_log(ix)
    END DO

  END SUBROUTINE update_coulomb_logarithm



  SUBROUTINE update_resistivity

    ! Loops over the resistivity grid and updates each point based on the
    ! resistivity model in that cell

    INTEGER :: ix

    DO ix = 1-ng, nx+ng
      SELECT CASE (resistivity_model(ix))
      CASE(c_resist_vacuum)
        resistivity(ix) = 0.0_num
      CASE(c_resist_milchberg)
        resistivity(ix) = calc_resistivity_milchberg(hy_Te(ix))
      CASE(c_resist_plastic)
        resistivity(ix) = calc_resistivity_plastic(hy_Te(ix))
      CASE(c_resist_rlm)
        resistivity(ix) = calc_resistivity_rlm(ix)
      END SELECT
    END DO

  END SUBROUTINE update_resistivity



  FUNCTION calc_resistivity_milchberg(Te)

    ! Calculates the resistivity for the background temperature Te, using the
    ! data from H. M. Milchberg, et al, 1988. Phys. Rev. Lett., 61(20), p.2364.
    !
    ! This applies to Al only. The fit to this data comes from Davies, et al,
    ! (2002). Phys. Rev. E, 65(2), 026407

    REAL(num), INTENT(IN) :: Te
    REAL(num) :: Te_eV
    REAL(num) :: calc_resistivity_milchberg

    ! Convert Te from Kelvin to eV
    Te_eV = Kelvin_to_eV * Te

    ! Calculate resistivity
    calc_resistivity_milchberg = Te_eV &
        / (5.0e6_num + 170.0_num*Te_eV**2.5_num + 3.0e5_num*Te_eV)

  END FUNCTION calc_resistivity_milchberg



  FUNCTION calc_resistivity_plastic(Te)

    ! Calculates the resistivity for the background temperature Te, using the
    ! heuristic model for plastic from Davies, et al, (1999). Phys. Rev. E,
    ! 59(5), 6032.

    REAL(num), INTENT(IN) :: Te
    REAL(num) :: Te_eV
    REAL(num) :: calc_resistivity_plastic

    ! Convert Te from Kelvin to eV
    Te_eV = Kelvin_to_eV * Te

    ! Calculate resistivity
    calc_resistivity_plastic = 1.0_num / (4.3e5_num + 1.3e3_num*Te_eV**1.5_num)

  END FUNCTION calc_resistivity_plastic



  FUNCTION calc_resistivity_rlm(ix)

    ! Calculates the resistivity for the background temperature Te, using a
    ! reduced form of the Lee-More resistivity. The approximations used here are
    ! similar to those used in the hybrid PIC code Zephyros.
    !
    ! The full Lee-More model involves calculating the chemical potential, and
    ! can be found in Lee, More (1984). Phys. Fluids. 27(5). Note that this
    ! paper is written in Gaussian CGS units, whereas EPOCH is in SI

    INTEGER, INTENT(IN) :: ix
    REAL(num) :: rlm_Z_star, rlm_cou_log, rlm_Te, rlm_ni
    REAL(num) :: fermi_temp, ve, ion_term, idebye_huckel2
    REAL(num) :: ion_sphere_rad, debye_huckel, class_impact, uncert_lim
    REAL(num) :: b_min, b_max, cou_log
    REAL(num) :: tau_cold, tau_hot, tau
    REAL(num) :: calc_resistivity_rlm

    ! Copy out array variables
    rlm_Z_star = ion_charge(ix)
    rlm_cou_log = ion_cou_log(ix)
    rlm_Te = hy_Te(ix)
    rlm_ni = ion_ni(ix)

    ! Ion sphere radius (interatomic distance)
    ion_sphere_rad = (rlm_c_sphere*rlm_ni)**(-1.0_num/3.0_num)

    ! Thermal electron velocity (see between Lee-More eq. 20 & 21)
    ve = SQRT(rlm_c_ve * rlm_Te)

    ! Collision time just above melting point (Lee-More section 3 end), with fit
    ! parameter rlm_1
    tau_cold = ion_sphere_rad / ve * rlm_1

    ! Collision time in the high temperature limit (Lee-More eq. 27 & 28a)
    ! This is tau*A^alpha, as used in Lee-More 23a
    tau_hot = rlm_c_hot * (kb*rlm_Te)**1.5_num &
        / (rlm_Z_star**2 * rlm_ni * rlm_cou_log)

    ! Approximation: as temperature falls, tau_hot becomes unphysically small,
    ! so switch to tau_cold
    tau = MAX(tau_cold, tau_hot)

    ! Inverted conductivity eq. 23a. Added a global scaling fit parameter
    calc_resistivity_rlm = m0/(rlm_Z_star * rlm_ni * tau * q0**2) * rlm_2

  END FUNCTION calc_resistivity_rlm



  FUNCTION thomas_fermi_ionisation(ix)

    ! Calculates the average ionisation state of the background ions using the
    ! algorithm outlined in table IV, More, R. M. (1985). "Pressure ionization,
    ! resonances, and the continuity of bound and free states". In Advances in
    ! atomic and molecular physics (Vol. 21, pp. 305-356)
    !
    ! For ease of comparison, we use their variable names with the prefix tf_
    !
    ! This subroutine uses an average global background, where solid_array
    ! properties have been averaged

    INTEGER, INTENT(IN) :: ix
    REAL(num) :: tf_T0, tf_Tf
    REAL(num) :: tf_A, tf_B, tf_C
    REAL(num) :: tf_Q1, tf_Q
    REAL(num) :: tf_x
    REAL(num) :: thomas_fermi_ionisation

    tf_T0 = hy_Te(ix) * Kelvin_to_eV * ion_Z_avg(ix)**(-4.0_num/3.0_num)
    tf_Tf = tf_T0/(1.0_num + tf_T0)

    tf_A = 0.003323_num*tf_T0**0.9718_num + 9.26148e-5_num*tf_T0**3.10165_num
    tf_B = -EXP(-1.7630_num + 1.43175_num*tf_Tf + 0.31546_num*tf_Tf**7)
    tf_C = -0.366667_num*tf_Tf + 0.983333_num

    tf_Q1 = tf_A * ion_reduced_density(ix)**tf_B
    tf_Q = (ion_reduced_density(ix)**tf_C + tf_Q1**tf_C)**(1.0_num/tf_C)

    tf_x = 14.3139_num * tf_Q**0.6624_num
    thomas_fermi_ionisation = ion_Z_avg(ix)*tf_x &
        / (1.0_num + tf_x + SQRT(1.0_num + 2.0_num*tf_x))

  END FUNCTION thomas_fermi_ionisation



  FUNCTION lee_more_coulomb_log(ix)

    ! Calculates the coulomb logarithm using the method outlined in Lee, More
    ! (1984). Phys. Fluids. 27(5). Note that this paper is written in Gaussian
    ! CGS units, whereas EPOCH is in SI

    INTEGER, INTENT(IN) :: ix
    REAL(num) :: rlm_Z_star, rlm_Ti, rlm_Te, rlm_ni
    REAL(num) :: fermi_temp, ve, ion_term, idebye_huckel2
    REAL(num) :: ion_sphere_rad, debye_huckel, class_impact, uncert_lim
    REAL(num) :: b_min, b_max
    REAL(num) :: lee_more_coulomb_log

    ! Copy out array variables
    rlm_Z_star = ion_charge(ix)
    rlm_Te = hy_Te(ix)
    rlm_ni = ion_ni(ix)

    IF (use_ion_temp) THEN
      rlm_Ti = hy_Ti(ix)
    ELSE
      rlm_Ti = rlm_Te
    END IF

    ! Ion sphere radius (interatomic distance)
    ion_sphere_rad = (4.0_num*pi*rlm_ni/3.0_num)**(-1.0_num/3.0_num)

    ! Fermi temperature
    fermi_temp = h_bar**2/(2.0_num*m0*q0*kb) &
        * (3.0_num*pi*rlm_Z_star*rlm_ni)**(2.0_num/3.0_num)

    ! Debye-Huckel screening length (Lee-More eq. 19)
    IF (rlm_Ti < c_tiny) THEN
      ion_term = 0.0_num
    ELSE
      ion_term = rlm_Z_star/rlm_Ti
    END IF
    idebye_huckel2 = rlm_Z_star*rlm_ni*q0**2/(epsilon0*kb) &
        * (1.0_num/(kb*SQRT(rlm_Te**2 + fermi_temp**2)) + ion_term)
    debye_huckel = 1.0_num/SQRT(idebye_huckel2)

    ! Thermal electron velocity (see between Lee-More eq. 20 & 21)
    ve = SQRT(3.0_num*kb*rlm_Te/m0)

    ! Classical distance of closest approach (Lee-More eq. 20)
    class_impact = rlm_Z_star*q0**2/(4.0_num*pi*epsilon0*m0*ve**2)

    ! High energy b_min limit is uncertainty principle (Lee-More eq. 21)
    uncert_lim = h_planck/(2.0_num*m0*ve)

    ! Impact parameters
    ! Minimum (Lee-More eq. 22)
    b_min = MAX(class_impact, uncert_lim)
    ! Maximum (see between Lee-More eq. 19 and 20 - Debye-Huckel screening
    ! breaks down when less than the interatomic distance)
    b_max = MAX(debye_huckel, ion_sphere_rad)

    ! Coulomb logarithm (Lee-More eq. 17 and following discussion)
    lee_more_coulomb_log = MAX(0.5_num*LOG(1.0_num + (b_max/b_min)**2), 2.0_num)

  END FUNCTION lee_more_coulomb_log

#endif
END MODULE hy_resistivity
