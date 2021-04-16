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
! hy_elastic_davies.F90
!
! This module performs the elastic scatter process using the Davies algorithm.

MODULE hy_elastic_davies
#ifdef HYBRID

  USE particles
  USE hy_shared

  IMPLICIT NONE

  ! Precalculate repeated variables
  REAL(num), PRIVATE, PARAMETER :: ipart_mc  = 1.0_num / mc0

CONTAINS

  SUBROUTINE Davies_elastic_scatter

    ! Calculates collisional scattering using the method discussed in  Davies,
    ! et al, (2002). Phys. Rev. E, 65(2), 026407

    INTEGER :: ispecies
    INTEGER(i8) :: ipart
    REAL(num) :: sum_dtheta, part_ne
    REAL(num) :: p, gamma, v
    REAL(num) :: rand_scatter, delta_theta, delta_phi
    REAL(num) :: part_x, part_y, part_z
    INTEGER :: i_sol
    TYPE(particle), POINTER :: current, next

    ! Generate a normally distributed random number for the scatter angle.
    ! The argument here refers to the standard deviation. The mean is
    ! assumed zero.
    rand_scatter = random_box_muller(1.0_num)

    ! Loop over all electron species
    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      IF (.NOT. species_list(ispecies)%species_type == c_species_id_electron) &
          CYCLE

      ! Loop over all particles in the species
      DO ipart = 1, species_list(ispecies)%attached_list%count
        next => current%next

        ! Extract particle variables
        part_x = current%part_pos(1) - x_grid_min_local
        part_y = current%part_pos(2) - y_grid_min_local
        part_z = current%part_pos(3) - z_grid_min_local
        p = SQRT(current%part_p(1)**2 + current%part_p(2)**2 &
            + current%part_p(3)**2)
        ! No scatter for immobile electrons
        IF (p < c_tiny) THEN
          current => next
          CYCLE
        END IF
        gamma = SQRT((p * ipart_mc)**2 + 1.0_num)
        v = p * ipart_mc * c / gamma

        ! Extract solid variables averaged over all solids present in this cell
        sum_dtheta = 0.0_num
        DO i_sol = 1, solid_count
          CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, &
              part_ne, solid_array(i_sol)%el_density)

          ! If ne is zero, then the solid is not in this cell, so don't
          ! calculate the constants
          IF (part_ne < TINY(1.0_num)) CYCLE

          ! Very low energy particles (~50 eV) will make ln_lambda_S go negative
          ! causing a NaN in delta_theta. The hybrid mode isn't designed for low
          ! energy electrons, so delta_theta will be ignored in this case
          sum_dtheta = sum_dtheta + solid_array(i_sol)%theta_fac * part_ne &
              * LOG(MAX(1.0_num, solid_array(i_sol)%ln_s * p)) * dt / v
        END DO

        ! Collisional changes to direction
        delta_theta = SQRT(sum_dtheta) / p * rand_scatter
        delta_phi = 2.0_num * pi * random()

        ! Apply rotation
        CALL rotate_p(current, COS(delta_theta), delta_phi, p)

        current => next
      END DO
    END DO

  END SUBROUTINE Davies_elastic_scatter

#endif
END MODULE hy_elastic_davies
