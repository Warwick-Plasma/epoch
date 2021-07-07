! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2012      Martin Ramsay <M.G.Ramsay@warwick.ac.uk>
! Copyright (C) 2009      Chris Brady <C.S.Brady@warwick.ac.uk>
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
! hy_laser.F90
!
! This module holds the hybrid particle injector. This allows the user to inject
! electrons into the simulation based on expected distributions from laser
! parameters.
!
! The structure is a mix between the physics_packages/injectors.F90 file, and
! the laser.f90 file.

MODULE hy_laser

  USE partlist
  USE evaluator
  USE particles

  IMPLICIT NONE

  REAL(num), PARAMETER, PRIVATE :: c_a0_term = q0 /(m0 * c) &
      * SQRT(2.0_num/(epsilon0 * c))

  REAL(num), PRIVATE :: particle_energy

CONTAINS

  SUBROUTINE init_hy_laser(boundary, hy_laser)

    INTEGER, INTENT(IN) :: boundary
    TYPE(hy_laser_block), INTENT(INOUT) :: hy_laser

    hy_laser%boundary = boundary
    hy_laser%ppc = -1
    hy_laser%species = -1
    hy_laser%omega_func_type = -1
    hy_laser%use_time_function = .FALSE.
    hy_laser%use_omega_function = .FALSE.
    hy_laser%intensity = 0.0_num
    hy_laser%efficiency = -1.0_num
    hy_laser%omega = -1.0_num
    hy_laser%t_start = 0.0_num
    hy_laser%t_end = t_end
    hy_laser%has_t_end = .FALSE.
    NULLIFY(hy_laser%next)

    ! Physical model choices for injection
    hy_laser%mean = -1
    hy_laser%e_dist = -1
    hy_laser%ang_dist = -1

    ! User defined energies and weights (if both are set, ignore I, omega, eff)
    hy_laser%user_mean_KE = -1.0_num
    hy_laser%user_weight = -1.0_num
    hy_laser%ignore_las = .FALSE.

    ! Angular distribution variables
    hy_laser%user_theta_max = pi
    hy_laser%cos_n_power = 0.0_num
    hy_laser%top_hat_L = 0.5_num
    hy_laser%mean_mult = 10.0_num
    hy_laser%sheng_angle = 0.0_num
    hy_laser%theta_mean = 0.0_num
    hy_laser%phi_mean = 0.0_num
    hy_laser%use_moore_max = .FALSE.
    hy_laser%use_sheng_dir = .FALSE.

  END SUBROUTINE init_hy_laser



  SUBROUTINE deallocate_hy_laser(hy_laser)

    TYPE(hy_laser_block), POINTER :: hy_laser

    IF (hy_laser%use_time_function) &
        CALL deallocate_stack(hy_laser%time_function)
    IF (hy_laser%use_omega_function) &
        CALL deallocate_stack(hy_laser%omega_function)
    DEALLOCATE(hy_laser)

  END SUBROUTINE deallocate_hy_laser



  SUBROUTINE deallocate_hy_lasers

    TYPE(hy_laser_block), POINTER :: current, next

    current => hy_lasers
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_hy_laser(current)
      current => next
    END DO

  END SUBROUTINE deallocate_hy_lasers



  SUBROUTINE attach_hy_laser(hy_laser)

    ! Subroutine to attach a created hy_laser object to the correct boundary

    INTEGER :: boundary
    TYPE(hy_laser_block), POINTER :: hy_laser, current

    boundary = hy_laser%boundary

    n_lasers(boundary) = n_lasers(boundary) + 1

    IF (ASSOCIATED(hy_lasers)) THEN
      ! If we already have hybrid lasers in hy_lasers, cycle to the end of the
      ! list and add the current hy_laser to the end
      current => hy_lasers
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      END DO
      current%next => hy_laser

    ELSE
      ! The hy_lasers list is empty, our current hy_laser is the start of the
      ! list
      hy_lasers => hy_laser
    END IF

  END SUBROUTINE attach_hy_laser



  SUBROUTINE populate_pack_from_hy_laser(hy_laser, parameters)

    ! This routine populates the constant elements of a parameter pack
    ! from a hy_laser. When evaluating spatially varying functions, we need to
    ! feed in an ix,iy pair - this sorts out the spatial coordinate associated
    ! with which boundary the laser is on (position along the boundary is
    ! separate)

    TYPE(hy_laser_block), POINTER :: hy_laser
    TYPE(parameter_pack), INTENT(INOUT) :: parameters

    parameters%pack_ix = 0

    SELECT CASE(hy_laser%boundary)
      CASE(c_bd_x_min)
        parameters%pack_ix = 0
      CASE(c_bd_x_max)
        parameters%pack_ix = nx
    END SELECT

  END SUBROUTINE populate_pack_from_hy_laser



  FUNCTION hy_laser_time_profile(hy_laser)

    ! Return the temporal laser intensity weighting

    TYPE(hy_laser_block), POINTER :: hy_laser
    REAL(num) :: hy_laser_time_profile
    INTEGER :: err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_hy_laser(hy_laser, parameters)
    IF (hy_laser%use_time_function) THEN
      hy_laser_time_profile = evaluate_with_parameters(hy_laser%time_function, &
          parameters, err)
      RETURN
    END IF

  END FUNCTION hy_laser_time_profile



  SUBROUTINE hy_laser_update_omega(hy_laser)

    ! Returns the angular frequency of the laser at this timestep (not spatially
    ! varying)

    TYPE(hy_laser_block), POINTER :: hy_laser
    INTEGER :: err
    TYPE(parameter_pack) :: parameters

    ! If not time-varying, use_omega_function is set to false to prevent
    ! recalculating the constant angular frequency
    IF (.NOT. hy_laser%use_omega_function) RETURN

    err = 0
    CALL populate_pack_from_hy_laser(hy_laser, parameters)
    hy_laser%omega = &
        evaluate_with_parameters(hy_laser%omega_function, parameters, err)
    ! User has specified frequency, not angular frequency
    IF (hy_laser%omega_func_type == c_of_freq) &
        hy_laser%omega = 2.0_num * pi * hy_laser%omega
    ! Use has specified wavelength, not angular frequency
    IF (hy_laser%omega_func_type == c_of_lambda) &
        hy_laser%omega = 2.0_num * pi * c / hy_laser%omega

  END SUBROUTINE hy_laser_update_omega



  SUBROUTINE update_hy_laser_omegas

    ! Loop through every laser and update the omega values

    TYPE(hy_laser_block), POINTER :: current

    current => hy_lasers
    DO WHILE(ASSOCIATED(current))
      IF (current%use_omega_function) CALL hy_laser_update_omega(current)
      current => current%next
    END DO

  END SUBROUTINE update_hy_laser_omegas



  SUBROUTINE run_hybrid_lasers

    ! Cycle through the hy_lasers list containing all hy_laser injectors, and
    ! run each hy_laser one at a time

    TYPE(hy_laser_block), POINTER :: current

    CALL update_hy_laser_omegas

    current => hy_lasers
    DO WHILE(ASSOCIATED(current))
      CALL run_single_laser(current)
      current => current%next
    END DO

  END SUBROUTINE run_hybrid_lasers



  SUBROUTINE run_single_laser(laser)

    ! Injects electrons into the simulation due to a given laser, based on the
    ! current intensity, wavelength of this laser, through the specified
    ! boundary.

    TYPE(hy_laser_block), POINTER :: laser
    INTEGER :: boundary
    REAL(num) :: bdy_pos, bdy_space
    TYPE(particle), POINTER :: new, weight_part
    TYPE(particle_list) :: plist
    TYPE(parameter_pack) :: parameters
    REAL(num), DIMENSION(3) :: p_dir
    INTEGER :: ipart
    REAL(num) :: vol
    REAL(num) :: time_frac, intens_cell
    REAL(num) :: mean_energy, energy_cell
    REAL(num) :: el_number, el_weight
    REAL(num) :: p_mag, v_bdy
    LOGICAL :: got_list
    REAL(num) :: sum_weight, weight_factor

    ! No injections if the laser is off
    IF (time < laser%t_start .OR. time > laser%t_end) RETURN

    ! If you have a moving window that has started moving then unless you
    ! EXPLICITLY give a t_end value to the injector stop the injector
    IF (move_window .AND. window_started .AND. .NOT. laser%has_t_end) &
        RETURN

    ! These lines check if the injection boundary is on the current rank
    boundary = laser%boundary
    IF (.NOT. is_boundary(boundary)) RETURN

    ! Precalculate boundary dependent parameters (repeated each boundary)
    IF (boundary == c_bd_x_min) THEN
      parameters%pack_ix = 0                ! Index coordinate of boundary
      bdy_pos = x_min                       ! Real coordinate of boundary
      bdy_space = -dx                       ! Transformation to behind boundary
      p_dir = (/1.0_num, 0.0_num, 0.0_num/) ! Initial momentum direction
    ELSE IF (boundary == c_bd_x_max) THEN
      parameters%pack_ix = nx
      bdy_pos = x_max
      bdy_space = dx
      p_dir = (/-1.0_num, 0.0_num, 0.0_num/)
    ELSE
      RETURN
    END IF

    ! Cell volume (dy and dz = 1m in epoch1d)
    vol = ABS(bdy_space)

    ! This distribution needs to change the weights of particles. This flag
    ! tells us when weight_part points to the first particle
    IF (laser%e_dist == e_dist_exp_weight) THEN
      got_list = .FALSE.
    END IF

    ! Evaluate functions for the current time
    IF (.NOT. laser%ignore_las) CALL hy_laser_update_omega(laser)
    time_frac = MAX(hy_laser_time_profile(laser), 0.0_num)
    IF (time_frac < c_tiny) RETURN

    CALL create_empty_partlist(plist)

    ! Don't inject low weight particles
    parameters%use_grid_position = .TRUE.

    ! Mean electron kinetic energy for each cell
    IF (laser%mean == c_mean_E_val) THEN
      mean_energy = laser%user_mean_KE
    ELSE
      intens_cell = laser%intensity * time_frac
      mean_energy = get_mean_energy(intens_cell, laser%omega, laser%mean)
    END IF

    ! Corresponding number of real particles
    IF (laser%e_dist == e_dist_mono_weight) THEN
      ! User defined weight, ignore laser parameters
      el_weight = laser%user_weight * time_frac
    ELSE
      energy_cell = intens_cell * dt
      el_number = energy_cell * laser%efficiency / mean_energy
      IF (laser%e_dist == e_dist_exp_weight) THEN
        ! Track total weight for special weight treatment
        sum_weight = 0.0_num
      ELSE
        ! All other distributions have constant particle weight
        el_weight = el_number/REAL(laser%ppc, num)
      END IF
    END IF

    ! Create particles
    DO ipart = 1, laser%ppc
      CALL create_particle(new)

      ! Sample momentum of new e-
      p_mag = sample_p_mag(mean_energy, laser)
      new%part_p = p_dir * p_mag

      ! Rotate mean e- direction
      new%particle_energy = SQRT((p_mag*c)**2 + m0c2**2)
      CALL sample_rotation(new, laser, p_mag)

      ! Push particle back a random fraction of the distance this particle
      ! would travel in the boundary direction in one timestep
      v_bdy = new%part_p(1) * c**2 / new%particle_energy
      new%part_pos = bdy_pos - random() * v_bdy * dt

      IF (laser%e_dist == e_dist_exp_weight) THEN
        ! Exponentially distribute the weights
        new%weight = EXP(-(new%particle_energy-m0c2)/mean_energy)
        sum_weight = sum_weight + new%weight
      ELSE
        ! Uniform weight
        new%weight = el_weight
      END IF

#ifdef PER_PARTICLE_CHARGE_MASS
      new%charge = species_list(laser%species)%charge
      new%mass = species_list(laser%species)%mass
#endif
      CALL add_particle_to_partlist(plist, new)
    END DO

    ! Final normalisation to distributed weights, such that the sum is
    ! el_number
    IF (laser%e_dist == e_dist_exp_weight) THEN
      weight_factor = el_number/sum_weight
      DO ipart = 1, laser%ppc
        ! Point to the next particle in the list
        IF (.NOT. got_list) THEN
          weight_part=>plist%head
          got_list = .TRUE.
        ELSE
          weight_part=>weight_part%next
        END IF
        ! Scale weights to sum to el_number
        weight_part%weight = weight_part%weight*weight_factor
      END DO
    END IF

    CALL append_partlist(species_list(laser%species)%attached_list, plist)

  END SUBROUTINE run_single_laser



  FUNCTION get_mean_energy(intensity, omega, model)

    ! Returns the mean electron kinetic energy corresponding to the requested
    ! model

    REAL(num), INTENT(IN) :: intensity, omega
    INTEGER, INTENT(IN) :: model
    REAL(num) :: get_mean_energy, wavelength

    IF (model == c_mean_a0) THEN
      ! KE = a0*m*c²
      get_mean_energy = c_a0_term * SQRT(intensity) / omega * m0c2
    ELSE IF (model == c_mean_wilks) THEN
      ! Wilks scaling: KE = m*c²(SQRT(1 + I*lambda²/(1.37e18W/cm²*um²))-1)
      ! Wilks, S.C. et al (1992). Phys. Rev. Lett., 69(9), 1383.
      wavelength = 2.0_num * pi * c / omega
      get_mean_energy = m0c2 &
          * (SQRT(1.0_num + intensity*wavelength**2/1.37e10_num) - 1.0_num)
    END IF

  END FUNCTION get_mean_energy



  FUNCTION sample_p_mag(KE_mean, laser)

    ! Samples a momentum drawn from the energy distribution specified by the
    ! user

    REAL(num), INTENT(IN) :: KE_mean
    TYPE(hy_laser_block), POINTER, INTENT(IN) :: laser
    REAL(num) :: sample_p_mag, sample_KE

    IF (laser%e_dist == e_dist_mono .OR. laser%e_dist == e_dist_mono_weight) &
        THEN
      ! Mono-energetic distribution, return momentum corresponding to KE_mean
      sample_KE = KE_mean

    ELSE IF (laser%e_dist == e_dist_tophat) THEN
      ! Top-hat distribution, specified with some fractional width L
      sample_KE = KE_mean &
          * (laser%top_hat_L * (2.0_num*random() - 1.0_num) + 1.0_num)

    ELSE IF (laser%e_dist == e_dist_exp) THEN
      ! Exponential distribution, with scale length KE_mean
      sample_KE = -KE_mean * LOG(1.0_num - random())

    ELSE IF (laser%e_dist == e_dist_exp_weight) THEN
      ! Uniform distribution from 0 to KE_Mean * mean_mult (weights will be exp)
      sample_KE = random() * laser%mean_mult * KE_mean

    ELSE IF (laser%e_dist == e_dist_mono_las_weight) THEN
      ! Ignore KE_mean (used for weight calc) and use a user-defined KE
      sample_KE = laser%las_weight_KE

    END IF

    sample_p_mag = SQRT((sample_KE+m0c2)**2 - m0c2**2)/c

  END FUNCTION sample_p_mag



  SUBROUTINE sample_rotation(part, laser, part_p)

    ! Rotate the injected particle. Firstly, sample a theta/phi pair from the
    ! chosen angular distribution, then apply a global transformation to theta
    ! and phi (user input) to model electron beams which aren't injected
    ! perpendicular to the boundary

    TYPE(particle), POINTER :: part
    TYPE(hy_laser_block), POINTER :: laser
    REAL(num) :: part_p
    REAL(num) :: part_gamma, sh1, sh2, theta, phi
    REAL(num) :: theta_max, costheta, npow

    ! Get electron gamma factor (particle_energy is defined in run_single_laser
    ! and has scope with all subroutines within this module)
    IF (laser%use_moore_max .OR. laser%use_sheng_dir) THEN
      part_gamma = particle_energy / m0c2
    END IF

    ! From Sheng, Z.M., et al. (2000). Phys. Rev. Lett., 85(25), 5340. This
    ! paper gives an analytic form of the injection direction of electrons from
    ! a laser coming in at angle laser%sheng_angle to the surface normal. The
    ! form appearing in this code is from eq. 2, the case where delta*Phi = 0.
    ! Apply distribution rotations on top of this Sheng direction if requested
    IF (laser%use_sheng_dir) THEN
      sh1 = 1.0_num/(part_gamma - 1.0_num)
      sh2 = 2.0_num + (part_gamma + 1.0_num)/TAN(laser%sheng_angle)**2
      theta = ATAN(1.0_num/SQRT(sh1 * sh2))
      phi = 0.0_num
      CALL rotate_p(part, COS(theta), phi, part_p)
    END IF

    ! User can define their own theta and phi transformations. Apply
    ! distribution rotations on top of these
    CALL rotate_p(part, COS(laser%theta_mean), laser%phi_mean, part_p)

    ! From Moore, C.I., et al. (1999). Phys. Rev. Lett., 82(8), 1688. This
    ! uses an energy-dependent maximum theta for the angular distribution of
    ! TAN(theta) = SQRT(2/(gamma+1)). If this is larger than the user-defined
    ! theta cut-off, use the user value
    IF (laser%use_moore_max) THEN
      theta_max = &
          MIN(ATAN(SQRT(2.0_num/(part_gamma - 1.0_num))), laser%user_theta_max)
    ELSE
      theta_max = laser%user_theta_max
    END IF

    ! Sample theta from distribution
    IF (laser%ang_dist == c_ang_uniform) THEN
      ! Uniformly distributed angular scatter between 0 and theta_max
      costheta = 1.0_num - random()*(1.0_num - COS(theta_max))

    ELSE IF (laser%ang_dist == c_ang_cos) THEN
      ! PDF: f(Omega) = A * (COS(theta))**npow
      ! Omega: solid angle [sr]
      ! npow: User defined number
      npow = laser%cos_n_power
      costheta = (1.0_num - (1.0_num - COS(theta_max)**(npow + 1.0_num)) &
          * random())**(1.0_num/(npow+1.0_num))

    ELSE IF (laser%ang_dist == c_ang_beam) THEN
      ! No further scatter
      RETURN
    END IF

    ! Apply rotation from distribution angles
    phi = 2.0_num * pi * random()
    CALL rotate_p(part, costheta, phi, part_p)

  END SUBROUTINE sample_rotation

END MODULE hy_laser
