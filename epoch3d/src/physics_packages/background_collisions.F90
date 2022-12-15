! Copyright (C) 2009-2019 University of Warwick
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
! background_collisions module
!
! To explain the purpose of this module, it helps to explain the other
! EPOCH collisions module. The Nanbu-Perez model used for the standard EPOCH
! collision routines comes from work in the papers:
!   Nanbu (1997) - "Theory of cumulative small-angle collisions in plasmas"
!   Perez (2012) - "Improved modelling of relativistic collisions and
!                   collisional ionisation in particle-in-cell codes"
! where the Nanbu paper describes how much a particle is scattered after
! traversing a certain distance through a scattering medium, and Perez discusses
! a relativistic extension and an implementation to PIC codes.
!
! The collisions module in EPOCH is slow as the Nanbu scatter angle depends on
! the relative velocity between spatially local particles. Hence, particles must
! be attributed to a local cell, and paired to estimate a range of relative
! velocities. However, when the speed of one particle is far greater than the
! second, then the relative velocity is approximately the fast particle speed,
! and pairing is not required.
!
! In this module, the scattering of particles in a fast species is calculated
! from a slow background species. The species to be considered "background" are
! user-defined in the input deck.

MODULE background_collisions

  USE partlist
  USE calc_df
  USE setup

  IMPLICIT NONE

  INTEGER, ALLOCATABLE :: back_id(:), back_id_all(:), id2array_all(:)
  INTEGER :: back_count_all
  LOGICAL :: use_lab
  REAL(num), ALLOCATABLE :: p1_lab(:,:)
  REAL(num), ALLOCATABLE :: back_ni_all(:,:,:,:), back_cou_log(:,:,:)

  ! Collisions only calculated once every k time-steps: dt_coll = k * dt
  REAL(num) :: dt_coll

  ! s12_a =  dt / (4 * pi * eps0^2 *c^4)
  REAL(num) :: s12_a

  ! s12_b = s12_a * q1^2 * q2^2 / (m1 * m2)
  REAL(num), ALLOCATABLE :: s12_b(:)

  ! Variables from run_background_collisions needed by calc_s12
  REAL(num) :: m1, m2, gamma1, v1, p1, n2, cou_log_12

  ! Variables from calc_s12 needed in run_background_collisions
  REAL(num) :: p1_com, v_com, gamma_com, gamma1_com

CONTAINS

  SUBROUTINE setup_background_collisions

    INTEGER :: ispecies, i_back, i1, i2, i_bi
    REAL(num) :: ke1_low, ke1_high, p1_low, p1_high, p1_mid, tol

    ! Are any background species present?
    use_background_collisions = .FALSE.
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%coll_background) THEN
        use_background_collisions = .TRUE.
        EXIT
      END IF
    END DO

    IF (use_background_collisions) THEN
      ! Count total number of background species
      back_count_all = 0
      DO ispecies = 1, n_species
        IF (species_list(ispecies)%coll_background) THEN
          back_count_all = back_count_all + 1
        END IF
      END DO

      ! Create two arrays:
      ! back_id_all: Species ID of all species which act as a background
      ! id2array_all: Reverse array, giving the back_id_all index for a given
      !               species ID
      ALLOCATE(back_id_all(back_count_all))
      ALLOCATE(id2array_all(n_species))
      i_back = 1
      DO ispecies = 1, n_species
        IF (species_list(ispecies)%coll_background) THEN
          back_id_all(i_back) = ispecies
          id2array_all(ispecies) = i_back
          i_back = i_back + 1
        ELSE
          id2array_all(ispecies) = 0
        END IF
      END DO

      ! Time interval between collision calculations
      dt_coll = REAL(coll_n_step, num) * dt

      ! Number of steps between background variable update
      IF (coll_subcycle_back) THEN
        back_n_step = MAX(NINT(back_update_dt / dt), 1)
      ELSE
        ! Recalculate every step if back_update_dt is not given
        back_n_step = 1
      END IF

      ! Need to calculate background variables on the first run
      coll_back_recalc = .TRUE.

      ! Pre-define constant factor s12_a, used for s12 calculation
      s12_a = dt_coll / (4.0_num * pi * epsilon0**2 * c**4)

      ! For each fast-background pair, calculate the fast particle momentum, p1
      ! below which there is negligible difference between the lab and CoM (*)
      ! frames
      ! This is defined such that |p1 - p1*|/p1 < rel_cutoff
      ALLOCATE(p1_lab(n_species, n_species))
      p1_lab = 0.0_num
      IF (use_rel_cutoff) THEN
        DO i2 = 1, n_species
          DO i1 = 1, n_species
            ! Ensure i1 is the fast species and i2 is the background
            IF (coll_pairs_state(i1, i2) == c_coll_background_2nd) THEN

              ! Mass variables temporarily saved to module variables
              m1 = species_list(i1)%mass
              m2 = species_list(i2)%mass

              ! Use bi-section to find the critical momentum, starting between
              ! 1 eV and 100 GeV kinetic energy
              ke1_low = 1.0_num * q0
              ke1_high = 100.0e9_num * q0

              ! Bi-section algorithm runs until fractional difference between
              ! |p1 - p1*|/p1 and rel_cutoff is less than 0.001
              p1_low = SQRT(((ke1_low + m1 * c**2) / c)**2 - m1**2 * c**2)
              p1_high = SQRT(((ke1_high + m1 * c**2) / c)**2 - m1**2 * c**2)
              i_bi = 1
              CALL p1_lab_bisection(p1_low, p1_high, p1_mid, tol)
              DO WHILE (tol > 1.0e-3_num)
                CALL p1_lab_bisection(p1_low, p1_high, p1_mid, tol)
                ! Stop if not converging
                i_bi = i_bi + 1
                IF (i_bi == 50) THEN
                  ! If convergence has failed, assume CoM-frame always important
                  p1_mid = 0.0_num
                  EXIT
                END IF
              END DO

              ! Save momentum, below which collision can be treated in lab frame
              p1_lab(i1,i2) = p1_mid
            END IF
          END DO
        END DO
      END IF
    END IF

  END SUBROUTINE setup_background_collisions



  SUBROUTINE p1_lab_bisection(p1_low, p1_high, p1_mid, tol)

    ! A bisection algorithm used to determine p1, such that:
    !   |p1 - p1*|/p1 = rel_cutoff
    ! Where p1 and p1* are the fast particle momenta in the lab and CoM frame
    ! respectively. Masses m1 and m2 are module variables.
    !
    ! The subroutine reads in upper and lower p1 bounds (p1_low, p1_high), and
    ! calculates the fractional difference between |p1 - p1*|/p1 and rel_cutoff,
    ! in the variable tol. The subroutine also returns p1 halfway between p1_low
    ! and p1_high, which is used as the p1 value if tol is low enough

    REAL(num), INTENT(INOUT) :: p1_low, p1_high
    REAL(num), INTENT(OUT) :: p1_mid, tol
    REAL(num) :: mass_fac, frac

    p1_mid = 0.5_num * (p1_low + p1_high)

    ! Calculate p1* as done in calc_s12() for p1_mid
    gamma1 = SQRT(1.0_num + (p1_mid/(m1*c))**2)
    v1 = p1_mid / (gamma1 * m1)
    mass_fac = (gamma1 * m1 + m2)
    v_com = p1_mid / mass_fac
    gamma_com = 1.0_num / SQRT(1.0_num - (v_com / c)**2)
    gamma1_com = (1.0_num - v_com * v1 / c**2) * gamma_com * gamma1
    p1_com = p1_mid &
        + ((gamma_com - 1.0_num) * v1 - gamma_com * v_com) * m1 * gamma1

    ! Fractional difference between p1_low in lab and CoM frames
    frac = ABS(p1_mid - p1_com) / p1_mid

    ! Fractional difference between frac and rel_cutoff
    tol = ABS(frac - rel_cutoff) / rel_cutoff

    ! High p1 gives high fractional difference, and low p1 has low difference
    IF (frac > rel_cutoff) THEN
      p1_high = p1_mid
    ELSE
      p1_low = p1_mid
    END IF

  END SUBROUTINE p1_lab_bisection



  SUBROUTINE run_background_collisions

    ! Subroutine called by the main EPOCH PIC loop. This loops over all fast
    ! particle species, and applies collisions with slow background particles,
    ! based off user-defined definitions in the input deck. Collisions proceed
    ! according to Nanbu and Perez (see header), but modified to assume the slow
    ! species is stationary.

    INTEGER :: i_fast, i_back, j_back, back_count, spec2
    INTEGER :: cell_x, cell_y, cell_z
    TYPE(particle), POINTER :: fast
    REAL(num) :: q1, q2, m1c2_sq, im1c2
    REAL(num) :: p1_sq, p1_dir(3), energy1
    REAL(num) :: s12, cos_theta, phi, p1_vec_com(3)

    ! Calculate number densities for background species. This may be sub-cycled,
    ! but must be recalculated after a load-balance (coll_back_recalc)
    IF (MODULO(step, back_n_step) == 0 .OR. coll_back_recalc) THEN
      CALL update_background_variables
      IF (coll_back_recalc) coll_back_recalc = .FALSE.
    END IF

    ! Only run collisions if the current step is a multiple of the variable
    ! coll_n_step (default has coll_n_step = 1, collision every step)
    IF (.NOT. MODULO(step, coll_n_step) == coll_n_step - 1) RETURN

    ! Loop over fast particle species
    DO i_fast = 1, n_species

      ! Ignore species which aren't the "fast" in a fast-background pair
      IF (.NOT. species_list(i_fast)%coll_fast) CYCLE

      ! Count background collision partner species
      back_count = 0
      DO i_back = 1, n_species
        IF (.NOT. coll_pairs_state(i_fast, i_back) == c_coll_background_2nd) &
            CYCLE
        back_count = back_count + 1
      END DO

      ! Store species ID of all background species for this current fast species
      ALLOCATE(back_id(back_count))
      j_back = 1
      DO i_back = 1, n_species
        IF (.NOT. coll_pairs_state(i_fast, i_back) == c_coll_background_2nd) &
            CYCLE
        back_id(j_back) = i_back
        j_back = j_back + 1
      END DO

      ! Extract variables related to the fast particle species as a whole
      ! Variable names match those present in Perez et al when possible
      m1 = species_list(i_fast)%mass
      q1 = species_list(i_fast)%charge
      m1c2_sq = (m1 * c**2)**2
      im1c2 = 1.0_num / (m1 * c**2)

      ! Pre-calculate species dependent pre-factor (for speed)
      ALLOCATE(s12_b(back_count))
      DO i_back = 1, back_count
        spec2 = back_id(i_back)
        m2 = species_list(spec2)%mass
        q2 = species_list(spec2)%charge
        s12_b(i_back) = s12_a * (q1*q2)**2 / (m1 * m2)
      END DO

      ! Apply collisions to all particles in the current fast species
      fast => species_list(i_fast)%attached_list%head
      DO WHILE(ASSOCIATED(fast))

        ! Calculate grid index for number density calculation
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x = FLOOR((fast%part_pos(1) - x_grid_min_local) / dx) + 1
        cell_y = FLOOR((fast%part_pos(2) - y_grid_min_local) / dy) + 1
        cell_z = FLOOR((fast%part_pos(3) - z_grid_min_local) / dz) + 1
#else
        cell_x = FLOOR((fast%part_pos(1) - x_grid_min_local) / dx + 1.5_num)
        cell_y = FLOOR((fast%part_pos(2) - y_grid_min_local) / dy + 1.5_num)
        cell_z = FLOOR((fast%part_pos(3) - z_grid_min_local) / dz + 1.5_num)
#endif

        ! Calculate background-averaged Coulomb logarithm in this cell
        IF (coulomb_log_auto) THEN
          cou_log_12 = back_cou_log(cell_x, cell_y, cell_z)
        ELSE
          cou_log_12 = coulomb_log
        END IF

        ! Loop over background species, apply collisions from each
        DO i_back = 1, back_count

          ! Species ID for current background species
          spec2 = back_id(i_back)

          ! Calculate fast particle kinematic variables (will change after coll)
          p1_sq = fast%part_p(1)**2 + fast%part_p(2)**2 + fast%part_p(3)**2
          p1 = SQRT(p1_sq)
          p1_dir = fast%part_p / p1
          energy1 = SQRT(p1_sq * c**2 + m1c2_sq)
          gamma1 = energy1 * im1c2
          v1 = p1 * c**2 / energy1

          ! If there is little difference between lab and CoM frames, use lab
          ! frame variables
          IF (use_rel_cutoff .AND. p1 < p1_lab(i_fast, spec2)) THEN
            use_lab = .TRUE.
          ELSE
            use_lab = .FALSE.
          END IF

          ! Calculate properties of the current background
          n2 = back_ni_all(id2array_all(spec2), cell_x, cell_y, cell_z)
          m2 = species_list(spec2)%mass
          q2 = species_list(spec2)%charge

          ! Calculate s12 contribution from this fast/background pair
          ! CoM variables are saved to module variables during this subroutine
          CALL calc_s12(i_back, s12)

          ! Calculate scatter angle in the CoM frame
          CALL calc_cos_theta(s12, cos_theta)

          ! Rotate momentum in the CoM frame
          p1_vec_com = p1_dir * p1_com
          phi = 2.0_num * pi * random()
          CALL rotate_p_vec(p1_vec_com, cos_theta, phi, p1_com)

          ! Transform momentum back to lab frame if needed (Perez (2012) eq. 13)
          ! Using: vc parallel to p1*, so DOT(vc,p1f*) = vc*(p1f*)*cos(theta)
          IF (.NOT. use_lab) THEN
            fast%part_p = p1_vec_com &
                + p1_dir * ((gamma_com - 1.0_num) * p1_com * cos_theta &
                + gamma1_com * gamma_com * m1 * v_com)
          ELSE
            fast%part_p = p1_vec_com
          END IF
        END DO

        ! Consider collisions for the next particle in the list
        fast => fast%next
      END DO

      DEALLOCATE(back_id, s12_b)
    END DO

  END SUBROUTINE run_background_collisions



  SUBROUTINE update_background_variables

    ! Creates arrays of arrays for variables associated with background species,
    ! which includes:
    !   back_ni_all(i, ix, iy, iz): Number density for background species i at
    !                               cell (ix,iy)
    !   back_cou_log(ix, iy, iz):  The Coulomb logarithm in cell (ix,iy,iz)
    !
    ! Approximation: Coulomb logarithm according to Perez (2012) has a p1 and n2
    ! dependence, so each fast particle should have its own cou_log for each
    ! background species. This is too expensive. Instead, scan through particles
    ! of all fast species and find the average p of the fast particles in each
    ! cell. Do similar averaging for the background species

    INTEGER :: ispecies, i_back, ix, iy, iz, spec2
    REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: num_dens, temp, fast_num_dens
    REAL(num), ALLOCATABLE :: pair_num_dens(:,:,:), pair_num_dens_tot(:,:,:)
    REAL(num), ALLOCATABLE :: q1q2(:,:,:), current_ekbar(:,:,:)
    REAL(num), ALLOCATABLE :: fast_ekbar(:,:,:), fast_p(:,:,:)
    REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: debye_length, b0, delta_x, r_min
    REAL(num) :: qi, q2, mi, ni, ti, bmax, bmin, n2_tot

    IF (ALLOCATED(back_ni_all)) DEALLOCATE(back_ni_all)
    IF (ALLOCATED(back_cou_log)) DEALLOCATE(back_cou_log)

    ALLOCATE(back_ni_all(back_count_all,1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(back_cou_log(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(num_dens(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    IF (coulomb_log_auto) THEN
      ALLOCATE(temp(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(debye_length(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(b0(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(delta_x(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(r_min(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(fast_num_dens(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(current_ekbar(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(fast_ekbar(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(fast_p(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(q1q2(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(pair_num_dens(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(pair_num_dens_tot(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      fast_num_dens = 0.0_num
      fast_ekbar = 0.0_num
      debye_length = 0.0_num
      q1q2 = 0.0_num
      pair_num_dens_tot = 0.0_num
    END IF

    ! Loop over all species present in the simulation. Start summing varialbes
    ! which will later be used for mean calculations
    i_back = 1
    DO ispecies = 1, n_species

      ! Need number density of background species, or all species for Coulomb
      ! logarithm calculation
      IF (coulomb_log_auto .OR. ispecies == back_id_all(i_back)) THEN
        CALL calc_number_density(num_dens, ispecies)
      END IF

      ! If calculating the logarithm, sum (Debye length)^(-2) contributions
      IF (coulomb_log_auto) THEN
        CALL calc_temperature(temp, ispecies)
        qi = species_list(ispecies)%charge
        ! Sum cell by cell, preventing 0 Debye length (stops 1/0 errors)
        DO iz = 1-ng,nz+ng
          DO iy = 1-ng,ny+ng
            DO ix = 1-ng,nx+ng
              ni = MAX(num_dens(ix,iy,iz), 1.0_num)
              ti = MAX(temp(ix,iy,iz), 1.0_num)
              debye_length(ix,iy,iz) = debye_length(ix,iy,iz) &
                  + (qi**2 * ni) / (epsilon0 * kb * ti)
            END DO
          END DO
        END DO
      END IF

      ! Save number density if this is a background species
      IF (ispecies == back_id_all(i_back)) THEN
        back_ni_all(i_back,:,:,:) = num_dens
        i_back = i_back + 1
      END IF

      ! Inspect fast species if calculating the Coulomb logarithm
      IF (coulomb_log_auto .AND. species_list(ispecies)%coll_fast) THEN
        ! Sum num_dens, and <KE>*num_dens for all fast particle species
        fast_num_dens = fast_num_dens + num_dens
        CALL calc_ekbar(current_ekbar, ispecies)
        fast_ekbar = fast_ekbar + current_ekbar * num_dens

        ! Also sum p*num_dens for this species
        mi = species_list(ispecies)%mass
        fast_p = fast_p &
            + SQRT((current_ekbar / (mi * c**2) + 1.0_num)**2 - 1.0_num) &
            * mi * c * num_dens

        ! Sum q1*q2*(n1 + n2) for all fast/background pairs
        DO spec2 = 1, n_species
          IF (coll_pairs_state(ispecies, spec2) == c_coll_background_2nd) THEN
            q2 = species_list(spec2)%charge
            CALL calc_number_density(pair_num_dens, spec2)
            pair_num_dens = pair_num_dens + num_dens
            q1q2 = q1q2 + ABS(qi * q2) * pair_num_dens
            pair_num_dens_tot = pair_num_dens_tot + pair_num_dens
          END IF
        END DO
      END IF
    END DO

    ! Calculate the Coulomb logarithm
    IF (coulomb_log_auto) THEN

      ! Calculate averaged quantities in each cell, with lower limits on
      ! variables to prevent 1/0 errors
      DO iz = 1-ng,nz+ng
        DO iy = 1-ng,ny+ng
          DO ix = 1-ng,nx+ng
            ! <KE> of all fast partciles in a cell
            fast_ekbar(ix,iy,iz) = MAX(fast_ekbar(ix,iy,iz), q0) &
                / MAX(fast_num_dens(ix,iy,iz), 1.0_num)
            ! <momentum> of all fast particles in a cell
            fast_p(ix,iy,iz) = MAX(fast_p(ix,iy,iz), q0 / c) &
                / MAX(fast_num_dens(ix,iy,iz), 1.0_num)
            ! Heisenburg uncertainty length for fast particles in each cell
            delta_x(ix,iy,iz) = h_planck / (2.0_num * fast_p(ix,iy,iz))
            ! <q1*q2> averaged over all fast/background pairs
            q1q2(ix,iy,iz) = q1q2(ix,iy,iz) &
                / MAX(pair_num_dens_tot(ix,iy,iz), 1.0_num)
          END DO
        END DO
      END DO

      ! Invert (debye_length)^(-2) array to return the debye length
      debye_length = 1.0_num / SQRT(debye_length)

      ! Impact parameter for a 2-rad scatter angle (using p*v ~= 2*KE)
      b0 = q1q2 / (8.0_num * pi * epsilon0 * fast_ekbar)

      ! Calculate interatomic distance for background particles
      r_min = 0.0_num
      IF (use_cold_correction) THEN
        ! Sum all background ni
        DO i_back = 1, back_count_all
          r_min = r_min + back_ni_all(i_back, :, :, :)
        END DO

        ! Calculate r_min
        DO iz = 1-ng,nz+ng
          DO iy = 1-ng,ny+ng
            DO ix = 1-ng,nx+ng
              n2_tot = MAX(r_min(ix,iy,iz), 1.0_num)
              r_min(ix,iy,iz) = &
                  (4.0_num * pi * n2_tot / 3.0_num)**(-1.0_num / 3.0_num)
            END DO
          END DO
        END DO
      END IF

      ! Coulomb logarithm
      DO iz = 1-ng,nz+ng
        DO iy = 1-ng,ny+ng
          DO ix = 1-ng,nx+ng
            bmax = MAX(debye_length(ix,iy,iz), r_min(ix,iy,iz))
            bmin = MAX(delta_x(ix,iy,iz), b0(ix,iy,iz))
            back_cou_log(ix,iy,iz) = MAX(1.0_num, LOG(bmax / bmin))
          END DO
        END DO
      END DO

      DEALLOCATE(debye_length, b0, delta_x, r_min)
      DEALLOCATE(temp, fast_num_dens, current_ekbar, fast_ekbar, fast_p)
      DEALLOCATE(q1q2, pair_num_dens, pair_num_dens_tot)
    END IF

    DEALLOCATE(num_dens)

  END SUBROUTINE update_background_variables



  SUBROUTINE calc_s12(i_back, s12)

    ! Subroutine used to calculate the isotropy parameter s12. This variable is
    ! related to the size of the "step" taken by the particle, weighted by the
    ! average single-collision angle. It is given by eq. (9) in Perez (2012). A
    ! low temperature correction is discussed in Section II.C.
    !
    ! This background model assumes species 2 is at rest, using approximations:
    !   1. gamma_2 = 1
    !   2. gamma_2* = gamma_c
    !   3. p2 = 0
    !   4. v_c is parallel to v_1, implied by Perez (2012) eq. 1.
    !   5. Not particle pairing, looping over species 1 only, ignore n12
    !      correction discussed in Section II.B. and used in Section II.C.
    !   6. Ignore mean-free-path of species 2, use lambda1 >= (4/3 pi n2)^(-1/3)
    !
    ! i_back is the index of the background particle. The remaining variables
    ! are module variables which this function has access to.

    INTEGER, INTENT(IN) :: i_back
    REAL(num), INTENT(OUT) :: s12
    REAL(num) :: mass_fac, p1_com_sq, p1_com_cu, s12_c, s12_d
    REAL(num) :: ilambda1, s12_cold

    mass_fac = (gamma1 * m1 + m2)

    ! Calculate variables in the centre-of-mass frame
    IF (.NOT. use_lab) THEN

      ! Centre of mass velocity (Perez (2012), eq. 1)
      v_com = p1 / mass_fac

      ! Centre of mass gamma factor
      gamma_com = 1.0_num / SQRT(1.0_num - (v_com / c)**2)

      ! Gamma factor of fast particle in the centre of mass frame (gamma_1*)
      gamma1_com = (1.0_num - v_com * v1 / c**2) * gamma_com * gamma1

      ! Fast particle momentum in centre of mass frame  (Perez (2012), eq. 2)
      p1_com = p1 + ((gamma_com - 1.0_num) * v1 - gamma_com * v_com) &
          * m1 * gamma1

    ELSE
      ! If p1 is similar to p1*, then just use lab frame variables
      gamma_com = 1.0_num
      gamma1_com = gamma1
      p1_com = p1
    END IF

    p1_com_sq = p1_com**2
    p1_com_cu = p1_com_sq * p1_com

    ! s12 (Perez (2012) eq. 9)
    s12_c = gamma_com * n2 * cou_log_12 / (gamma1 * p1_com_cu * mass_fac)
    s12_d = (gamma1_com * gamma_com * m1 * m2 * c**2 + p1_com_sq)**2
    s12 = s12_b(i_back) * s12_c * s12_d

    ! Cold particle correction (Perez (2012) eq. 20)
    ! MATLAB proto-typing on Shock Ignition plasma suggests this correction is
    ! important for electrons under 3 eV - for a keV plasma near 0.25 critical
    ! density, this can be ignored
    IF (use_cold_correction) THEN
      ilambda1 = (4.0_num * pi * n2 / 3.0_num)**(1.0_num/3.0_num)
      s12_cold = (m1 + m2) * v1 * ilambda1 * dt_coll / m2
      s12 = MIN(s12, s12_cold)
    END IF

  END SUBROUTINE calc_s12



  SUBROUTINE calc_cos_theta(s12, cos_theta)

    ! For a given isotropy parameter s12, this subroutine calculates the cosine
    ! of the cumulative scatter angle theta in the CoM frame. This is done
    ! using eq. 10 and eq. 11 in Perez (2012)

    REAL(num), INTENT(IN) :: s12
    REAL(num), INTENT(OUT) :: cos_theta
    REAL(num) :: a_perez, ia, s12_2, s12_3, s12_4, s12_5

    ! Sampling in Perez (2012), eq. 11
    IF (s12 < 0.1_num) THEN
      ! 5e-9 limit prevents |cos(theta)| > 1
      cos_theta = 1.0_num + s12 * LOG(MAX(random(), 5e-9_num))
      RETURN
    ELSE IF (s12 < 3.0_num) THEN
      s12_2 = s12**2
      s12_3 = s12_2 * s12
      s12_4 = s12_3 * s12
      s12_5 = s12_4 * s12
      ia = 0.0056958_num + 0.9560202_num * s12 - 0.508139_num * s12_2 &
          + 0.47913906_num * s12_3 - 0.12788975_num * s12_4 &
          + 0.02389567_num * s12_5
      a_perez = 1.0_num / ia
    ELSE IF (s12 < 6.0_num) THEN
      a_perez = 3.0_num * EXP(-s12)
      ia = 1.0_num / a_perez
    ELSE
      ! Corrected from Perez, which writes 2r + 1 (cos(theta) between 2 and 3)
      cos_theta = 2.0_num * random() - 1.0_num
      RETURN
    END IF

    ! Use A to get cos(theta), through Perez (2012), eq. 10
    cos_theta = ia * LOG(EXP(-a_perez) + 2.0_num * random() * SINH(a_perez))
    cos_theta = MAX(MIN(cos_theta, 1.0_num), -1.0_num)

  END SUBROUTINE calc_cos_theta



  SUBROUTINE rotate_p_vec(p_vec, cos_theta, phi, p_mag)

    ! Let the polar axis be defined as the direction of the initial momentum
    ! p_vec. This subroutine rotates p_vec by the polar angle theta (we read in
    ! cos(theta)), and the azimuthal angle phi, without changing the magnitude.
    !
    ! If we have already calculated the magnitude of the particle's momentum,
    ! p_mag, this can also be fed into the subroutine to speed up the
    ! calculation
    !
    ! Rotation scheme covered by Peplow (1999) - "Direction cosines and
    ! polarization vectors for monte carlo photon scattering"

    REAL(num), INTENT(INOUT) :: p_vec(3)
    REAL(num), INTENT(IN) :: cos_theta, phi
    REAL(num), OPTIONAL :: p_mag
    REAL(num) :: p, frac_p, pcos_theta, sin_theta, psin_theta
    REAL(num) :: ux, uy, uz, pfrac_uz, cos_phi, term_1, term_2

    ! Extract particle momentum
    IF (PRESENT(p_mag)) THEN
      p = p_mag
    ELSE
      p = SQRT(p_vec(1)**2 + p_vec(2)**2 + p_vec(3)**2)
    END IF
    frac_p = 1.0_num / p

    ! Precalculate repeated terms
    sin_theta = SQRT(1.0_num - cos_theta**2)
    pcos_theta = p * cos_theta
    psin_theta = p * sin_theta
    uz = p_vec(3) * frac_p

    IF (ABS(1.0_num - uz) < 1.0e-10_num) THEN
      ! Special case if the polar direction points along z
      p_vec(1) = psin_theta * COS(phi)
      p_vec(2) = psin_theta * SIN(phi)
      p_vec(3) = pcos_theta * SIGN(1.0_num, uz)
    ELSE
      ! Precalculate repeated terms
      ux = p_vec(1) * frac_p
      uy = p_vec(2) * frac_p
      pfrac_uz = p / SQRT(1.0_num - uz**2)
      cos_phi = COS(phi)
      term_1 = sin_theta * cos_phi * uz * pfrac_uz + pcos_theta
      term_2 = sin_theta * SIN(phi) * pfrac_uz

      p_vec(1) = ux * term_1 - uy * term_2
      p_vec(2) = uy * term_1 + ux * term_2
      p_vec(3) = uz * pcos_theta &
          + cos_phi * sin_theta * (uz**2 - 1.0_num) * pfrac_uz
    END IF

  END SUBROUTINE rotate_p_vec



  SUBROUTINE deallocate_background_collisions

    IF (use_background_collisions) THEN
      DEALLOCATE(back_id_all, id2array_all, p1_lab)
    END IF

  END SUBROUTINE deallocate_background_collisions

END MODULE background_collisions
