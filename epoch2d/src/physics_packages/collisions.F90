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

! Relativistic binary collision module
! written by M. G. Ramsay and H. Schmitz
! based on algorithm by Sentoku & Kemp [J Comput Phys, 227, 6846 (2008)]

MODULE collisions

  USE calc_df
#ifdef PREFETCH
  USE prefetch
#endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: particle_collisions, setup_collisions, collisional_ionisation
  PUBLIC :: deallocate_collisions
#ifdef COLLISIONS_TEST
  PUBLIC :: test_collisions
#endif

  ABSTRACT INTERFACE
    SUBROUTINE scatter_proto(current, impact, mass1, mass2, charge1, charge2, &
        weight1, weight2, idens, jdens, log_lambda, factor, np)

      IMPORT particle, num

      TYPE(particle), POINTER :: current, impact
      REAL(num), INTENT(IN) :: mass1, mass2
      REAL(num), INTENT(IN) :: charge1, charge2
      REAL(num), INTENT(IN) :: weight1, weight2
      REAL(num), INTENT(IN) :: idens, jdens, log_lambda
      REAL(num), INTENT(IN) :: factor, np
    END SUBROUTINE scatter_proto


    SUBROUTINE intra_collisions_proto(p_list, mass, charge, weight, &
        dens, log_lambda, user_factor)

      IMPORT particle_list, num

      TYPE(particle_list), INTENT(INOUT) :: p_list
      REAL(num), INTENT(IN) :: mass, charge, weight
      REAL(num), INTENT(IN) :: user_factor
      REAL(num), INTENT(IN) :: dens, log_lambda

    END SUBROUTINE intra_collisions_proto



    SUBROUTINE inter_collisions_proto(p_list1, p_list2, mass1, mass2, &
      charge1, charge2, weight1, weight2, &
      idens, jdens, log_lambda, user_factor )

      IMPORT particle_list, num

      TYPE(particle_list), INTENT(INOUT) :: p_list1
      TYPE(particle_list), INTENT(INOUT) :: p_list2

      REAL(num), INTENT(IN) :: mass1, charge1, weight1
      REAL(num), INTENT(IN) :: mass2, charge2, weight2

      REAL(num), INTENT(IN) :: idens, jdens
      REAL(num), INTENT(IN) :: log_lambda
      REAL(num), INTENT(IN) :: user_factor

    END SUBROUTINE inter_collisions_proto

  END INTERFACE

  PROCEDURE(scatter_proto), POINTER, SAVE :: scatter_fn => NULL()
  PROCEDURE(intra_collisions_proto), POINTER, SAVE :: intra_coll_fn => NULL()
  PROCEDURE(inter_collisions_proto), POINTER, SAVE :: inter_coll_fn => NULL()

  REAL(num), PARAMETER :: eps = EPSILON(1.0_num)
  REAL(num), PARAMETER :: one_m_2eps = 1.0_num - 2.0_num * eps
  REAL(num), PARAMETER :: one_p_2eps = 1.0_num + 2.0_num * eps

  REAL(num), PARAMETER :: e_rest = m0 * c**2
  REAL(num), PARAMETER :: e_rest_ev = e_rest / ev
  REAL(num), PARAMETER :: mrbeb_const = 2.0_num * pi * a0**2 * alpha**4
  REAL(num), PARAMETER :: cc = c**2

#ifndef PER_SPECIES_WEIGHT
  REAL(num), DIMENSION(3,0:2), PARAMETER :: a_bell = RESHAPE( &
      (/ 0.5250_num, 0.5300_num, 0.1300_num, &
         0.0000_num, 0.6000_num, 0.3880_num, &
         0.0000_num, 0.0000_num, 0.3500_num /) * 1e-13_num, (/3,3/) )

  REAL(num), DIMENSION(3,0:2,7), PARAMETER :: b_bell = RESHAPE( &
      (/-0.5100_num, -0.4100_num,  0.2500_num, &
         0.0000_num, -0.4000_num, -0.2000_num, &
         0.0000_num,  0.0000_num,  1.6000_num, &
         0.2000_num,  0.1500_num, -1.5000_num, &
         0.0000_num, -0.7100_num, -0.2356_num, &
         0.0000_num,  0.0000_num, -3.0000_num, &
         0.0500_num,  0.1500_num,  2.4000_num, &
         0.0000_num,  0.6550_num,  0.5355_num, &
         0.0000_num,  0.0000_num,  4.0000_num, &
        -0.0250_num, -0.2000_num,  3.2200_num, &
         0.0000_num,  0.4250_num,  3.1500_num, &
         0.0000_num,  0.0000_num,  2.0000_num, &
        -0.1000_num, -0.1500_num, -3.6670_num, &
         0.0000_num, -0.7500_num, -8.5000_num, &
         0.0000_num,  0.0000_num, -5.0000_num, &
         0.0000_num,  0.0000_num,  0.0000_num, &
         0.0000_num,  0.0000_num,  5.0500_num, &
         0.0000_num,  0.0000_num, -1.5000_num, &
         0.0000_num,  0.0000_num,  0.0000_num, &
         0.0000_num,  0.0000_num,  0.3700_num, &
         0.0000_num,  0.0000_num,  3.5000_num /) * 1e-13_num, (/3,3,7/) )

  REAL(num), DIMENSION(0:2), PARAMETER :: &
      l_bell = (/ 1.27_num, 0.542_num, 0.95_num /) * 1e-13_num
#endif

  REAL(num), DIMENSION(:,:), ALLOCATABLE :: meanx, meany, meanz, part_count

CONTAINS

  SUBROUTINE particle_collisions

    INTEGER :: ispecies, jspecies
    INTEGER(i8) :: ix, iy
    TYPE(particle_list), POINTER :: p_list1
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: idens, jdens
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: itemp, jtemp, log_lambda
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: iekbar
    REAL(num) :: user_factor, q1, q2, m1, m2, w1, w2
    LOGICAL :: collide_species

    ALLOCATE(idens(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(jdens(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(itemp(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(jtemp(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(log_lambda(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meanx(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meany(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meanz(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(part_count(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(iekbar(1-ng:nx+ng,1-ng:ny+ng))

    DO ispecies = 1, n_species
      ! Currently no support for photon collisions so just cycle round
      IF (species_list(ispecies)%species_type == c_species_id_photon) &
          CYCLE
      ! Currently no support for collisions involving chargeless particles
      IF (ABS(species_list(ispecies)%charge) <= c_tiny) &
          CYCLE

      collide_species = .FALSE.
      DO jspecies = ispecies, n_species
        user_factor = coll_pairs(ispecies, jspecies)
        IF (user_factor > 0) THEN
          collide_species = .TRUE.
          EXIT
        END IF
      END DO

      IF (.NOT.collide_species) CYCLE

      CALL calc_coll_number_density(idens, ispecies)

      m1 = species_list(ispecies)%mass
      q1 = species_list(ispecies)%charge
      w1 = species_list(ispecies)%weight

      DO iy = 1, ny
      DO ix = 1, nx
        p_list1 => species_list(ispecies)%secondary_list(ix,iy)
        CALL shuffle_particle_list_random(p_list1)
      END DO ! ix
      END DO ! iy

      DO jspecies = ispecies, n_species
        ! Currently no support for photon collisions so just cycle round
        IF (species_list(jspecies)%species_type == c_species_id_photon) &
            CYCLE
        ! Currently no support for collisions involving chargeless particles
        IF (ABS(species_list(jspecies)%charge) <= c_tiny) &
            CYCLE
        user_factor = coll_pairs(ispecies, jspecies)
        IF (user_factor <= 0) CYCLE

        IF (ispecies /= jspecies) THEN
          CALL calc_coll_number_density(jdens, jspecies)
        END IF

        m2 = species_list(jspecies)%mass
        q2 = species_list(jspecies)%charge
        w2 = species_list(jspecies)%weight

        IF (coulomb_log_auto) THEN
          CALL calc_coll_ekbar(iekbar, ispecies)
          CALL calc_coll_temperature_ev(itemp, ispecies)
          IF (ispecies == jspecies) THEN
            log_lambda = calc_coulomb_log(iekbar, itemp, idens, idens, &
                q1, q1, m1)
          ELSE
            CALL calc_coll_temperature_ev(jtemp, jspecies)
            log_lambda = calc_coulomb_log(iekbar, jtemp, idens, jdens, &
                q1, q2, m1)
          END IF
        ELSE
          log_lambda = coulomb_log
        END IF

        DO iy = 1, ny
        DO ix = 1, nx
          IF (ispecies == jspecies) THEN
            CALL intra_coll_fn( &
                species_list(ispecies)%secondary_list(ix,iy), &
                m1, q1, w1, idens(ix,iy), &
                log_lambda(ix,iy), user_factor)
          ELSE
            CALL inter_species_collisions( &
                species_list(ispecies)%secondary_list(ix,iy), &
                species_list(jspecies)%secondary_list(ix,iy), &
                m1, m2, q1, q2, w1, w2, idens(ix,iy), jdens(ix,iy), &
                log_lambda(ix,iy), user_factor)
          END IF
        END DO ! ix
        END DO ! iy
      END DO ! jspecies
    END DO ! ispecies

    DEALLOCATE(idens, jdens, itemp, jtemp, log_lambda)
    DEALLOCATE(meanx, meany, meanz, part_count)
    DEALLOCATE(iekbar)

  END SUBROUTINE particle_collisions



  SUBROUTINE collisional_ionisation

#ifndef PER_SPECIES_WEIGHT
    INTEGER :: ispecies, jspecies, ion_species, e_species, n1, n2, l
    INTEGER(i8) :: ix, iy
    TYPE(particle_list), POINTER :: p_list1
    TYPE(particle_list) :: ionising_e, ejected_e
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: idens, jdens, e_dens
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: itemp, jtemp, e_temp
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: log_lambda, e_log_lambda
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: iekbar, e_ekbar
    REAL(num) :: user_factor, e_user_factor, q1, q2, m1, m2, w1, w2
    REAL(num) :: q_e, m_e, w_e, q_full, ionisation_energy
    LOGICAL :: use_coulomb_log_auto_i, use_coulomb_log_auto

    DO iy = 1, ny
    DO ix = 1, nx
      DO ispecies = 1, n_species
        p_list1 => species_list(ispecies)%secondary_list(ix,iy)
        CALL shuffle_particle_list_random(p_list1)
      END DO
    END DO ! ix
    END DO ! iy

    ALLOCATE(idens(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(jdens(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(e_dens(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(itemp(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(jtemp(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(e_temp(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(log_lambda(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(e_log_lambda(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meanx(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meany(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(meanz(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(part_count(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(iekbar(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(e_ekbar(1-ng:nx+ng,1-ng:ny+ng))

    CALL create_empty_partlist(ionising_e)
    CALL create_empty_partlist(ejected_e)

    DO ispecies = 1, n_species
      ! Currently no support for photon collisions so just cycle round
      IF (species_list(ispecies)%species_type == c_species_id_photon) &
          CYCLE
      ! Currently no support for collisions involving chargeless particles
      ! unless ionisation occurs
      use_coulomb_log_auto_i = .TRUE.
      IF (ABS(species_list(ispecies)%charge) <= c_tiny) THEN
        IF (.NOT. species_list(ispecies)%ionise) CYCLE
        use_coulomb_log_auto_i = .FALSE.
      END IF
      CALL calc_coll_number_density(idens, ispecies)
      CALL calc_coll_temperature_ev(itemp, ispecies)
      CALL calc_coll_ekbar(iekbar, ispecies)

      m1 = species_list(ispecies)%mass
      q1 = species_list(ispecies)%charge
      w1 = species_list(ispecies)%weight

      IF (species_list(ispecies)%ionise) THEN
        e_species = species_list(ispecies)%release_species
        CALL calc_coll_number_density(e_dens, e_species)
        CALL calc_coll_temperature_ev(e_temp, e_species)
        m_e = species_list(e_species)%mass
        q_e = species_list(e_species)%charge
        w_e = species_list(e_species)%weight
        n1 = species_list(ispecies)%n
        l = species_list(ispecies)%l
        ionisation_energy = species_list(ispecies)%ionisation_energy / ev
        ion_species = species_list(ispecies)%ionise_to_species
        n2 = species_list(ispecies)%n
        DO WHILE(species_list(ion_species)%ionise)
          ion_species = species_list(ion_species)%ionise_to_species
        END DO
        q_full = species_list(ion_species)%charge
        ion_species = species_list(ispecies)%ionise_to_species
      END IF

      DO jspecies = ispecies, n_species
        ! Currently no support for photon collisions so just cycle round
        IF (species_list(jspecies)%species_type == c_species_id_photon) &
            CYCLE
        ! Currently no support for collisions involving chargeless particles
        ! unless ionisation occurs
        use_coulomb_log_auto = use_coulomb_log_auto_i
        IF (ABS(species_list(jspecies)%charge) <= c_tiny) THEN
          IF (.NOT. (species_list(ispecies)%electron &
              .AND. species_list(jspecies)%ionise)) CYCLE
          use_coulomb_log_auto = .FALSE.
        END IF
        user_factor = coll_pairs(ispecies, jspecies)
        IF (user_factor <= 0) CYCLE

        CALL calc_coll_number_density(jdens, jspecies)
        CALL calc_coll_temperature_ev(jtemp, jspecies)

        m2 = species_list(jspecies)%mass
        q2 = species_list(jspecies)%charge
        w2 = species_list(jspecies)%weight

        IF (species_list(ispecies)%electron &
            .AND. species_list(jspecies)%ionise) THEN
          e_species = species_list(jspecies)%release_species
          CALL calc_coll_number_density(e_dens, e_species)
          CALL calc_coll_temperature_ev(e_temp, e_species)
          m_e = species_list(e_species)%mass
          q_e = species_list(e_species)%charge
          w_e = species_list(e_species)%weight
          n1 = species_list(jspecies)%n
          l = species_list(jspecies)%l
          ionisation_energy = species_list(jspecies)%ionisation_energy / ev
          ion_species = species_list(jspecies)%ionise_to_species
          n2 = species_list(ion_species)%n
          DO WHILE(species_list(ion_species)%ionise)
            ion_species = species_list(ion_species)%ionise_to_species
          END DO
          q_full = species_list(ion_species)%charge
          ion_species = species_list(jspecies)%ionise_to_species
        END IF

        IF (coulomb_log_auto) THEN
          IF (use_coulomb_log_auto) THEN
            log_lambda = calc_coulomb_log(iekbar, jtemp, idens, jdens, &
                q1, q2, m1)
          ELSE
            log_lambda = 0
          END IF
          IF (species_list(ispecies)%electron &
              .AND. species_list(jspecies)%ionise) THEN
            e_log_lambda = calc_coulomb_log(iekbar, e_temp, idens, &
                e_dens, q1, q_e, m1)
            e_user_factor = coll_pairs(ispecies, ion_species)
          ELSE IF (species_list(ispecies)%ionise &
              .AND. species_list(jspecies)%electron) THEN
            CALL calc_coll_ekbar(e_ekbar, e_species)
            e_log_lambda = calc_coulomb_log(e_ekbar, jtemp, e_dens, &
                jdens, q_e, q2, m_e)
            e_user_factor = coll_pairs(ion_species, jspecies)
          END IF
        ELSE
          log_lambda = coulomb_log
          e_log_lambda = coulomb_log
        END IF

        IF (ispecies == jspecies) THEN
          DO iy = 1, ny
          DO ix = 1, nx
            CALL intra_coll_fn( &
                species_list(ispecies)%secondary_list(ix,iy), &
                m1, q1, w1, idens(ix,iy), log_lambda(ix,iy), user_factor)
          END DO ! ix
          END DO ! iy
        ELSE IF (species_list(ispecies)%ionise &
            .AND. species_list(jspecies)%electron) THEN
          DO iy = 1, ny
          DO ix = 1, nx
            ! Perform collisional ionisation before calculating scatter
            CALL preionise(species_list(jspecies)%secondary_list(ix,iy), &
                species_list(ispecies)%secondary_list(ix,iy), &
                species_list(ion_species)%secondary_list(ix,iy), &
                ionising_e, ejected_e, m2, m1, q2, q1, jdens(ix,iy), &
                q_full, ionisation_energy, n1, n2, l)
            ! Scatter ionising impact electrons off of ejected target electrons
            ! unless specified otherwise in input deck
            IF (e_user_factor > 0.0_num) THEN
              CALL inter_species_collisions(ejected_e, ionising_e, &
                  m_e, m2, q_e, q2, w_e, w2, &
                  e_dens(ix,iy), jdens(ix,iy), &
                  e_log_lambda(ix,iy), e_user_factor)
            END IF
            ! Scatter non-ionising impact electrons off of remaining unionised
            ! targets provided target has charge
            IF (ABS(q1) > c_tiny) THEN
              CALL inter_species_collisions( &
                  species_list(ispecies)%secondary_list(ix,iy), &
                  species_list(jspecies)%secondary_list(ix,iy), &
                  m1, m2, q1, q2, w1, w2, idens(ix,iy), jdens(ix,iy), &
                  log_lambda(ix,iy), user_factor)
            END IF
            ! Put ions and electrons into respective lists
            CALL append_partlist( &
                species_list(jspecies)%secondary_list(ix,iy), ionising_e)
            CALL append_partlist( &
                species_list(species_list(ispecies)%release_species)&
                %secondary_list(ix,iy), ejected_e)
          END DO ! ix
          END DO ! iy
        ELSE IF (species_list(ispecies)%electron &
            .AND. species_list(jspecies)%ionise) THEN
          DO iy = 1, ny
          DO ix = 1, nx
            ! Perform collisional ionisation before calculating scatter
            CALL preionise(species_list(ispecies)%secondary_list(ix,iy), &
                species_list(jspecies)%secondary_list(ix,iy), &
                species_list(ion_species)%secondary_list(ix,iy), &
                ionising_e, ejected_e, m1, m2, q1, q2, idens(ix,iy), &
                q_full, ionisation_energy, n1, n2, l)
            ! Scatter ionising impact electrons off of ejected target electrons
            ! unless specified otherwise in input deck
            IF (e_user_factor > 0.0_num) THEN
              CALL inter_species_collisions(ejected_e, ionising_e, &
                  m1, m_e, q1, q_e, w1, w_e, &
                  idens(ix,iy), e_dens(ix,iy), &
                  e_log_lambda(ix,iy), e_user_factor)
            END IF
            ! Scatter non-ionising impact electrons off of remaining unionised
            ! targets provided target has charge
            IF (ABS(q2) > c_tiny) THEN
              CALL inter_species_collisions( &
                  species_list(ispecies)%secondary_list(ix,iy), &
                  species_list(jspecies)%secondary_list(ix,iy), &
                  m1, m2, q1, q2, w1, w2, idens(ix,iy), jdens(ix,iy), &
                  log_lambda(ix,iy), user_factor)
            END IF
            ! Put electrons into respective lists
            CALL append_partlist( &
                species_list(ispecies)%secondary_list(ix,iy), ionising_e)
            CALL append_partlist( &
                species_list(species_list(jspecies)%release_species)&
                %secondary_list(ix,iy), ejected_e)
          END DO ! ix
          END DO ! iy
        ELSE
          DO iy = 1, ny
          DO ix = 1, nx
            CALL inter_species_collisions( &
                species_list(ispecies)%secondary_list(ix,iy), &
                species_list(jspecies)%secondary_list(ix,iy), &
                m1, m2, q1, q2, w1, w2, idens(ix,iy), jdens(ix,iy), &
                log_lambda(ix,iy), user_factor)
          END DO ! ix
          END DO ! iy
        END IF
      END DO ! jspecies
    END DO ! ispecies

    DEALLOCATE(idens, jdens, itemp, jtemp, log_lambda)
    DEALLOCATE(meanx, meany, meanz, part_count)
    DEALLOCATE(e_dens, e_temp, e_log_lambda)
    DEALLOCATE(iekbar, e_ekbar)
#endif

  END SUBROUTINE collisional_ionisation



#ifndef PER_SPECIES_WEIGHT
  SUBROUTINE preionise(electrons, ions, ionised, ionising_e, &
      ejected_e, e_mass, ion_mass, e_charge, ion_charge, e_dens, &
      full_ion_charge, ionisation_energy, n1, n2, l)

    TYPE(particle_list), INTENT(INOUT) :: electrons, ions, ionised
    TYPE(particle_list), INTENT(INOUT) :: ionising_e, ejected_e

    REAL(num), INTENT(IN) :: e_mass, e_charge, ion_mass, ion_charge, e_dens
    REAL(num), INTENT(IN) :: ionisation_energy, full_ion_charge

    INTEGER, INTENT(IN) :: n1, n2, l

    TYPE(particle), POINTER :: electron, ion, ejected_electron, next_ion, next_e

    LOGICAL, DIMENSION(:), ALLOCATABLE :: lost_ke, was_ionised

    REAL(num) :: factor, np, eiics, fion, gr, red_inc, red_ion,&
        i_p2, rot_y, rot_z, e_p2_i, e_e, beta_i, gamma_i, e_p_rot(3), e_ke_i, &
        e_v_i, mrbeb_c, t, tp, bp, bt2, bb2
    REAL(num) :: ionisation_energy_inv, red_ion_inv, prob_factor, denominator
    INTEGER(KIND=8) :: e_count, ion_count, pcount, i, k

    ! Inter-species collisions
    e_count = electrons%count
    ion_count = ions%count

    ! If there aren't enough particles to collide, then don't bother
    IF (e_count == 0 .OR. ion_count == 0) RETURN

    pcount = MAX(e_count, ion_count)
    factor = 0.0_num
    np = 0.0_num

    ALLOCATE(lost_ke(e_count), was_ionised(ion_count))
    lost_ke = .FALSE.
    was_ionised = .FALSE.

    ! temporarily join tail to the head of the lists to make them circular
    electrons%tail%next => electrons%head
    ions%tail%next => ions%head
    electron => electrons%head
    ion => ions%head

    DO k = 1, pcount
      np = np + electron%weight
      factor = factor + MIN(electron%weight, ion%weight)
      electron => electron%next
      ion => ion%next
    END DO

    electron => electrons%head
    ion => ions%head

    ionisation_energy_inv = 1.0_num / ionisation_energy
    red_ion = e_rest_ev * ionisation_energy_inv
    red_ion_inv = 1.0_num / red_ion
    ! Area must be multiplied by 1e-4 to convert from cm^2 to m^2
    prob_factor = -e_dens * np / factor * dt * 1e-4_num

    DO k = 1, pcount
      i_p2 = DOT_PRODUCT(ion%part_p, ion%part_p)
      ! Angles for rotation such that ion velocity |v| = v_x
      IF (i_p2 > 0.0_num) THEN
        IF (ABS(ion%part_p(1)) > c_tiny) THEN
          rot_y = DATAN(ion%part_p(3) / ion%part_p(1))
        ELSE
          rot_y = pi / 2.0_num
        END IF
        denominator = ion%part_p(1) * DCOS(rot_y) + ion%part_p(3) * DSIN(rot_y)
        IF (ABS(denominator) > c_tiny) THEN
          rot_z = DATAN(-ion%part_p(2) / denominator)
        ELSE
          rot_z = pi / 2.0_num
        END IF
        ! Rotate electron momentum into ion frame to simplify Lorentz transform
        e_p_rot = (/ (electron%part_p(1) * DCOS(rot_y) + electron%part_p(3) &
            * DSIN(rot_y)) * DCOS(rot_z) - electron%part_p(2) * DSIN(rot_z), &
            (electron%part_p(1) * DCOS(rot_y) + electron%part_p(3) &
            * DSIN(rot_y)) * DSIN(rot_z) + electron%part_p(2) * DCOS(rot_z), &
            electron%part_p(3) * DCOS(rot_y) - electron%part_p(1) &
            * DSIN(rot_y) /)
        rot_y = -rot_y
        rot_z = -rot_z
        e_p_rot = (/ (electron%part_p(1) * DCOS(rot_y) + electron%part_p(3) &
            * DSIN(rot_y)) * DCOS(rot_z) - electron%part_p(2) * DSIN(rot_z), &
            (electron%part_p(1) * DCOS(rot_y) + electron%part_p(3) &
            * DSIN(rot_y)) * DSIN(rot_z) + electron%part_p(2) * DCOS(rot_z), &
            electron%part_p(3) * DCOS(rot_y) - electron%part_p(1) &
            * DSIN(rot_y) /)
        ! Lorentz transform relativistic electron kinetic energy to ion frame
        gamma_i = SQRT(i_p2 / (ion_mass * c)**2 + 1.0_num)
        beta_i = SQRT(1.0_num - 1.0_num / gamma_i**2)
        e_e = c * SQRT(DOT_PRODUCT(e_p_rot, e_p_rot) + e_mass * e_rest)
        e_ke_i = (gamma_i * (e_e - beta_i * e_p_rot(1) * c) - e_rest) / ev
        ! Lorentz transform electron momentum
        e_p_rot(1) = gamma_i * (e_p_rot(1) - beta_i * e_e / c)
        ! Find electron velocity in ion frame
        e_p2_i = DOT_PRODUCT(e_p_rot, e_p_rot)
        e_v_i = SQRT(e_p2_i / (e_mass**2 + e_p2_i / cc))
      ELSE
        e_p2_i = DOT_PRODUCT(electron%part_p, electron%part_p)
        e_ke_i = c * (SQRT(e_p2_i + e_mass * e_rest) - e_mass * c) / ev
        e_v_i = SQRT(e_p2_i / (e_mass**2 + e_p2_i / cc))
      END IF
      ! Must enforce that electrons with insufficient kinetic energies cannot
      ! cause ionisation, as all cross sectional models used show massively
      ! increasing electron impact ionisation cross section as kinetic energy
      ! tends to zero
      IF (e_ke_i >= ion%weight / electron%weight * ionisation_energy &
          .AND. .NOT. was_ionised(MOD(k - 1, ion_count) + 1)) THEN
        ! Find cross section
        red_inc = e_ke_i * ionisation_energy_inv
        ! Use MBELL model for atomic number < 36
        IF (n1 < 4 .AND. l < 3) THEN
          ! Relativistic correction for high energy incident electrons
          gr = (1.0_num + 2.0_num * red_ion) / (red_inc + 2.0_num * red_ion) &
              * ((red_inc + red_ion) / (1.0_num + red_ion))**2 &
              * ((1.0_num + red_inc) * (red_inc + 2.0_num * red_ion) &
              * (1.0_num + red_ion)**2 &
              / (red_ion**2 * (1.0_num + 2.0_num * red_ion) &
              + red_inc * (red_inc + 2.0_num * red_ion) &
              * (1.0_num + red_ion)**2))**1.5_num

          ! Ionic correction for effect of charge of ion target upon incident
          ! electron
          fion = 1.0_num + 3.0_num * (ion_charge / (full_ion_charge &
              * red_inc))**l_bell(l)
          eiics = 0.0_num
          DO i = 1, 7
            eiics = eiics + b_bell(n1,l,i) * (1.0_num - 1.0_num/red_inc)**i
          END DO
          ! BELL cross section (cm^2)
          eiics = (a_bell(n1,l) * LOG(red_inc) + eiics) / (e_ke_i &
              * ionisation_energy)
          ! Relativisitic MBELL cross section (cm^2)
          eiics = fion * gr * eiics
        ! Use MRBEB model everywhere else
        ELSE
          t  = red_inc
          tp = e_ke_i / e_rest_ev
          bp = red_ion_inv
          bt2 = 1.0_num - 1.0_num / (1.0_num + tp)**2
          bb2 = 1.0_num - 1.0_num / (1.0_num + bp)**2

          mrbeb_c = hartree / ionisation_energy / 2.0_num &
              * (0.3_num * (ion_charge / q0 / n1)**2 &
              + 0.7_num * ((ion_charge / q0 + 1.0_num) / n2)**2)

          ! MRBEB cross section (cm^2)
          eiics = mrbeb_const / (bt2 + mrbeb_c * bb2) / bp &
              * (0.5_num * (LOG(bt2 / (1.0_num - bt2)) - bt2 &
              - LOG(2.0_num * bp)) * (1.0_num - 1.0_num / t**2) &
              + 1.0_num - 1.0_num / t &
              - LOG(t) / (t + 1.0_num) * (1.0_num + 2.0_num * tp) &
              / (1.0_num + 0.5_num * tp)**2 &
              + bp**2 / (1.0_num + 0.5_num * tp)**2 * (t - 1.0_num) / 2.0_num)
        END IF
        IF (random() < 1.0_num - EXP(prob_factor * eiics * e_v_i)) THEN
          ! Mark ionisation as occurring
          was_ionised(MOD(k - 1, ion_count) + 1) = .TRUE.
          lost_ke(MOD(k - 1, e_count) + 1) = .TRUE.
          IF (i_p2 > 0.0_num) THEN
            ! Reduce electron momentum by ionisation energy
            e_p_rot = SQRT(((ev / c * (e_ke_i - ion%weight / electron%weight &
                * ionisation_energy + e_rest_ev))**2 - e_mass * e_rest) &
                / e_p2_i) * e_p_rot
            ! Inverse Lorentz transform electron momentum back into simulation
            ! frame
            e_p_rot(1) = gamma_i * (e_p_rot(1) + beta_i &
                * SQRT(DOT_PRODUCT(e_p_rot, e_p_rot) + e_mass * e_rest))
            ! Undo frame rotation
            rot_y = -rot_y
            rot_z = -rot_z
            e_p_rot = (/ (e_p_rot(1) * DCOS(rot_y) + e_p_rot(3) &
                * DSIN(rot_y)) * DCOS(rot_z) - e_p_rot(2) * DSIN(rot_z), &
                (e_p_rot(1) * DCOS(rot_y) + e_p_rot(3) * DSIN(rot_y)) &
                * DSIN(rot_z) + e_p_rot(2) * DCOS(rot_z), e_p_rot(3) &
                * DCOS(rot_y) - e_p_rot(1) * DSIN(rot_y) /)
            ! If numerical error causes the electron to gain energy we catch it
            ! and apply the non-relativistic ionisation energy correction
            IF (DOT_PRODUCT(electron%part_p, electron%part_p) &
                < DOT_PRODUCT(e_p_rot, e_p_rot)) THEN
              electron%part_p = SQRT(((ev / c * (e_ke_i - ion%weight &
                / electron%weight * ionisation_energy + e_rest_ev))**2 &
                - e_mass * e_rest) / e_p2_i) * electron%part_p
            ELSE
              electron%part_p = e_p_rot
            END IF
          ELSE
            electron%part_p = SQRT(((ev / c * (e_ke_i - ion%weight &
                / electron%weight * ionisation_energy + e_rest_ev))**2 &
                - e_mass * e_rest) / e_p2_i) * electron%part_p
          END IF
        END IF
      END IF
      ion => ion%next
      electron => electron%next
    END DO

    ! restore the tail of the lists
    NULLIFY(electrons%tail%next)
    NULLIFY(ions%tail%next)

    electron => electrons%head
    ion => ions%head

    DO k = 1, MAX(e_count, ion_count)
      IF (k <= e_count) THEN
        next_e => electron%next
        IF (lost_ke(k)) THEN
          CALL remove_particle_from_partlist(electrons, electron)
          CALL add_particle_to_partlist(ionising_e, electron)
        END IF
      END IF
      IF (k <= ion_count) THEN
        next_ion => ion%next
        IF (was_ionised(k)) THEN
          CALL create_particle(ejected_electron)
          ejected_electron%weight = ion%weight
          ejected_electron%part_pos = ion%part_pos
          ! Ionise whilst conserving momentum
          ejected_electron%part_p = e_mass/ion_mass*ion%part_p
          ion%part_p = ion%part_p - ejected_electron%part_p
#ifdef PER_PARTICLE_CHARGE_MASS
          ejected_electron%charge = e_charge
          ejected_electron%mass = e_mass
          ion%charge = ion%charge - ejected_electron%charge
          ion%mass = ion%mass - ejected_electron%mass
#endif
#ifdef PARTICLE_DEBUG
          ejected_electron%processor = rank
          ejected_electron%processor_at_t0 = rank
#endif
          CALL add_particle_to_partlist(ejected_e, ejected_electron)
          CALL remove_particle_from_partlist(ions, ion)
          CALL add_particle_to_partlist(ionised, ion)
        END IF
      END IF
      electron => next_e
      ion => next_ion
    END DO

    DEALLOCATE(lost_ke, was_ionised)

  END SUBROUTINE preionise
#endif



  SUBROUTINE intra_collisions_sk(p_list, mass, charge, weight, &
      dens, log_lambda, user_factor)
    ! Perform collisions between particles of the same species.

    TYPE(particle_list), INTENT(INOUT) :: p_list
    REAL(num), INTENT(IN) :: mass, charge, weight
    REAL(num), INTENT(IN) :: user_factor
    REAL(num), INTENT(IN) :: dens, log_lambda
    TYPE(particle), POINTER :: current, impact
    REAL(num) :: factor, np
    INTEGER(i8) :: icount, k, pcount
    REAL(num), DIMENSION(3) :: p1, p2 ! Pre-collision momenta
    REAL(num), DIMENSION(3) :: p1_norm, p2_norm ! Normalised momenta
    REAL(num), DIMENSION(3) :: vc ! Velocity of COM frame wrt lab frame
    REAL(num), DIMENSION(3) :: p3, p4
    REAL(num), DIMENSION(3) :: p5, p6
    REAL(num), DIMENSION(3) :: v3, v4
    REAL(num), DIMENSION(3) :: vr, vcr
    REAL(num), DIMENSION(3) :: c1, c2, c3
    REAL(num) :: m1, m2, q1, q2, w1, w2 ! Masses and charges
    REAL(num) :: e1, e2 ! Pre-collision energies
    REAL(num) :: e3, e4, e5, e6
    REAL(num) :: gamma_rel, gamma_rel2, gamma_rel_m1, gamma_rel_r
    REAL(num) :: tvar ! Dummy variable for temporarily storing values
    REAL(num) :: vc_sq, vc_sq_cc, p1_vc, p2_vc, p3_mag
    REAL(num) :: delta, sin_theta, cos_theta, tan_theta_cm, tan_theta_cm2
    REAL(num) :: vrabs, denominator, wr
    REAL(num) :: nu, ran1, ran2

    factor = 0.0_num
    np = 0.0_num

    ! Intra-species collisions
    icount = p_list%count

    ! If there aren't enough particles to collide, then don't bother
    IF (icount <= 1) RETURN

    ! Number of collisions
    pcount = icount / 2 + MOD(icount, 2_i8)

#ifdef PER_SPECIES_WEIGHT
    np = icount * weight
    ! Factor of 2 due to intra species collisions
    ! See Section 4.1 of Nanbu
    factor = user_factor * pcount * weight * 2.0_num
#else
    ! temporarily join tail to the head of the list to make it circular
    p_list%tail%next => p_list%head

    current => p_list%head
    impact => current%next
    DO k = 1, pcount
      np = np + current%weight + impact%weight
      factor = factor + MIN(current%weight, impact%weight)
      current => impact%next
      impact => current%next
#ifdef PREFETCH
      CALL prefetch_particle(current)
      CALL prefetch_particle(impact)
#endif
    END DO
    factor = user_factor / factor / 2.0_num
#endif

    ! If possible, use per-species properties
    m1 = mass
    m2 = mass
    q1 = charge
    q2 = charge
    w1 = weight
    w2 = weight

    current => p_list%head
    impact => current%next
    DO k = 1, pcount
#ifdef PER_PARTICLE_CHARGE_MASS
      m1 = current%mass
      m2 = current%mass
      q1 = current%charge
      q2 = current%charge
#endif
#ifndef PER_SPECIES_WEIGHT
      w1 = current%weight
      w2 = current%weight
#endif
      CALL scatter_fn(current, impact, m1, m2, q1, q2, &
          w1, w2, dens, dens, log_lambda, factor, np)
      current => impact%next
      impact => current%next
#ifdef PREFETCH
      CALL prefetch_particle(current)
      CALL prefetch_particle(impact)
#endif
    END DO

    ! restore the tail of the list
    NULLIFY(p_list%tail%next)

  END SUBROUTINE intra_collisions_sk



  SUBROUTINE intra_collisions_np(p_list, mass, charge, weight, &
      dens, log_lambda, user_factor)
    ! Perform collisions between particles of the same species.

    TYPE(particle_list), INTENT(INOUT) :: p_list
    REAL(num), INTENT(IN) :: mass, charge, weight
    REAL(num), INTENT(IN) :: user_factor
    REAL(num), INTENT(IN) :: dens, log_lambda
    TYPE(particle), POINTER :: current, impact
    REAL(num) :: factor, np
    INTEGER(i8) :: icount, k, pcount
    REAL(num) :: ran1, ran2, s12, cosp, sinp, s_fac, v_rel
    REAL(num) :: sinp_cos, sinp_sin, s_prime
    REAL(num) :: a, a_inv, p_perp, p_tot, v_sq, gamma_rel_inv
    REAL(num) :: p_perp2, p_perp_inv
    REAL(num), DIMENSION(3) :: p1, p2, p3, p4, vc, v1, v2, p5, p6
    REAL(num), DIMENSION(3) :: p1_norm, p2_norm
    REAL(num), DIMENSION(3,3) :: mat
    REAL(num) :: p_mag, p_mag2, fac, gc, vc_sq, wr
    REAL(num) :: gm1, gm2, gm3, gm4, gm, gc_m1_vc
    REAL(num) :: m1, m2, q1, q2, w1, w2, e1, e5, e2, e6
    REAL(num), PARAMETER :: pi4_eps2_c4 = 4.0_num * pi * epsilon0**2 * c**4
    REAL(num), PARAMETER :: two_thirds = 2.0_num / 3.0_num
    REAL(num), PARAMETER :: pi_fac = &
                                (4.0_num * pi / 3.0_num)**(1.0_num / 3.0_num)

    factor = 0.0_num
    np = 0.0_num

    ! Intra-species collisions
    icount = p_list%count

    ! If there aren't enough particles to collide, then don't bother
    IF (icount <= 1) RETURN

    ! Number of collisions
    pcount = icount / 2 + MOD(icount, 2_i8)

#ifdef PER_SPECIES_WEIGHT
    np = icount * weight
    ! Factor of 2 due to intra species collisions
    ! See Section 4.1 of Nanbu
    factor = user_factor * pcount * weight * 2.0_num
#else
    ! temporarily join tail to the head of the list to make it circular
    p_list%tail%next => p_list%head

    current => p_list%head
    impact => current%next
    DO k = 1, pcount
      np = np + current%weight + impact%weight
      factor = factor + MIN(current%weight, impact%weight)
      current => impact%next
      impact => current%next
#ifdef PREFETCH
      CALL prefetch_particle(current)
      CALL prefetch_particle(impact)
#endif
    END DO
    factor = user_factor / factor / 2.0_num
#endif

    ! If possible, use per-species properties
    m1 = mass
    m2 = mass
    q1 = charge
    q2 = charge
    w1 = weight
    w2 = weight

    current => p_list%head
    impact => current%next
    DO k = 1, pcount
#ifdef PER_PARTICLE_CHARGE_MASS
      m1 = current%mass
      m2 = current%mass
      q1 = current%charge
      q2 = current%charge
#endif
#ifndef PER_SPECIES_WEIGHT
      w1 = current%weight
      w2 = current%weight
#endif
      CALL scatter_fn(current, impact, m1, m2, q1, q2, &
          w1, w2, dens, dens, log_lambda, factor, np)
      current => impact%next
      impact => current%next
#ifdef PREFETCH
      CALL prefetch_particle(current)
      CALL prefetch_particle(impact)
#endif
    END DO

    ! restore the tail of the list
    NULLIFY(p_list%tail%next)

  END SUBROUTINE intra_collisions_np



  SUBROUTINE inter_species_collisions(p_list1, p_list2, mass1, mass2, &
      charge1, charge2, weight1, weight2, &
      idens, jdens, log_lambda, user_factor )

    TYPE(particle_list), INTENT(INOUT) :: p_list1
    TYPE(particle_list), INTENT(INOUT) :: p_list2

    REAL(num), INTENT(IN) :: mass1, charge1, weight1
    REAL(num), INTENT(IN) :: mass2, charge2, weight2

    REAL(num), INTENT(IN) :: idens, jdens
    REAL(num), INTENT(IN) :: log_lambda
    REAL(num), INTENT(IN) :: user_factor

    TYPE(particle), POINTER :: current, impact

    REAL(num) :: factor, np
    INTEGER(i8) :: icount, jcount, pcount, k

    REAL(num) :: m1, m2, q1 q2, w1, w2

    factor = 0.0_num
    np = 0.0_num

    ! Inter-species collisions
    icount = p_list1%count
    jcount = p_list2%count
    pcount = MAX(icount, jcount)

    IF (icount > 0 .AND. jcount > 0) THEN
      ! temporarily join tail to the head of the lists to make them circular
      p_list1%tail%next => p_list1%head
      p_list2%tail%next => p_list2%head

#ifdef PER_SPECIES_WEIGHT
      np = icount * weight1
      factor = pcount * MIN(weight1, weight2)
#else
      current => p_list1%head
      impact => p_list2%head

      IF (icount >= jcount) THEN
        DO k = 1, icount
          np = np + current%weight
          current => current%next
        END DO
      ELSE
        DO k = 1, jcount
          np = np + impact%weight
          impact => impact%next
        END DO
      END IF

      current => p_list1%head
      impact => p_list2%head

      DO k = 1, pcount
        factor = factor + MIN(current%weight, impact%weight)
        current => current%next
        impact => impact%next
#ifdef PREFETCH
        CALL prefetch_particle(current)
        CALL prefetch_particle(impact)
#endif
      END DO
#endif

      ! If possible, use per-species properties
      m1 = mass1
      m2 = mass2
      q1 = charge1
      q2 = charge2
      w1 = weight1
      w2 = weight2

      current => p_list1%head
      impact => p_list2%head
      DO k = 1, pcount
#ifdef PER_PARTICLE_CHARGE_MASS
        m1 = current%mass
        m2 = impact%mass
        q1 = current%charge
        q2 = impact%charge
#endif
#ifndef PER_SPECIES_WEIGHT
        w1 = current%weight
        w2 = impact%weight
#endif
        CALL scatter_fn(current, impact, m1, m2, q1, q2, &
            w1, w2, idens, jdens, log_lambda, &
            user_factor / factor, np)
        current => current%next
        impact => impact%next
#ifdef PREFETCH
        CALL prefetch_particle(current)
        CALL prefetch_particle(impact)
#endif
      END DO

      ! restore the tail of the lists
      NULLIFY(p_list1%tail%next)
      NULLIFY(p_list2%tail%next)
    END IF

  END SUBROUTINE inter_species_collisions



  ! Binary collision scattering operator based jointly on:
  ! Perez et al. PHYSICS OF PLASMAS 19, 083104 (2012), and
  ! K. Nanbu and S. Yonemura, J. Comput. Phys. 145, 639 (1998)
  SUBROUTINE scatter_np(current, impact, m1, m2, q1, q2, &
      w1, w2, idens, jdens, log_lambda, factor, np)

    TYPE(particle), POINTER :: current, impact
    REAL(num), INTENT(IN) :: m1, m2
    REAL(num), INTENT(IN) :: q1, q2
    REAL(num), INTENT(IN) :: w1, w2
    REAL(num), INTENT(IN) :: idens, jdens, log_lambda
    REAL(num), INTENT(IN) :: factor, np
    REAL(num) :: ran1, ran2, s12, cosp, sinp, s_fac, v_rel
    REAL(num) :: sinp_cos, sinp_sin, s_prime
    REAL(num) :: a, a_inv, p_perp, p_tot, v_sq, gamma_rel_inv
    REAL(num) :: p_perp2, p_perp_inv
    REAL(num), DIMENSION(3) :: p1, p2, p3, p4, vc, v1, v2, p5, p6
    REAL(num), DIMENSION(3) :: p1_norm, p2_norm
    REAL(num), DIMENSION(3,3) :: mat
    REAL(num) :: p_mag, p_mag2, fac, gc, vc_sq, wr
    REAL(num) :: gm1, gm2, gm3, gm4, gm, gc_m1_vc
    REAL(num) :: e1, e5, e2, e6
    REAL(num), PARAMETER :: pi4_eps2_c4 = 4.0_num * pi * epsilon0**2 * c**4
    REAL(num), PARAMETER :: two_thirds = 2.0_num / 3.0_num
    REAL(num), PARAMETER :: pi_fac = &
                                (4.0_num * pi / 3.0_num)**(1.0_num / 3.0_num)

    p1 = current%part_p / c
    p2 = impact%part_p / c

    p1_norm = p1 / m0
    p2_norm = p2 / m0

    ! Two stationary particles can't collide, so don't try
    IF (DOT_PRODUCT(p1_norm, p1_norm) < eps &
        .AND. DOT_PRODUCT(p2_norm, p2_norm) < eps) RETURN

    ! Ditto for two particles with the same momentum
    vc = (p1_norm - p2_norm)
    IF (DOT_PRODUCT(vc, vc) < eps) RETURN

    p1_norm = p1 / m1
    gm1 = SQRT(DOT_PRODUCT(p1_norm, p1_norm) + 1.0_num) * m1

    p2_norm = p2 / m2
    gm2 = SQRT(DOT_PRODUCT(p2_norm, p2_norm) + 1.0_num) * m2

    gm = gm1 + gm2

    ! Pre-collision velocities
    v1 = p1 / gm1
    v2 = p2 / gm2

    ! Pre-collision energies
    e1 = c * SQRT(DOT_PRODUCT(p1, p1) + (m1 * c)**2)
    e2 = c * SQRT(DOT_PRODUCT(p2, p2) + (m2 * c)**2)

    ! Velocity of centre-of-momentum (COM) reference frame
    vc = (p1 + p2) / gm
    vc_sq = DOT_PRODUCT(vc, vc)

    gamma_rel_inv = SQRT(1.0_num - vc_sq)
    gc = 1.0_num / gamma_rel_inv

    gc_m1_vc = (gc - 1.0_num) / vc_sq

    p3 = p1 + (gc_m1_vc * DOT_PRODUCT(vc, v1) - gc) * gm1 * vc

    v_sq = DOT_PRODUCT(vc, v1)
    gm3 = (1.0_num - v_sq) * gc * gm1
    v_sq = DOT_PRODUCT(vc, v2)
    gm4 = (1.0_num - v_sq) * gc * gm2

    p_mag2 = DOT_PRODUCT(p3, p3)
    p_mag = SQRT(p_mag2)

    s_fac = idens * jdens * dt * factor * dx * dy
    fac = (q1 * q2)**2 * log_lambda * s_fac / (pi4_eps2_c4 * gm1 * gm2)
    s12 = fac * gc * p_mag * c / gm * (gm3 * gm4 / p_mag2 + 1.0_num)**2

    ! Cold plasma upper limit for s12
    v_rel = gm * p_mag * c / (gm3 * gm4 * gc)
    s_prime = pi_fac * s_fac * (m1 + m2) * v_rel &
        / MAX(m1 * idens**two_thirds, m2 * jdens**two_thirds)

    s12 = MIN(s12, s_prime)

    ran1 = random()
    ran2 = random() * 2.0_num * pi

    ! Inversion from Perez et al. PHYSICS OF PLASMAS 19, 083104 (2012)
    IF (s12 < 0.1_num) THEN
      cosp = 1.0_num + s12 * LOG(MAX(ran1, 5e-9_num))
    ELSE IF (s12 >= 0.1_num .AND. s12 < 3.0_num) THEN
      a_inv = 0.0056958_num + (0.9560202_num + (-0.508139_num &
           + (0.47913906_num + (-0.12788975_num + 0.02389567_num &
           * s12) * s12) * s12) * s12) * s12
      a = 1.0_num / a_inv
      cosp = a_inv * LOG(EXP(-a) + 2.0_num * ran1 * SINH(a))
    ELSE IF (s12 >= 3.0_num .AND. s12 < 6.0_num) THEN
      a = 3.0_num * EXP(-s12)
      cosp = LOG(EXP(-a) + 2.0_num * ran1 * SINH(a)) / a
    ELSE
      cosp = 2.0_num * ran1 - 1.0_num
    END IF

    sinp = SIN(ACOS(cosp))

    ! Calculate new momenta according to rotation by angle p
    p_perp2 = p3(1)**2 + p3(2)**2
    p_perp = SQRT(p_perp2)
    p_tot = SQRT(p_perp2 + p3(3)**2)
    p_perp_inv = 1.0_num / (p_perp + c_tiny)

    mat(1,1) =  p3(1) * p3(3) * p_perp_inv
    mat(1,2) = -p3(2) * p_tot * p_perp_inv
    mat(1,3) =  p3(1)
    mat(2,1) =  p3(2) * p3(3) * p_perp_inv
    mat(2,2) =  p3(1) * p_tot * p_perp_inv
    mat(2,3) =  p3(2)
    mat(3,1) = -p_perp
    mat(3,2) =  0.0_num
    mat(3,3) =  p3(3)

    sinp_cos = sinp * COS(ran2)
    sinp_sin = sinp * SIN(ran2)

    p3(1) = mat(1,1) * sinp_cos + mat(1,2) * sinp_sin + mat(1,3) * cosp
    p3(2) = mat(2,1) * sinp_cos + mat(2,2) * sinp_sin + mat(2,3) * cosp
    p3(3) = mat(3,1) * sinp_cos + mat(3,2) * sinp_sin + mat(3,3) * cosp

    p4 = -p3

    p5 = (p3 + (gc_m1_vc * DOT_PRODUCT(vc, p3) + gm3 * gc) * vc) * c
    p6 = (p4 + (gc_m1_vc * DOT_PRODUCT(vc, p4) + gm4 * gc) * vc) * c

    e5 = c * SQRT(DOT_PRODUCT(p5, p5) + (m1 * c)**2)
    e6 = c * SQRT(DOT_PRODUCT(p6, p6) + (m2 * c)**2)

    wr = w1 / w2
    IF (wr > one_p_2eps) THEN
      CALL weighted_particles_correction(w2 / w1, p1, p5, e1, e5, m1)
    ELSE IF (wr < one_m_2eps) THEN
      CALL weighted_particles_correction(w1 / w2, p2, p6, e2, e6, m2)
    END IF

    ! Update particle properties
    current%part_p = p5
    impact%part_p = p6

  END SUBROUTINE scatter_np




  ! Binary collision scattering operator based on that of
  ! Y. Sentoku and A. J. Kemp. [J Comput Phys, 227, 6846 (2008)]
  SUBROUTINE scatter_sk(current, impact, m1, m2, q1, q2, &
      w1, w2, idens, jdens, log_lambda, factor, np)

    ! Here the Coulomb collisions are performed by rotating the momentum
    ! vector of one of the particles in the centre of momentum reference
    ! frame. Then, as the collisions are assumed to be elastic, the momentum
    ! of the second particle is simply the opposite of the first.

    TYPE(particle), POINTER :: current, impact
    REAL(num), INTENT(IN) :: m1, m2
    REAL(num), INTENT(IN) :: q1, q2
    REAL(num), INTENT(IN) :: w1, w2
    REAL(num), INTENT(IN) :: idens, jdens, log_lambda
    REAL(num), INTENT(IN) :: factor, np
    REAL(num), DIMENSION(3) :: p1, p2 ! Pre-collision momenta
    REAL(num), DIMENSION(3) :: p1_norm, p2_norm ! Normalised momenta
    REAL(num), DIMENSION(3) :: vc ! Velocity of COM frame wrt lab frame
    REAL(num), DIMENSION(3) :: p3, p4
    REAL(num), DIMENSION(3) :: p5, p6
    REAL(num), DIMENSION(3) :: v3, v4
    REAL(num), DIMENSION(3) :: vr, vcr
    REAL(num), DIMENSION(3) :: c1, c2, c3
    REAL(num) :: e1, e2 ! Pre-collision energies
    REAL(num) :: e3, e4, e5, e6
    REAL(num) :: gamma_rel, gamma_rel2, gamma_rel_m1, gamma_rel_r
    REAL(num) :: tvar ! Dummy variable for temporarily storing values
    REAL(num) :: vc_sq, vc_sq_cc, p1_vc, p2_vc, p3_mag
    REAL(num) :: delta, sin_theta, cos_theta, tan_theta_cm, tan_theta_cm2
    REAL(num) :: vrabs, denominator, wr
    REAL(num) :: nu, ran1, ran2

    ! Copy all of the necessary particle data into variables with easier to
    ! read names
    p1 = current%part_p
    p2 = impact%part_p
    p1_norm = p1 / mc0
    p2_norm = p2 / mc0

    ! Two stationary particles can't collide, so don't try
    IF (DOT_PRODUCT(p1_norm, p1_norm) < eps &
        .AND. DOT_PRODUCT(p2_norm, p2_norm) < eps) RETURN

    ! Ditto for two particles with the same momentum
    vc = (p1_norm - p2_norm)
    IF (DOT_PRODUCT(vc, vc) < eps) RETURN

    ! Pre-collision energies
    e1 = c * SQRT(DOT_PRODUCT(p1, p1) + (m1 * c)**2)
    e2 = c * SQRT(DOT_PRODUCT(p2, p2) + (m2 * c)**2)

    ! Velocity of centre-of-momentum (COM) reference frame
    vc = (p1 + p2) * cc / (e1 + e2)
    vc_sq = DOT_PRODUCT(vc, vc)
    vc_sq_cc = vc_sq / cc

    gamma_rel2 = 1.0_num / (1.0_num - vc_sq_cc)
    gamma_rel = SQRT(gamma_rel2)
    gamma_rel_m1 = gamma_rel2 * vc_sq_cc / (gamma_rel + 1.0_num)

    ! Lorentz momentum transform to get into COM frame
    p1_vc = DOT_PRODUCT(p1, vc)
    p2_vc = DOT_PRODUCT(p2, vc)
    tvar = p1_vc * gamma_rel_m1 / (vc_sq + c_tiny)
    p3 = p1 + vc * (tvar - gamma_rel * e1 / cc)
    tvar = p2_vc * gamma_rel_m1 / (vc_sq + c_tiny)
    p4 = p2 + vc * (tvar - gamma_rel * e2 / cc)

    p3_mag = SQRT(DOT_PRODUCT(p3, p3))

    ! Lorentz energy transform
    e3 = gamma_rel * (e1 - p1_vc)
    e4 = gamma_rel * (e2 - p2_vc)
    ! Pre-collision velocities in COM frame
    v3 = p3 * cc / e3
    v4 = p4 * cc / e4

    ! Relative velocity
    tvar = 1.0_num - (DOT_PRODUCT(v3, v4) / cc)
    vr = (v3 - v4) / tvar
    vrabs = SQRT(DOT_PRODUCT(vr, vr))

    ! Collision frequency
    nu = coll_freq(vrabs, log_lambda, m1, m2, q1, q2, MIN(idens, jdens))
    nu = MIN(nu * factor * np * dt, 0.02_num)

    ! NOTE: nu is now the number of collisions per timestep, NOT collision
    ! frequency

    ! New coordinate system to simplify scattering.
    CALL new_coords(vr, c1, c2, c3)

    ! this is to ensure that 0 < ran1 < 1
    ! ran1=0 gives NaN in logarithm
    ! ran1=1 could give positive logarithm due to rounding errors
    ran1 = (1.0_num - 1.0e-10_num) * random() + 0.5e-10_num
    ran2 = 2.0_num * pi * random()

    ! Box Muller method for random from Gaussian distribution,
    ! mean 0, variance of nu
    ! Possible place for speed up by caching the second Box Muller number
    ! and using it later
    ! SQRT(-2.0_num * nu * LOG(ran1)) * COS(ran2)
    delta = SQRT(-2.0_num * nu * LOG(ran1)) * SIN(ran2)
    ran2 = 2.0_num * pi * random()

    ! angle theta in the One Particle at Rest frame
    sin_theta = 2.0_num * delta / (1.0_num + delta**2)
    cos_theta = (1.0_num - delta**2) / (1.0_num + delta**2)

    ! Transform angles from particle j's rest frame to COM frame
    ! Note azimuthal angle (ran2) is invariant under this transformation
    IF (m1 > m2) THEN
      vcr = v3
    ELSE
      vcr = v4
    END IF
    gamma_rel_r = 1.0_num / SQRT(1.0_num - (DOT_PRODUCT(vcr, vcr) / cc))

    denominator = gamma_rel_r * (cos_theta - SQRT(DOT_PRODUCT(vcr, vcr)) &
        / MAX(vrabs, c_tiny))
    IF (ABS(denominator) > SQRT(c_tiny)) THEN
      tan_theta_cm = sin_theta / denominator
      tan_theta_cm2 = tan_theta_cm**2
    ELSE
      tan_theta_cm = c_largest_number
      tan_theta_cm2 = c_largest_number
    END IF

    sin_theta = tan_theta_cm / SQRT(1.0_num + tan_theta_cm2)
    cos_theta = 1.0_num / SQRT(1.0_num + tan_theta_cm2)

    ! Post-collision momenta in COM frame
    p3 = p3_mag * (c1 * cos_theta + c2 * sin_theta * COS(ran2) &
        + c3 * sin_theta * SIN(ran2))
    p4 = -p3

    ! Lorentz momentum transform to get back to lab frame
    tvar = DOT_PRODUCT(p3, vc) * gamma_rel_m1 / vc_sq
    p5 = p3 + vc * (tvar + gamma_rel * e3 / cc)
    tvar = DOT_PRODUCT(p4, vc) * gamma_rel_m1 / vc_sq
    p6 = p4 + vc * (tvar + gamma_rel * e4 / cc)

    e5 = c * SQRT(DOT_PRODUCT(p5, p5) + (m1 * c)**2)
    e6 = c * SQRT(DOT_PRODUCT(p6, p6) + (m2 * c)**2)

    wr = w1 / w2
    IF (wr > one_p_2eps) THEN
      CALL weighted_particles_correction(w2 / w1, p1, p5, e1, e5, m1)
    ELSE IF (wr < one_m_2eps) THEN
      CALL weighted_particles_correction(w1 / w2, p2, p6, e2, e6, m2)
    END IF

    ! Update particle properties
    current%part_p = p5
    impact%part_p = p6

  END SUBROUTINE scatter_sk



  PURE FUNCTION coll_freq(vrabs, log_lambda, m1, m2, q1, q2, jdens)

    REAL(num), INTENT(IN) :: vrabs, log_lambda, m1, m2, q1, q2
    REAL(num), INTENT(IN) :: jdens
    REAL(num) :: mu, coll_freq, numerator, denominator
    REAL(num), PARAMETER :: fac = 4.0_num * pi * epsilon0**2

    mu = (m1 * m2) / (m1 + m2)

    IF (vrabs > 0.0_num) THEN
      numerator = (q1 * q2)**2 * jdens * log_lambda
      denominator = fac * mu**2 * vrabs**3
      IF (denominator <= 0.0_num &
          .OR. EXPONENT(numerator) - EXPONENT(denominator) &
          >= c_maxexponent) THEN
        coll_freq = 0.0_num
      ELSE
        coll_freq = numerator / denominator
      END IF
    ELSE
      coll_freq = 0.0_num
    END IF

  END FUNCTION coll_freq



  SUBROUTINE weighted_particles_correction(wtr, p, p_scat, en, en_scat, mass)

    ! This is the correction to the particle according to
    ! Sentoku and Kemp (2008) formulas 21 to 26.
    REAL(num), INTENT(INOUT) :: p_scat(3)
    REAL(num), INTENT(IN) :: p(3)
    REAL(num), INTENT(IN) :: wtr
    REAL(num), INTENT(IN) :: en, en_scat
    REAL(num), INTENT(IN) :: mass

    REAL(num) :: p_after(3)
    REAL(num) :: en_after
    REAL(num) :: gamma_en, gamma_p
    REAL(num) :: delta_p, p_mag, p_trans_mag
    REAL(num) :: c1(3), c2(3), c3(3)
    REAL(num) :: phi

    en_after = (1.0_num - wtr) * en + wtr * en_scat
    p_after  = (1.0_num - wtr) * p  + wtr * p_scat
    p_mag = SQRT(DOT_PRODUCT(p_after, p_after))
    gamma_en = en_after / (mass * cc)
    gamma_p = SQRT(1.0_num + (p_mag / mass / c)**2)

    ! This if-statement is just to take care of possible rounding errors
    ! gamma_p should always be smaller than gamma_en
    IF (gamma_p < gamma_en) THEN
      ! magnitude of the momentum correction
      delta_p = mass * c * SQRT(gamma_en**2 - gamma_p**2)
      p_trans_mag = SQRT(p_after(2)**2 + p_after(3)**2)

      CALL new_coords(p_after, c1, c2, c3)

      phi = 2.0_num * pi * random()

      ! Correcting for the loss in energy by adding a perpendicular
      ! momentum correction
      p_scat = p_after + delta_p * (c2 * COS(phi) + c3 * SIN(phi))
    END IF

  END SUBROUTINE weighted_particles_correction



  SUBROUTINE new_coords(vector, c1, c2, c3)

    ! New orthonormal basis vectors for calculating scattering angles:
    ! c1 = v / |v|
    ! c2 = (v x e1) / |v x e1|
    ! c3 = ((v x e1) x v) / |(v x e1) x v|
    ! where e1 is a unit vector parallel to x-axis. I.e., e1 = (1,0,0)
    ! New x-axis, c1, is parallel to the input vector
    ! New y-axis, c2, is perpendicular to c1 and the x-axis
    ! New z-axis ,c3, is perpendicular to c1 and c2
    REAL(num), DIMENSION(3), INTENT(IN) :: vector
    REAL(num), DIMENSION(3), INTENT(OUT) :: c1, c2, c3
    REAL(num) :: vtrans, vmag

    vmag = SQRT(DOT_PRODUCT(vector, vector))
    vtrans = SQRT(vector(2)**2 + vector(3)**2)

    IF (vtrans > c_tiny) THEN
      c1 = vector / vmag
      c2 = (/ 0.0_num, vector(3), -vector(2) /)
      c2 = c2 / vtrans
      c3 = (/ vtrans**2, &
          -(vector(1) * vector(2)), &
          -(vector(1) * vector(3)) /)
      c3 = c3 / (vmag * vtrans)
    ELSE
      c1 = (/ 1.0_num, 0.0_num, 0.0_num /)
      c2 = (/ 0.0_num, 1.0_num, 0.0_num /)
      c3 = (/ 0.0_num, 0.0_num, 1.0_num /)
    END IF

  END SUBROUTINE new_coords



  SUBROUTINE shuffle_particle_list_random(p_list)

    TYPE(particle_list), INTENT(INOUT) :: p_list
    TYPE(particle), POINTER :: particle1, particle2

    INTEGER :: i, idx, swap_idx
    INTEGER :: p_num

    p_num = INT(p_list%count)

    ! Nothing to be done
    IF (p_num <= 2) RETURN

    ! First make sure that the sorting array is large enough
    ! This should happen infrequently
    IF (p_num > coll_sort_array_size) THEN
      DEALLOCATE(coll_sort_array)

      ! make the sort array somewhat larger to avoid frequent deallocation
      ! and reallocation
      coll_sort_array_size = (11 * p_num) / 10 + 10
      ALLOCATE(coll_sort_array(coll_sort_array_size))
    END IF

    ! Copy all the particle pointers into the array and create random
    ! sort indices
    particle1 => p_list%head
    DO i = 1,p_num
      coll_sort_array(i)%particle => particle1
      particle1 => particle1%next
    END DO

    ! Shuffle particles using Durstenfeld's algorithm
    DO idx = p_num,2,-1
      swap_idx = FLOOR(idx * random()) + 1
      particle1 => coll_sort_array(idx)%particle
      coll_sort_array(idx)%particle => coll_sort_array(swap_idx)%particle
      coll_sort_array(swap_idx)%particle => particle1
    END DO

    ! Finally we have to copy back to the list
    ! Do head first
    particle1 => coll_sort_array(1)%particle
    p_list%head => particle1
    NULLIFY(particle1%prev)

    ! Then do all the particle links between head and tail.
    DO i = 2, p_num
      particle2 => particle1
      particle1 => coll_sort_array(i)%particle

      particle1%prev => particle2
      particle2%next => particle1
    END DO

    ! Finally set the tail (at the end of the loop, particle is pointing to
    ! the tail)
    p_list%tail => particle1
    NULLIFY(particle1%next)

  END SUBROUTINE shuffle_particle_list_random



  PURE FUNCTION calc_coulomb_log(ekbar1, temp2, dens1, dens2, q1, q2, m1)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: ekbar1, temp2
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: dens1, dens2
    REAL(num), INTENT(IN) :: q1, q2, m1
    REAL(num), DIMENSION(1-ng:nx+ng,1-ng:ny+ng) :: calc_coulomb_log
    REAL(num) :: b0, dB, bmin, bmax
    REAL(num) :: local_ekbar1, local_temp2, gamm
    INTEGER :: i, j

    calc_coulomb_log = 0.0_num
    DO j = 1-ng, ny+ng
    DO i = 1-ng, nx+ng
      local_ekbar1 = MAX(ekbar1(i,j), 100.0_num * q0)
      local_temp2 = MAX(temp2(i,j), 100.0_num)
      IF (dens1(i,j) <= 1.0_num .OR. dens2(i,j) <= 1.0_num) THEN
        calc_coulomb_log(i,j) = 1.0_num
      ELSE
        bmax = SQRT(epsilon0 * q0 * local_temp2 / (ABS(q2) * q0 * dens2(i,j)))
        b0 = ABS(q1 * q2) / (8.0_num * pi * epsilon0 * local_ekbar1)
        gamm = (local_ekbar1 / (m1 * cc)) + 1.0_num
        dB = 2.0_num * pi * h_bar / (SQRT(gamm**2 - 1.0_num) * m1 * c)
        bmin = MAX(b0, dB)
        calc_coulomb_log(i,j) = MAX(1.0_num, LOG(bmax / bmin))
      END IF
    END DO
    END DO

  END FUNCTION calc_coulomb_log



  SUBROUTINE calc_coll_number_density(data_array, ispecies)

    ! This routine calculates an approximate number density
    ! It assumes each particle only contributes to a single cell
    ! Returns zero if species_type == c_species_id_photon

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: ispecies
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: idx
    INTEGER :: jx, jy
    TYPE(particle), POINTER :: current
#include "particle_head.inc"

    data_array = 0.0_num

    IF (species_list(ispecies)%species_type == c_species_id_photon) RETURN

    idx = 1.0_num / dx / dy

    wdata = species_list(ispecies)%weight

    DO jy = 1, ny
    DO jx = 1, nx
      current => species_list(ispecies)%secondary_list(jx,jy)%head

      DO WHILE (ASSOCIATED(current))
#ifndef PER_SPECIES_WEIGHT
        wdata = current%weight
#endif

#include "particle_to_grid.inc"
        data_array(cell_x, cell_y) = &
            data_array(cell_x, cell_y) + wdata

        current => current%next
      END DO
    END DO ! jx
    END DO ! jy

    data_array = data_array * idx

  END SUBROUTINE calc_coll_number_density



  SUBROUTINE calc_coll_temperature_ev(sigma, ispecies)

    ! This subroutine calculates the grid-based temperature of a given
    ! particle species.
    ! It is almost identical to the calc_temperature subroutine in calc_df,
    ! except it uses the secondary_list rather than the attached_list.

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: sigma
    INTEGER, INTENT(IN) :: ispecies
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_pmx, part_pmy, part_pmz, sqrt_part_m
    ! The weight of a particle
    REAL(num) :: part_w
    REAL(num) :: gf
    INTEGER :: ix, iy
    INTEGER :: jx, jy
    TYPE(particle), POINTER :: current
#include "particle_head.inc"

    meanx = 0.0_num
    meany = 0.0_num
    meanz = 0.0_num
    part_count = 0.0_num
    sigma = 0.0_num

    sqrt_part_m  = SQRT(species_list(ispecies)%mass)
    part_w = species_list(ispecies)%weight

    DO jy = 1, ny
    DO jx = 1, nx
      current => species_list(ispecies)%secondary_list(jx,jy)%head

      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          gf = gx(ix) * gy(iy) * part_w
          meanx(cell_x+ix, cell_y+iy) = &
              meanx(cell_x+ix, cell_y+iy) + gf * part_pmx
          meany(cell_x+ix, cell_y+iy) = &
              meany(cell_x+ix, cell_y+iy) + gf * part_pmy
          meanz(cell_x+ix, cell_y+iy) = &
              meanz(cell_x+ix, cell_y+iy) + gf * part_pmz
          part_count(cell_x+ix, cell_y+iy) = &
              part_count(cell_x+ix, cell_y+iy) + gf
        END DO
        END DO
        current => current%next
      END DO
    END DO ! jx
    END DO ! jy

    CALL calc_boundary(meanx)
    CALL calc_boundary(meany)
    CALL calc_boundary(meanz)
    CALL calc_boundary(part_count)

    part_count = MAX(part_count, 1.e-6_num)

    meanx = meanx / part_count
    meany = meany / part_count
    meanz = meanz / part_count

    sqrt_part_m  = SQRT(species_list(ispecies)%mass)

    part_count = 0.0_num
    DO jy = 1, ny
    DO jx = 1, nx
      current => species_list(ispecies)%secondary_list(jx,jy)%head

      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          gf = gx(ix) * gy(iy)
          sigma(cell_x+ix, cell_y+iy) = &
              sigma(cell_x+ix, cell_y+iy) + gf &
              * ((part_pmx - meanx(cell_x+ix, cell_y+iy))**2 &
              + (part_pmy - meany(cell_x+ix, cell_y+iy))**2 &
              + (part_pmz - meanz(cell_x+ix, cell_y+iy))**2)
          part_count(cell_x+ix, cell_y+iy) = &
              part_count(cell_x+ix, cell_y+iy) + gf
        END DO
        END DO
        current => current%next
      END DO
    END DO ! jx
    END DO ! jy

    CALL calc_boundary(sigma)
    CALL calc_boundary(part_count)

    ! 3/2 kT = <p^2>/(2m)
    sigma = sigma / MAX(part_count, 1.e-6_num) / q0 / 3.0_num

  END SUBROUTINE calc_coll_temperature_ev



  SUBROUTINE calc_coll_ekbar(data_array, ispecies)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: ispecies
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc, part_u2
    ! The weight of a particle
    REAL(num) :: part_w
    REAL(num) :: gf
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, gamma_rel, gamma_rel_m1
    INTEGER :: ix, iy
    INTEGER :: jx, jy
    TYPE(particle), POINTER :: current
#include "particle_head.inc"

    data_array = 0.0_num
    part_count = 0.0_num
    part_mc  = 1.0_num
    part_w = 1.0_num

    part_mc = c * species_list(ispecies)%mass
    part_w = species_list(ispecies)%weight
    fac = part_mc * part_w * c

    DO jy = 1, ny
    DO jx = 1, nx
      current => species_list(ispecies)%secondary_list(jx,jy)%head

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
#endif
        fac = part_mc * part_w * c
#else
#ifndef PER_SPECIES_WEIGHT
        part_w = current%weight
        fac = part_mc * part_w * c
#endif
#endif

        IF (species_list(ispecies)%species_type /= c_species_id_photon) THEN
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc

          part_u2 = part_ux**2 + part_uy**2 + part_uz**2
          gamma_rel = SQRT(part_u2 + 1.0_num)
          gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)

          wdata = gamma_rel_m1 * fac
        ELSE
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
          wdata = current%particle_energy * part_w
#else
          wdata = 0.0_num
#endif
        END IF

#include "particle_to_grid.inc"

        DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          gf = gx(ix) * gy(iy)
          data_array(cell_x+ix, cell_y+iy) = &
              data_array(cell_x+ix, cell_y+iy) + gf * wdata
          part_count(cell_x+ix, cell_y+iy) = &
              part_count(cell_x+ix, cell_y+iy) + gf * part_w
        END DO
        END DO

        current => current%next
      END DO
    END DO ! jx
    END DO ! jy

    CALL calc_boundary(data_array)
    CALL calc_boundary(part_count)

    data_array = data_array / MAX(part_count, c_tiny)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    END DO

  END SUBROUTINE calc_coll_ekbar



  SUBROUTINE setup_collisions

    ALLOCATE(coll_pairs(n_species, n_species))
    coll_pairs = 1.0_num
    coll_sort_array_size = 1
    ALLOCATE(coll_sort_array(coll_sort_array_size))

    IF (use_nanbu) THEN
      scatter_fn => scatter_np
      intra_coll_fn => intra_collisions_np
    ELSE
      scatter_fn => scatter_sk
      intra_coll_fn => intra_collisions_sk
    END IF

  END SUBROUTINE setup_collisions



  SUBROUTINE deallocate_collisions

    INTEGER :: stat

    DEALLOCATE(coll_pairs, coll_sort_array, STAT=stat)

  END SUBROUTINE deallocate_collisions



#ifdef COLLISIONS_TEST
  SUBROUTINE test_collisions

    CALL setup_collisions
    CALL test_scatter
    CALL test_inter_species
    CALL test_intra_species
    CALL test_shuffle

  END SUBROUTINE test_collisions



  SUBROUTINE test_scatter

    TYPE(particle), POINTER :: part1, part2
    REAL(num) :: p1(3), p2(3)  ! momenta of the collision partners
    REAL(num) :: mass1, mass2     ! masses of the collision partners
    REAL(num) :: charge1, charge2 ! charges of the collision partners
    REAL(num) :: wt1, wt2         ! weights of the collision partners
    REAL(num) :: density          ! particle density
    REAL(num) :: t_factor         ! time step correction factor

    REAL(num) :: en1_before, en2_before
    REAL(num) :: en1_after, en2_after
    REAL(num) :: en_error
    REAL(num) :: en_sqr
    REAL(num) :: p_error(3)
    REAL(num) :: p_sqr

    INTEGER :: i, N

    CALL create_particle(part1)
    CALL create_particle(part2)

    N = 1000000

    dt = 1e-8
    density = 1e22
    t_factor = 1.0_num

    !==============================================================
    ! electron-electron collisions (equal weighting)
    !==============================================================

    WRITE(*,*) '==================================================='
    WRITE(*,*) '============   SCATTER TEST   ====================='
    WRITE(*,*) '==================================================='
    WRITE(*,*)
    WRITE(*,*) 'Testing electron-electron collisions (equal weighting)'

    mass1 = m0
    mass2 = m0
    charge1 = -q0
    charge2 = -q0
    wt1 = 1.0_num
    wt2 = 1.0_num

    en_sqr = 0.0_num
    p_sqr = 0.0_num

    en_error = 0.0_num
    p_error = 0.0_num

    DO i = 1, N
      p1(1) = 5 * mass1 * c * random()
      p1(2) = 5 * mass1 * c * random()
      p1(3) = 5 * mass1 * c * random()
      p2(1) = 5 * mass2 * c * random()
      p2(2) = 5 * mass2 * c * random()
      p2(3) = 5 * mass2 * c * random()

      part1%part_p = p1
      part2%part_p = p2
      part1%weight = wt1
      part2%weight = wt2

      en1_before = cc * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_before = cc * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error + p1 + p2
      p_sqr = p_sqr + SUM((p1 + p2)**2)
      en_error = en_error + en1_before + en2_before
      en_sqr = en_sqr + (en1_before + en2_before)**2

      CALL scatter(part1, part2, mass1, mass2, charge1, charge2, wt1, wt2, &
          density, density, 1.0e4_num, 1.0e4_num, 10.0_num, t_factor)

      p1 = part1%part_p
      p2 = part2%part_p

      en1_after = cc * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_after = cc * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error - p1 - p2
      en_error = en_error - en1_after - en2_after
    END DO

    WRITE(*,'(''  Errors after '',I10,'' iterations'')') N
    WRITE(*,'(''    p: '',ES15.8)') SQRT(SUM(p_error**2) / p_sqr)
    WRITE(*,'(''    E: '',ES15.8)') SQRT(en_error**2 / en_sqr)

    !==============================================================
    ! electron-ion collisions (equal weighting)
    !==============================================================

    WRITE(*,*) 'Testing electron-ion collisions (equal weighting)'

    mass1 = m0
    mass2 = 1836.2_num * 105.0_num * m0
    charge1 = -q0
    charge2 = -22.0_num * q0
    wt1 = 1.0_num
    wt2 = 1.0_num

    en_sqr = 0.0_num
    p_sqr = 0.0_num

    en_error = 0.0_num
    p_error = 0.0_num

    DO i = 1, N
      p1(1) = 5 * mass1 * c * random()
      p1(2) = 5 * mass1 * c * random()
      p1(3) = 5 * mass1 * c * random()
      p2(1) = 5 * mass2 * c * random()
      p2(2) = 5 * mass2 * c * random()
      p2(3) = 5 * mass2 * c * random()

      part1%part_p = p1
      part2%part_p = p2
      part1%weight = wt1
      part2%weight = wt2

      en1_before = cc * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_before = cc * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error + p1 + p2
      p_sqr = p_sqr + SUM((p1 + p2)**2)
      en_error = en_error + en1_before + en2_before
      en_sqr = en_sqr + (en1_before + en2_before)**2

      ! alternate particle1 and particle2
      IF (MOD(i, 2) == 0) THEN
        CALL scatter(part1, part2, mass1, mass2, charge1, charge2, wt1, wt2, &
            density, density, 1.0e4_num, 1.0e4_num, 5.0_num, t_factor)
      ELSE
        CALL scatter(part2, part1, mass2, mass1, charge2, charge1, wt2, wt1, &
            density, density, 1.0e4_num, 1.0e4_num, 5.0_num, t_factor)
      END IF

      p1 = part1%part_p
      p2 = part2%part_p

      en1_after = cc * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_after = cc * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error - p1 - p2
      en_error = en_error - en1_after - en2_after
    END DO

    WRITE(*,'(''  Errors after '',I10,'' iterations'')') N
    WRITE(*,'(''    p: '',ES15.8)') SQRT(SUM(p_error**2) / p_sqr)
    WRITE(*,'(''    E: '',ES15.8)') SQRT(en_error**2 / en_sqr)

    !==============================================================
    ! electron-electron collisions (random weighting)
    !==============================================================

    WRITE(*,*) 'Testing electron-electron collisions (random weighting)'

    mass1 = m0
    mass2 = m0
    charge1 = -q0
    charge2 = -q0

    en_sqr = 0.0_num
    p_sqr = 0.0_num

    en_error = 0.0_num
    p_error = 0.0_num

    DO i = 1, N
      p1(1) = 5 * mass1 * c * random()
      p1(2) = 5 * mass1 * c * random()
      p1(3) = 5 * mass1 * c * random()
      p2(1) = 5 * mass2 * c * random()
      p2(2) = 5 * mass2 * c * random()
      p2(3) = 5 * mass2 * c * random()
      wt1 = random() + 1e-10_num
      wt2 = random() + 1e-10_num

      part1%part_p = p1
      part2%part_p = p2
      part1%weight = wt1
      part2%weight = wt2

      en1_before = cc * wt1 * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_before = cc * wt2 * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error + wt1 * p1 + wt2 * p2
      p_sqr = p_sqr + SUM((wt1 * p1 + wt2 * p2)**2)
      en_error = en_error + en1_before + en2_before
      en_sqr = en_sqr + (en1_before + en2_before)**2

      CALL scatter(part1, part2, mass1, mass2, charge1, charge2, wt1, wt2, &
            density, density, 1.0e4_num, 1.0e4_num, 5.0_num, t_factor)

      p1 = part1%part_p
      p2 = part2%part_p

      en1_after = cc * wt1 * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_after = cc * wt2 * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error - wt1 * p1 - wt2 * p2
      en_error = en_error - en1_after - en2_after
    END DO

    WRITE(*,'(''  Errors after '',I10,'' iterations'')') N
    WRITE(*,'(''    p: '',ES15.8)') SQRT(SUM(p_error**2) / p_sqr)
    WRITE(*,'(''    E: '',ES15.8)') SQRT(en_error**2 / en_sqr)

    !==============================================================
    ! electron-ion collisions (random weighting)
    !==============================================================

    WRITE(*,*) 'Testing electron-ion collisions (random weighting)'

    mass1 = m0
    mass2 = 1836.2_num * 105.0_num * m0
    charge1 = -q0
    charge2 = -22.0_num * q0

    en_sqr = 0.0_num
    p_sqr = 0.0_num

    en_error = 0.0_num
    p_error = 0.0_num

    DO i = 1, N
      p1(1) = 5 * mass1 * c * random()
      p1(2) = 5 * mass1 * c * random()
      p1(3) = 5 * mass1 * c * random()
      p2(1) = 5 * mass2 * c * random()
      p2(2) = 5 * mass2 * c * random()
      p2(3) = 5 * mass2 * c * random()
      wt1 = random() + 1e-10_num
      wt2 = random() + 1e-10_num

      part1%part_p = p1
      part2%part_p = p2
      part1%weight = wt1
      part2%weight = wt2

      en1_before = cc * wt1 * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_before = cc * wt2 * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error + wt1 * p1 + wt2 * p2
      p_sqr = p_sqr + SUM((wt1 * p1 + wt2 * p2)**2)
      en_error = en_error + en1_before + en2_before
      en_sqr = en_sqr + (en1_before + en2_before)**2

      ! alternate particle1 and particle2
      IF (MOD(i, 2) == 0) THEN
        CALL scatter(part1, part2, mass1, mass2, charge1, charge2, wt1, wt2, &
            density, density, 1.0e4_num, 1.0e4_num, 5.0_num, t_factor)
      ELSE
        CALL scatter(part2, part1, mass2, mass1, charge2, charge1, wt2, wt1, &
            density, density, 1.0e4_num, 1.0e4_num, 5.0_num, t_factor)
      END IF

      p1 = part1%part_p
      p2 = part2%part_p

      en1_after = cc * wt1 * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_after = cc * wt2 * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error - wt1 * p1 - wt2 * p2
      en_error = en_error - en1_after - en2_after
    END DO

    WRITE(*,'(''  Errors after '',I10,'' iterations'')') N
    WRITE(*,'(''    p: '',ES15.8)') SQRT(SUM(p_error**2) / p_sqr)
    WRITE(*,'(''    E: '',ES15.8)') SQRT(en_error**2 / en_sqr)

  END SUBROUTINE test_scatter



  SUBROUTINE test_inter_species

    TYPE(particle_list) :: partlist1
    TYPE(particle_list) :: partlist2
    TYPE(particle), POINTER :: part
    REAL(num) :: mass1, mass2     ! masses of the collision partners
    REAL(num) :: charge1, charge2 ! charges of the collision partners
    REAL(num) :: plist1_length, plist2_length

    INTEGER(i8) :: N, max_num, i, j
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: histo1, histo2
    INTEGER :: histo1max, histo2max
    INTEGER(i8) :: cnt1, cnt2
    INTEGER(i8) :: a, b
    INTEGER :: error

    dt = 1.0e-8_num

    N = 10000
    max_num = 200

    plist1_length = 0.0_num
    plist2_length = 0.0_num

    ALLOCATE(histo1(0:2*max_num))
    ALLOCATE(histo2(0:2*max_num))

    WRITE(*,*) '==================================================='
    WRITE(*,*) '======   INTER SPECIES COLLISION TEST   ==========='
    WRITE(*,*) '==================================================='
    WRITE(*,*)
    WRITE(*,*) 'Testing electron-electron collisions (random weighting)'

    mass1 = m0
    mass2 = m0
    charge1 = -q0
    charge2 = -q0

    DO i = 1, N
      ! allocate particles
      CALL create_allocated_partlist(partlist1, INT(max_num*random(), KIND=i8))
      CALL create_allocated_partlist(partlist2, INT(max_num*random(), KIND=i8))

      ! fill particle values

      part => partlist1%head
      DO WHILE (ASSOCIATED(part))
        part%part_p(1) = 5 * mass1 * c * random()
        part%part_p(2) = 5 * mass1 * c * random()
        part%part_p(3) = 5 * mass1 * c * random()
        part%coll_count = 0

        part%weight = random() + 1e-10_num

        part => part%next
      END DO

      part => partlist2%head
      DO WHILE (ASSOCIATED(part))
        part%part_p(1) = 5 * mass2 * c * random()
        part%part_p(2) = 5 * mass2 * c * random()
        part%part_p(3) = 5 * mass2 * c * random()
        part%coll_count = 0

        part%weight = random() + 1e-10_num

        part => part%next
      END DO

      ! call scattering routine
      CALL inter_species_collisions(partlist1, partlist2, &
          mass1, mass2, charge1, charge2, 1.0_num, 1.0_num, &
          1.0e15_num, 1.0e15_num, 1.0e4_num, 1.0e4_num, 5.0_num, 1.0_num)

      ! diagnostics
      histo1 = 0
      histo2 = 0
      histo1max = 0
      histo2max = 0

!      WRITE(*,'(''  Particle List 1: '',I10)') partlist1%count
      cnt1 = partlist1%count
      cnt2 = partlist2%count
      plist1_length = plist1_length + cnt1
      plist2_length = plist2_length + cnt2

      ! creating histograms

      part => partlist1%head
      DO j = 1, partlist1%count
!        WRITE(*,'(''    '',I4,'': '',I10)') j, part%coll_count
        IF (part%coll_count > histo1max) histo1max = part%coll_count
        histo1(part%coll_count) = histo1(part%coll_count) + 1
        part => part%next
      END DO

!      WRITE(*,'(''  Particle List 2: '',I10)') partlist2%count

      part => partlist2%head
      DO WHILE (ASSOCIATED(part))
!        WRITE(*,'(''    '',I4,'': '',I10)') j, part%coll_count
        IF (part%coll_count > histo2max) histo2max = part%coll_count
        histo2(part%coll_count) = histo2(part%coll_count) + 1
        part => part%next
      END DO

! only need to output if something is wrong
!      WRITE(*,*) '  Histogram 1'
!      DO j = 1, histo1max
!        WRITE(*,'(''    '',I4,'': '',I10)') j, histo1(j)
!      END DO
!
!      WRITE(*,*) '  Histogram 2'
!      DO j = 1, histo2max
!        WRITE(*,'(''    '',I4,'': '',I10)') j, histo2(j)
!      END DO

      ! performing check on histograms

      ! first, subtract expected values from the array
      IF ((cnt1 > 0) .AND. (cnt2 > 0)) THEN
        IF (cnt1 > cnt2) THEN
          histo1(2) = histo1(2) - cnt1
          a = cnt1 / cnt2
          b = cnt1 - cnt2 * a
          histo2(2*a) = histo2(2*a) - (cnt2 - b)
          histo2(2*a+2) = histo2(2*a+2) - b
        ELSE IF (cnt1 < cnt2) THEN
          histo2(2) = histo2(2) - cnt2
          a = cnt2 / cnt1
          b = cnt2 - cnt1 * a
          histo1(2*a) = histo1(2*a) - (cnt1 - b)
          histo1(2*a+2) = histo1(2*a+2) - b
        ELSE
          histo1(2) = histo1(2) - cnt1
          histo2(2) = histo2(2) - cnt1
        END IF
      END IF

      ! now check that both arrays are zero
      error = 0
      DO j = 1, histo1max
        IF (histo1(j) /= 0) error = error + 1
      END DO
      DO j = 1, histo2max
        IF (histo2(j) /= 0) error = error + 1
      END DO

      IF (error > 0) THEN
        WRITE(*,*) '  Error in inter species collisions'
        WRITE(*,'(''    Iteration:   '',I10)') i
        WRITE(*,'(''    List counts: '',I10,'', '',I10)') cnt1, cnt2
        STOP
      END IF

      ! free particle list
      CALL destroy_partlist(partlist1)
      CALL destroy_partlist(partlist2)
    END DO

    WRITE(*,*) '  SUCCESS!'
    WRITE(*,'(''    Number of iterations:   '',I10)') N
    WRITE(*,'(''    Average list lengths:   '',F10.5,'', '',F10.5)') &
        plist1_length / N, plist2_length / N

    DEALLOCATE(histo1)
    DEALLOCATE(histo2)

  END SUBROUTINE test_inter_species



  SUBROUTINE test_intra_species

    TYPE(particle_list) :: partlist
    TYPE(particle), POINTER :: part
    REAL(num) :: mass, charge ! mass and charg of the collision partners

    INTEGER(i8) :: N, max_num, i, j
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: histo
    INTEGER :: histo_max
    INTEGER :: error
    REAL(num) :: plist_length

    N = 10000
    max_num = 200

    ALLOCATE(histo(0:2*max_num))

    WRITE(*,*) '==================================================='
    WRITE(*,*) '======   INTRA SPECIES COLLISION TEST   ==========='
    WRITE(*,*) '==================================================='
    WRITE(*,*)
    WRITE(*,*) 'Testing collisions (random weighting)'

    mass = m0
    charge = -q0
    plist_length = 0

    DO i = 1, N
      ! allocate particles
      CALL create_allocated_partlist(partlist, INT(max_num*random(), KIND=i8))

      ! fill particle values

      part => partlist%head
      DO WHILE (ASSOCIATED(part))
        part%part_p(1) = 5 * mass * c * random()
        part%part_p(2) = 5 * mass * c * random()
        part%part_p(3) = 5 * mass * c * random()
        part%coll_count = 0

        part%weight = random() + 1e-10_num

        part => part%next
      END DO

      ! call scattering routine
      CALL intra_coll_fn(partlist, mass, charge, &
          1.0_num, 1.0e15_num, 1.0e4_num, 5.0_num, 1.0_num)

      ! diagnostics
      histo = 0
      histo_max = 0
      plist_length = plist_length + partlist%count

!      WRITE(*,'(''  Particle List: '',I10)') partlist%count

      part => partlist%head
      DO j = 1, partlist%count
        ! WRITE(*,'(''    '',I4,'': '',I10)') j, part%coll_count
        IF (part%coll_count > histo_max) histo_max = part%coll_count
        histo(part%coll_count) = histo(part%coll_count) + 1
        part => part%next
      END DO

! only need to output if something is wrong
!      WRITE(*,*) '  Histogram'
!      DO j = 1, histo_max
!        WRITE(*,'(''    '',I4,'': '',I10)') j, histo(j)
!      END DO

      ! performing check on histograms

      ! first, subtract expected value from the array
      histo(2) = histo(2) - partlist%count

      ! now check that array is zero
      error = 0
      DO j = 1, histo_max
        IF (histo(j) /= 0) error = error + 1
      END DO

      IF (error > 0) THEN
        WRITE(*,*) '  Error in intra species collisions'
        WRITE(*,'(''    Iteration:   '',I10)') i
        WRITE(*,'(''    List count: '',I10)') partlist%count
        STOP
      END IF

      ! free particle list
      CALL destroy_partlist(partlist)
    END DO

    WRITE(*,*) '  SUCCESS!'
    WRITE(*,'(''    Number of iterations:   '',I10)') N
    WRITE(*,'(''    Average list length:   '',F10.5)') plist_length / N

    DEALLOCATE(histo)

  END SUBROUTINE test_intra_species



  SUBROUTINE test_shuffle

    TYPE(particle_list) :: partlist
    TYPE(particle), POINTER :: part
    INTEGER :: N, max_num, min_num, k
    INTEGER(i8) :: i, j, plist_length
    INTEGER :: iterations
    REAL(num), DIMENSION(:), ALLOCATABLE :: histo
    INTEGER, DIMENSION(:), ALLOCATABLE :: minp, maxp
    REAL(num), DIMENSION(:), ALLOCATABLE :: std_dev

    N = 10000
    max_num = 200
    min_num = 1
    iterations = 10

    ALLOCATE(histo(0:max_num))
    ALLOCATE(minp(0:max_num))
    ALLOCATE(maxp(0:max_num))
    ALLOCATE(std_dev(0:max_num))

    WRITE(*,*) '==================================================='
    WRITE(*,*) '============   RANDOM SHUFFLE TEST   =============='
    WRITE(*,*) '==================================================='
    WRITE(*,*)

    DO k = 1, iterations
      plist_length = INT((max_num - min_num) * random() + min_num, i8)
      histo = 0.0_num
      minp = max_num
      maxp = 0
      std_dev = 0.0_num

      WRITE(*,'(''  List length: '',I10)') plist_length

      DO i = 1, N
        ! allocate particles
        CALL create_allocated_partlist(partlist, plist_length)

        ! fill particle values
        !   note: we only need particle position which is stored in coll_count
        part => partlist%head
        DO j = 1, plist_length
          IF (.NOT. ASSOCIATED(part)) WRITE(*,*) '    !!not associated!!'
          part%coll_count = INT(j)
          part => part%next
        END DO

        ! now shuffle
        CALL shuffle_particle_list_random(partlist)

        ! perform statistics
        part => partlist%head
        DO j = 1, plist_length
          histo(j) = histo(j) + part%coll_count
          IF (minp(j) > part%coll_count) minp(j) = part%coll_count
          IF (maxp(j) < part%coll_count) maxp(j) = part%coll_count
          std_dev(j) = std_dev(j) + part%coll_count**2
          part => part%next
        END DO

        CALL destroy_partlist(partlist)
      END DO

      WRITE(*,'(''   Statistics ('',I10,'' runs)'')') N
      WRITE(*,*) '    avg        std_dev       min        max'
      DO i = 1, plist_length
        WRITE(*,'(''    '',F10.5,''    '',F10.5,''    '',I10,''    '',I10)') &
            histo(i) / N, SQRT(std_dev(i) / N - (histo(i) / N)**2), &
            minp(i), maxp(i)
      END DO

    END DO

    DEALLOCATE(histo)

  END SUBROUTINE test_shuffle

#endif

END MODULE collisions
