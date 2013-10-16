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
#ifdef COLLISIONS_TEST
  PUBLIC :: test_collisions
#endif

  REAL(num) :: collision_count, large_angle_collision
  REAL(num), PARAMETER :: eps = EPSILON(1.0_num)

  REAL(num) :: nu_avg
  REAL(num), PARAMETER :: e_rest = m0 * c**2
  REAL(num), PARAMETER :: e_rest_ev = e_rest / ev
  REAL(num), PARAMETER :: mrbeb_const = 2.0_num * pi * a0**2 * alpha**4
  INTEGER :: nu_count

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

  REAL(num), DIMENSION(:), ALLOCATABLE :: meanx, meany, meanz, part_count

CONTAINS

  SUBROUTINE particle_collisions

    INTEGER :: ispecies, jspecies
    INTEGER(i8) :: ix
    TYPE(particle_list), POINTER :: p_list1
    REAL(num), DIMENSION(:), ALLOCATABLE :: idens, jdens
    REAL(num), DIMENSION(:), ALLOCATABLE :: itemp, jtemp, log_lambda
    REAL(num) :: user_factor, q1, q2, m1, m2, w1, w2
    LOGICAL :: collide_species

    ALLOCATE(idens(-2:nx+3))
    ALLOCATE(jdens(-2:nx+3))
    ALLOCATE(itemp(-2:nx+3))
    ALLOCATE(jtemp(-2:nx+3))
    ALLOCATE(log_lambda(-2:nx+3))
    ALLOCATE(meanx(-2:nx+3))
    ALLOCATE(meany(-2:nx+3))
    ALLOCATE(meanz(-2:nx+3))
    ALLOCATE(part_count(-2:nx+3))

    DO ispecies = 1, n_species
      ! Currently no support for photon collisions so just cycle round
      IF (species_list(ispecies)%species_type .EQ. c_species_id_photon) &
          CYCLE
      ! Currently no support for collisions involving chargeless particles
      IF (ABS(species_list(ispecies)%charge) .LE. c_tiny) &
          CYCLE

      collide_species = .FALSE.
      DO jspecies = ispecies, n_species
        user_factor = coll_pairs(ispecies, jspecies)
        IF (user_factor .GT. 0) THEN
          collide_species = .TRUE.
          EXIT
        ENDIF
      ENDDO

      IF (.NOT.collide_species) CYCLE

      CALL calc_coll_number_density(idens, ispecies)
      CALL calc_coll_temperature(itemp, ispecies)

      m1 = species_list(ispecies)%mass
      q1 = species_list(ispecies)%charge
      w1 = species_list(ispecies)%weight
      itemp = itemp * kb / q0

      DO ix = 1, nx
        p_list1 => species_list(ispecies)%secondary_list(ix)
        CALL shuffle_particle_list_random(p_list1)
      ENDDO ! ix

      DO jspecies = ispecies, n_species
        ! Currently no support for photon collisions so just cycle round
        IF (species_list(jspecies)%species_type .EQ. c_species_id_photon) &
            CYCLE
        ! Currently no support for collisions involving chargeless particles
        IF (ABS(species_list(jspecies)%charge) .LE. c_tiny) &
            CYCLE
        user_factor = coll_pairs(ispecies, jspecies)
        IF (user_factor .LE. 0) CYCLE

        CALL calc_coll_number_density(jdens, jspecies)
        CALL calc_coll_temperature(jtemp, jspecies)

        m2 = species_list(jspecies)%mass
        q2 = species_list(jspecies)%charge
        w2 = species_list(jspecies)%weight
        jtemp = jtemp * kb / q0

        IF (coulomb_log_auto) THEN
          log_lambda = calc_coulomb_log(itemp, jdens, q1, q2)
        ELSE
          log_lambda = coulomb_log
        ENDIF

        DO ix = 1, nx
          IF (ispecies .EQ. jspecies) THEN
            CALL intra_species_collisions( &
                species_list(ispecies)%secondary_list(ix), &
                m1, q1, w1, idens(ix), itemp(ix), &
                log_lambda(ix), user_factor)
          ELSE
            CALL inter_species_collisions( &
                species_list(ispecies)%secondary_list(ix), &
                species_list(jspecies)%secondary_list(ix), &
                m1, m2, q1, q2, w1, w2, idens(ix), jdens(ix), &
                itemp(ix), jtemp(ix), log_lambda(ix), &
                user_factor)
          ENDIF
        ENDDO ! ix
      ENDDO ! jspecies
    ENDDO ! ispecies

    DEALLOCATE(idens, jdens, itemp, jtemp, log_lambda)
    DEALLOCATE(meanx, meany, meanz, part_count)

  END SUBROUTINE particle_collisions



#ifdef PER_PARTICLE_WEIGHT
  SUBROUTINE collisional_ionisation

    INTEGER :: ispecies, jspecies, ion_species, e_species, n1, n2, l
    INTEGER(i8) :: ix
    TYPE(particle_list), POINTER :: p_list1
    TYPE(particle_list) :: ionising_e, ejected_e
    REAL(num), DIMENSION(:), ALLOCATABLE :: idens, jdens, e_dens
    REAL(num), DIMENSION(:), ALLOCATABLE :: itemp, jtemp, e_temp
    REAL(num), DIMENSION(:), ALLOCATABLE :: log_lambda, e_log_lambda
    REAL(num) :: user_factor, e_user_factor, q1, q2, m1, m2, w1, w2
    REAL(num) :: q_e, m_e, w_e, q_full, ionisation_energy
    LOGICAL :: use_coulomb_log_auto_i, use_coulomb_log_auto

    DO ix = 1, nx
      DO ispecies = 1, n_species
        p_list1 => species_list(ispecies)%secondary_list(ix)
        CALL shuffle_particle_list_random(p_list1)
      ENDDO
    ENDDO ! ix

    ALLOCATE(idens(-2:nx+3))
    ALLOCATE(jdens(-2:nx+3))
    ALLOCATE(e_dens(-2:nx+3))
    ALLOCATE(itemp(-2:nx+3))
    ALLOCATE(jtemp(-2:nx+3))
    ALLOCATE(e_temp(-2:nx+3))
    ALLOCATE(log_lambda(-2:nx+3))
    ALLOCATE(e_log_lambda(-2:nx+3))
    ALLOCATE(meanx(-2:nx+3))
    ALLOCATE(meany(-2:nx+3))
    ALLOCATE(meanz(-2:nx+3))
    ALLOCATE(part_count(-2:nx+3))

    CALL create_empty_partlist(ionising_e)
    CALL create_empty_partlist(ejected_e)

    DO ispecies = 1, n_species
      ! Currently no support for photon collisions so just cycle round
      IF (species_list(ispecies)%species_type .EQ. c_species_id_photon) &
          CYCLE
      ! Currently no support for collisions involving chargeless particles
      ! unless ionisation occurs
      use_coulomb_log_auto_i = .TRUE.
      IF (species_list(ispecies)%charge .EQ. 0.0_num) THEN
        IF (.NOT. species_list(ispecies)%ionise) CYCLE
        use_coulomb_log_auto_i = .FALSE.
      ENDIF
      CALL calc_coll_number_density(idens, ispecies)
      CALL calc_coll_temperature(itemp, ispecies)

      m1 = species_list(ispecies)%mass
      q1 = species_list(ispecies)%charge
      w1 = species_list(ispecies)%weight
      itemp = itemp * kb / q0

      IF (species_list(ispecies)%ionise) THEN
        e_species = species_list(ispecies)%release_species
        CALL calc_coll_number_density(e_dens, e_species)
        CALL calc_coll_temperature(e_temp, e_species)
        m_e = species_list(e_species)%mass
        q_e = species_list(e_species)%charge
        w_e = species_list(e_species)%weight
        e_temp = e_temp * kb / q0
        n1 = species_list(ispecies)%n
        l = species_list(ispecies)%l
        ionisation_energy = species_list(ispecies)%ionisation_energy / ev
        ion_species = species_list(ispecies)%ionise_to_species
        n2 = species_list(ispecies)%n
        DO WHILE(species_list(ion_species)%ionise)
          ion_species = species_list(ion_species)%ionise_to_species
        ENDDO
        q_full = species_list(ion_species)%charge
        ion_species = species_list(ispecies)%ionise_to_species
      ENDIF

      DO jspecies = ispecies, n_species
        ! Currently no support for photon collisions so just cycle round
        IF (species_list(jspecies)%species_type .EQ. c_species_id_photon) &
            CYCLE
        ! Currently no support for collisions involving chargeless particles
        ! unless ionisation occurs
        use_coulomb_log_auto = use_coulomb_log_auto_i
        IF (species_list(jspecies)%charge .EQ. 0.0_num) THEN
          IF (.NOT. (species_list(ispecies)%electron &
              .AND. species_list(jspecies)%ionise)) CYCLE
          use_coulomb_log_auto = .FALSE.
        ENDIF
        user_factor = coll_pairs(ispecies, jspecies)
        IF (user_factor .LE. 0) CYCLE

        CALL calc_coll_number_density(jdens, jspecies)
        CALL calc_coll_temperature(jtemp, jspecies)

        m2 = species_list(jspecies)%mass
        q2 = species_list(jspecies)%charge
        w2 = species_list(jspecies)%weight
        jtemp = jtemp * kb / q0

        IF (species_list(ispecies)%electron &
            .AND. species_list(jspecies)%ionise) THEN
          e_species = species_list(jspecies)%release_species
          CALL calc_coll_number_density(e_dens, e_species)
          CALL calc_coll_temperature(e_temp, e_species)
          m_e = species_list(e_species)%mass
          q_e = species_list(e_species)%charge
          w_e = species_list(e_species)%weight
          e_temp = e_temp * kb / q0
          n1 = species_list(jspecies)%n
          l = species_list(jspecies)%l
          ionisation_energy = species_list(jspecies)%ionisation_energy / ev
          ion_species = species_list(jspecies)%ionise_to_species
          n2 = species_list(ion_species)%n
          DO WHILE(species_list(ion_species)%ionise)
            ion_species = species_list(ion_species)%ionise_to_species
          ENDDO
          q_full = species_list(ion_species)%charge
          ion_species = species_list(jspecies)%ionise_to_species
        ENDIF

        IF (coulomb_log_auto) THEN
          IF (use_coulomb_log_auto) THEN
            log_lambda = calc_coulomb_log(itemp, jdens, q1, q2)
          ELSE
            log_lambda = 0
          ENDIF
          IF (species_list(ispecies)%electron &
              .AND. species_list(jspecies)%ionise) THEN
            e_log_lambda = calc_coulomb_log(itemp, e_dens, q1, q_e)
            e_user_factor = coll_pairs(ispecies, ion_species)
          ELSE IF (species_list(ispecies)%ionise &
              .AND. species_list(jspecies)%electron) THEN
            e_log_lambda = calc_coulomb_log(e_temp, jdens, q_e, q2)
            e_user_factor = coll_pairs(ion_species, jspecies)
          ENDIF
        ELSE
          log_lambda = coulomb_log
          e_log_lambda = coulomb_log
        ENDIF

        IF (ispecies .EQ. jspecies) THEN
          DO ix = 1, nx
            CALL intra_species_collisions( &
                species_list(ispecies)%secondary_list(ix), &
                m1, q1, w1, idens(ix), itemp(ix), &
                log_lambda(ix), user_factor)
          ENDDO ! ix
        ELSE IF (species_list(ispecies)%ionise &
            .AND. species_list(jspecies)%electron) THEN
          DO ix = 1, nx
            ! Perform collisional ionisation before calculating scatter
            CALL preionise(species_list(jspecies)%secondary_list(ix), &
                species_list(ispecies)%secondary_list(ix), &
                species_list(ion_species)%secondary_list(ix), &
                ionising_e, ejected_e, m2, m1, q2, q1, jdens(ix), &
                q_full, ionisation_energy, n1, n2, l)
            ! Scatter ionising impact electrons off of ejected target electrons
            ! unless specified otherwise in input deck
            IF (e_user_factor .GT. 0.0_num) THEN
              CALL inter_species_collisions( &
                  ejected_e, ionising_e, &
                  m_e, m2, q_e, q2, w_e, w2, e_dens(ix), jdens(ix), &
                  e_temp(ix), jtemp(ix), e_log_lambda(ix), &
                  e_user_factor)
            ENDIF
            ! Scatter non-ionising impact electrons off of remaining unionised
            ! targets provided target has charge
            IF (q1 .NE. 0.0_num) THEN
              CALL inter_species_collisions( &
                  species_list(ispecies)%secondary_list(ix), &
                  species_list(jspecies)%secondary_list(ix), &
                  m1, m2, q1, q2, w1, w2, idens(ix), jdens(ix), &
                  itemp(ix), jtemp(ix), log_lambda(ix), &
                  user_factor)
            ENDIF
            ! Put ions and electrons into respective lists
            CALL append_partlist( &
                species_list(jspecies)%secondary_list(ix), ionising_e)
            CALL append_partlist( &
                species_list(jspecies)%secondary_list(ix), ejected_e)
          ENDDO ! ix
        ELSE IF (species_list(ispecies)%electron &
            .AND. species_list(jspecies)%ionise) THEN
          DO ix = 1, nx
            ! Perform collisional ionisation before calculating scatter
            CALL preionise(species_list(ispecies)%secondary_list(ix), &
                species_list(jspecies)%secondary_list(ix), &
                species_list(ion_species)%secondary_list(ix), &
                ionising_e, ejected_e, m1, m2, q1, q2, idens(ix), &
                q_full, ionisation_energy, n1, n2, l)
            ! Scatter ionising impact electrons off of ejected target electrons
            ! unless specified otherwise in input deck
            IF (e_user_factor .GT. 0.0_num) THEN
              CALL inter_species_collisions( &
                  ejected_e, ionising_e, &
                  m1, m_e, q1, q_e, w1, w_e, idens(ix), e_dens(ix), &
                  itemp(ix), e_temp(ix), e_log_lambda(ix), &
                  e_user_factor)
            ENDIF
            ! Scatter non-ionising impact electrons off of remaining unionised
            ! targets provided target has charge
            IF (q2 .NE. 0.0_num) THEN
              CALL inter_species_collisions( &
                  species_list(ispecies)%secondary_list(ix), &
                  species_list(jspecies)%secondary_list(ix), &
                  m1, m2, q1, q2, w1, w2, idens(ix), jdens(ix), &
                  itemp(ix), jtemp(ix), log_lambda(ix), &
                  user_factor)
            ENDIF
            ! Put electrons into respective lists
            CALL append_partlist( &
                species_list(ispecies)%secondary_list(ix), ionising_e)
            CALL append_partlist( &
                species_list(ispecies)%secondary_list(ix), ejected_e)
          ENDDO ! ix
        ELSE
          DO ix = 1, nx
            CALL inter_species_collisions( &
                species_list(ispecies)%secondary_list(ix), &
                species_list(jspecies)%secondary_list(ix), &
                m1, m2, q1, q2, w1, w2, idens(ix), jdens(ix), &
                itemp(ix), jtemp(ix), log_lambda(ix), &
                user_factor)
          ENDDO ! ix
        ENDIF
      ENDDO ! jspecies
    ENDDO ! ispecies

    DEALLOCATE(idens, jdens, itemp, jtemp, log_lambda)
    DEALLOCATE(meanx, meany, meanz, part_count)
    DEALLOCATE(e_dens, e_temp, e_log_lambda)

  END SUBROUTINE collisional_ionisation
#endif



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

    REAL(num) :: factor, np, eiics, reduced_energy, fion, gr, red_inc, red_ion,&
        i_p2, rot_y, rot_z, e_p2_i, e_e, beta_i, gamma_i, e_p_rot(3), e_ke_i, &
        e_v_i, gamma_e_i, mrbeb_c, t, tp, bp, bt2, bb2
    REAL(num) :: ionisation_energy_inv, red_ion_inv, prob_factor
    INTEGER(KIND=8) :: e_count, ion_count, pcount, i, k

    ! Inter-species collisions
    e_count = electrons%count
    ion_count = ions%count

    ! If there aren't enough particles to collide, then don't bother
    IF(e_count .EQ. 0 .OR. ion_count .EQ. 0) RETURN

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
    ENDDO

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
      IF(i_p2 .GT. 0.0_num) THEN
        IF(ion%part_p(1) .NE. 0.0_num) THEN
          rot_y = DATAN(ion%part_p(3) / ion%part_p(1))
        ELSE
          rot_y = pi / 2.0_num
        ENDIF
        IF(ion%part_p(1) * DCOS(rot_y) + ion%part_p(3) * DSIN(rot_y) &
            .NE. 0.0_num) THEN
          rot_z = DATAN(-ion%part_p(2) / (ion%part_p(1) * DCOS(rot_y) &
              + ion%part_p(3) * DSIN(rot_y)))
        ELSE
          rot_z = pi / 2.0_num
        ENDIF
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
        e_v_i = SQRT(e_p2_i / (e_mass**2 + e_p2_i / c**2))
      ELSE
        e_p2_i = DOT_PRODUCT(electron%part_p, electron%part_p)
        e_ke_i = c * (SQRT(e_p2_i + e_mass * e_rest) - e_mass * c) / ev
        e_v_i = SQRT(e_p2_i / (e_mass**2 + e_p2_i / c**2))
      ENDIF
      ! Must enforce that electrons with insufficient kinetic energies cannot
      ! cause ionisation, as all cross sectional models used show massively
      ! increasing electron impact ionisation cross section as kinetic energy
      ! tends to zero
      IF (e_ke_i .GE. ion%weight / electron%weight * ionisation_energy &
          .AND. .NOT. was_ionised(MOD(k - 1, ion_count) + 1)) THEN
        ! Find cross section
        red_inc = e_ke_i * ionisation_energy_inv
        ! Use MBELL model for atomic number < 36
        IF (n1 .LT. 4 .AND. l .LT. 3) THEN
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
          ENDDO
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
              * 0.7_num * ((ion_charge / q0 + 1.0_num) / n2)**2)

          ! MRBEB cross section (cm^2)
          eiics = mrbeb_const / (bt2 + mrbeb_c * bb2) / bp &
              * (0.5_num * (LOG(bt2 / (1.0_num - bt2)) - bt2 &
              - LOG(2.0_num * bp)) * (1.0_num - 1.0_num / t**2) &
              + 1.0_num - 1.0_num / t &
              - LOG(t) / (t + 1.0_num) * (1.0_num + 2.0_num * tp) &
              / (1.0_num + 0.5_num * tp)**2 &
              + bp**2 / (1.0_num + 0.5_num * tp)**2 * (t - 1.0_num) / 2.0_num)
        END IF
        IF (random() .LT. 1.0_num - EXP(prob_factor * eiics * e_v_i)) THEN
          ! Mark ionisation as occurring
          was_ionised(MOD(k - 1, ion_count) + 1) = .TRUE.
          lost_ke(MOD(k - 1, e_count) + 1) = .TRUE.
          IF(i_p2 .GT. 0.0_num) THEN
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
            IF(DOT_PRODUCT(electron%part_p, electron%part_p) .LT. &
                DOT_PRODUCT(e_p_rot, e_p_rot)) THEN
              electron%part_p = SQRT(((ev / c * (e_ke_i - ion%weight &
                / electron%weight * ionisation_energy + e_rest_ev))**2 &
                - e_mass * e_rest) / e_p2_i) * electron%part_p
            ELSE
              electron%part_p = e_p_rot
            ENDIF
          ELSE
            electron%part_p = SQRT(((ev / c * (e_ke_i - ion%weight &
                / electron%weight * ionisation_energy + e_rest_ev))**2 &
                - e_mass * e_rest) / e_p2_i) * electron%part_p
          ENDIF
        ENDIF
      ENDIF
      ion => ion%next
      electron => electron%next
    ENDDO

    ! restore the tail of the lists
    NULLIFY(electrons%tail%next)
    NULLIFY(ions%tail%next)

    electron => electrons%head
    ion => ions%head

    DO k = 1, MAX(e_count, ion_count)
      IF (k .LE. e_count) THEN
        next_e => electron%next
        IF (lost_ke(k)) THEN
          CALL remove_particle_from_partlist(electrons, electron)
          CALL add_particle_to_partlist(ionising_e, electron)
        ENDIF
      ENDIF
      IF (k .LE. ion_count) THEN
        next_ion => ion%next
        IF (was_ionised(k)) THEN
          ALLOCATE(ejected_electron)
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
        ENDIF
      ENDIF
      electron => next_e
      ion => next_ion
    ENDDO

    DEALLOCATE(lost_ke, was_ionised)

  END SUBROUTINE preionise



  SUBROUTINE intra_species_collisions(p_list, mass, charge, weight, &
      dens, temp, log_lambda, user_factor)
    ! Perform collisions between particles of the same species.

    TYPE(particle_list), INTENT(INOUT) :: p_list
    REAL(num), INTENT(IN) :: mass, charge, weight
    REAL(num), INTENT(IN) :: user_factor
    TYPE(particle), POINTER :: current, impact
    REAL(num) :: factor, np
    REAL(num) :: dens, temp, log_lambda
    INTEGER(i8) :: icount, k

    factor = 0.0_num
    np = 0.0_num

    ! Intra-species collisions
    icount = p_list%count

    ! If there aren't enough particles to collide, then don't bother
    IF (icount .LE. 1) RETURN

    ! No collisions in cold plasma so return
    IF (temp .LE. c_tiny) RETURN

#ifndef PER_PARTICLE_WEIGHT
    np = icount * weight
    factor = user_factor
#else
    current => p_list%head
    impact => current%next
    DO k = 2, icount-2, 2
      np = np + current%weight + impact%weight
      factor = factor + MIN(current%weight, impact%weight)
      current => impact%next
      impact => current%next
#ifdef PREFETCH
      CALL prefetch_particle(current)
      CALL prefetch_particle(impact)
#endif
    ENDDO
    np = np + current%weight + impact%weight
    factor = factor + MIN(current%weight, impact%weight)

    IF (MOD(icount, 2_i8) .NE. 0) THEN
      np = np + impact%next%weight
      factor = factor + MIN(current%weight, impact%next%weight)
      factor = factor + MIN(impact%weight, impact%next%weight)
    ENDIF

    factor = user_factor * np / factor
#endif

    current => p_list%head
    impact => current%next
    DO k = 2, icount-2, 2
      CALL scatter(current, impact, mass, mass, charge, charge, &
          weight, weight, dens, dens, temp, temp, log_lambda, factor)
      current => impact%next
      impact => current%next
#ifdef PREFETCH
      CALL prefetch_particle(current)
      CALL prefetch_particle(impact)
#endif
    ENDDO

    IF (MOD(icount, 2_i8) .EQ. 0) THEN
      CALL scatter(current, impact, mass, mass, charge, charge, &
          weight, weight, dens, dens, temp, temp, log_lambda, factor)
    ELSE
      CALL scatter(current, impact, mass, mass, charge, charge, &
          weight, weight, dens, dens, temp, temp, log_lambda, 0.5_num*factor)
      current => impact%next
      impact => current%prev%prev
      CALL scatter(current, impact, mass, mass, charge, charge, &
          weight, weight, dens, dens, temp, temp, log_lambda, 0.5_num*factor)
      current => current%prev
      impact => current%next
      CALL scatter(current, impact, mass, mass, charge, charge, &
          weight, weight, dens, dens, temp, temp, log_lambda, 0.5_num*factor)
    ENDIF

  END SUBROUTINE intra_species_collisions



  SUBROUTINE inter_species_collisions(p_list1, p_list2, mass1, mass2, &
      charge1, charge2, weight1, weight2, &
      idens, jdens, itemp, jtemp, log_lambda, user_factor )

    TYPE(particle_list), INTENT(INOUT) :: p_list1
    TYPE(particle_list), INTENT(INOUT) :: p_list2

    REAL(num), INTENT(IN) :: mass1, charge1, weight1
    REAL(num), INTENT(IN) :: mass2, charge2, weight2

    REAL(num), INTENT(IN) :: idens, jdens
    REAL(num), INTENT(IN) :: itemp, jtemp, log_lambda
    REAL(num), INTENT(IN) :: user_factor

    TYPE(particle), POINTER :: current, impact

    REAL(num) :: factor, np
    INTEGER(i8) :: icount, jcount, pcount, k

    factor = 0.0_num
    np = 0.0_num

    ! No collisions in cold plasma so return
    IF (itemp .LE. c_tiny .AND. jtemp .LE. c_tiny) RETURN

    ! Inter-species collisions
    icount = p_list1%count
    jcount = p_list2%count
    pcount = MAX(icount, jcount)

    IF (icount .GT. 0 .AND. jcount .GT. 0) THEN
      ! temporarily join tail to the head of the lists to make them circular
      p_list1%tail%next => p_list1%head
      p_list2%tail%next => p_list2%head

#ifndef PER_PARTICLE_WEIGHT
      np = pcount * weight1
      factor = pcount * MIN(weight1, weight2)
#else
      current => p_list1%head
      impact => p_list2%head

      DO k = 1, pcount
        np = np + current%weight
        factor = factor + MIN(current%weight, impact%weight)
        current => current%next
        impact => impact%next
#ifdef PREFETCH
        CALL prefetch_particle(current)
        CALL prefetch_particle(impact)
#endif
      ENDDO
#endif

      current => p_list1%head
      impact => p_list2%head
      DO k = 1, pcount
        CALL scatter(current, impact, mass1, mass2, charge1, charge2, &
            weight1, weight2, idens, jdens, itemp, jtemp, &
            log_lambda, user_factor * np / factor)
        current => current%next
        impact => impact%next
#ifdef PREFETCH
        CALL prefetch_particle(current)
        CALL prefetch_particle(impact)
#endif
      ENDDO

      ! restore the tail of the lists
      NULLIFY(p_list1%tail%next)
      NULLIFY(p_list2%tail%next)
    ENDIF

  END SUBROUTINE inter_species_collisions



  SUBROUTINE scatter(current, impact, mass1, mass2, charge1, charge2, &
      weight1, weight2, idens, jdens, itemp, jtemp, log_lambda, factor)

    ! Here the Coulomb collisions are performed by rotating the momentum
    ! vector of one of the particles in the centre of momentum reference
    ! frame. Then, as the collisions are assumed to be elastic, the momentum
    ! of the second particle is simply the opposite of the first.
    ! This algorithm is based on that of Y. Sentoku and A. J. Kemp.
    ! [J Comput Phys, 227, 6846 (2008)]

    TYPE(particle), POINTER :: current, impact
    REAL(num), INTENT(IN) :: mass1, mass2
    REAL(num), INTENT(IN) :: charge1, charge2
    REAL(num), INTENT(IN) :: weight1, weight2
    REAL(num), INTENT(IN) :: idens, jdens, itemp, jtemp, log_lambda
    REAL(num), INTENT(IN) :: factor
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
    REAL(num) :: gc, gcr
    REAL(num) :: tvar ! Dummy variable for temporarily storing values
    REAL(num) :: vc_sq, vc_mag, p1_vc, p2_vc, p3_mag
    REAL(num) :: delta, sin_theta, cos_theta, tan_theta_cm, tan_theta_cm2
    REAL(num) :: vrabs
    REAL(num) :: nu, ran1, ran2
    !REAL(num) :: m_red

    ! Copy all of the necessary particle data into variables with easier to
    ! read names
    p1 = current%part_p
    p2 = impact%part_p
    p1_norm = p1 / mc0
    p2_norm = p2 / mc0

    ! Two stationary particles can't collide, so don't try
    IF (DOT_PRODUCT(p1_norm, p1_norm) .LT. eps &
        .AND. DOT_PRODUCT(p2_norm, p2_norm) .LT. eps) RETURN

    ! Ditto for two particles with the same momentum
    vc = (p1_norm - p2_norm)
    IF (DOT_PRODUCT(vc, vc) .LT. eps) RETURN

#ifdef PER_PARTICLE_CHARGE_MASS
    m1 = current%mass
    m2 = impact%mass
    q1 = current%charge
    q2 = impact%charge
#else
    m1 = mass1
    m2 = mass2
    q1 = charge1
    q2 = charge2
#endif

    ! Pre-collision energies
    e1 = c * SQRT(DOT_PRODUCT(p1, p1) + (m1 * c)**2)
    e2 = c * SQRT(DOT_PRODUCT(p2, p2) + (m2 * c)**2)

    ! Velocity of centre-of-momentum (COM) reference frame
    vc = (p1 + p2) * c**2 / (e1 + e2)
    vc_sq = DOT_PRODUCT(vc, vc)
    vc_mag = SQRT(vc_sq)
    gc = 1.0_num / SQRT(1.0_num - vc_sq / c**2)

    ! Lorentz momentum transform to get into COM frame
    p1_vc = DOT_PRODUCT(p1, vc)
    p2_vc = DOT_PRODUCT(p2, vc)
    tvar = p1_vc * (gc - 1.0_num) / (vc_sq + c_tiny)
    p3 = p1 + vc * (tvar - gc * e1 / c**2)
    tvar = p2_vc * (gc - 1.0_num) / (vc_sq + c_tiny)
    p4 = p2 + vc * (tvar - gc * e2 / c**2)

    p3_mag = SQRT(DOT_PRODUCT(p3, p3))

    ! Lorentz energy transform
    e3 = gc * (e1 - p1_vc)
    e4 = gc * (e2 - p2_vc)
    ! Pre-collision velocities in COM frame
    v3 = p3 * c**2 / e3
    v4 = p4 * c**2 / e4

    ! Relative velocity
    tvar = 1.0_num - (DOT_PRODUCT(v3, v4) / c**2)
    vr = (v3 - v4) / tvar
    vrabs = SQRT(DOT_PRODUCT(vr, vr))

    ! Collision frequency
    nu = coll_freq(vrabs, log_lambda, m1, m2, q1, q2, itemp, jtemp, jdens)
    nu = nu * factor * dt

!    m_red = mass1 * mass2 / (mass1 + mass2)
!    nu = ((idens * (charge1 * charge2)**2 * log_lambda) &
!        / (8.0_num * pi * (epsilon0**2) * (m_red**2) * (vrabs**3))) &
!        * gc * dt * factor

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

    ! angle theta in the One Particle at Rest frame
    sin_theta = (2.0_num * delta) / (1.0_num + delta**2)
    cos_theta = (1.0_num - delta**2) / (1.0_num + delta**2)

    ! Transform angles from particle j's rest frame to COM frame
    ! Note azimuthal angle (ran2) is invariant under this transformation
    vcr = -v4
    gcr = 1.0_num / SQRT(1.0_num - (DOT_PRODUCT(vcr, vcr) / c**2))

    tan_theta_cm = sin_theta &
        / (gcr * (cos_theta - SQRT(DOT_PRODUCT(vcr, vcr)) / vrabs))
    tan_theta_cm2 = tan_theta_cm**2

    sin_theta = SQRT(tan_theta_cm2 / (1 + tan_theta_cm2))
    cos_theta = SQRT(1.0_num / (1.0_num + tan_theta_cm2))

    ! Post-collision momenta in COM frame
    p3 = p3_mag * (c1 * cos_theta + c2 * sin_theta * COS(ran2) &
        + c3 * sin_theta * SIN(ran2))
    p4 = -p3

    ! Lorentz momentum transform to get back to lab frame
    tvar = DOT_PRODUCT(p3, vc) * (gc - 1.0_num) / vc_sq
    p5 = p3 + vc * (tvar + gc * e3 / c**2)
    tvar = DOT_PRODUCT(p4, vc) * (gc - 1.0_num) / vc_sq
    p6 = p4 + vc * (tvar + gc * e4 / c**2)

    e5 = c * SQRT(DOT_PRODUCT(p5, p5) + (m1 * c)**2)
    e6 = c * SQRT(DOT_PRODUCT(p6, p6) + (m2 * c)**2)

#ifdef PER_PARTICLE_WEIGHT
    w1 = current%weight
    w2 = impact%weight
#else
    w1 = weight1
    w2 = weight2
#endif

    IF (w1 .GT. w2) THEN
      CALL weighted_particles_correction(w2 / w1, p1, p5, e1, e5, m1)
    ELSEIF (w2 .GT. w1) THEN
      CALL weighted_particles_correction(w1 / w2, p2, p6, e2, e6, m2)
    ENDIF

    ! Update particle properties
    current%part_p = p5
    impact%part_p = p6

  END SUBROUTINE scatter



  PURE FUNCTION coll_freq(vrabs, log_lambda, m1, m2, q1, q2, itemp, jtemp, &
      jdens)

    REAL(num), INTENT(IN) :: vrabs, log_lambda, m1, m2, q1, q2
    REAL(num), INTENT(IN) :: itemp, jtemp, jdens
    REAL(num) :: mu, coll_freq

    mu = (m1 * m2) / (m1 + m2)
    coll_freq = velocity_collisions(vrabs, log_lambda, mu, q1, q2, jdens)
!    coll_freq = temperature_collisions(itemp, log_lambda, mu, q1, q2, jdens)
!    coll_freq = manheimer_collisions(vrabs, log_lambda, m1, m2, q1, q2, &
!        jtemp, jdens)
!    coll_freq = MAX(coll_freq, vrabs / (jdens**(1.0_num / 3.0_num)))

  END FUNCTION



  PURE FUNCTION velocity_collisions(vrabs, log_lambda, mu, q1, q2, jdens)

    REAL(num), INTENT(IN) :: vrabs, log_lambda, mu, q1, q2, jdens
    REAL(num), PARAMETER :: fac = 4.0_num * pi * epsilon0**2
    REAL(num) :: numerator, denominator
    REAL(num) :: velocity_collisions

    IF (vrabs .GT. 0.0_num) THEN
      numerator = (q1 * q2)**2 * jdens * log_lambda
      denominator = fac * mu**2 * vrabs**3
      IF (denominator .LE. 0.0_num &
          .OR. EXPONENT(numerator) - EXPONENT(denominator) &
          .GE. c_maxexponent) THEN
        velocity_collisions = 0.0_num
      ELSE
        velocity_collisions = numerator / denominator
      ENDIF
    ELSE
      velocity_collisions = 0.0_num
    ENDIF

  END FUNCTION velocity_collisions



  PURE FUNCTION temperature_collisions(itemp, log_lambda, mu, q1, q2, jdens)

    REAL(num), INTENT(IN) :: itemp, log_lambda, mu, q1, q2, jdens
    REAL(num) :: temperature_collisions

    IF (itemp .GT. c_tiny) THEN
      temperature_collisions = ((q1 * q2)**2 * jdens * log_lambda) &
          / (3.0_num * epsilon0**2 * SQRT(mu) &
          * (2.0_num * pi * q0 * itemp)**1.5_num)
    ELSE
      temperature_collisions = 0.0_num
    ENDIF

  END FUNCTION temperature_collisions



  PURE FUNCTION manheimer_collisions(vrabs, log_lambda, m1, m2, q1, q2, &
      jtemp, jdens)

    REAL(num), INTENT(IN) :: vrabs, log_lambda, m1, m2, q1, q2, jtemp, jdens
    REAL(num) :: gr, mu, ek, slow, fast
    REAL(num) :: manheimer_collisions

    ! Manheimer-like collision operator
    ! Valid for e-i and e-e collisions
    gr = 1.0_num / SQRT(1.0_num - (vrabs / c)**2)
    mu = m2 / 1.6726d-27
    ek = (gr - 1.0_num) * m1 * c**2 / q0

    IF (jtemp .LE. 0.0_num) THEN
      IF (ek .LE. 0.0_num) THEN
        manheimer_collisions = 0.0_num
      ELSE
        manheimer_collisions = 3.9d-6 / (SQRT(ek**3) + c_tiny)
      ENDIF
    ELSE
      IF (ek .LE. 0.0_num) THEN
        manheimer_collisions = 0.23_num * SQRT((mu / (jtemp + c_tiny))**3)
      ELSE
        slow = 0.23_num * SQRT((mu / (jtemp + c_tiny))**3)
        fast = 3.9d-6 / (SQRT(ek**3) + c_tiny)
        manheimer_collisions = slow / (1.0_num + slow / fast)
      ENDIF
    ENDIF
    manheimer_collisions = manheimer_collisions * jdens * log_lambda &
        * (q2 / q0)**2 * 1.0d-6

  END FUNCTION manheimer_collisions



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

    en_after = (1 - wtr) * en + wtr * en_scat
    p_after  = (1 - wtr) * p  + wtr * p_scat
    p_mag = SQRT(DOT_PRODUCT(p_after, p_after))
!    gamma_en = 1 + en_after / (mass * c**2)
    gamma_en = en_after / (mass * c**2)
    gamma_p = SQRT(1 + (p_mag / mass / c)**2)

    ! This if-statement is just to take care of possible rounding errors
    ! gamma_p should always be smaller than gamma_en
    IF (gamma_p .LT. gamma_en) THEN
      ! magnitude of the momentum correction
      delta_p = mass * c * SQRT(gamma_en**2 - gamma_p**2)
      p_trans_mag = SQRT(p_after(2)**2 + p_after(3)**2)

      CALL new_coords(p_after, c1, c2, c3)

      phi = 2.0_num * pi * random()

      ! Correcting for the loss in energy by adding a perpendicular
      ! momentum correction
      p_scat = p_after + delta_p * (c2 * COS(phi) + c3 * SIN(phi))
    ENDIF

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

    IF (vtrans .GT. c_tiny) THEN
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
    ENDIF

  END SUBROUTINE new_coords



  SUBROUTINE shuffle_particle_list_random(p_list)

    TYPE(particle_list), INTENT(INOUT) :: p_list
    TYPE(particle), POINTER :: particle1, particle2

    INTEGER :: i, idx, swap_idx
    INTEGER :: p_num

    p_num = INT(p_list%count)

    ! Nothing to be done
    IF (p_num .LE. 2) RETURN

    ! First make sure that the sorting array is large enough
    ! This should happen infrequently
    IF (p_num .GT. coll_sort_array_size) THEN
      DEALLOCATE(coll_sort_array)

      ! make the sort array somewhat larger to avoid frequent deallocation
      ! and reallocation
      coll_sort_array_size = (11 * p_num) / 10 + 10
      ALLOCATE(coll_sort_array(coll_sort_array_size))
    ENDIF

    ! Copy all the particle pointers into the array and create random
    ! sort indices
    particle1 => p_list%head
    DO i = 1,p_num
      coll_sort_array(i)%particle => particle1
      particle1 => particle1%next
    ENDDO

    ! Shuffle particles using Durstenfeld's algorithm
    DO idx = p_num,2,-1
      swap_idx = FLOOR(idx * random()) + 1
      particle1 => coll_sort_array(idx)%particle
      coll_sort_array(idx)%particle => coll_sort_array(swap_idx)%particle
      coll_sort_array(swap_idx)%particle => particle1
    ENDDO

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
    ENDDO

    ! Finally set the tail (at the end of the loop, particle is pointing to
    ! the tail)
    p_list%tail => particle1
    NULLIFY(particle1%next)

  END SUBROUTINE shuffle_particle_list_random



  PURE FUNCTION calc_coulomb_log(temp, dens, q1, q2)

    REAL(num), DIMENSION(-2:), INTENT(IN) :: temp, dens
    REAL(num), INTENT(IN) :: q1, q2
    REAL(num), DIMENSION(-2:nx+3) :: calc_coulomb_log
    REAL(num), PARAMETER :: cfac = 4.13d6 * SQRT((q0 / kb)**3)
    REAL(num), PARAMETER :: exp1 = 2.7182818284590452353602874713526625_num
    REAL(num) :: fac, efac, lfac, temp3, den, ratio
    INTEGER :: i

    fac  = cfac * (q0 / q1)**2 * ABS(q0 / q2)
    efac = (exp1 / fac)**2
    lfac = LOG(fac)

    DO i = -2, nx+3
      temp3 = temp(i)**3
      den = dens(i)
      IF (den .LE. 0.0_num &
          .OR. EXPONENT(temp3) - EXPONENT(den) .GE. c_maxexponent) THEN
        calc_coulomb_log(i) = 1.0_num
      ELSE
        ratio = temp3 / den
        IF (ratio .LE. efac) THEN
          calc_coulomb_log(i) = 1.0_num
        ELSE
          calc_coulomb_log(i) = lfac + 0.5_num * LOG(ratio)
        ENDIF
      ENDIF
    ENDDO

  END FUNCTION



  SUBROUTINE calc_coll_number_density(data_array, ispecies)

    ! This subroutine calculates the grid-based number density of a given
    ! particle species.
    ! It is almost identical to the calc_number_density subroutine in calc_df,
    ! except it uses the secondary_list rather than the attached_list.

    REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: ispecies
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: idx
    INTEGER :: ix
    INTEGER :: jx
    TYPE(particle), POINTER :: current
#include "particle_head.inc"

    data_array = 0.0_num

    idx   = 1.0_num / dx

#ifndef PER_PARTICLE_WEIGHT
    wdata = species_list(ispecies)%weight * idx
#endif
    DO jx = 1, nx
      current => species_list(ispecies)%secondary_list(jx)%head
      DO WHILE (ASSOCIATED(current))
#ifdef PER_PARTICLE_WEIGHT
        wdata = current%weight * idx
#endif

#include "particle_to_grid.inc"

        DO ix = sf_min, sf_max
          data_array(cell_x+ix) = data_array(cell_x+ix) + gx(ix) * wdata
        ENDDO ! ix

        current => current%next
      ENDDO
    ENDDO ! jx

    CALL calc_boundary(data_array)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_coll_number_density



  SUBROUTINE calc_coll_temperature(sigma, ispecies)

    ! This subroutine calculates the grid-based temperature of a given
    ! particle species.
    ! It is almost identical to the calc_temperature subroutine in calc_df,
    ! except it uses the secondary_list rather than the attached_list.

    REAL(num), DIMENSION(-2:), INTENT(OUT) :: sigma
    INTEGER, INTENT(IN) :: ispecies
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_pmx, part_pmy, part_pmz, sqrt_part_m
    ! The weight of a particle
    REAL(num) :: l_weight
    REAL(num) :: gf
    INTEGER :: ix
    INTEGER :: jx
    TYPE(particle), POINTER :: current
#include "particle_head.inc"

    meanx = 0.0_num
    meany = 0.0_num
    meanz = 0.0_num
    part_count = 0.0_num
    sigma = 0.0_num

#ifndef PER_PARTICLE_CHARGE_MASS
    sqrt_part_m  = SQRT(species_list(ispecies)%mass)
#endif
#ifndef PER_PARTICLE_WEIGHT
    l_weight = species_list(ispecies)%weight
#endif
    DO jx = 1, nx
      current => species_list(ispecies)%secondary_list(jx)%head
      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

        DO ix = sf_min, sf_max
          gf = gx(ix) * l_weight
          meanx(cell_x+ix) = meanx(cell_x+ix) + gf * part_pmx
          meany(cell_x+ix) = meany(cell_x+ix) + gf * part_pmy
          meanz(cell_x+ix) = meanz(cell_x+ix) + gf * part_pmz
          part_count(cell_x+ix) = part_count(cell_x+ix) + gf
        ENDDO ! ix
        current => current%next
      ENDDO
    ENDDO ! jx

    CALL calc_boundary(meanx)
    CALL calc_boundary(meany)
    CALL calc_boundary(meanz)
    CALL calc_boundary(part_count)

    part_count = MAX(part_count, 1.e-6_num)

    meanx = meanx / part_count
    meany = meany / part_count
    meanz = meanz / part_count

    part_count = 0.0_num
    DO jx = 1, nx
      current => species_list(ispecies)%secondary_list(jx)%head
      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

        DO ix = sf_min, sf_max
          gf = gx(ix)
          sigma(cell_x+ix) = sigma(cell_x+ix) + gf &
              * ((part_pmx - meanx(cell_x+ix))**2 &
              + (part_pmy - meany(cell_x+ix))**2 &
              + (part_pmz - meanz(cell_x+ix))**2)
          part_count(cell_x+ix) = part_count(cell_x+ix) + gf
        ENDDO ! ix
        current => current%next
      ENDDO
    ENDDO ! jx

    CALL calc_boundary(sigma)
    CALL calc_boundary(part_count)

    ! 3/2 kT = <p^2>/(2m)
    sigma = sigma / MAX(part_count, 1.e-6_num) / kb / 3.0_num

  END SUBROUTINE calc_coll_temperature



  SUBROUTINE setup_collisions

    ALLOCATE(coll_pairs(n_species, n_species))
    coll_pairs = 1.0_num
    coll_sort_array_size = 1
    ALLOCATE(coll_sort_array(coll_sort_array_size))

  END SUBROUTINE setup_collisions



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

    ALLOCATE(part1)
    ALLOCATE(part2)

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

      en1_before = c**2 * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_before = c**2 * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error + p1 + p2
      p_sqr = p_sqr + SUM((p1 + p2)**2)
      en_error = en_error + en1_before + en2_before
      en_sqr = en_sqr + (en1_before + en2_before)**2

      CALL scatter(part1, part2, mass1, mass2, charge1, charge2, wt1, wt2, &
          density, density, 1.0e4_num, 1.0e4_num, 10.0_num, t_factor)

      p1 = part1%part_p
      p2 = part2%part_p

      en1_after = c**2 * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_after = c**2 * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error - p1 - p2
      en_error = en_error - en1_after - en2_after
    ENDDO

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

      en1_before = c**2 * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_before = c**2 * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error + p1 + p2
      p_sqr = p_sqr + SUM((p1 + p2)**2)
      en_error = en_error + en1_before + en2_before
      en_sqr = en_sqr + (en1_before + en2_before)**2

      ! alternate particle1 and particle2
      IF (MOD(i, 2) .EQ. 0) THEN
        CALL scatter(part1, part2, mass1, mass2, charge1, charge2, wt1, wt2, &
            density, density, 1.0e4_num, 1.0e4_num, 5.0_num, t_factor)
      ELSE
        CALL scatter(part2, part1, mass2, mass1, charge2, charge1, wt2, wt1, &
            density, density, 1.0e4_num, 1.0e4_num, 5.0_num, t_factor)
      ENDIF

      p1 = part1%part_p
      p2 = part2%part_p

      en1_after = c**2 * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_after = c**2 * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error - p1 - p2
      en_error = en_error - en1_after - en2_after
    ENDDO

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


      en1_before = c**2 * wt1 * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_before = c**2 * wt2 * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error + wt1 * p1 + wt2 * p2
      p_sqr = p_sqr + SUM((wt1 * p1 + wt2 * p2)**2)
      en_error = en_error + en1_before + en2_before
      en_sqr = en_sqr + (en1_before + en2_before)**2

      CALL scatter(part1, part2, mass1, mass2, charge1, charge2, wt1, wt2, &
            density, density, 1.0e4_num, 1.0e4_num, 5.0_num, t_factor)

      p1 = part1%part_p
      p2 = part2%part_p

      en1_after = c**2 * wt1 * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_after = c**2 * wt2 * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error - wt1 * p1 - wt2 * p2
      en_error = en_error - en1_after - en2_after
    ENDDO

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

      en1_before = c**2 * wt1 * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_before = c**2 * wt2 * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error + wt1 * p1 + wt2 * p2
      p_sqr = p_sqr + SUM((wt1 * p1 + wt2 * p2)**2)
      en_error = en_error + en1_before + en2_before
      en_sqr = en_sqr + (en1_before + en2_before)**2

      ! alternate particle1 and particle2
      IF (MOD(i, 2) .EQ. 0) THEN
        CALL scatter(part1, part2, mass1, mass2, charge1, charge2, wt1, wt2, &
            density, density, 1.0e4_num, 1.0e4_num, 5.0_num, t_factor)
      ELSE
        CALL scatter(part2, part1, mass2, mass1, charge2, charge1, wt2, wt1, &
            density, density, 1.0e4_num, 1.0e4_num, 5.0_num, t_factor)
      ENDIF

      p1 = part1%part_p
      p2 = part2%part_p

      en1_after = c**2 * wt1 * SQRT(DOT_PRODUCT(p1, p1) + (mass1 * c)**2)
      en2_after = c**2 * wt2 * SQRT(DOT_PRODUCT(p2, p2) + (mass2 * c)**2)

      p_error = p_error - wt1 * p1 - wt2 * p2
      en_error = en_error - en1_after - en2_after
    ENDDO

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

    INTEGER :: N, max_num, i, j
    INTEGER, DIMENSION(:), ALLOCATABLE :: histo1, histo2
    INTEGER :: histo1max, histo2max
    INTEGER :: cnt1, cnt2
    INTEGER :: a, b
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
      ENDDO

      part => partlist2%head
      DO WHILE (ASSOCIATED(part))
        part%part_p(1) = 5 * mass2 * c * random()
        part%part_p(2) = 5 * mass2 * c * random()
        part%part_p(3) = 5 * mass2 * c * random()
        part%coll_count = 0

        part%weight = random() + 1e-10_num

        part => part%next
      ENDDO

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
        IF (part%coll_count .GT. histo1max) histo1max = part%coll_count
        histo1(part%coll_count) = histo1(part%coll_count) + 1
        part => part%next
      ENDDO

!      WRITE(*,'(''  Particle List 2: '',I10)') partlist2%count

      part => partlist2%head
      DO WHILE (ASSOCIATED(part))
!        WRITE(*,'(''    '',I4,'': '',I10)') j, part%coll_count
        IF (part%coll_count .GT. histo2max) histo2max = part%coll_count
        histo2(part%coll_count) = histo2(part%coll_count) + 1
        part => part%next
      ENDDO

! only need to output if something is wrong
!      WRITE(*,*) '  Histogram 1'
!      DO j = 1, histo1max
!        WRITE(*,'(''    '',I4,'': '',I10)') j, histo1(j)
!      ENDDO
!
!      WRITE(*,*) '  Histogram 2'
!      DO j = 1, histo2max
!        WRITE(*,'(''    '',I4,'': '',I10)') j, histo2(j)
!      ENDDO

      ! performing check on histograms

      ! first, subtract expected values from the array
      IF ((cnt1 .GT. 0) .AND. (cnt2 .GT. 0)) THEN
        IF (cnt1 .GT. cnt2) THEN
          histo1(2) = histo1(2) - cnt1
          a = cnt1 / cnt2
          b = cnt1 - cnt2 * a
          histo2(2*a) = histo2(2*a) - (cnt2 - b)
          histo2(2*a+2) = histo2(2*a+2) - b
        ELSE IF (cnt1 .LT. cnt2) THEN
          histo2(2) = histo2(2) - cnt2
          a = cnt2 / cnt1
          b = cnt2 - cnt1 * a
          histo1(2*a) = histo1(2*a) - (cnt1 - b)
          histo1(2*a+2) = histo1(2*a+2) - b
        ELSE
          histo1(2) = histo1(2) - cnt1
          histo2(2) = histo2(2) - cnt1
        ENDIF
      ENDIF

      ! now check that both arrays are zero
      error = 0
      DO j = 1, histo1max
        IF (histo1(j) .NE. 0) error = error + 1
      ENDDO
      DO j = 1, histo2max
        IF (histo2(j) .NE. 0) error = error + 1
      ENDDO

      IF (error .GT. 0) THEN
        WRITE(*,*) '  Error in inter species collisions'
        WRITE(*,'(''    Iteration:   '',I10)') i
        WRITE(*,'(''    List counts: '',I10,'', '',I10)') cnt1, cnt2
        STOP
      ENDIF

      ! free particle list
      CALL destroy_partlist(partlist1)
      CALL destroy_partlist(partlist2)
    ENDDO

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

    INTEGER :: N, max_num, i, j
    INTEGER, DIMENSION(:), ALLOCATABLE :: histo
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
      ENDDO

      ! call scattering routine
      CALL intra_species_collisions(partlist, mass, charge, &
          1.0_num, 1.0e15_num, 1.0e4_num, 5.0_num, 1.0_num)

      ! diagnostics
      histo = 0
      histo_max = 0
      plist_length = plist_length + partlist%count

!      WRITE(*,'(''  Particle List: '',I10)') partlist%count

      part => partlist%head
      DO j = 1, partlist%count
        ! WRITE(*,'(''    '',I4,'': '',I10)') j, part%coll_count
        IF (part%coll_count .GT. histo_max) histo_max = part%coll_count
        histo(part%coll_count) = histo(part%coll_count) + 1
        part => part%next
      ENDDO

! only need to output if something is wrong
!      WRITE(*,*) '  Histogram'
!      DO j = 1, histo_max
!        WRITE(*,'(''    '',I4,'': '',I10)') j, histo(j)
!      ENDDO

      ! performing check on histograms

      ! first, subtract expected value from the array
      histo(2) = histo(2) - partlist%count

      ! now check that array is zero
      error = 0
      DO j = 1, histo_max
        IF (histo(j) .NE. 0) error = error + 1
      ENDDO

      IF (error .GT. 0) THEN
        WRITE(*,*) '  Error in intra species collisions'
        WRITE(*,'(''    Iteration:   '',I10)') i
        WRITE(*,'(''    List count: '',I10)') partlist%count
        STOP
      ENDIF

      ! free particle list
      CALL destroy_partlist(partlist)
    ENDDO

    WRITE(*,*) '  SUCCESS!'
    WRITE(*,'(''    Number of iterations:   '',I10)') N
    WRITE(*,'(''    Average list length:   '',F10.5)') plist_length / N

    DEALLOCATE(histo)

  END SUBROUTINE test_intra_species



  SUBROUTINE scatter_count(particle1, particle2, full)

    TYPE(particle), INTENT(INOUT) :: particle1, particle2
    LOGICAL, INTENT(IN) :: full
    INTEGER :: coll_num

    IF (full) THEN
      coll_num = 2
    ELSE
      coll_num = 1
    ENDIF

    particle1%coll_count = particle1%coll_count + coll_num
    particle2%coll_count = particle2%coll_count + coll_num

  END SUBROUTINE scatter_count



  SUBROUTINE test_shuffle

    TYPE(particle_list) :: partlist
    TYPE(particle), POINTER :: part
    INTEGER :: N, max_num, min_num, i, j, k
    INTEGER(i8) :: plist_length
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
      plist_length = (max_num - min_num) * random() + min_num
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
          part%coll_count = j
          part => part%next
        ENDDO

        ! now shuffle
        CALL shuffle_particle_list_random(partlist)

        ! perform statistics
        part => partlist%head
        DO j = 1, plist_length
          histo(j) = histo(j) + part%coll_count
          if (minp(j) .GT. part%coll_count) minp(j) = part%coll_count
          if (maxp(j) .LT. part%coll_count) maxp(j) = part%coll_count
          std_dev(j) = std_dev(j) + part%coll_count**2
          part => part%next
        ENDDO

        CALL destroy_partlist(partlist)
      ENDDO

      WRITE(*,'(''   Statistics ('',I10,'' runs)'')') N
      WRITE(*,*) '    avg        std_dev       min        max'
      DO i = 1, plist_length
        WRITE(*,'(''    '',F10.5,''    '',F10.5,''    '',I10,''    '',I10)') &
            histo(i) / N, SQRT(std_dev(i) / N - (histo(i) / N)**2), &
            minp(i), maxp(i)
      ENDDO

    ENDDO

    DEALLOCATE(histo)

  END SUBROUTINE test_shuffle



  SUBROUTINE check_particle_data

    INTEGER :: ispecies
    INTEGER :: ipart
    REAL(num) :: part_x
    REAL(num) :: part_px, part_py, part_pz
    LOGICAL :: haveNaN
    TYPE(particle), POINTER :: current

    haveNaN = .FALSE.

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO ipart = 1, species_list(ispecies)%attached_list%count
        part_x  = current%part_pos - x_grid_min_local
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)

        IF (part_x .NE. part_x) THEN
          WRITE (*,*) 'WARNING x = NaN on node ', rank, ipart
          haveNaN = .TRUE.
        ENDIF
        IF (part_px .NE. part_px) THEN
          WRITE (*,*) 'WARNING px = NaN on node ', rank, ipart
          haveNaN = .TRUE.
        ENDIF
        IF (part_py .NE. part_py) THEN
          WRITE (*,*) 'WARNING py = NaN on node ', rank, ipart
          haveNaN = .TRUE.
        ENDIF
        IF (part_pz .NE. part_pz) THEN
          WRITE (*,*) 'WARNING pz = NaN on node ', rank, ipart
          haveNaN = .TRUE.
        ENDIF

        IF (haveNaN) THEN
          STOP
        ENDIF

        current => current%next
      ENDDO
    ENDDO

  END SUBROUTINE check_particle_data

#endif

END MODULE collisions
