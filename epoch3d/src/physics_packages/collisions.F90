! Relativistic binary collision module
! written by M. G. Ramsay and H. Schmitz
! based on algorithm by Sentoku & Kemp [J Comput Phys, 227, 6846 (2008)]

MODULE collisions

  USE random_generator
  USE boundary
  USE calc_df

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: particle_collisions, setup_collisions

  REAL(num) :: collision_count, large_angle_collision

  REAL(num) :: nu_avg
  INTEGER :: nu_count

  REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: meanx, meany, meanz, part_count

CONTAINS

  SUBROUTINE particle_collisions

    INTEGER :: ispecies, jspecies
    INTEGER(i8) :: ix, iy, iz
    TYPE(particle_list), POINTER :: p_list1
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: idens, jdens
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: itemp, jtemp, log_lambda
    REAL(num) :: user_factor, q1, q2, m1, m2, w1, w2

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          DO ispecies = 1, n_species
            p_list1 => species_list(ispecies)%secondary_list(ix,iy,iz)
            CALL shuffle_particle_list_random(p_list1)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ALLOCATE(idens(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(jdens(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(itemp(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(jtemp(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(log_lambda(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(meanx(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(meany(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(meanz(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(part_count(-2:nx+3,-2:ny+3,-2:nz+3))

    DO ispecies = 1, n_species
      ! Currently no support for photon collisions so just cycle round
      IF (species_list(ispecies)%species_type .EQ. c_species_id_photon) &
          CYCLE
      CALL calc_coll_number_density(idens, ispecies)
      CALL calc_coll_temperature(itemp, ispecies)

      m1 = species_list(ispecies)%mass
      q1 = species_list(ispecies)%charge
      w1 = species_list(ispecies)%weight
      itemp = itemp * kb / q0

      DO jspecies = ispecies, n_species
        ! Currently no support for photon collisions so just cycle round
        IF (species_list(jspecies)%species_type .EQ. c_species_id_photon) &
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

        DO iz = 1, nz
        DO iy = 1, ny
        DO ix = 1, nx
          IF (ispecies .EQ. jspecies) THEN
            CALL intra_species_collisions( &
                species_list(ispecies)%secondary_list(ix,iy,iz), &
                m1, q1, w1, idens(ix,iy,iz), itemp(ix,iy,iz), &
                log_lambda(ix,iy,iz), user_factor)
          ELSE
            CALL inter_species_collisions( &
                species_list(ispecies)%secondary_list(ix,iy,iz), &
                species_list(jspecies)%secondary_list(ix,iy,iz), &
                m1, m2, q1, q2, w1, w2, idens(ix,iy,iz), jdens(ix,iy,iz), &
                itemp(ix,iy,iz), jtemp(ix,iy,iz), log_lambda(ix,iy,iz), &
                user_factor)
          ENDIF
        ENDDO ! ix
        ENDDO ! iy
        ENDDO ! iz
      ENDDO ! jspecies
    ENDDO ! ispecies

    DEALLOCATE(idens, jdens, itemp, jtemp, log_lambda)
    DEALLOCATE(meanx, meany, meanz, part_count)

  END SUBROUTINE particle_collisions



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
    tvar = p1_vc * (gc - 1.0_num) / vc_sq
    p3 = p1 + vc * (tvar - gc * e1 / c**2)
    tvar = p2_vc * (gc - 1.0_num) / vc_sq
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
    cos_theta = 1.0_num - (2.0_num * delta**2) / (1.0_num + delta**2)

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
    REAL(num) :: gr, mu, ek, slow, fast
    REAL(num) :: coll_freq

    ! Manheimer-like collision operator
    ! Valid for e-i and e-e collisions
    gr = 1.0_num / SQRT(1.0_num - (vrabs / c)**2)
    mu = m2 / 1.6726d-27
    ek = (gr - 1.0_num) * m1 * c**2 / q0
    slow = 0.23_num * (mu / jtemp)**1.5_num
    fast = 3.9d-6 / ek**1.5_num
    coll_freq = slow / (1.0_num + slow / fast)
    IF (jtemp .LE. 0.0_num) coll_freq = fast
    IF (ek .LE. 0.0_num) coll_freq = 0.0_num
    IF (coll_freq .GT. 0.0_num) THEN
      coll_freq = coll_freq * jdens * log_lambda * (q2 / q0)**2 / 1.0d6
    ELSE
      coll_freq = 0.0_num
    ENDIF

    ! Velocity-dependent collision operator
    !mu = (m1 * m2) / (m1 + m2)
    !coll_freq = ((q1 * q2)**2 * jdens * log_lambda) &
    !    / (4.0_num * pi * (epsilon0 * mu)**2 * vrabs**3)
    !coll_freq = ((q1 * q2)**2 * jdens * log_lambda) &
    !    / (3.0_num * epsilon0**2 * SQRT(mu) &
    !    * (2.0_num * pi * q0 * itemp)**1.5_num)

    !coll_freq = MAX(coll_freq, vrabs / (jdens**(1.0_num / 3.0_num)))

  END FUNCTION



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

    IF (vtrans .NE. 0.0_num) THEN
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



  ! Swaps to entries in coll_sort_array
  ! Hopefully the routine will be inlined by the compiler
  SUBROUTINE swap_coll_sort_elements(i, j)

    TYPE(particle_sort_element) :: swap
    INTEGER :: i, j

    swap = coll_sort_array(i)
    coll_sort_array(i) = coll_sort_array(j)
    coll_sort_array(j) = swap

  END SUBROUTINE swap_coll_sort_elements



  ! Attaches random numbers to each particle and performs a QuickSort
  ! using the coll_sort_array.
  ! Particle numbers per cell aren't too large, so this should not result
  ! in a serious memory overhead
  SUBROUTINE shuffle_particle_list_random(p_list)

    TYPE(particle_list), INTENT(INOUT) :: p_list
    TYPE(particle), POINTER :: particle1, patricle2

    INTEGER :: i
    INTEGER :: p_num

    ! stack for holding block information
    ! large enough to sort 2^100 elements
    INTEGER, DIMENSION(100) :: sblock_start, sblock_end
    INTEGER :: sblocks

    INTEGER :: b_start, b_end ! abbreviations

    INTEGER :: pivot, store_index
    REAL(num) :: pivot_val

    p_num = INT(p_list%count)

    ! Nothing to be done
    IF (p_num .LE. 2) RETURN

    ! First make sure that the sorting array is large enough
    ! This should not happen
    IF (p_num .GT. coll_sort_array_size) THEN
      IF (ASSOCIATED(coll_sort_array)) DEALLOCATE(coll_sort_array)

      ! make the sort array somewhat larger to avoid frequent deallocation
      ! and reallocation
      coll_sort_array_size = (11 * INT(p_list%count)) / 10 + 10
      ALLOCATE(coll_sort_array(coll_sort_array_size))
    ENDIF

    ! Copy all the particle pointers into the array and create random
    ! sort indices
    i = 1
    particle1 => p_list%head
    DO WHILE (ASSOCIATED(particle1))
      coll_sort_array(i)%particle => particle1
      coll_sort_array(i)%sort_index = random()
      particle1 => particle1%next
      i = i + 1
    ENDDO

    sblock_start(1) = 1
    sblock_end(1) = p_num
    sblocks = 1

    DO WHILE(sblocks .GT. 0)
      b_start = sblock_start(sblocks)
      b_end = sblock_end(sblocks)
      sblocks = sblocks - 1

      pivot = (b_start + b_end) / 2
      pivot_val = coll_sort_array(pivot)%sort_index

      CALL swap_coll_sort_elements(pivot, b_end) ! Move pivot to end
      store_index = b_start
      DO i = b_start, b_end - 1
        IF (coll_sort_array(i)%sort_index .LE. pivot_val) THEN
           CALL swap_coll_sort_elements(i, store_index)
           store_index = store_index + 1
        ENDIF
      ENDDO
      ! Move pivot to its final place
      CALL swap_coll_sort_elements(store_index, b_end)

      ! Create two sub-blocks if they contain at least two elements
      ! Remember, the pivot element is in its right place
      IF (b_start .LT. store_index - 1) THEN
        sblocks = sblocks + 1
        sblock_start(sblocks) = b_start
        sblock_end(sblocks) = store_index - 1
      ENDIF
      IF (store_index + 1 .LT. b_end) THEN
        sblocks = sblocks + 1
        sblock_start(sblocks) = store_index + 1
        sblock_end(sblocks) = b_end
      ENDIF

    ENDDO

    ! Finally we have to copy back to the list
    ! Do head first
    particle1 => coll_sort_array(1)%particle
    p_list%head => particle1
    NULLIFY(particle1%prev)

    ! Then do all the particle links between head and tail.
    DO i = 2, p_num
      patricle2 => particle1
      particle1 => coll_sort_array(i)%particle

      particle1%prev => patricle2
      patricle2%next => particle1
    ENDDO

    ! Finally set the tail (at the end of the loop, particle is pointing to
    ! the tail)
    p_list%tail => particle1
    NULLIFY(particle1%next)

  END SUBROUTINE shuffle_particle_list_random



  PURE FUNCTION calc_coulomb_log(temp, dens, q1, q2)

    REAL(num), DIMENSION(-2:nx+3,-2:ny+3,-2:nz+3), INTENT(IN) :: temp, dens
    REAL(num), INTENT(IN) :: q1, q2
    REAL(num), DIMENSION(-2:nx+3,-2:ny+3,-2:nz+3) :: calc_coulomb_log

    calc_coulomb_log = LOG(4.13d6 * q0**3 / ABS(q1**2 * q2) &
        * (temp / kb * q0)**1.5_num / SQRT(dens))
    WHERE (calc_coulomb_log .LT. 1.0_num) calc_coulomb_log = 1.0_num
    WHERE (dens .EQ. 0.0_num) calc_coulomb_log = 1.0_num

  END FUNCTION



  SUBROUTINE calc_coll_number_density(data_array, ispecies)

    ! This subroutine calculates the grid-based number density of a given
    ! particle species.
    ! It is almost identical to the calc_number_density subroutine in calc_df,
    ! except it uses the secondary_list rather than the attached_list.

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: ispecies
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: idx
    INTEGER :: ix, iy, iz
    INTEGER :: jx, jy, jz
    TYPE(particle), POINTER :: current
#include "particle_head.inc"

    data_array = 0.0_num

    idx   = 1.0_num / dx / dy / dz

#ifndef PER_PARTICLE_WEIGHT
    wdata = species_list(ispecies)%weight * idx
#endif
    DO jz = 1, nz
    DO jy = 1, ny
    DO jx = 1, nx
      current => species_list(ispecies)%secondary_list(jx,jy,jz)%head
      DO WHILE (ASSOCIATED(current))
#ifdef PER_PARTICLE_WEIGHT
        wdata = current%weight * idx
#endif

#include "particle_to_grid.inc"

        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata
            ENDDO
          ENDDO
        ENDDO

        current => current%next
      ENDDO
    ENDDO
    ENDDO
    ENDDO

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

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: sigma
    INTEGER, INTENT(IN) :: ispecies
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_pmx, part_pmy, part_pmz, sqrt_part_m
    ! The weight of a particle
    REAL(num) :: l_weight
    REAL(num) :: gf
    INTEGER :: ix, iy, iz
    INTEGER :: jx, jy, jz
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
    DO jz = 1, nz
    DO jy = 1, ny
    DO jx = 1, nx
      current => species_list(ispecies)%secondary_list(jx,jy,jz)%head
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

        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * gz(iz) * l_weight
              meanx(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  meanx(cell_x+ix, cell_y+iy, cell_z+iz) + gf * part_pmx
              meany(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  meany(cell_x+ix, cell_y+iy, cell_z+iz) + gf * part_pmy
              meanz(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  meanz(cell_x+ix, cell_y+iy, cell_z+iz) + gf * part_pmz
              part_count(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  part_count(cell_x+ix, cell_y+iy, cell_z+iz) + gf
            ENDDO
          ENDDO
        ENDDO
        current => current%next
      ENDDO
    ENDDO
    ENDDO
    ENDDO

    CALL calc_boundary(meanx)
    CALL calc_boundary(meany)
    CALL calc_boundary(meanz)
    CALL calc_boundary(part_count)

    part_count = MAX(part_count, 1.e-6_num)

    meanx = meanx / part_count
    meany = meany / part_count
    meanz = meanz / part_count

    part_count = 0.0_num
    DO jz = 1, nz
    DO jy = 1, ny
    DO jx = 1, nx
      current => species_list(ispecies)%secondary_list(jx,jy,jz)%head
      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * gz(iz)
              sigma(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  sigma(cell_x+ix, cell_y+iy, cell_z+iz) + gf &
                  * ((part_pmx - meanx(cell_x+ix, cell_y+iy, cell_z+iz))**2 &
                  + (part_pmy - meany(cell_x+ix, cell_y+iy, cell_z+iz))**2 &
                  + (part_pmz - meanz(cell_x+ix, cell_y+iy, cell_z+iz))**2)
              part_count(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  part_count(cell_x+ix, cell_y+iy, cell_z+iz) + gf
            ENDDO
          ENDDO
        ENDDO
        current => current%next
      ENDDO
    ENDDO
    ENDDO
    ENDDO

    CALL calc_boundary(sigma)
    CALL calc_boundary(part_count)

    ! N/2 kT = <p^2>/(2m), where N is the number of degrees of freedom
    sigma = sigma / MAX(part_count, 1.e-6_num) / kb / REAL(dof)

  END SUBROUTINE calc_coll_temperature



  SUBROUTINE setup_collisions

    ALLOCATE(coll_pairs(1:n_species, 1:n_species))
    coll_pairs = 1.0_num

  END SUBROUTINE setup_collisions



#ifdef COLLISIONS_TEST
  SUBROUTINE test_collisions

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
    INTEGER :: histo_max
    INTEGER :: error

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
    REAL(num) :: part_x, part_y, part_z
    REAL(num) :: part_px, part_py, part_pz
    LOGICAL :: haveNaN
    TYPE(particle), POINTER :: current

    haveNaN = .FALSE.

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO ipart = 1, species_list(ispecies)%attached_list%count
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local
        part_z  = current%part_pos(3) - z_min_local
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)

        IF (part_x .NE. part_x) THEN
          WRITE (*,*) 'WARNING x = NaN on node ', rank, ipart
          haveNaN = .TRUE.
        ENDIF
        IF (part_y .NE. part_y) THEN
          WRITE (*,*) 'WARNING y = NaN on node ', rank, ipart
          haveNaN = .TRUE.
        ENDIF
        IF (part_z .NE. part_z) THEN
          WRITE (*,*) 'WARNING z = NaN on node ', rank, ipart
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
