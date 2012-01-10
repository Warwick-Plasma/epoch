MODULE iterators

  USE particle_pointer_advance

  IMPLICIT NONE

  SAVE

  TYPE(particle_species), POINTER :: current_species

CONTAINS

  ! iterator for particle positions
  FUNCTION iterate_particles(array, n_points, start, direction)

    REAL(num) :: iterate_particles
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN) :: direction
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
        part_count = part_count + 1
        array(part_count) = cur%part_pos(direction) - window_shift(direction)
        cur=>cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    n_points = part_count

    iterate_particles = 0

  END FUNCTION iterate_particles



#ifdef PER_PARTICLE_WEIGHT
  ! iterator for particle weight
  ! Only present if you are using the PER_PARTICLE_WEIGHT
  ! Precompiler option
  FUNCTION iterate_weight(array, n_points, start)

    REAL(num) :: iterate_weight
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
        part_count = part_count + 1
        array(part_count) = cur%weight
        cur=>cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    n_points = part_count

    iterate_weight = 0

  END FUNCTION iterate_weight
#endif



  ! iterator for particle momenta
  FUNCTION iterate_px(array, n_points, start)

    REAL(num) :: iterate_px
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
        part_count = part_count + 1
        array(part_count) = cur%part_p(1)
        cur=>cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    n_points = part_count

    iterate_px = 0

  END FUNCTION iterate_px



  FUNCTION iterate_py(array, n_points, start)

    REAL(num) :: iterate_py
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
        part_count = part_count + 1
        array(part_count) = cur%part_p(2)
        cur=>cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    n_points = part_count

    iterate_py = 0

  END FUNCTION iterate_py



  FUNCTION iterate_pz(array, n_points, start)

    REAL(num) :: iterate_pz
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
        part_count = part_count + 1
        array(part_count) = cur%part_p(3)
        cur=>cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    n_points = part_count

    iterate_pz = 0

  END FUNCTION iterate_pz



  ! iterator for particle velocities
  FUNCTION iterate_vx(array, n_points, start)

    REAL(num) :: iterate_vx
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count
    REAL(num) :: part_mc2, gamma_mass

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
#ifndef PER_PARTICLE_CHARGE_MASS
    part_mc2 = (current_species%mass * c)**2
#endif
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
        part_count = part_count + 1
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc2 = (cur%mass * c)**2
#endif
        gamma_mass = SQRT(SUM(cur%part_p**2) + part_mc2) / c
        array(part_count) = cur%part_p(1) / gamma_mass
        cur=>cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    n_points = part_count

    iterate_vx = 0

  END FUNCTION iterate_vx



  FUNCTION iterate_vy(array, n_points, start)

    REAL(num) :: iterate_vy
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count
    REAL(num) :: part_mc2, gamma_mass

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
#ifndef PER_PARTICLE_CHARGE_MASS
    part_mc2 = (current_species%mass * c)**2
#endif
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
        part_count = part_count + 1
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc2 = (cur%mass * c)**2
#endif
        gamma_mass = SQRT(SUM(cur%part_p**2) + part_mc2) / c
        array(part_count) = cur%part_p(2) / gamma_mass
        cur=>cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    n_points = part_count

    iterate_vy = 0

  END FUNCTION iterate_vy



  FUNCTION iterate_vz(array, n_points, start)

    REAL(num) :: iterate_vz
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count
    REAL(num) :: part_mc2, gamma_mass

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
#ifndef PER_PARTICLE_CHARGE_MASS
    part_mc2 = (current_species%mass * c)**2
#endif
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
        part_count = part_count + 1
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc2 = (cur%mass * c)**2
#endif
        gamma_mass = SQRT(SUM(cur%part_p**2) + part_mc2) / c
        array(part_count) = cur%part_p(3) / gamma_mass
        cur=>cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    n_points = part_count

    iterate_vz = 0

  END FUNCTION iterate_vz



  ! iterator for particle charge
  FUNCTION iterate_charge(array, n_points, start)

    REAL(num) :: iterate_charge
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
        part_count = part_count + 1
#ifdef PER_PARTICLE_CHARGE_MASS
        array(part_count) = cur%charge
#else
        array(part_count) = current_species%charge
#endif
        cur=>cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    n_points = part_count

    iterate_charge = 0

  END FUNCTION iterate_charge



  ! iterator for particle mass
  FUNCTION iterate_mass(array, n_points, start)

    REAL(num) :: iterate_mass
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
        part_count = part_count + 1
#ifdef PER_PARTICLE_CHARGE_MASS
        array(part_count) = cur%mass
#else
        array(part_count) = current_species%mass
#endif
        cur=>cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    n_points = part_count

    iterate_mass = 0

  END FUNCTION iterate_mass



#ifdef PARTICLE_DEBUG
  ! iterator for particle processor
  FUNCTION iterate_processor(array, n_points, start)

    REAL(num) :: iterate_processor
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
        part_count = part_count + 1
        array(part_count) = REAL(cur%processor, num)
        IF (cur%processor .GE. nproc) PRINT *, "Bad Processor"
        cur=>cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    n_points = part_count

    iterate_processor = 0

  END FUNCTION iterate_processor



  FUNCTION iterate_processor0(array, n_points, start)

    REAL(num) :: iterate_processor0
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
        part_count = part_count + 1
        array(part_count) = REAL(cur%processor_at_t0, num)
        IF (cur%processor .GE. nproc) PRINT *, "Bad Processor"
        cur=>cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    n_points = part_count

    iterate_processor0 = 0

  END FUNCTION iterate_processor0
#endif



#if PARTICLE_ID || PARTICLE_ID4
  ! iterator for particle id
  FUNCTION iterate_id(array, n_points, start)

    REAL(num) :: iterate_id
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
        part_count = part_count + 1
        array(part_count) = REAL(cur%id, num)
        cur=>cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    n_points = part_count

    iterate_id = 0

  END FUNCTION iterate_id
#endif

END MODULE iterators
