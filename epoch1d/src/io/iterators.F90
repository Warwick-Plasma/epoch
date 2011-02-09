MODULE iterators

  USE particle_pointer_advance

  IMPLICIT NONE

  TYPE :: setting_block
    LOGICAL :: restart
  END TYPE setting_block

  SAVE

  TYPE(setting_block) :: iterator_settings

CONTAINS

  ! iterator for particle positions
  SUBROUTINE iterate_particles(data, n_points, direction, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    INTEGER, INTENT(IN) :: direction
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      IF ((.NOT. current_family%dump) &
          .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)

      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = cur%part_pos - window_shift
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      IF (.NOT. ASSOCIATED(cur)) &
          CALL advance_particle_family(current_family, current_list, cur)
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_particles



  ! iterator for particle species
  SUBROUTINE iterate_species(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      IF ((.NOT. current_family%dump) &
          .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)

      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = REAL(current_family%id, num)
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      IF (.NOT. ASSOCIATED(cur)) &
          CALL advance_particle_family(current_family, current_list, cur)
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_species



#ifdef PER_PARTICLE_WEIGHT
  ! iterator for particle weight
  ! Only present if you are using the PER_PARTICLE_WEIGHT
  ! Precompiler option
  SUBROUTINE iterate_weight(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      IF ((.NOT. current_family%dump) &
          .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)

      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = cur%weight
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      IF (.NOT. ASSOCIATED(cur)) &
          CALL advance_particle_family(current_family, current_list, cur)
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_weight
#endif



  ! iterator for particle momenta
  SUBROUTINE iterate_px(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      IF ((.NOT. current_family%dump) &
          .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)

      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = cur%part_p(1)
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      IF (.NOT. ASSOCIATED(cur)) &
          CALL advance_particle_family(current_family, current_list, cur)
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_px



  SUBROUTINE iterate_py(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      IF ((.NOT. current_family%dump) &
          .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)

      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = cur%part_p(2)
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      IF (.NOT. ASSOCIATED(cur)) &
          CALL advance_particle_family(current_family, current_list, cur)
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_py



  SUBROUTINE iterate_pz(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      IF ((.NOT. current_family%dump) &
          .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)

      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = cur%part_p(3)
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      IF (.NOT. ASSOCIATED(cur)) &
          CALL advance_particle_family(current_family, current_list, cur)
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_pz



  ! iterator for particle velocities
  SUBROUTINE iterate_vx(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count
    REAL(num) :: part_mc2, gamma_mass

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      IF ((.NOT. current_family%dump) &
          .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)

#ifndef PER_PARTICLE_CHARGE_MASS
      part_mc2 = (current_family%mass * c)**2
#endif
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGE_MASS
          part_mc2 = (cur%mass * c)**2
#endif
          gamma_mass = SQRT(SUM(cur%part_p**2) + part_mc2) / c
          data(part_count) = cur%part_p(1) / gamma_mass
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      IF (.NOT. ASSOCIATED(cur)) &
          CALL advance_particle_family(current_family, current_list, cur)
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_vx



  SUBROUTINE iterate_vy(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count
    REAL(num) :: part_mc2, gamma_mass

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      IF ((.NOT. current_family%dump) &
          .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)

#ifndef PER_PARTICLE_CHARGE_MASS
      part_mc2 = (current_family%mass * c)**2
#endif
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGE_MASS
          part_mc2 = (cur%mass * c)**2
#endif
          gamma_mass = SQRT(SUM(cur%part_p**2) + part_mc2) / c
          data(part_count) = cur%part_p(2) / gamma_mass
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      IF (.NOT. ASSOCIATED(cur)) &
          CALL advance_particle_family(current_family, current_list, cur)
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_vy



  SUBROUTINE iterate_vz(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count
    REAL(num) :: part_mc2, gamma_mass

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      IF ((.NOT. current_family%dump) &
          .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)

#ifndef PER_PARTICLE_CHARGE_MASS
      part_mc2 = (current_family%mass * c)**2
#endif
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGE_MASS
          part_mc2 = (cur%mass * c)**2
#endif
          gamma_mass = SQRT(SUM(cur%part_p**2) + part_mc2) / c
          data(part_count) = cur%part_p(3) / gamma_mass
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      IF (.NOT. ASSOCIATED(cur)) &
          CALL advance_particle_family(current_family, current_list, cur)
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_vz



  ! iterator for particle charge
  SUBROUTINE iterate_charge(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      IF ((.NOT. current_family%dump) &
          .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)

      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGE_MASS
          data(part_count) = cur%charge
#else
          data(part_count) = current_family%charge
#endif
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      IF (.NOT. ASSOCIATED(cur)) &
          CALL advance_particle_family(current_family, current_list, cur)
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_charge



  ! iterator for particle mass
  SUBROUTINE iterate_mass(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      IF ((.NOT. current_family%dump) &
          .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)

      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGE_MASS
          data(part_count) = cur%mass
#else
          data(part_count) = current_family%mass
#endif
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      IF (.NOT. ASSOCIATED(cur)) &
          CALL advance_particle_family(current_family, current_list, cur)
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_mass



#ifdef PARTICLE_DEBUG
  ! iterator for particle processor
  SUBROUTINE iterate_processor(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      IF ((.NOT. current_family%dump) &
          .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)

      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = REAL(cur%processor, num)
          IF (cur%processor .GE. nproc) PRINT *, "Bad Processor"
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      IF (.NOT. ASSOCIATED(cur)) &
          CALL advance_particle_family(current_family, current_list, cur)
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_processor



  SUBROUTINE iterate_processor0(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      IF ((.NOT. current_family%dump) &
          .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)

      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = REAL(cur%processor_at_t0, num)
          IF (cur%processor .GE. nproc) PRINT *, "Bad Processor"
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      IF (.NOT. ASSOCIATED(cur)) &
          CALL advance_particle_family(current_family, current_list, cur)
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_processor0
#endif

END MODULE iterators
