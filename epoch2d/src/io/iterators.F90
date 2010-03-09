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
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN) :: direction
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family
    INTEGER(8) :: part_count

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
          data(part_count) = cur%part_pos(direction)-window_shift(direction)
          cur=>cur%next
        ENDDO
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        ! If this species isn't to be dumped and this isn't a restart dump then
        IF ((.NOT. current_family%dump) &
            .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_particles



  ! iterator for particle charge
  SUBROUTINE iterate_charge(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGEMASS
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
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_charge



#ifdef PER_PARTICLE_WEIGHT

  ! iterator for particle weight
  ! Only present if you are using the PER_PARTICLE_WEIGHT
  ! Precompiler option
  SUBROUTINE iterate_weight(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
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
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF ((.NOT. current_family%dump) &
            .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_weight
#endif



  ! iterator for particle mass
  SUBROUTINE iterate_mass(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGEMASS
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
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_mass



#ifdef PART_DEBUG
  ! iterator for particle processor
  SUBROUTINE iterate_processor(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
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
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_processor



  SUBROUTINE iterate_processor0(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
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
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_processor0
#endif



  ! iterator for particle processor
  SUBROUTINE iterate_species(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
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
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF ((.NOT. current_family%dump) &
            .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_species



  ! iterator for particle velocities
  SUBROUTINE iterate_vx(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count
    REAL(num) :: root, part_m

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGEMASS
          part_m = cur%mass
#else
          part_m = current_family%mass
#endif
#ifdef NEWTONIAN
          root = part_m
#else
          root = SQRT(part_m**2 &
              + (cur%part_p(1)**2 + cur%part_p(2)**2 + cur%part_p(3)**2)/c**2)
#endif
          IF (root .NE. 0.0_num) root = 1.0_num/root
          data(part_count) = cur%part_p(1) * root
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_vx



  SUBROUTINE iterate_vy(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count
    REAL(num) :: root, part_m

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGEMASS
          part_m = cur%mass
#else
          part_m = current_family%mass
#endif
#ifdef NEWTONIAN
          root = part_m
#else
          root = SQRT(part_m**2 &
              + (cur%part_p(1)**2 + cur%part_p(2)**2 + cur%part_p(3)**2)/c**2)
#endif
          IF (root .NE. 0.0_num) root = 1.0_num/root
          data(part_count) = cur%part_p(2) * root
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_vy



  SUBROUTINE iterate_vz(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count
    REAL(num) :: root, part_m
    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
      DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. n_points))
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
          part_count = part_count+1
#ifdef PER_PARTICLE_CHARGEMASS
          part_m = cur%mass
#else
          part_m = current_family%mass
#endif
#ifdef NEWTONIAN
          root = part_m
#else
          root = SQRT(part_m**2 &
              + (cur%part_p(1)**2 + cur%part_p(2)**2 + cur%part_p(3)**2)/c**2)
#endif
          IF (root .NE. 0.0_num) root = 1.0_num/root
          data(part_count) = cur%part_p(3) * root
          cur=>cur%next
        ENDDO
        ! If the current partlist is exhausted, switch to the next one
        IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
      ENDDO

      ! If the current particle_family is exhausted, then switch to the next one
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF (.NOT. current_family%dump) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_vz



  ! iterator for particle momenta
  SUBROUTINE iterate_px(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
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
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF ((.NOT. current_family%dump) &
            .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_px



  SUBROUTINE iterate_py(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
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
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF ((.NOT. current_family%dump) &
            .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_py



  SUBROUTINE iterate_pz(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    TYPE(particle_list), POINTER, SAVE :: current_list
    TYPE(particle_family), POINTER, SAVE :: current_family

    IF (start)  THEN
      CALL start_particle_family(current_family, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_family) .AND. (part_count .LT. n_points))
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
      DO WHILE (.NOT. ASSOCIATED(cur))
        CALL advance_particle_family(current_family, current_list, cur)
        IF (.NOT. ASSOCIATED(current_family)) EXIT
        IF ((.NOT. current_family%dump) &
            .AND. (.NOT. iterator_settings%restart)) NULLIFY(cur)
      ENDDO
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_pz

END MODULE iterators
