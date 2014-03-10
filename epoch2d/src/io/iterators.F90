MODULE iterators

  USE particle_pointer_advance
#if PARTICLE_ID || PARTICLE_ID4
  USE partlist
#endif

  IMPLICIT NONE

  SAVE

  TYPE(particle_species), POINTER :: current_species

CONTAINS

  ! iterator for particle positions
  FUNCTION it_output_position(array, npoint_it, start, direction, param)

    REAL(num) :: it_output_position
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: npoint_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN) :: direction
    INTEGER, INTENT(IN), OPTIONAL :: param
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. npoint_it))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
        part_count = part_count + 1
        array(part_count) = cur%part_pos(direction) - window_shift(direction)
        cur => cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    npoint_it = part_count

    it_output_position = 0

  END FUNCTION it_output_position



  FUNCTION it_output_real(array, npoint_it, start, param)

    REAL(num) :: it_output_real
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: npoint_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count, ndim
    REAL(num) :: part_m, part_mc, part_mcc, part_mc2, gamma_mass, csqr, charge

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
    ENDIF

    part_count = 0
    csqr = c**2
    charge = current_species%charge
    part_m   = current_species%mass
    part_mc  = part_m * c
    part_mc2 = part_mc**2
    part_mcc = part_m * c**2

    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. npoint_it))
      SELECT CASE (param)
#ifdef PER_PARTICLE_WEIGHT
      CASE (c_dump_part_weight) ! particle weight
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
          part_count = part_count + 1
          array(part_count) = cur%weight
          cur => cur%next
        ENDDO
#endif

      CASE (c_dump_part_px)
        ndim = 1
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
          part_count = part_count + 1
          array(part_count) = cur%part_p(ndim)
          cur => cur%next
        ENDDO

      CASE (c_dump_part_py)
        ndim = 2
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
          part_count = part_count + 1
          array(part_count) = cur%part_p(ndim)
          cur => cur%next
        ENDDO

      CASE (c_dump_part_pz)
        ndim = 3
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
          part_count = part_count + 1
          array(part_count) = cur%part_p(ndim)
          cur => cur%next
        ENDDO

      CASE (c_dump_part_vx)
        ndim = 1
#ifdef PHOTONS
        IF (current_species%species_type .NE. c_species_id_photon) THEN
#endif
          DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
            part_count = part_count + 1
#ifdef PER_PARTICLE_CHARGE_MASS
            part_mc2 = (cur%mass * c)**2
#endif
            gamma_mass = SQRT(SUM(cur%part_p**2) + part_mc2) / c
            array(part_count) = cur%part_p(ndim) / gamma_mass
            cur => cur%next
          ENDDO
#ifdef PHOTONS
        ELSE
          DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
            part_count = part_count + 1
            array(part_count) = cur%part_p(ndim) * csqr / cur%particle_energy
            cur => cur%next
          ENDDO
        ENDIF
#endif

      CASE (c_dump_part_vy)
        ndim = 2
#ifdef PHOTONS
        IF (current_species%species_type .NE. c_species_id_photon) THEN
#endif
          DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
            part_count = part_count + 1
#ifdef PER_PARTICLE_CHARGE_MASS
            part_mc2 = (cur%mass * c)**2
#endif
            gamma_mass = SQRT(SUM(cur%part_p**2) + part_mc2) / c
            array(part_count) = cur%part_p(ndim) / gamma_mass
            cur => cur%next
          ENDDO
#ifdef PHOTONS
        ELSE
          DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
            part_count = part_count + 1
            array(part_count) = cur%part_p(ndim) * csqr / cur%particle_energy
            cur => cur%next
          ENDDO
        ENDIF
#endif

      CASE (c_dump_part_vz)
        ndim = 3
#ifdef PHOTONS
        IF (current_species%species_type .NE. c_species_id_photon) THEN
#endif
          DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
            part_count = part_count + 1
#ifdef PER_PARTICLE_CHARGE_MASS
            part_mc2 = (cur%mass * c)**2
#endif
            gamma_mass = SQRT(SUM(cur%part_p**2) + part_mc2) / c
            array(part_count) = cur%part_p(ndim) / gamma_mass
            cur => cur%next
          ENDDO
#ifdef PHOTONS
        ELSE
          DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
            part_count = part_count + 1
            array(part_count) = cur%part_p(ndim) * csqr / cur%particle_energy
            cur => cur%next
          ENDDO
        ENDIF
#endif

      CASE (c_dump_part_charge)
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
          part_count = part_count + 1
#ifdef PER_PARTICLE_CHARGE_MASS
          array(part_count) = cur%charge
#else
          array(part_count) = charge
#endif
          cur => cur%next
        ENDDO

      CASE (c_dump_part_mass)
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
          part_count = part_count + 1
#ifdef PER_PARTICLE_CHARGE_MASS
          array(part_count) = cur%mass
#else
          array(part_count) = part_m
#endif
          cur => cur%next
        ENDDO

      CASE (c_dump_part_ek)
        IF (current_species%species_type .NE. c_species_id_photon) THEN
          DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
            part_count = part_count + 1
#ifdef PER_PARTICLE_CHARGE_MASS
            part_m   = cur%mass
            part_mcc = part_m * c**2
            part_mc2 = (part_m * c)**2
#endif
            array(part_count) = &
                c * SQRT(SUM(cur%part_p**2) + part_mc2) - part_mcc
            cur => cur%next
          ENDDO
#ifdef PHOTONS
        ELSE
          DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
            part_count = part_count + 1
            array(part_count) = cur%particle_energy
            cur => cur%next
          ENDDO
#endif
        ENDIF

      CASE (c_dump_part_gamma)
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
          part_count = part_count + 1
#ifdef PER_PARTICLE_CHARGE_MASS
          part_mc = cur%mass * c
#endif
          array(part_count) = SQRT(SUM((cur%part_p/part_mc)**2) + 1.0_num)
          cur => cur%next
        ENDDO

      CASE (c_dump_part_rel_mass)
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
          part_count = part_count + 1
#ifdef PER_PARTICLE_CHARGE_MASS
          part_m  = cur%mass
          part_mc = part_m * c
#endif
          array(part_count) = &
              part_m * SQRT(SUM((cur%part_p/part_mc)**2) + 1.0_num)
          cur => cur%next
        ENDDO

#ifdef PHOTONS
      CASE (c_dump_part_opdepth)
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
          part_count = part_count + 1
          array(part_count) = cur%optical_depth
          cur => cur%next
        ENDDO

      CASE (c_dump_part_qed_energy)
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
          part_count = part_count + 1
          array(part_count) = cur%particle_energy
          cur => cur%next
        ENDDO

#ifdef TRIDENT_PHOTONS
      CASE (c_dump_part_opdepth_tri)
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
          part_count = part_count + 1
          array(part_count) = cur%optical_depth_tri
          cur => cur%next
        ENDDO
#endif
#endif
      END SELECT
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) CALL advance_particle_list(current_list, cur)
    ENDDO
    npoint_it = part_count

    it_output_real = 0

  END FUNCTION it_output_real



#if defined(PARTICLE_ID4) || defined(PARTICLE_DEBUG)
  FUNCTION it_output_integer4(array, npoint_it, start, param)

    INTEGER(i4) :: it_output_integer4
    INTEGER(i4), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: npoint_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
#if defined(PARTICLE_ID4)
      IF (param .EQ. c_dump_part_id) THEN
        CALL generate_particle_ids(current_list)
      ENDIF
#endif
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. npoint_it))
      SELECT CASE (param)
#ifdef PARTICLE_DEBUG
      CASE (c_dump_part_proc)
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
          part_count = part_count + 1
          array(part_count) = cur%processor
          IF (cur%processor .GE. nproc) PRINT *, 'Bad Processor'
          cur => cur%next
        ENDDO

      CASE (c_dump_part_proc0)
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
          part_count = part_count + 1
          array(part_count) = cur%processor_at_t0
          IF (cur%processor .GE. nproc) PRINT *, 'Bad Processor'
          cur => cur%next
        ENDDO
#endif
#if defined(PARTICLE_ID4)
      CASE (c_dump_part_id)
        DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
          part_count = part_count + 1
          array(part_count) = cur%id
          cur => cur%next
        ENDDO
#endif
      END SELECT
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) THEN
        CALL advance_particle_list(current_list, cur)
#if defined(PARTICLE_ID4)
        IF (param .EQ. c_dump_part_id .AND. ASSOCIATED(current_list)) THEN
          CALL generate_particle_ids(current_list)
        ENDIF
#endif
      ENDIF
    ENDDO
    npoint_it = part_count

    it_output_integer4 = 0

  END FUNCTION it_output_integer4
#endif



#if defined(PARTICLE_ID)
  FUNCTION it_output_integer8(array, npoint_it, start, param)

    INTEGER(i8) :: it_output_integer8
    INTEGER(i8), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: npoint_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    TYPE(particle), POINTER, SAVE :: cur
    TYPE(particle_list), POINTER, SAVE :: current_list
    INTEGER :: part_count

    IF (start)  THEN
      CALL start_particle_list(current_species, current_list, cur)
      CALL generate_particle_ids(current_list)
    ENDIF

    part_count = 0
    DO WHILE (ASSOCIATED(current_list) .AND. (part_count .LT. npoint_it))
      DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. npoint_it))
        part_count = part_count + 1
        array(part_count) = cur%id
        cur => cur%next
      ENDDO
      ! If the current partlist is exhausted, switch to the next one
      IF (.NOT. ASSOCIATED(cur)) THEN
        CALL advance_particle_list(current_list, cur)
        IF (ASSOCIATED(current_list)) THEN
          CALL generate_particle_ids(current_list)
        ENDIF
      ENDIF
    ENDDO
    npoint_it = part_count

    it_output_integer8 = 0

  END FUNCTION it_output_integer8
#endif

END MODULE iterators
