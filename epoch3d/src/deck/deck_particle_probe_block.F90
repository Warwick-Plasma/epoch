MODULE deck_particle_probe_block

  USE probes
  USE strings_advanced

#ifndef PARTICLE_PROBES
CONTAINS

  SUBROUTINE probe_deck_dummy

  END SUBROUTINE probe_deck_dummy

#else
  SAVE
  TYPE(particle_probe), POINTER :: working_probe
  REAL(num) :: point2(c_ndims), point3(c_ndims)
  LOGICAL :: got_point, got_normal, got_x

CONTAINS

  SUBROUTINE probe_deck_initialise

  END SUBROUTINE probe_deck_initialise



  SUBROUTINE probe_deck_finalise

  END SUBROUTINE probe_deck_finalise



  SUBROUTINE probe_block_start

    IF (deck_state .EQ. c_ds_first) RETURN

    ALLOCATE(working_probe)
    CALL init_probe(working_probe)
    got_point = .FALSE.
    got_normal = .FALSE.
    got_x = .FALSE.

  END SUBROUTINE probe_block_start



  SUBROUTINE probe_block_end

    LOGICAL :: discard
    REAL(num), DIMENSION(c_ndims) :: r1, r2

    IF (deck_state .EQ. c_ds_first) RETURN

    discard = .FALSE.
    IF (got_point) THEN
      IF (rank .EQ. 0) THEN
        IF (got_x) THEN
          DO io = stdout,du,du-stdout ! Print to stdout and to file
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Both "x1", etc. and "point" were used in probe block.'
            WRITE(io,*) 'Only "point" and "normal" will be used.'
          ENDDO
        ENDIF
      ENDIF
      IF (.NOT. got_normal) discard = .TRUE.
    ELSE
      IF (.NOT.got_x) THEN
        discard = .TRUE.
      ELSE
        ! Old style configuration supplied. Need to calculate the normal.
        ! The probe calculates the signed distance from a point to a plane
        ! using Hessian normal form
        r1 = point2 - working_probe%point
        r2 = point3 - working_probe%point
        ! r1 (cross) r2
        working_probe%normal = (/r1(2)*r2(3) - r1(3)*r2(2), &
            r1(3)*r2(1) - r1(1)*r2(3), r1(1)*r2(2) - r1(2)*r2(1)/)

        IF (SUM(ABS(working_probe%normal)) .EQ. 0) discard = .TRUE.
      ENDIF
    ENDIF

    IF (discard) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout,du,du-stdout ! Print to stdout and to file
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Position not fully specified for distribution ', &
            'function. It will be discarded.'
        ENDDO
      ENDIF

      DEALLOCATE(working_probe)
      NULLIFY(working_probe)
    ELSE
      ! Normalise the normal. Not really necessary but doesn't hurt.
      working_probe%normal = &
          working_probe%normal / SQRT(SUM(working_probe%normal**2))

      CALL attach_probe(working_probe)
    ENDIF

  END SUBROUTINE probe_block_end



  FUNCTION probe_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode, ispecies, io

    errcode = c_err_none
    IF (deck_state .EQ. c_ds_first) RETURN
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    ! get particle probe diagnostics (rolling total of all particles which
    ! pass through a given region of real space (defined by a point on a plane
    ! and the normal to that plane.
    IF (str_cmp(element, "dump")) THEN
      working_probe%dump = as_integer(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "point")) THEN
      got_point = .TRUE.
      CALL get_vector(value, working_probe%point, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "normal")) THEN
      got_normal = .TRUE.
      CALL get_vector(value, working_probe%normal, errcode)
      RETURN
    ENDIF

    ! Top left
    IF (str_cmp(element, "x_tl")) THEN
      got_x = .TRUE.
      working_probe%point(1) = as_real(value, errcode)
      RETURN
    ENDIF
    IF (str_cmp(element, "y_tl")) THEN
      got_x = .TRUE.
      working_probe%point(2) = as_real(value, errcode)
      RETURN
    ENDIF
    IF (str_cmp(element, "z_tl")) THEN
      got_x = .TRUE.
      working_probe%point(3) = as_real(value, errcode)
      RETURN
    ENDIF

    ! Bottom right
    IF (str_cmp(element, "x_br")) THEN
      got_x = .TRUE.
      point2(1) = as_real(value, errcode)
      RETURN
    ENDIF
    IF (str_cmp(element, "y_br")) THEN
      got_x = .TRUE.
      point2(2) = as_real(value, errcode)
      RETURN
    ENDIF
    IF (str_cmp(element, "z_br")) THEN
      got_x = .TRUE.
      point2(3) = as_real(value, errcode)
      RETURN
    ENDIF

    ! Top right
    IF (str_cmp(element, "x_tr")) THEN
      got_x = .TRUE.
      point3(1) = as_real(value, errcode)
      RETURN
    ENDIF
    IF (str_cmp(element, "y_tr")) THEN
      got_x = .TRUE.
      point3(2) = as_real(value, errcode)
      RETURN
    ENDIF
    IF (str_cmp(element, "z_tr")) THEN
      got_x = .TRUE.
      point3(3) = as_real(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "probe_species")) THEN
      ispecies = as_integer(value, errcode)
      IF (errcode .EQ. c_err_none) THEN
        IF (ispecies .GT. 0 .AND. ispecies .LE. n_species) THEN
          working_probe%probe_species=>species_list(ispecies)
        ELSE
          IF (rank .EQ. 0) THEN
            DO io = stdout, du, du - stdout ! Print to stdout and to file
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Unable to attach probe to non existant species ', &
                  ispecies
            ENDDO
          ENDIF
          errcode = c_err_bad_value
        ENDIF
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "ek_min")) THEN
      working_probe%ek_min = as_real(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "ek_max")) THEN
      working_probe%ek_max = as_real(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "name")) THEN
      working_probe%name = TRIM(value)
      RETURN
    ENDIF

    errcode = c_err_unknown_element

  END FUNCTION probe_block_handle_element



  FUNCTION probe_block_check() RESULT(errcode)

    INTEGER :: errcode
    errcode = c_err_none

  END FUNCTION probe_block_check
#endif

END MODULE deck_particle_probe_block
