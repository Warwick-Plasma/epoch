MODULE deck_eio_particle_probe_block

  USE probes
  USE strings_advanced

#ifndef PARTICLE_PROBES
CONTAINS

  SUBROUTINE probe_deck_dummy

  END SUBROUTINE probe_deck_dummy

#else
  SAVE
  TYPE(particle_probe), POINTER :: working_probe
  LOGICAL :: got_point, got_normal

CONTAINS

  SUBROUTINE probe_block_start

    ALLOCATE(working_probe)
    CALL init_probe(working_probe)
    got_point = .FALSE.
    got_normal = .FALSE.

  END SUBROUTINE probe_block_start



  SUBROUTINE probe_block_end

    LOGICAL :: discard

    discard = .FALSE.
    IF (.NOT. got_point) discard = .TRUE.
    IF (.NOT. got_normal) discard = .TRUE.

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
      working_probe%normal = SIGN(1.0_num, working_probe%normal)

      CALL attach_probe(working_probe)
    ENDIF

  END SUBROUTINE probe_block_end



  FUNCTION handle_probe_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_probe_deck, ispecies, io

    handle_probe_deck = c_err_none

    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    ! get particle probe diagnostics (rolling total of all particles which
    ! pass through a given region of real space (defined by a point on a plane
    ! and the normal to that plane.
    IF (str_cmp(element, "dump")) THEN
      working_probe%dump = as_integer(value, handle_probe_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "point") .OR. str_cmp(element, "probe_point")) THEN
      got_point = .TRUE.
      working_probe%point = as_real(value, handle_probe_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "normal")) THEN
      got_normal = .TRUE.
      working_probe%normal = as_real(value, handle_probe_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "left_to_right")) THEN
      got_normal = as_logical(value, handle_probe_deck)
      IF (got_normal) THEN
        working_probe%normal = 1.0_num
      ELSE
        working_probe%normal = -1.0_num
      ENDIF
      got_normal = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, "probe_species")) THEN
      ispecies = as_integer(value, handle_probe_deck)
      IF (handle_probe_deck .EQ. c_err_none) THEN
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
          handle_probe_deck = c_err_bad_value
        ENDIF
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "ek_min")) THEN
      working_probe%ek_min = as_real(value, handle_probe_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "ek_max")) THEN
      working_probe%ek_max = as_real(value, handle_probe_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "name")) THEN
      working_probe%name = TRIM(value)
      RETURN
    ENDIF

    handle_probe_deck = c_err_unknown_element

  END FUNCTION handle_probe_deck
#endif

END MODULE deck_eio_particle_probe_block
