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

CONTAINS

  SUBROUTINE probe_block_start

    ALLOCATE(working_probe)
    CALL init_probe(working_probe)

  END SUBROUTINE probe_block_start



  SUBROUTINE probe_block_end

    ! Check whether or not the probe is valid
    IF (working_probe%vertex_bottom(2) .EQ. working_probe%vertex_top(2)) THEN
      IF (rank .EQ. 0) &
          PRINT*, "Probe y1 and y2 must be different. probe ", &
              TRIM(working_probe%name), " abandoned."
      DEALLOCATE(working_probe)
      NULLIFY(working_probe)
    ELSE
      CALL attach_probe(working_probe)
    ENDIF

  END SUBROUTINE probe_block_end



  FUNCTION handle_probe_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_probe_deck, ispecies

    handle_probe_deck = c_err_none

    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    ! get particle probe diagnostics (rolling total of all particles which
    ! pass through a given region of real space (defined by the line between
    ! two points in 2D).
    IF (str_cmp(element, "dump")) THEN
      working_probe%dump = as_integer(value, handle_probe_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "x1")) THEN
      working_probe%vertex_bottom(1) = as_real(value, handle_probe_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "y1")) THEN
      working_probe%vertex_bottom(2) = as_real(value, handle_probe_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "x2")) THEN
      working_probe%vertex_top(1) = as_real(value, handle_probe_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "y2")) THEN
      working_probe%vertex_top(2) = as_real(value, handle_probe_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "probe_species")) THEN
      ispecies = as_integer(value, handle_probe_deck)
      IF (handle_probe_deck .EQ. c_err_none) THEN
        IF (ispecies .GT. 0 .AND. ispecies .LE. n_species) THEN
          working_probe%probe_species=>particle_species(ispecies)
        ELSE
          IF (rank .EQ. 0) &
              PRINT *, "Unable to attach probe to non existant species ", &
                  ispecies
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

  END FUNCTION handle_probe_deck
#endif

END MODULE deck_eio_particle_probe_block
