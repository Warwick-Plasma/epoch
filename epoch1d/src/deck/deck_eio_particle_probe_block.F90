MODULE deck_eio_particle_probe_block

  USE probes
  USE shared_data
  USE strings_advanced
#ifndef PARTICLE_PROBES
CONTAINS

  SUBROUTINE probe_deck_dummy
  END SUBROUTINE probe_deck_dummy
#else

  SAVE
  TYPE(particle_probe),POINTER :: working_probe

CONTAINS

  SUBROUTINE probe_block_start

    ALLOCATE(working_probe)
    CALL init_probe(working_probe)

  END SUBROUTINE probe_block_start

  SUBROUTINE probe_block_end

    CALL attach_probe(working_probe)

  END SUBROUTINE probe_block_end

  FUNCTION handle_probe_deck(element,value)
    CHARACTER(*),INTENT(IN) :: element,value
    INTEGER :: handle_probe_deck, ispecies

    handle_probe_deck=c_err_none

    IF (element .EQ. blank .OR. value .EQ. blank) RETURN
    IF (str_cmp(element,"dump")) THEN
       working_probe%dump=as_integer(value,handle_probe_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"probe_point")) THEN
       working_probe%probe_point=as_real(value,handle_probe_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"left_to_right")) THEN
       working_probe%left_to_right=as_logical(value,handle_probe_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"probe_species")) THEN
       ispecies=as_integer(value,handle_probe_deck)
       IF (handle_probe_deck .EQ. c_err_none) THEN
          IF (ispecies .GT. 0 .AND. ispecies .LE. n_species) THEN
             working_probe%probe_species=>particle_species(ispecies)
          ELSE
             IF (rank .EQ. 0) PRINT *,"Unable to attach probe to non existant species ",ispecies
             handle_probe_deck=c_err_bad_value
          ENDIF
       ENDIF
       RETURN
    ENDIF

    IF (str_cmp(element,"ek_min")) THEN
       working_probe%ek_min=as_real(value,handle_probe_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"ek_max")) THEN
       working_probe%ek_max=as_real(value,handle_probe_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"name")) THEN
       working_probe%name=TRIM(value)
       RETURN
    ENDIF
  END FUNCTION handle_probe_deck

#endif

END MODULE deck_eio_particle_probe_block
