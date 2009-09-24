MODULE deck_eio_particle_probe_block

  USE probes
  USE shared_data
  USE strings_advanced
#ifndef PARTICLE_PROBES
CONTAINS

  SUBROUTINE Probe_Dummy
  END SUBROUTINE Probe_Dummy
#else

  SAVE
  TYPE(Particle_Probe),POINTER :: Working_Probe

CONTAINS

  SUBROUTINE Probe_Block_Start

    ALLOCATE(Working_Probe)
    CALL Init_Probe(Working_Probe)

  END SUBROUTINE Probe_Block_Start

  SUBROUTINE Probe_Block_End

    CALL Attach_Probe(Working_Probe)

  END SUBROUTINE Probe_Block_End

  FUNCTION HandleProbeDeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER :: HandleProbeDeck, iSpecies

    HandleProbeDeck=ERR_NONE

    IF (Element .EQ. blank .OR. Value .EQ. blank) RETURN
    IF (StrCmp(Element,"dump")) THEN
       working_probe%dump=AsInteger(Value,HandleProbeDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"probe_point")) THEN
       working_probe%probe_point=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"left_to_right")) THEN
       working_probe%LeftToRight=AsLogical(Value,HandleProbeDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"probe_species")) THEN
       iSpecies=AsInteger(Value,HandleProbeDeck)
       IF (HandleProbeDeck .EQ. ERR_NONE) THEN
          IF (iSpecies .GT. 0 .AND. iSpecies .LE. nSpecies) THEN
             working_probe%probe_species=>ParticleSpecies(iSpecies)
          ELSE
             IF (rank .EQ. 0) PRINT *,"Unable to attach probe to non existant species ",iSpecies
             HandleProbeDeck=ERR_BAD_VALUE
          ENDIF
       ENDIF
       RETURN
    ENDIF

    IF (StrCmp(Element,"ek_min")) THEN
       working_probe%ek_min=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"ek_max")) THEN
       working_probe%ek_max=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"name")) THEN
       working_probe%name=TRIM(Value)
       RETURN
    ENDIF
  END FUNCTION HandleProbeDeck

#endif

END MODULE deck_eio_particle_probe_block
