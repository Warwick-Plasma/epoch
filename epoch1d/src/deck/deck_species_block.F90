MODULE deck_species_block

  USE shared_data
  USE strings
  USE strings_advanced
  USE setup

  IMPLICIT NONE

  SAVE

  INTEGER :: SpeciesLoaded

CONTAINS



  FUNCTION HandleSpeciesDeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    CHARACTER(LEN=EntryLength) :: Part1
    INTEGER :: Part2
    INTEGER :: HandleSpeciesDeck
    INTEGER :: partswitch
    LOGICAL :: Handled

    HandleSpeciesDeck=ERR_NONE 
    IF (Value .EQ. blank) RETURN
    HandleSpeciesDeck=ERR_UNKNOWN_ELEMENT

    Handled=.FALSE.
    IF (StrCmp(Element,"n_species")) THEN
       nspecies=AsInteger(Value,HandleSpeciesDeck)
       IF (nspecies .GT. 0) THEN
          IF (Rank .EQ. 0) PRINT '("Code running with ",i2," species")',nspecies
          CALL Setup_Species
       ENDIF
       HandleSpeciesDeck=ERR_NONE
       Handled=.TRUE.
    ENDIF

    IF (nspecies .LE. 0) THEN
       IF (rank .EQ. 0) THEN
          PRINT *,"Invalid number of species specified"
       ENDIF
       HandleSpeciesDeck=ERR_BAD_VALUE
    ENDIF

    IF (.NOT. ASSOCIATED(ParticleSpecies)) THEN
       Extended_Error_String="n_species"
       IF (rank .EQ. 0) THEN
          PRINT *,"Attempting to set species data before setting n_species"
       ENDIF
       HandleSpeciesDeck=ERR_REQUIRED_ELEMENT_NOT_SET
       RETURN
    ENDIF
    IF (Handled) RETURN

    CALL SplitOffInt(Element,Part1,Part2,HandleSpeciesDeck)

    IF (Part2 .LT. 1 .OR. Part2 .GT. nspecies) THEN
       IF (Rank .EQ. 0) PRINT '("Ignoring attempt to set property ",a," for non existent species ",i2)',TRIM(ADJUSTL(Part1)),Part2
       HandleSpeciesDeck=ERR_NONE
       RETURN
    ENDIF

    partswitch=0
    IF (StrCmp(Part1,"name")) THEN
       ParticleSpecies(Part2)%Name = TRIM(Value)
       IF (rank .EQ. 0) PRINT '("Name of species ",i2," is ",a)',Part2,TRIM(Value)
       HandleSpeciesDeck=ERR_NONE
       RETURN
    ENDIF
    IF (StrCmp(Part1,"mass")) THEN
       HandleSpeciesDeck=ERR_NONE
       ParticleSpecies(Part2)%Mass=AsReal(Value,HandleSpeciesDeck)*M0
       RETURN
    ENDIF
    IF (StrCmp(Part1,"charge")) THEN
       HandleSpeciesDeck=ERR_NONE
       ParticleSpecies(Part2)%Charge=AsReal(Value,HandleSpeciesDeck)*Q0
       RETURN
    ENDIF
    IF (StrCmp(Part1,"frac") .OR. StrCmp(Part1,"fraction")) THEN
       HandleSpeciesDeck=ERR_NONE
       IF (npart_global .GE. 0) THEN
          ParticleSpecies(Part2)%Count=AsReal(Value,HandleSpeciesDeck)*npart_global
       ELSE
          Extended_Error_String="npart"
          HandleSpeciesDeck=ERR_REQUIRED_ELEMENT_NOT_SET
       ENDIF
       RETURN
    ENDIF
    IF (StrCmp(Part1,"npart")) THEN
       HandleSpeciesDeck=ERR_NONE
       ParticleSpecies(Part2)%Count=AsLongInteger(Value,HandleSpeciesDeck)
       RETURN
    ENDIF
    IF (StrCmp(Part1,"dump")) THEN
       HandleSpeciesDeck=ERR_NONE
       ParticleSpecies(Part2)%Dump=AsLogical(Value,HandleSpeciesDeck)
       RETURN
    ENDIF
    !*************************************************************
    !This section sets properties for tracer particles
    !*************************************************************
    IF(StrCmp(Part1,"tracer")) THEN
       HandleSpeciesDeck=ERR_NONE
#ifdef TRACER_PARTICLES
       ParticleSpecies(Part2)%Tracer=AsLogical(Value,HandleSpeciesDeck)
#else
       HandleSpeciesDeck=ERR_PP_OPTIONS_WRONG
       Extended_Error_String="-DTRACER_PARTICLES"
#endif
    ENDIF
    !*************************************************************
    !This section sets properties for particle splitting
    !*************************************************************
    IF (StrCmp(Part1,"split")) THEN
       HandleSpeciesDeck=ERR_NONE
#ifdef SPLIT_PARTICLES_AFTER_PUSH
       ParticleSpecies(Part2)%Split=AsLogical(Value,HandleSpeciesDeck)
#else
       HandleSpeciesDeck=ERR_PP_OPTIONS_WRONG
       Extended_Error_String="-DSPLIT_PARTICLES_AFTER_PUSH"
#endif
       RETURN
    ENDIF
    IF (StrCmp(Part1,"npart_max")) THEN
       HandleSpeciesDeck=ERR_NONE
#ifdef SPLIT_PARTICLES_AFTER_PUSH
       ParticleSpecies(Part2)%nPart_Max=AsLongInteger(Value,HandleSpeciesDeck)
#else
       HandleSpeciesDeck=ERR_PP_OPTIONS_WRONG
       Extended_Error_String="-DSPLIT_PARTICLES_AFTER_PUSH"
#endif
       RETURN
    ENDIF


    !*************************************************************
    !This section sets properties for ionisation
    !*************************************************************
    IF (StrCmp(Part1,"ionise")) THEN
       HandleSpeciesDeck=ERR_NONE
#ifdef PART_IONISE
       ParticleSpecies(Part2)%ionise=AsLogical(Value,HandleSpeciesDeck)
#else
       HandleSpeciesDeck=ERR_PP_OPTIONS_WRONG
       Extended_Error_String="-DPART_IONISE"
#endif
       RETURN
    ENDIF
    IF (StrCmp(Part1,"ionise_to_species")) THEN
       HandleSpeciesDeck=ERR_NONE
#ifdef PART_IONISE
       ParticleSpecies(Part2)%ionise_to_species=AsInteger(Value,HandleSpeciesDeck)
#else
       HandleSpeciesDeck=ERR_PP_OPTIONS_WRONG
       Extended_Error_String="-DPART_IONISE"
#endif
       RETURN
    ENDIF
    IF (StrCmp(Part1,"release_species_on_ionise")) THEN
       HandleSpeciesDeck=ERR_NONE
#ifdef PART_IONISE
       ParticleSpecies(Part2)%release_species=AsInteger(Value,HandleSpeciesDeck)
#else
       HandleSpeciesDeck=ERR_PP_OPTIONS_WRONG
       Extended_Error_String="-DPART_IONISE"
#endif
       RETURN
    ENDIF
    IF (StrCmp(Part1,"ionisation_energy")) THEN
       HandleSpeciesDeck=ERR_NONE
#ifdef PART_IONISE
       ParticleSpecies(Part2)%ionisation_energy=AsReal(Value,HandleSpeciesDeck)
#else
       HandleSpeciesDeck=ERR_PP_OPTIONS_WRONG
       Extended_Error_String="-DPART_IONISE"
#endif
       RETURN
    ENDIF




  END FUNCTION HandleSpeciesDeck

  FUNCTION CheckSpeciesBlock()

    INTEGER :: CheckSpeciesBlock

!!$    DO iSpecies=1,nSpecies
!!$       IF (StrCmp(ParticleSpecies(iSpecies)%name,blank)) THEN
!!$          CheckSpeciesBlock=ERR_MISSING_ELEMENTS
!!$       ENDIF
!!$    ENDDO

    !Should do error checking but can't be bothered at the moment
    CheckSpeciesBlock=ERR_NONE

  END FUNCTION CheckSpeciesBlock


END MODULE deck_species_block
