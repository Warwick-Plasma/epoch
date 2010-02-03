MODULE deck_species_block

  USE shared_data
  USE strings
  USE strings_advanced

  IMPLICIT NONE

  SAVE

  INTEGER :: SpeciesLoaded

CONTAINS



  FUNCTION HandleSpeciesDeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    CHARACTER(30) :: Part1
    INTEGER :: Part2
    INTEGER :: HandleSpeciesDeck
    INTEGER :: partswitch
    LOGICAL :: Handled
    CHARACTER(LEN=9) :: string

    HandleSpeciesDeck=ERR_NONE 
    IF (Value .EQ. blank) RETURN

    !    PRINT *,"IN ",Element,Value

    HandleSpeciesDeck=ERR_UNKNOWN_ELEMENT

    Handled=.FALSE.
    IF (StrCmp(Element,"n_species")) THEN
       nspecies=AsInteger(Value,HandleSpeciesDeck)
       IF (nspecies .GT. 0) THEN
          IF (Rank .EQ. 0) THEN
             CALL IntegerAsString(nspecies,string)
             PRINT *,"Code running with ",TRIM(ADJUSTL(string))," species"
          ENDIF
          ALLOCATE(ParticleSpecies(1:nspecies))
          !ALLOCATE(Species_Name(1:nspecies))
       ENDIF
       HandleSpeciesDeck=ERR_NONE
       Handled=.TRUE.
    ENDIF

    IF (nspecies .LE. 0) THEN
       IF (rank .EQ. 0) THEN
          PRINT *,"Either invalid number of species specified or attempting to set species data before setting n_species"
       ENDIF
       HandleSpeciesDeck=ERR_MISSING_ELEMENTS
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
       IF (rank .EQ. 0) THEN
          CALL IntegerAsString(part2,string)
          PRINT *,"Name of species ",TRIM(ADJUSTL(string))," is ",TRIM(Value)
       ENDIF
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
    IF (StrCmp(Part1,"split")) THEN
       HandleSpeciesDeck=ERR_NONE
#ifdef SPLIT_PARTICLES_AFTER_PUSH
       ParticleSpecies(Part2)%Split=AsLogical(Value,HandleSpeciesDeck)
#endif
       RETURN
    ENDIF
    IF (StrCmp(Part1,"npart_max")) THEN
       HandleSpeciesDeck=ERR_NONE
#ifdef SPLIT_PARTICLES_AFTER_PUSH
       ParticleSpecies(Part2)%nPart_Max=AsLongInteger(Value,HandleSpeciesDeck)
#endif
       RETURN
    ENDIF



  END FUNCTION HandleSpeciesDeck

  FUNCTION CheckSpeciesBlock()

    INTEGER :: CheckSpeciesBlock

    !Should do error checking but can't be bothered at the moment
    CheckSpeciesBlock=ERR_NONE

  END FUNCTION CheckSpeciesBlock

END MODULE deck_species_block
