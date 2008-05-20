MODULE deck_species_block

  USE shared_data
  USE strings

  IMPLICIT NONE

  SAVE

  INTEGER :: SpeciesLoaded

CONTAINS

  FUNCTION HandleSpeciesDeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    CHARACTER(30) :: Part1
    INTEGER :: Part2
    INTEGER :: HandleSpeciesDeck
    INTEGER :: loop,elementselected,partswitch
    LOGICAL :: Handled
    REAL(num) :: conv_val

    HandleSpeciesDeck=ERR_NONE 

    IF (ICHAR(Element(1:1)) == 0) RETURN

    HandleSpeciesDeck=ERR_UNKNOWN_ELEMENT

    Handled=.FALSE.
    IF (StrCmp(Element,"n_species")) THEN
       nspecies=AsInteger(Value,HandleSpeciesDeck)
       IF (nspecies .GT. 0) THEN
          ALLOCATE(species(1:nspecies,1:3))
          ALLOCATE(Species_Name(1:nspecies))
          IF (rank .EQ. 0) PRINT '("Running code with ",i2," species of particle")',nspecies
       ENDIF
       HandleSpeciesDeck=ERR_NONE
       Handled=.TRUE.
    ENDIF

    IF (.NOT. ALLOCATED(species)) THEN
       IF (rank .EQ. 0) THEN
          PRINT *,"Either invalid number of species specified or attempting to set species data before setting n_species"
       ENDIF
       HandleSpeciesDeck=ERR_BAD_VALUE
       RETURN
    ENDIF
    IF (Handled) RETURN

    CALL SplitOffInt(Element,Part1,Part2,HandleSpeciesDeck)

    IF (Part2 .LT. 1 .OR. Part2 .GT. nspecies) THEN
       IF (rank .EQ. 0) THEN
          PRINT *,"Attempting to set non_existant particle species, ",part2,", ignoring"
       ENDIF
       HandleSpeciesDeck=ERR_NONE
       RETURN
    ENDIF

    IF (StrCmp(Part1,"name")) THEN
       Species_Name(Part2)%Value = TRIM(Value)
       IF (rank .EQ. 0) PRINT '("Name of species ",i2," is ",a)',Part2,TRIM(Value)
       HandleSpeciesDeck=ERR_NONE
       RETURN
    ENDIF
    IF (StrCmp(Part1,"charge")) THEN
       HandleSpeciesDeck=ERR_NONE
       species(Part2,1)=AsReal(Value,HandleSpeciesDeck) * Q0
    ENDIF
    IF (StrCmp(Part1,"mass")) THEN
       HandleSpeciesDeck=ERR_NONE
       species(Part2,2)=AsReal(Value,HandleSpeciesDeck) * M0
    ENDIF
    IF (StrCmp(Part1,"fraction")) THEN
       HandleSpeciesDeck=ERR_NONE
       species(Part2,3)=AsReal(Value,HandleSpeciesDeck)
    ENDIF


  END FUNCTION HandleSpeciesDeck

  FUNCTION CheckSpeciesBlock

    INTEGER :: CheckSpeciesBlock

    !Should do error checking but can't be bothered at the moment
    CheckSpeciesBlock=ERR_NONE

  END FUNCTION CheckSpeciesBlock

END MODULE deck_species_block
