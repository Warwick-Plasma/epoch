MODULE deck_ic_species_block

  USE shared_data
  USE strings_advanced
  USE shared_parser_data
  USE strings
  USE shunt
  USE evaluator

  IMPLICIT NONE

  SAVE

  INTEGER :: SpeciesLoaded

CONTAINS

  FUNCTION HandleICSpeciesDeck(Species_ID,Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER,INTENT(IN) :: Species_ID
    CHARACTER(30) :: Part1
    INTEGER :: Part2
    INTEGER :: HandleICSpeciesDeck
    INTEGER :: loop,elementselected,partswitch
    LOGICAL :: Handled,Temp
    REAL(num) :: conv_val
    TYPE(primitivestack) :: output
    INTEGER :: ix,iy,ERR

    HandleICSpeciesDeck=ERR_NONE
    IF (Element .EQ. blank .OR. Value .EQ. blank) RETURN

    IF (Species_ID .LT. 0 .OR. Species_ID .GT. nspecies) THEN
       IF (rank .EQ. 0) PRINT *,"Attempting to set non-existent species initial conditions. Ignoring."
       RETURN
    ENDIF

    IF (StrCmp(Element,"minrho")) THEN
       InitialConditions(Species_ID)%MinRho=AsReal(Value,HandleICSpeciesDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"maxrho")) THEN
       InitialConditions(Species_ID)%MaxRho=AsReal(Value,HandleICSpeciesDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"rho") .OR. StrCmp(Element,"number_density")) THEN
       CALL EvaluateStringInSpace(Value,InitialConditions(Species_ID)%rho(-2:nx+3,-2:ny+3,-2:nz+3),(/-2,nx+3/),(/-2,ny+3/),(/-2,nz+3/),HandleICSpeciesDeck)
       RETURN
    ENDIF

	IF (StrCmp(Element,"mass_density")) THEN
		CALL EvaluateStringInSpace(Value,InitialConditions(Species_ID)%rho(-2:nx+3,-2:ny+3,-2:nz+3),(/-2,nx+3/),(/-2,ny+3/),(/-2,nz+3/),HandleICSpeciesDeck)
		InitialConditions(Species_ID)%rho=InitialConditions(Species_ID)%rho/ParticleSpecies(Species_ID)%mass
		RETURN
	ENDIF
	
	IF (StrCmp(Element,"drift_x")) THEN
		InitialConditions(Species_ID)%drift(1)=AsReal(Value,HandleICSpeciesDeck)
		RETURN
	ENDIF

	IF (StrCmp(Element,"drift_y")) THEN
		InitialConditions(Species_ID)%drift(2)=AsReal(Value,HandleICSpeciesDeck)
		RETURN
	ENDIF

	IF (StrCmp(Element,"drift_z")) THEN
		InitialConditions(Species_ID)%drift(3)=AsReal(Value,HandleICSpeciesDeck)
		RETURN
	ENDIF

    IF (StrCmp(Element,"temp")) THEN
       CALL EvaluateStringInSpace(Value,InitialConditions(Species_ID)%temp(-1:nx+2,-1:ny+2,-1:nz+2,1),(/-1,nx+2/),(/-1,ny+2/),(/-1,nz+2/),HandleICSpeciesDeck)
       InitialConditions(Species_ID)%temp(-1:nx+2,-1:ny+2,-1:nz+2,2)=InitialConditions(Species_ID)%temp(-1:nx+2,-1:ny+2,-1:nz+2,1)
       InitialConditions(Species_ID)%temp(-1:nx+2,-1:ny+2,-1:nz+2,3)=InitialConditions(Species_ID)%temp(-1:nx+2,-1:ny+2,-1:nz+2,1)
       RETURN
    ENDIF


    IF (StrCmp(Element,"temp_x")) THEN
       CALL EvaluateStringInSpace(Value,InitialConditions(Species_ID)%temp(-1:nx+2,-1:ny+2,-1:nz+2,1),(/-1,nx+2/),(/-1,ny+2/),(/-1,nz+2/),HandleICSpeciesDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"temp_y")) THEN
       CALL EvaluateStringInSpace(Value,InitialConditions(Species_ID)%temp(-1:nx+2,-1:ny+2,-1:nz+2,2),(/-1,nx+2/),(/-1,ny+2/),(/-1,nz+2/),HandleICSpeciesDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"temp_z")) THEN
       CALL EvaluateStringInSpace(Value,InitialConditions(Species_ID)%temp(-1:nx+2,-1:ny+2,-1:nz+2,3),(/-1,nx+2/),(/-1,ny+2/),(/-1,nz+2/),HandleICSpeciesDeck)
       RETURN
    ENDIF

  END FUNCTION HandleICSpeciesDeck

  FUNCTION CheckICSpeciesBlock()

    INTEGER :: CheckICSpeciesBlock

    !Should do error checking but can't be bothered at the moment
    CheckICSpeciesBlock=ERR_NONE

  END FUNCTION CheckICSpeciesBlock

END MODULE deck_ic_species_block
