MODULE deck

  USE shared_data
  USE strings
  USE initial_conditions
  USE deck_control_block
  USE deck_boundaries_block
  USE deck_species_block
  USE deck_io_block

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Read_Deck

  SAVE
  CHARACTER(len=30) :: Current_Block_Name,Blank
  LOGICAL :: InvalidBlock
CONTAINS

  SUBROUTINE Read_Deck

    character :: u1
    integer :: pos=1,flip=1,s,f,elements=0,LUN
    LOGICAL :: IsComment
    TYPE(Entry), DIMENSION(2) :: DeckValues
    CHARACTER(len=14+Data_Dir_Max_Length)::syscmd
    LOGICAL :: Terminate =.FALSE.
    INTEGER :: errcode_deck

    errcode_deck=ERR_NONE

    LUN=5

    ControlBlockDone=.FALSE.
    BoundaryBlockDone=.FALSE.

    IsComment=.FALSE.

    IF (rank .EQ. 0) THEN
       LUN=10
       OPEN(unit=10,file="./input.deck")
       OPEN(unit=40,file="./deck.status")
       DeckValues(1)%Value=""
       DeckValues(2)%Value=""
       !Use non-advancing IO to pop characters off the deck file one at a time
       !Use basic token parsing to split into two substrings across an "=" or ":" symbol
       DO
          Errcode_Deck=ERR_NONE
          READ(LUN,"(A1)",advance='no',size=s,iostat=f),u1
          IF (u1 .EQ. '#') IsComment=.TRUE.
          IF (.NOT. IsComment .AND. u1 .NE. '#') THEN
             IF (u1 .NE. '=' .AND. u1 .NE. " " .AND. u1 .NE. char(32) .AND. u1 .NE. char(9) .AND. u1 .NE. ':' .AND. f .EQ. 0) THEN
                DeckValues(flip)%Value(pos:pos)=u1
                pos=pos+1
             ENDIF
             IF (u1 .EQ. '=' .OR. u1 .EQ. ':') THEN
                flip=2
                pos=1
             ENDIF
          ENDIF
          IF (f .EQ. -2) IsComment=.FALSE.
          IF (f .EQ. -2 .AND. pos .GT. 1) THEN
             elements=elements+1
             flip=1
             pos=1
             DeckValues(1)%value=TRIM(ADJUSTL(DeckValues(1)%value))
             DeckValues(2)%value=TRIM(ADJUSTL(DeckValues(2)%value))
             CALL MPI_BCAST(1,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
             CALL MPI_BCAST(DeckValues(1)%value,30,MPI_CHARACTER,0,MPI_COMM_WORLD,errcode)
             CALL MPI_BCAST(DeckValues(2)%value,30,MPI_CHARACTER,0,MPI_COMM_WORLD,errcode)             
             CALL HandleDeckElement(DeckValues(1)%value,DeckValues(2)%value,Errcode_deck)
             DeckValues(1)%value=""
             DeckValues(2)%value=""
             IsComment=.FALSE.
          ENDIF
          IF (f .EQ. -1) THEN
             CALL MPI_BCAST(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)   
             IF (deckfile) CLOSE(LUN)
             EXIT
          ENDIF
          Terminate=Terminate .OR. IAND(errcode_deck,ERR_TERMINATE) .NE. 0
       ENDDO
    ELSE
       DO
          Errcode_deck=ERR_NONE
          CALL MPI_BCAST(f,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
          IF (f .EQ. 0) EXIT
          CALL MPI_BCAST(DeckValues(1)%value,30,MPI_CHARACTER,0,MPI_COMM_WORLD,errcode)
          CALL MPI_BCAST(DeckValues(2)%value,30,MPI_CHARACTER,0,MPI_COMM_WORLD,errcode)             
          CALL HandleDeckElement(DeckValues(1)%value,DeckValues(2)%value,Errcode_deck)
          DeckValues(1)%value=""
          DeckValues(2)%value=""
          Terminate=Terminate .OR. IAND(errcode_deck,ERR_TERMINATE) .NE. 0
       ENDDO
    ENDIF

    CALL MPI_BARRIER(MPI_COMM_WORLD,errcode)

    !Don't check compulsory blocks if going to bomb anyway, just stinks up the output file
    IF (.NOT. Terminate) CALL CheckCompulsoryBlocks(errcode_deck)
    Terminate=Terminate .OR. IAND(errcode_deck,ERR_TERMINATE) .NE. 0
    !Fatal error, cause code to bomb
    IF (Terminate .AND. RANK .EQ. 0) THEN
       PRINT *,""
       WRITE(40,*) ""
       PRINT *,'***FATAL ERROR*** The code cannot parse the input deck sufficiently to run. Please read the output file "deck.status" in the current directory to check for errors.'
       WRITE(40,*) "***FATAL ERROR*** The code cannot parse the input deck sufficiently to run. Please read this file and correct any errors mentioned"
       PRINT *,""
       PRINT *,""
       PRINT *,""
       PRINT *,""
       PRINT *,""
       PRINT *,""
    ENDIF
    CLOSE(40)



    IF (Terminate) CALL MPI_ABORT(MPI_COMM_WORLD,errcode)

    IF (rank .EQ. 0) THEN
       !Create the data directory
       syscmd=" "
       WRITE(syscmd,"('mkdir 'a' >& /dev/null')"),TRIM(data_dir)
       CALL System(Trim(syscmd))

       !Copy the input deck into the data directory
       syscmd=" "
       WRITE(syscmd,"('cp input.deck 'a'/.')"),TRIM(data_dir)
       CALL System(TRIM(syscmd))

       !Move the deck diagnostics into the data directory
       syscmd=" "
       WRITE(syscmd,"('mv deck.status 'a'/.')"),TRIM(data_dir)
       CALL System(TRIM(syscmd))
    ENDIF


  END SUBROUTINE Read_Deck

  SUBROUTINE HandleDeckElement(Element,Value,Errcode_deck)

    CHARACTER(*),INTENT(IN) :: Element
    CHARACTER(*),INTENT(IN) :: Value
    INTEGER,INTENT(INOUT) :: Errcode_deck
    INTEGER :: Handled
    INTEGER :: State,rankcheck
    INTEGER, SAVE :: errcount

    rankcheck=0

    State=0
    IF (StrCmp(Element,"begin")) THEN
       Errcode_deck=HandleBlock(Value,Blank,Blank)
       InvalidBlock=IAND(Errcode_deck,ERR_UNKNOWN_BLOCK) .NE. 0
       IF(InvalidBlock) THEN
          IF (RANK .EQ. rankcheck) THEN
             PRINT *,char(9),"Unknown block ",TRIM(Value)," in input deck, ignoring"
          ENDIF
       ENDIF
       errcount=0
       Current_Block_Name=Value
       IF (rank .EQ. rankcheck) THEN
          WRITE(40,*),"Beginning ", TRIM(ADJUSTL(Value)), " block"
          WRITE(40,*),""
       ENDIF
       !Reset errcode_deck here because reporting ERR_UNKNOWN_ELEMENT is OK
       ERRCODE_DECK=ERR_NONE
       RETURN
    ENDIF
    IF (StrCmp(Element,"end")) THEN
       InvalidBlock=.TRUE.
       IF (rank .EQ. rankcheck) THEN
          WRITE(40,*),""
          WRITE(40,*),"Ending ",TRIM(ADJUSTL(Value)), " block"
          WRITE(40,*),""
          IF (errcount .NE. 0) THEN 
             WRITE(40,*) "***WARNING*** block ",TRIM(ADJUSTL(Value))," contains errors"
             WRITE(40,*) ""
          ENDIF
       ENDIF
       RETURN
    ENDIF

    !Check invalid block to limit amount of rubbish that appears
    !If the input deck is invalid
    IF (.NOT. InvalidBlock) THEN
       Errcode_deck=HandleBlock(Current_Block_Name,Element,Value)
    ELSE
       RETURN
    ENDIF

    IF (Errcode_deck==ERR_NONE) THEN
       IF (rank .EQ. rankcheck) WRITE(40,*),char(9),"Element ", TRIM(ADJUSTL(Element))," = ",TRIM(ADJUSTL(Value)), " handled OK"
       RETURN
    ENDIF
    !Test for error conditions
    !If an error is fatal then OR Errcode_deck with ERR_TERMINATE
    IF (IAND(Errcode_deck,ERR_UNKNOWN_ELEMENT) /= 0) THEN
       IF (rank .EQ. rankcheck) THEN
          WRITE(40,*) ""
          PRINT *,"***WARNING*** Unrecognised element ",TRIM(Element), " in input deck. Code will continue to run, but behaviour is undefined"
          WRITE(40,*) "***WARNING*** Unrecognised element ",TRIM(Element), " in input deck. Code will continue to run, but behaviour is undefined"
          PRINT *,""
          WRITE(40,*) ""
       ENDIF
    ENDIF
    IF (IAND(Errcode_deck,ERR_PRESET_ELEMENT) /= 0) THEN
       IF (rank .EQ. rankcheck) THEN
          WRITE(40,*) ""
          PRINT *,"***WARNING*** Element ",TRIM(Element), " is set multiple times in this deck. Code will continue using first value in deck"
          WRITE(40,*) "***WARNING*** Element ",TRIM(Element), " is set multiple times in this deck. Code will continue using first value in deck"
          PRINT *,""
          WRITE(40,*) ""
       ENDIF
    ENDIF
    IF (IAND(Errcode_deck, ERR_PRESET_ELEMENT_USE_LATER) /= 0) THEN
       IF (rank .EQ. rankcheck) THEN
          WRITE(40,*) ""
          PRINT *,"***WARNING*** Element ",TRIM(Element), " is set multiple times in this deck. Code will continue using last value in deck"
          WRITE(40,*) "***WARNING*** Element ",TRIM(Element), " is set multiple times in this deck. Code will continue using last value in deck"
          PRINT *,""
          WRITE(40,*) ""
       ENDIF
    ENDIF
    IF (IAND(Errcode_deck, ERR_BAD_VALUE) /= 0) THEN
       IF (rank .EQ. rankcheck) THEN
          WRITE(40,*) ""
          PRINT *,"***ERROR*** Value ",TRIM(Value)," in element ",TRIM(Element)," is invalid or could not be parsed. Code will terminate."
          WRITE(40,*) "***ERROR*** Value ",TRIM(Value)," in element ",TRIM(Element)," is invalid or could not be parsed. Code will terminate."
          PRINT *,""
          WRITE(40,*) ""
          Errcode_deck=IOR(Errcode_deck,ERR_TERMINATE)
       ENDIF
    ENDIF
    IF (IAND(Errcode_deck , ERR_OTHER) /= 0) THEN
       IF (rank .EQ. rankcheck) THEN
          WRITE(40,*) ""
          PRINT *,"***ERROR*** You have managed to find an impossible situation in this code. Good for you. Just because of that, code will terminate."
          WRITE(40,*) "***ERROR*** You have managed to find an impossible situation in this code. Good for you. Just because of that, code will terminate."
          PRINT *,""
          WRITE(40,*) ""
          Errcode_deck=IOR(Errcode_deck,ERR_TERMINATE)
       ENDIF
    ENDIF

    ErrCount=ErrCount+1

  END SUBROUTINE HandleDeckElement

  FUNCTION HandleBlock(BlockName,BlockElement,BlockValue)

    CHARACTER(len=*),INTENT(IN) :: BlockName, BlockElement,BlockValue
    INTEGER :: HandleBlock

    HandleBlock=ERR_UNKNOWN_BLOCK
    !Test for known blocks
    IF (StrCmp(BlockName,"control"))  THEN
       HandleBlock=HandleControlDeck(BlockElement,BlockValue)
       RETURN
    ENDIF
    IF (StrCmp(BlockName,"boundaries")) THEN
       HandleBlock=HandleBoundaryDeck(BlockElement,BlockValue)
       RETURN
    ENDIF
    IF (StrCmp(BlockName,"species")) THEN
       HandleBlock=HandleSpeciesDeck(BlockElement,BlockValue)
       RETURN
    ENDIF
    IF (StrCmp(BlockName,"output")) THEN
       HandleBlock=HandleIODeck(BlockElement,BlockValue)
       RETURN
    ENDIF
    !Pass through to the custom block
    HandleBlock=HandleCustomBlock(BlockName,BlockElement,BlockValue)

  END FUNCTION HandleBlock



  !These subroutines are there to check for the basic minimal compulsory blocks are present
  !They're a bit ugly, but they seem to be the easiest way to do it without adding complexity to the code
  SUBROUTINE CheckCompulsoryBlocks(errcode_deck)

    INTEGER :: index
    LOGICAL :: Problem_Found
    INTEGER,INTENT(INOUT) :: errcode_deck

    Problem_Found=.FALSE.

    errcode_deck=CheckControlBlock()

    errcode_deck=IOR(errcode_deck,CheckBoundaryBlock())
    errcode_deck=IOR(errcode_deck,CheckSpeciesBlock())
    errcode_deck=IOR(errcode_deck,CheckIOBlock())
    errcode_deck=IOR(errcode_deck,CheckCustomBlocks())

    problem_found =(IAND(errcode_deck,ERR_MISSING_ELEMENTS) .NE. 0) 



    IF (Problem_Found) THEN
       errcode_deck=IOR(errcode_deck,ERR_TERMINATE)
       IF (rank .EQ. 0) THEN
          PRINT *,""
          PRINT *,"Not all required elements of input deck specified. Please fix input deck and rerun code"
          WRITE(40,*) ""
          WRITE(40,*) "Not all required elements of input deck specified. Please fix input deck and rerun code"
       ENDIF
    ELSE
       IF (rank .EQ. 0) THEN
          PRINT *,"Input deck complete and valid. Attempting to set up equilibrium"
          PRINT *,""
          WRITE(40,*) "Input deck complete and valid."
       ENDIF
    ENDIF

  END SUBROUTINE CheckCompulsoryBlocks

END MODULE deck
