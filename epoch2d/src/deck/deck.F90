MODULE deck

  !Basic operations
  USE shared_data
  USE strings
  USE strings_advanced
  !Deck internals
  USE deck_constant_block
  USE deck_deo_block
  !Deck Blocks
  USE deck_control_block
  USE deck_boundaries_block
  USE deck_species_block
  USE deck_io_block
  USE deck_window_block
  !Initial Condition Blocks
  USE deck_ic_laser_block
  USE deck_ic_species_block
  USE deck_ic_fields_block
  USE deck_ic_external_block
  !Extended IO Blocks
  USE deck_eio_dist_fn_block
#ifdef PARTICLE_PROBES
  USE deck_eio_particle_probe_block
#endif
  !Custom blocks
  USE custom_deck

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Read_Deck

  SAVE
  CHARACTER(len=EntryLength) :: Current_Block_Name
  LOGICAL :: InvalidBlock
CONTAINS

  !---------------------------------------------------------------------------------
  !These subroutines actually call the routines which read the deck blocks
  !---------------------------------------------------------------------------------

  !This subroutine is called when a new block is started
  !If a block NEEDS to do something when it starts, then
  !The revelevant subroutine should be called here
  SUBROUTINE StartBlock(BlockName)
    CHARACTER(LEN=*),INTENT(IN) :: BlockName

    IF (StrCmp(BlockName,"laser")) CALL Laser_Start
    IF (StrCmp(BlockName,"window")) CALL Window_Start
    IF (StrCmp(BlockName,"dist_fn")) CALL Dist_Fn_Start
#ifdef PARTICLE_PROBES
    IF (StrCmp(BlockName,"probe")) CALL Probe_Block_Start
#endif
    IF (StrCmp(BlockName,"species_external")) CALL Start_External
    IF (StrCmp(BlockName,"fields_external")) CALL Start_External

  END SUBROUTINE StartBlock

  !This subroutine is called when a new block is ended
  !If a block NEEDS to do something when it ends, then
  !The revelevant subroutine should be called here
  SUBROUTINE EndBlock(BlockName)
    CHARACTER(LEN=*),INTENT(IN) :: BlockName

    IF (StrCmp(BlockName,"laser")) CALL Laser_End
    IF (StrCmp(BlockName,"dist_fn")) CALL Dist_Fn_End
#ifdef PARTICLE_PROBES
    IF (StrCmp(BlockName,"probe")) CALL Probe_Block_End
#endif


  END SUBROUTINE EndBlock


  FUNCTION HandleBlock(BlockName,BlockElement,BlockValue)

    CHARACTER(len=*),INTENT(IN) :: BlockName, BlockElement,BlockValue
    CHARACTER(len=EntryLength) :: Part1
    INTEGER :: HandleBlock
    INTEGER :: Part2,Result

    HandleBlock=ERR_UNKNOWN_BLOCK
    !Constants can be defined in any deck state, so put them here
    IF (StrCmp(BlockName,"constant")) THEN
       HandleBlock=HandleConstantDeck(BlockElement,BlockValue)
       RETURN
    ENDIF
    IF (StrCmp(BlockName,"deo")) THEN
       HandleBlock=HandleDEODeck(BlockElement,BlockValue)
       RETURN
    ENDIF
    IF (Deck_State .EQ. DS_DECK) THEN
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
       IF (StrCmp(BlockName,"window")) THEN
          HandleBlock=HandleWindowDeck(BlockElement,BlockValue)
          RETURN
       ENDIF
    ELSE IF (Deck_State .EQ. DS_IC) THEN
       !Initial conditions blocks go here
       IF (StrCmp(BlockName,"fields")) THEN
          HandleBlock=HandleICFieldsDeck(BlockElement,BlockValue)
          RETURN
       ENDIF
       IF (StrCmp(BlockName,"fields_external")) THEN
          HandleBlock=HandleICFieldsDeck(BlockElement,BlockValue)
          RETURN
       ENDIF
       IF (StrCmp(BlockName,"laser")) THEN
          HandleBlock=HandleICLaserDeck(BlockElement,BlockValue)
          RETURN
       ENDIF
       Result=ERR_NONE
       CALL SplitOffInt(BlockName,Part1,Part2,Result)
       IF (Result .EQ. ERR_NONE) THEN
          IF (StrCmp(Part1,"species")) THEN
             HandleBlock=HandleICSpeciesDeck(part2,BlockElement,BlockValue)
             RETURN
          ENDIF
          IF (StrCmp(Part1,"species_external")) THEN
             HandleBlock=HandleICExternalSpeciesDeck(part2,BlockElement,BlockValue)
             RETURN
          ENDIF
       ENDIF
    ELSE IF (Deck_State .EQ. DS_EIO) THEN
       IF (StrCmp(BlockName,"dist_fn")) THEN
          HandleBlock=HandleEIODistFnDeck(BlockElement,BlockValue)
          RETURN
       ENDIF
       IF (StrCmp(BlockName,"probe")) THEN
#ifdef PARTICLE_PROBES
          HandleBlock=HandleProbeDeck(BlockElement,BlockValue)
          RETURN
#else
          HandleBlock=ERR_PP_OPTIONS_WRONG
          Extended_Error_String="-DPARTICLE_PROBES"
          RETURN
#endif
       ENDIF
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

    errcode_deck=ERR_NONE

    IF (Deck_State .EQ. DS_DECK) THEN
       errcode_deck=IOR(errcode_deck,CheckControlBlock())
       errcode_deck=IOR(errcode_deck,CheckBoundaryBlock())
       errcode_deck=IOR(errcode_deck,CheckSpeciesBlock())
       errcode_deck=IOR(errcode_deck,CheckIOBlock())
       errcode_deck=IOR(errcode_deck,CheckWindowBlock())
       errcode_deck=IOR(errcode_deck,CheckCustomBlocks())
    ELSE IF (Deck_State .EQ. DS_IC) THEN
       errcode_deck=IOR(errcode_deck,CheckICFieldsBlock())
       errcode_deck=IOR(errcode_deck,CheckICSpeciesBlock())
    ENDIF
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
          IF (Deck_State .EQ. DS_DECK) THEN
             PRINT *,"Input deck complete and valid. Attempting to set up equilibrium"
             PRINT *,""
             WRITE(40,*) "Input deck complete and valid."
          ELSE IF (Deck_State .EQ. DS_IC) THEN
             PRINT *,"Initial conditions complete and valid. Attempting to load particles"
             PRINT *,""
             WRITE(40,*) "Initial conditions complete and valid."
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE CheckCompulsoryBlocks

  !---------------------------------------------------------------------------------
  !These subroutines are the in depth detail of how the parser works
  !---------------------------------------------------------------------------------
  FUNCTION GetFreeLUN()

    !This subroutine simply cycles round until it finds a free lun between MinLun and MaxLun
    INTEGER :: GetFreeLUN
    INTEGER :: LUN
    INTEGER, PARAMETER :: MinLun=10, MaxLun=20
    LOGICAL :: Open

    Open=.TRUE.

    LUN=MinLun
    DO
       INQUIRE(UNIT=LUN,OPENED=Open)
       IF (.NOT. Open) EXIT
       LUN=LUN+1
       IF (LUN .GT. MaxLun) THEN
          IF (Rank .EQ. 0) THEN
             PRINT *,"***FATAL ERROR*** unable to open LUN for input deck read"
          ENDIF
          CALL MPI_ABORT(MPI_COMM_WORLD,errcode)
       ENDIF
    ENDDO

    GetFreeLun=LUN

  END FUNCTION GetFreeLUN

  RECURSIVE SUBROUTINE Read_Deck(filename,FirstCall)

    CHARACTER(len=*),INTENT(IN) :: filename
    LOGICAL, INTENT(IN) :: FirstCall
    character :: u1
    integer :: pos=1,flip=1,s,f,elements=0,LUN
    LOGICAL :: IsComment
    TYPE(Entry), DIMENSION(2) :: DeckValues
    CHARACTER(len=45+Data_Dir_Max_Length)::syscmd, DeckFilename, StatusFilename
    LOGICAL :: Terminate =.FALSE., Exists
    INTEGER :: errcode_deck
    LOGICAL :: White_Space_Over

    !No error yet
    errcode_deck=ERR_NONE
    !Characteristic string which represents a "blank" string
    Blank="BLANKBLANK"

    LUN=5

    !Make the whole filename by adding the data_dir to the filename
    DeckFileName=TRIM(ADJUSTL(data_dir))// '/' // TRIM(ADJUSTL(filename))

    !Deck_State tells the code whether it's parsing the normal input deck 
    !Or the initial conditions. You can add more states if you want.
    !Just search for Deck_State
    IF (Deck_State .EQ. DS_DECK) THEN
       StatusFileName=TRIM(ADJUSTL(data_dir))// '/' // "deck.status"
    ELSE IF (Deck_State .EQ. DS_IC) THEN
       StatusFileName=TRIM(ADJUSTL(data_dir))// '/' // "ic.status"
    ELSE IF(Deck_State .EQ. DS_EIO) THEN
       StatusFileName=TRIM(ADJUSTL(data_dir))// '/' // "eio.status"
    ENDIF

    !If this is the first time that this deck has been called then do some housekeeping
    !Put any initialisation code that is needed in here
    IF (FirstCall) THEN
       ControlBlockDone=.FALSE.
       BoundaryBlockDone=.FALSE.
    ENDIF

    !Is comment is a flag which tells the code when a # character has been found and everything beyond it
    !Is a comment
    IsComment=.FALSE.

    !Rank 0 reads the file and then passes it out to the other nodes using MPI_BCAST
    IF (rank .EQ. 0) THEN
       !Check whether or not the input deck file requested exists
       INQUIRE(File=DeckFilename,Exist=Exists)
       IF (.NOT. Exists) THEN
          PRINT *,"***ERROR*** Input deck file ",DeckFilename," does not exist. Create the file and rerun the code."
          CALL MPI_ABORT(MPI_COMM_WORLD,errcode) 
       ENDIF
       !Get a free LUN. Don't use a constant LUN to allow for recursion
       LUN=GetFreeLUN()
       OPEN(unit=LUN,file=TRIM(ADJUSTL(DeckFilename)))
       IF (FirstCall) OPEN(unit=40,file=StatusFilename)
       DeckValues(1)%Value=""
       DeckValues(2)%Value=""
       !Use non-advancing IO to pop characters off the deck file one at a time
       !Use basic token parsing to split into two substrings across an "=" or ":" symbol
       DO
          Errcode_Deck=ERR_NONE
          !Read a character
          !When you reach an EOL character iostat returns -2
          !When you reach an EOF iostat returns -1
          READ(LUN,"(A1)",advance='no',size=s,iostat=f),u1
          !If the character is a # then switch to comment mode
          IF (u1 .EQ. '#') IsComment=.TRUE.
          !If not in comment mode then use the character
          IF (.NOT. IsComment) THEN
             !If the current character isn't a special character then just stick it in the buffer
             !             IF (u1 .NE. '=' .AND. u1 .NE. char(32) .AND. u1 .NE. char(9) .AND. u1 .NE. ':' .AND. f .EQ. 0) THEN
             IF (u1 .NE. '=' .AND. u1 .NE. char(9) .AND. u1 .NE. ':' .AND. f .EQ. 0) THEN
                IF ((u1 .NE. ' ' .AND. u1 .NE. CHAR(32)).OR. White_Space_Over) THEN
                   DeckValues(flip)%Value(pos:pos)=u1
                   pos=pos+1
                   White_Space_Over=.TRUE.
                ENDIF
             ENDIF
             !If it's equals or : then you're parsing the other part of the expression
             IF (u1 .EQ. '=' .OR. u1 .EQ. ':') THEN
                flip=2
                pos=1
             ENDIF
          ENDIF
          !If f=-2 then you've reached the end of the line, so comment state is definitely false
          IF (f .EQ. -2) IsComment=.FALSE.
          !If you've not read a blank line then
          IF (f .EQ. -2 .AND. pos .GT. 1) THEN
             elements=elements+1
             flip=1
             pos=1
             DeckValues(1)%value=TRIM(ADJUSTL(DeckValues(1)%value))
             DeckValues(2)%value=TRIM(ADJUSTL(DeckValues(2)%value))
             CALL MPI_BCAST(1,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
             CALL MPI_BCAST(DeckValues(1)%value,EntryLength,MPI_CHARACTER,0,MPI_COMM_WORLD,errcode)
             CALL MPI_BCAST(DeckValues(2)%value,EntryLength,MPI_CHARACTER,0,MPI_COMM_WORLD,errcode)             
             CALL HandleDeckElement(DeckValues(1)%value,DeckValues(2)%value,Errcode_deck)
             DeckValues(1)%value=""
             DeckValues(2)%value=""
             IsComment=.FALSE.
             White_Space_Over=.FALSE.
          ENDIF
          IF (f .EQ. -1) THEN
             CALL MPI_BCAST(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)   
             CLOSE(LUN)
             EXIT
          ENDIF
          Terminate=Terminate .OR. IAND(errcode_deck,ERR_TERMINATE) .NE. 0
          IF (Terminate) EXIT
       ENDDO
    ELSE
       DO
          Errcode_deck=ERR_NONE
          CALL MPI_BCAST(f,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
          IF (f .EQ. 0) EXIT
          CALL MPI_BCAST(DeckValues(1)%value,EntryLength,MPI_CHARACTER,0,MPI_COMM_WORLD,errcode)
          CALL MPI_BCAST(DeckValues(2)%value,EntryLength,MPI_CHARACTER,0,MPI_COMM_WORLD,errcode)             
          CALL HandleDeckElement(DeckValues(1)%value,DeckValues(2)%value,Errcode_deck)
          DeckValues(1)%value=""
          DeckValues(2)%value=""
          Terminate=Terminate .OR. IAND(errcode_deck,ERR_TERMINATE) .NE. 0
          IF (Terminate) EXIT
       ENDDO
    ENDIF


    CALL MPI_BARRIER(MPI_COMM_WORLD,errcode)

    !Don't check compulsory blocks if going to bomb anyway, just stinks up the output file
    IF (.NOT. Terminate .AND. FirstCall) CALL CheckCompulsoryBlocks(errcode_deck)
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

    IF (FirstCall) CLOSE(40)

    IF (Terminate) CALL MPI_ABORT(MPI_COMM_WORLD,errcode)

    CALL MPI_BARRIER(MPI_COMM_WORLD,errcode)


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

    IF (StrCmp(Element,"import")) THEN
       InvalidBlock=.TRUE.
       IF (rank .EQ. rankcheck) THEN
          WRITE(40,*),""
          WRITE(40,*),"Importing ",TRIM(ADJUSTL(Value)), " file"
          WRITE(40,*),""
       ENDIF
       CALL Read_Deck(TRIM(ADJUSTL(Value)),.FALSE.)
       RETURN
    ENDIF

    IF (StrCmp(Element,"begin")) THEN
       Errcode_deck=HandleBlock(Value,Blank,Blank)
       InvalidBlock=IAND(Errcode_deck,ERR_UNKNOWN_BLOCK) .NE. 0
       IF(InvalidBlock) THEN
          IF (RANK .EQ. rankcheck) THEN
             PRINT *,char(9),"Unknown block ",TRIM(Value)," in input deck, ignoring"
          ENDIF
       ENDIF
       CALL StartBlock(Value)
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
       CALL EndBlock(Current_Block_Name)
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
    !If an error is fatal then set Terminate to .TRUE.
    IF (IAND(Errcode_deck,ERR_UNKNOWN_ELEMENT) /= 0) THEN
       IF (rank .EQ. rankcheck) THEN
          WRITE(40,*) ""
          PRINT *,""
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
          PRINT *,""
          PRINT *,"***WARNING*** Element ",TRIM(Element), " is set multiple times in this deck. Code will continue using last value in deck"
          WRITE(40,*) "***WARNING*** Element ",TRIM(Element), " is set multiple times in this deck. Code will continue using last value in deck"
          PRINT *,""
          WRITE(40,*) ""
       ENDIF
    ENDIF
    IF (IAND(Errcode_deck, ERR_BAD_VALUE) /= 0) THEN
       IF (rank .EQ. rankcheck) THEN
          WRITE(40,*) ""
          PRINT *,""
          PRINT *,"***ERROR*** Value ",TRIM(Value)," in element ",TRIM(Element)," is invalid or could not be parsed. Code will terminate."
          WRITE(40,*) "***ERROR*** Value ",TRIM(Value)," in element ",TRIM(Element)," is invalid or could not be parsed. Code will terminate."
          PRINT *,""
          WRITE(40,*) ""
          Errcode_deck=IOR(Errcode_deck,ERR_TERMINATE)
       ENDIF
    ENDIF
    IF (IAND(Errcode_deck, ERR_REQUIRED_ELEMENT_NOT_SET) /= 0) THEN
       IF (rank .EQ. rankcheck) THEN
          WRITE(40,*) ""
          PRINT *,""
          PRINT *,"***ERROR*** Value ",TRIM(Value)," in element ",TRIM(Element)," cannot be set because a prerequisite element, ", TRIM(Extended_Error_String),",has not been set. Code will terminate"
          WRITE(40,*) "***ERROR*** Value ",TRIM(Value)," in element ",TRIM(Element)," cannot be set because a prerequisite element, ", TRIM(Extended_Error_String),",has not been set. Code will terminate"
          PRINT *,""
          WRITE(40,*) ""
          Errcode_deck=IOR(Errcode_deck,ERR_TERMINATE)
       ENDIF
    ENDIF
    IF (IAND(Errcode_deck , ERR_PP_OPTIONS_WRONG) /= 0) THEN
       IF (rank .EQ. rankcheck) THEN
          WRITE(40,*) ""
          PRINT *,""
          PRINT *,"***ERROR*** The element ",TRIM(Element)," of block ",TRIM(Current_Block_Name)," cannot be set because the code has not been compiled with the correct preprocessor options.",&
               "Code will continue, but to use selected features, please recompile with ",TRIM(Extended_Error_String)," option"
          WRITE(40,*) "***ERROR*** The element ",TRIM(Element)," of block ",TRIM(Current_Block_Name)," cannot be set because the code has not been compiled with the correct preprocessor options.",&
               "Code will continue, but to use selected features, please recompile with ",TRIM(Extended_Error_String)," option"
          PRINT *,""
          WRITE(40,*) ""
       ENDIF
    ENDIF
    IF (IAND(Errcode_deck , ERR_OTHER) /= 0) THEN
       IF (rank .EQ. rankcheck) THEN
          WRITE(40,*) ""
          PRINT *,""
          PRINT *,"***ERROR*** You have managed to find an impossible situation in this code. Good for you. Just because of that, code will terminate."
          WRITE(40,*) "***ERROR*** You have managed to find an impossible situation in this code. Good for you. Just because of that, code will terminate."
          PRINT *,""
          WRITE(40,*) ""
          Errcode_deck=IOR(Errcode_deck,ERR_TERMINATE)
       ENDIF
    ENDIF

    ErrCount=ErrCount+1

  END SUBROUTINE HandleDeckElement

END MODULE deck
