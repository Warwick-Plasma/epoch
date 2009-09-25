MODULE deck_eio_dist_fn_block

  USE shared_data
  USE strings_advanced
  USE dist_fn

  IMPLICIT NONE

  SAVE

  TYPE(Distribution_Function_Block),POINTER :: WorkingBlock
CONTAINS

  FUNCTION HandleEIODistFnDeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER :: HandleEIODistFnDeck

    CHARACTER(LEN=EntryLength) :: Part1
    INTEGER :: Part2
    INTEGER :: work
    REAL(num) :: work1,work2

    HandleEIODistFnDeck=ERR_NONE
    IF (Element .EQ. blank .OR. Value .EQ. blank) RETURN

    IF (StrCmp(Element,"name")) THEN
       WorkingBlock%Name=Value
       RETURN
    ENDIF

    IF (StrCmp(Element,"ndims")) THEN
       work=AsInteger(Value,HandleEIODistFnDeck)
       IF (work .EQ. 2 .OR. work .EQ. 3) THEN
          WorkingBlock%nDims=work
       ELSE
          IF (rank .EQ. 0) PRINT *,"Distribution functions can only be 2D or 3D"
          HandleEIODistFnDeck=ERR_BAD_VALUE
       ENDIF
       RETURN
    ENDIF
    IF (WorkingBlock%nDims .EQ. -1) THEN
       IF (rank .EQ. 0) PRINT *,"Must set number of dimensions before setting other distribution function properties."
       Extended_Error_String="nDims"
       HandleEIODistFnDeck=ERR_REQUIRED_ELEMENT_NOT_SET
       RETURN
    ENDIF

    IF (StrCmp(Element,"dumpmask")) THEN
       WorkingBlock%DumpMask=AsInteger(Value,HandleEIODistFnDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"restrict_x")) THEN
       CALL SplitRange(Value,work1,work2,HandleEIODistFnDeck)
       IF (HandleEIODistFnDeck .NE. ERR_NONE) RETURN
       WorkingBlock%Use_Restrictions(1)=.TRUE.
       WorkingBlock%Restrictions(1,:)=(/work1,work2/)
    ENDIF
    IF (StrCmp(Element,"restrict_y")) THEN
       CALL SplitRange(Value,work1,work2,HandleEIODistFnDeck)
       IF (HandleEIODistFnDeck .NE. ERR_NONE) RETURN
       WorkingBlock%Use_Restrictions(2)=.TRUE.
       WorkingBlock%Restrictions(2,:)=(/work1,work2/)
    ENDIF
    IF (StrCmp(Element,"restrict_px")) THEN
       CALL SplitRange(Value,work1,work2,HandleEIODistFnDeck)
       IF (HandleEIODistFnDeck .NE. ERR_NONE) RETURN
       WorkingBlock%Use_Restrictions(3)=.TRUE.
       WorkingBlock%Restrictions(3,:)=(/work1,work2/)
    ENDIF
    IF (StrCmp(Element,"restrict_py")) THEN
       CALL SplitRange(Value,work1,work2,HandleEIODistFnDeck)
       IF (HandleEIODistFnDeck .NE. ERR_NONE) RETURN
       WorkingBlock%Use_Restrictions(4)=.TRUE.
       WorkingBlock%Restrictions(4,:)=(/work1,work2/)
    ENDIF
    IF (StrCmp(Element,"restrict_pz")) THEN
       CALL SplitRange(Value,work1,work2,HandleEIODistFnDeck)
       IF (HandleEIODistFnDeck .NE. ERR_NONE) RETURN
       WorkingBlock%Use_Restrictions(5)=.TRUE.
       WorkingBlock%Restrictions(5,:)=(/work1,work2/)
    ENDIF


    CALL SplitOffInt(Element,Part1,Part2,HandleEIODistFnDeck)
    IF (HandleEIODistFnDeck .NE. ERR_NONE) THEN
       HandleEIODistFnDeck = ERR_UNKNOWN_ELEMENT
       RETURN
    ENDIF
    IF (StrCmp(Part1,"direction")) THEN
       WorkingBlock%Directions(Part2)=AsReal(Value,HandleEIODistFnDeck)
       RETURN
    ENDIF
    IF (StrCmp(Part1,"range")) THEN
       CALL SplitRange(TRIM(Value),work1,work2,HandleEIODistFnDeck)
       IF (HandleEIODistFnDeck .NE. ERR_NONE) RETURN
       WorkingBlock%Ranges(Part2,1)=work1
       WorkingBlock%Ranges(Part2,2)=work2
       RETURN
    ENDIF
    IF (StrCmp(Part1,"resolution")) THEN
       WorkingBlock%Resolution(Part2)=AsInteger(Value,HandleEIODistFnDeck)
       RETURN
    ENDIF
    IF (StrCmp(Part1,"include_species_")) THEN
       IF (Part2 .LT. 1 .OR. Part2 .GT. nSpecies) THEN
          IF (rank .EQ. 0) PRINT *,"Species ",Part2," does not exist, ignoring attempt to set output state."
          HandleEIODistFnDeck=ERR_NONE
          RETURN
       ENDIF
       WorkingBlock%Use_Species(Part2)=AsLogical(Value,HandleEIODistFnDeck)
       RETURN
    ENDIF


    HandleEIODistFnDeck=ERR_UNKNOWN_ELEMENT

  END FUNCTION HandleEIODistFnDeck

  FUNCTION CheckEIODistFnBlock()

    INTEGER :: CheckEIODistFnBlock

    !Should do error checking but can't be bothered at the moment
    CheckEIODistFnBlock=ERR_NONE

  END FUNCTION CheckEIODistFnBlock

  SUBROUTINE Dist_Fn_Start
    !Every new laser uses the internal time function
    ALLOCATE(WorkingBlock)
    CALL Setup_Dist_Fn(WorkingBlock)

  END SUBROUTINE Dist_Fn_Start

  SUBROUTINE Dist_Fn_End

    CALL Attach_Dist_Fn(WorkingBlock)

  END SUBROUTINE Dist_Fn_End

END MODULE deck_eio_dist_fn_block
