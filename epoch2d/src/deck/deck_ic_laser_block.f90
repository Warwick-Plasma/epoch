MODULE deck_ic_laser_block

  USE shared_data
  USE strings_advanced
  USE shared_parser_data
  USE strings
  USE shunt
  USE evaluator
  USE laser

  IMPLICIT NONE

  SAVE

  TYPE(Laser_Block),POINTER :: WorkingLaser
  LOGICAL :: Direction_Set=.FALSE.
  INTEGER :: Direction

CONTAINS

  FUNCTION HandleICLaserDeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    CHARACTER(LEN=EntryLength) :: Part1
    INTEGER :: Part2
    INTEGER :: HandleICLaserDeck
    INTEGER :: loop,elementselected,partswitch
    LOGICAL :: Handled,Temp
    REAL(num) :: dummy
    TYPE(primitivestack) :: output
    INTEGER :: ix,iy,ERR

    HandleICLaserDeck=ERR_NONE
    IF (Element .EQ. blank .OR. Value .EQ. blank) RETURN

    IF(StrCmp(Element,"direction")) THEN
       !If the direction has already been set, simply ignore further calls to it
       IF (Direction_Set) RETURN
       Direction=AsDirection(Value,HandleICLaserDeck)
       Direction_Set=.TRUE.
       CALL Init_Laser(Direction,WorkingLaser)
       RETURN
    ENDIF

    IF (.NOT. Direction_Set) THEN
       IF (Rank .EQ. 0) THEN
          PRINT *,"***ERROR*** Cannot set laser properties before direction is set"
          WRITE(40,*) "***ERROR*** Cannot set laser properties before direction is set"
       ENDIF
       Extended_Error_String="direction"
       HandleICLaserDeck=ERR_REQUIRED_ELEMENT_NOT_SET
       RETURN
    ENDIF

    IF (StrCmp(Element,"amp")) THEN
       WorkingLaser%Amp=AsReal(Value,HandleICLaserDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"freq")) THEN
       WorkingLaser%Freq=AsReal(Value,HandleICLaserDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"profile")) THEN
       workinglaser%profile=0.0_num
       output%stackpoint=0
       CALL Tokenize(Value,output,HandleICLaserDeck)
       IF (WorkingLaser%Direction .EQ. BD_LEFT .OR. WorkingLaser%Direction .EQ. BD_RIGHT) THEN
          DO iy=1,ny
             workinglaser%profile(iy) = EvaluateAtPoint(output,0,iy,HandleICLaserDeck)
          ENDDO
       ELSE IF (WorkingLaser%Direction .EQ. BD_UP .OR. WorkingLaser%Direction .EQ. BD_DOWN) THEN
          DO ix=1,nx
             workinglaser%profile(ix) = EvaluateAtPoint(output,ix,0,HandleICLaserDeck)
          ENDDO
       ENDIF
       RETURN
    ENDIF
    IF (StrCmp(Element,"phase")) THEN
       workinglaser%phase=0.0_num
       output%stackpoint=0
       CALL Tokenize(Value,output,HandleICLaserDeck)
       IF (WorkingLaser%Direction .EQ. BD_LEFT .OR. WorkingLaser%Direction .EQ. BD_RIGHT) THEN
          DO iy=1,ny
             workinglaser%phase(iy) = EvaluateAtPoint(output,0,iy,HandleICLaserDeck)
          ENDDO
       ELSE IF (WorkingLaser%Direction .EQ. BD_UP .OR. WorkingLaser%Direction .EQ. BD_DOWN) THEN
          DO ix=1,nx
             workinglaser%phase(ix) = EvaluateAtPoint(output,ix,0,HandleICLaserDeck)
          ENDDO
       ENDIF
       RETURN
    ENDIF
    IF (StrCmp(Element,"t_start")) THEN
       workinglaser%t_start=AsTime(Value,HandleICLaserDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"t_end")) THEN
       workinglaser%t_end=AsTime(Value,HandleICLaserDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"t_profile")) THEN
       WorkingLaser%UseTimeFunction=.TRUE.
       WorkingLaser%TimeFunction%StackPoint=0
       CALL Tokenize(Value,WorkingLaser%TimeFunction,HandleICLaserDeck)
       !Evaluate it once to check that it's a valid block
       dummy=Evaluate(WorkingLaser%TimeFunction,HandleICLaserDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"pol")) THEN
       !Convert from degrees to radians
       WorkingLaser%pol=pi * AsReal(Value,HandleICLaserDeck)/180.0_num
       RETURN
    ENDIF
    IF (StrCmp(Element,"id")) THEN
       WorkingLaser%ID=AsInteger(Value,HandleICLaserDeck)
       RETURN
    ENDIF
    HandleICLaserDeck=ERR_UNKNOWN_ELEMENT

  END FUNCTION HandleICLaserDeck

  FUNCTION CheckICLaserBlock()

    INTEGER :: CheckICLaserBlock

    !Should do error checking but can't be bothered at the moment
    CheckICLaserBlock=ERR_NONE

  END FUNCTION CheckICLaserBlock

  SUBROUTINE Laser_Start
    !Every new laser uses the internal time function
    ALLOCATE(WorkingLaser)
    WorkingLaser%UseTimeFunction=.FALSE.
  END SUBROUTINE Laser_Start

  SUBROUTINE Laser_End
    CALL Attach_Laser(WorkingLaser)
    Direction_Set=.FALSE.
  END SUBROUTINE Laser_End

  FUNCTION AsTime(Value,ERR)

    CHARACTER(LEN=*),INTENT(IN) :: Value
    INTEGER,INTENT(INOUT) :: ERR
    REAL(num) :: AsTime

    IF (StrCmp(Value,"start")) THEN
       AsTime=0.0_num
       RETURN
    ENDIF

    IF (StrCmp(Value,"end")) THEN
       AsTime=t_end
       RETURN
    ENDIF

    AsTime=AsReal(Value,ERR)

  END FUNCTION AsTime

END MODULE deck_ic_laser_block
