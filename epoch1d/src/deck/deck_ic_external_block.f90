MODULE deck_ic_external_block

  USE shared_data
  USE strings_advanced
  USE simple_io

  IMPLICIT NONE

  SAVE

  INTEGER(KIND=MPI_OFFSET_KIND) :: offset
CONTAINS

  SUBROUTINE Start_External
    offset=0
  END SUBROUTINE Start_External

  FUNCTION HandleICExternalSpeciesDeck(Species_ID,Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER,INTENT(IN) :: Species_ID
    INTEGER :: HandleICExternalSpeciesDeck
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset =0 

    HandleICExternalSpeciesDeck=ERR_NONE
    IF (Element .EQ. blank .OR. Value .EQ. blank) RETURN

    IF (StrCmp(Element,"offset")) THEN
       offset=AsLongIntegerSimple(Value,HandleICExternalSpeciesDeck)
       RETURN
    ENDIF

    IF (Species_ID .LT. 0 .OR. Species_ID .GT. nspecies) THEN
       IF (rank .EQ. 0) PRINT *,"Attempting to set non-existent species initial conditions. Ignoring."
       RETURN
    ENDIF

    IF (StrCmp(Element,"minrho")) THEN
       InitialConditions(Species_ID)%MinRho=AsReal(Value,HandleICExternalSpeciesDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"maxrho")) THEN
       InitialConditions(Species_ID)%MaxRho=AsReal(Value,HandleICExternalSpeciesDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"rho")) THEN
       CALL Load_Single_Array_From_Data_File(Value,InitialConditions(Species_ID)%Rho,offset,HandleICExternalSpeciesDeck)  
       RETURN
    ENDIF

    IF (StrCmp(Element,"temp")) THEN
       CALL Load_Single_Array_From_Data_File(Value,InitialConditions(Species_ID)%Temp(:,1),offset&
            ,HandleICExternalSpeciesDeck)
       InitialConditions(Species_ID)%temp(-1:nx+2,2)=InitialConditions(Species_ID)%temp(-1:nx+2,1)
       InitialConditions(Species_ID)%temp(-1:nx+2,3)=InitialConditions(Species_ID)%temp(-1:nx+2,1)
       RETURN
    ENDIF


    IF (StrCmp(Element,"temp_x")) THEN
       CALL Load_Single_Array_From_Data_File(Value,InitialConditions(Species_ID)%Temp(:,1),offset&
            ,HandleICExternalSpeciesDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"temp_y")) THEN
       CALL Load_Single_Array_From_Data_File(Value,InitialConditions(Species_ID)%Temp(:,2),offset&
            ,HandleICExternalSpeciesDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"temp_z")) THEN
       CALL Load_Single_Array_From_Data_File(Value,InitialConditions(Species_ID)%Temp(:,3),offset&
            ,HandleICExternalSpeciesDeck)
       RETURN
    ENDIF
    HandleICExternalSpeciesDeck=ERR_UNKNOWN_ELEMENT

  END FUNCTION HandleICExternalSpeciesDeck

  FUNCTION HandleICExternalFieldsDeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER :: HandleICExternalFieldsDeck
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset =0 

    HandleICExternalFieldsDeck=ERR_NONE
    IF (Element .EQ. blank .OR. Value .EQ. blank) RETURN

    IF (StrCmp(Element,"ex")) THEN
       CALL Load_Single_Array_From_Data_File(Value,Ex,offset,HandleICExternalFieldsDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"ey")) THEN
       CALL Load_Single_Array_From_Data_File(Value,Ey,offset,HandleICExternalFieldsDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"ez")) THEN
       CALL Load_Single_Array_From_Data_File(Value,Ez,offset,HandleICExternalFieldsDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"bx")) THEN
       CALL Load_Single_Array_From_Data_File(Value,Bx,offset,HandleICExternalFieldsDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"by")) THEN
       CALL Load_Single_Array_From_Data_File(Value,By,offset,HandleICExternalFieldsDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"bz")) THEN
       CALL Load_Single_Array_From_Data_File(Value,Bz,offset,HandleICExternalFieldsDeck)
       RETURN
    ENDIF

  END FUNCTION HandleICExternalFieldsDeck

END MODULE deck_ic_external_block
