MODULE deck_ic_external_block

  USE shared_data
  USE strings_advanced
  USE mpi_subtype_control
  USE boundary

  IMPLICIT NONE

  SAVE

  INTEGER :: SpeciesLoaded

CONTAINS

  FUNCTION HandleICExternalSpeciesDeck(Species_ID,Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER,INTENT(IN) :: Species_ID
    CHARACTER(LEN=EntryLength) :: Part1
    INTEGER :: Part2
    INTEGER :: HandleICExternalSpeciesDeck
    INTEGER :: loop,elementselected,partswitch
    LOGICAL :: Handled,Temp
    REAL(num) :: conv_val
    TYPE(primitivestack) :: output
    INTEGER :: ix,iy,ERR

    HandleICExternalSpeciesDeck=ERR_NONE
    IF (Element .EQ. blank .OR. Value .EQ. blank) RETURN

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
       CALL LoadDataFile(Value,InitialConditions(Species_ID)%Rho,HandleICExternalSpeciesDeck)  
       RETURN
    ENDIF

    IF (StrCmp(Element,"temp")) THEN
       CALL LoadDataFile(Value,InitialConditions(Species_ID)%Temp(:,:,:,1),HandleICExternalSpeciesDeck)
       InitialConditions(Species_ID)%temp(-1:nx+2,-1:ny+2,-1:nz+2,2)=&
            InitialConditions(Species_ID)%temp(-1:nx+2,-1:ny+2,-1:nz+2,1)
       InitialConditions(Species_ID)%temp(-1:nx+2,-1:ny+2,-1:nz+2,3)=&
            InitialConditions(Species_ID)%temp(-1:nx+2,-1:ny+2,-1:nz+2,1)
       RETURN
    ENDIF


    IF (StrCmp(Element,"temp_x")) THEN
       CALL LoadDataFile(Value,InitialConditions(Species_ID)%Temp(:,:,:,1),HandleICExternalSpeciesDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"temp_y")) THEN
       CALL LoadDataFile(Value,InitialConditions(Species_ID)%Temp(:,:,:,2),HandleICExternalSpeciesDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"temp_z")) THEN
       CALL LoadDataFile(Value,InitialConditions(Species_ID)%Temp(:,:,:,3),HandleICExternalSpeciesDeck)
       RETURN
    ENDIF

  END FUNCTION HandleICExternalSpeciesDeck

  FUNCTION HandleICExternalFieldsDeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    CHARACTER(LEN=EntryLength) :: Part1
    INTEGER :: Part2
    INTEGER :: HandleICExternalFieldsDeck
    INTEGER :: loop,elementselected,partswitch
    LOGICAL :: Handled,Temp
    REAL(num) :: conv_val
    TYPE(primitivestack) :: output
    INTEGER :: ix,iy,ERR

    HandleICExternalFieldsDeck=ERR_NONE
    IF (Element .EQ. blank .OR. Value .EQ. blank) RETURN

    IF (StrCmp(Element,"ex")) THEN
       CALL LoadDataFile(Value,Ex,HandleICExternalFieldsDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"ey")) THEN
       CALL LoadDataFile(Value,Ey,HandleICExternalFieldsDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"ez")) THEN
       CALL LoadDataFile(Value,Ez,HandleICExternalFieldsDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"bx")) THEN
       CALL LoadDataFile(Value,Bx,HandleICExternalFieldsDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"by")) THEN
       CALL LoadDataFile(Value,By,HandleICExternalFieldsDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"bz")) THEN
       CALL LoadDataFile(Value,Bz,HandleICExternalFieldsDeck)
       RETURN
    ENDIF

  END FUNCTION HandleICExternalFieldsDeck

  SUBROUTINE LoadDataFile(FileName,Array,ERR)
    CHARACTER(LEN=*),INTENT(IN) :: FileName
    REAL(num),DIMENSION(-2:,-2:,-2:),INTENT(INOUT) :: Array
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER :: subtype, fh
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset=0

    CALL MPI_FILE_OPEN(comm,TRIM(Filename),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,errcode)
    IF (errcode .NE. 0) THEN
       IF (rank .EQ. 0) PRINT *,"File ",TRIM(FileName), " does not exist."
       ERR=IOR(ERR,ERR_BAD_VALUE)
       RETURN
    ENDIF
    subtype = Create_Current_Field_Subtype()
    CALL MPI_FILE_SET_VIEW(fh,offset,mpireal,subtype,"native",MPI_INFO_NULL,errcode)
    CALL MPI_FILE_READ_ALL(fh,Array(1:nx,1:ny,1:nz),nx*ny*nz,mpireal,status,errcode)
    CALL MPI_FILE_CLOSE(fh,errcode)
    CALL MPI_TYPE_FREE(subtype,errcode)

    CALL Field_BC(Array)
    CALL Field_Zero_Gradient(Array,.TRUE.)

  END SUBROUTINE LoadDataFile

END MODULE deck_ic_external_block
