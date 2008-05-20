MODULE input_arb

  USE shared_data
  USE iocommon
  USE inputfunctions
  USE output

  IMPLICIT NONE

CONTAINS

  !This subroutine is used to wrap a block containing program specific data
  !Which there is no general way of allowing other programs to read
  !It permits the use of a single string to idenitify the program that wrote it
  SUBROUTINE cfd_Get_Arb_Block(Reader)

    INTERFACE
       SUBROUTINE Reader(filehandle,current_displacement,Generator_Name)
         USE shared_data
         INTEGER,INTENT(IN) :: filehandle
         INTEGER(KIND=MPI_OFFSET_KIND),INTENT(IN) :: current_displacement
         CHARACTER(LEN=*),INTENT(IN) :: Generator_Name
       END SUBROUTINE Reader
    END INTERFACE

    CHARACTER(Len=MaxStringLen) :: Gen_Name

    CALL cfd_Skip_Block_Header
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, Gen_Name, MaxStringLen, MPI_CHARACTER, cfd_status, cfd_errcode)
    current_displacement=current_displacement+MaxStringLen
    CALL cfd_Skip_Block_Metadata
    CALL Reader(cfd_filehandle,current_displacement,Gen_Name)
    CALL cfd_Skip_Block


  END SUBROUTINE cfd_Get_Arb_Block

END MODULE input_arb
