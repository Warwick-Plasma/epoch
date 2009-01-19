MODULE output_arb

  USE shared_data
  USE iocommon
  USE output

  IMPLICIT NONE

CONTAINS

  !This subroutine is used to wrap a block containing program specific data
  !Which there is no general way of allowing other programs to read
  !It permits the use of a single string to idenitify the program that wrote it
  SUBROUTINE cfd_Write_Arb_Block(Name,Class,Generator_Desc,Data_Length,Writer)

    CHARACTER(Len=*),INTENT(IN) :: Name,Class,Generator_Desc
    INTEGER(8),INTENT(IN) :: Data_Length
    INTERFACE
       SUBROUTINE Writer(filehandle,current_displacement)
         USE shared_data
         INTEGER,INTENT(IN) :: filehandle
         INTEGER(KIND=MPI_OFFSET_KIND),INTENT(IN) :: current_displacement
       END SUBROUTINE Writer
    END INTERFACE

    INTEGER(KIND=8) :: mdlength,blocklength
    INTEGER(KIND=MPI_OFFSET_KIND) :: initial_displacement

    !Outputs general block header as described in cfd_Write_Block_Header and then a single
    !string

    mdlength=1*MaxStringLen
    blocklength=mdlength+Data_Length

    CALL cfd_Write_Block_Header(name,class,TYPE_ARB_DB,blocklength,mdlength,default_rank)
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
       CALL cfd_Safe_Write_String(Generator_Desc)
    ENDIF
    current_displacement=current_displacement+MaxStringLen

    CALL Writer(cfd_filehandle,current_displacement)
    current_displacement=current_displacement+Data_Length


  END SUBROUTINE Cfd_Write_Arb_Block

END MODULE output_arb
