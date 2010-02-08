MODULE input_arb

  USE shared_data
  USE iocommon
  USE input_functions
  USE output

  IMPLICIT NONE

CONTAINS

  !This subroutine is used to wrap a block containing program specific data
  !Which there is no general way of allowing other programs to read
  !It permits the use of a single string to idenitify the program that wrote it
  SUBROUTINE cfd_get_arb_block(reader)

    INTERFACE
       SUBROUTINE reader(filehandle,current_displacement,generator_name)
         USE shared_data
         INTEGER,INTENT(IN) :: filehandle
         INTEGER(KIND=MPI_OFFSET_KIND),INTENT(IN) :: current_displacement
         CHARACTER(len=*),INTENT(IN) :: generator_name
       END SUBROUTINE reader
    END INTERFACE

    CHARACTER(len=max_string_len) :: gen_name

    CALL cfd_skip_block_header
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, gen_name, max_string_len, MPI_CHARACTER, cfd_status, cfd_errcode)
    current_displacement=current_displacement+max_string_len
    CALL cfd_skip_block_metadata
    CALL reader(cfd_filehandle,current_displacement,gen_name)
    CALL cfd_skip_block


  END SUBROUTINE cfd_get_arb_block

END MODULE input_arb
