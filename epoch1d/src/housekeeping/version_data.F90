MODULE version_data

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: c_code_name = "EPOCH1D"
  INTEGER(4), PARAMETER :: c_version = 3, c_revision = 1
  INTEGER(4), PARAMETER :: c_code_io_version = 1
  CHARACTER(LEN=*), PARAMETER :: c_commit_id = _COMMIT
  CHARACTER(LEN=*), PARAMETER :: c_compile_machine = _MACHINE
  CHARACTER(LEN=*), PARAMETER :: c_compile_flags = "unknown"
  INTEGER(4), PARAMETER :: c_compile_date = _DATE
  CHARACTER(LEN=16) :: version_string
  CHARACTER(LEN=70) :: ascii_header

END MODULE version_data
