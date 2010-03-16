MODULE version_data

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: c_code_name = "EPOCH3D"
  INTEGER(4), PARAMETER :: c_version = 1, c_revision = 3
  INTEGER(4), PARAMETER :: c_code_io_version = 1
  CHARACTER(LEN=*), PARAMETER :: c_commit_id = _COMMIT
  CHARACTER(LEN=*), PARAMETER :: c_compile_machine = _MACHINE
  CHARACTER(LEN=*), PARAMETER :: c_compile_flags = _FLAGS
  INTEGER(4), PARAMETER :: c_compile_date = _DATE

END MODULE version_data
