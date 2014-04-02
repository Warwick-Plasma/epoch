MODULE sdf_source_info

IMPLICIT NONE

CHARACTER(LEN=1) :: sdf_bytes_git_version = ''
CHARACTER(LEN=1) :: sdf_bytes_compile_date_string = ''
CHARACTER(LEN=1) :: sdf_bytes_compile_machine_info = ''
CHARACTER(LEN=1) :: sdf_bytes_compiler_info = ''
CHARACTER(LEN=1) :: sdf_bytes_compiler_flags = ''
INTEGER, PARAMETER :: sdf_bytes_compile_date = 1
CHARACTER(LEN=1) :: sdf_bytes_checksum_type = ''
CHARACTER(LEN=1) :: sdf_bytes_checksum = ''
CHARACTER(LEN=1) :: sdf_bytes_mimetype = ''
INTEGER, PARAMETER :: sdf_bytes_padding = 0
INTEGER, PARAMETER :: sdf_bytes_len = 0
INTEGER(8) :: sdf_bytes(1)
CHARACTER(LEN=1) :: sdf_diff_bytes_checksum_type = ''
CHARACTER(LEN=1) :: sdf_diff_bytes_checksum = ''
CHARACTER(LEN=1) :: sdf_diff_bytes_mimetype = ''
INTEGER, PARAMETER :: sdf_diff_bytes_padding = 0
INTEGER, PARAMETER :: sdf_diff_bytes_len = 0
INTEGER(8) :: sdf_diff_bytes(1)

END MODULE sdf_source_info
