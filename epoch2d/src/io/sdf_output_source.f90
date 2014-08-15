MODULE sdf_output_source

  USE sdf_common
  USE sdf_source_info
  USE sdf_output

  IMPLICIT NONE

CONTAINS

  SUBROUTINE sdf_write_source_info(h)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=c_id_length) :: stitched_ids(3)
    CHARACTER(LEN=c_id_length) :: time_string
    CHARACTER(LEN=512) :: string_array(6)
    INTEGER :: n = 0

    IF (SIZE(sdf_bytes) .GT. 1 .OR. &
          TRIM(sdf_bytes_checksum_type) .NE. '') THEN
      n = n + 1
      CALL sdf_safe_copy_id(h, 'sdf_source/source', stitched_ids(n))
      CALL sdf_write_datablock(h, stitched_ids(n), &
          'SDF source code', sdf_bytes, &
          sdf_bytes_padding, sdf_bytes_mimetype, &
          sdf_bytes_checksum_type, sdf_bytes_checksum)
    ENDIF

    IF (SIZE(sdf_bytes) .EQ. 1 .AND. SIZE(sdf_diff_bytes) .GT. 1) THEN
      n = n + 1
      CALL sdf_safe_copy_id(h, 'sdf_source/diff', stitched_ids(n))
      CALL sdf_write_datablock(h, stitched_ids(n), &
          'SDF repository differences', sdf_diff_bytes, &
          sdf_diff_bytes_padding, sdf_diff_bytes_mimetype, &
          sdf_diff_bytes_checksum_type, sdf_diff_bytes_checksum)
    ENDIF

    n = n + 1
    CALL sdf_safe_copy_id(h, 'sdf_source/info', stitched_ids(n))
    WRITE(time_string, '(I20)') sdf_bytes_compile_date

    string_array(1) = TRIM(sdf_bytes_git_version)
    string_array(2) = TRIM(sdf_bytes_compile_date_string)
    string_array(3) = TRIM(ADJUSTL(time_string))
    string_array(4) = TRIM(sdf_bytes_compile_machine_info)
    string_array(5) = TRIM(sdf_bytes_compiler_info)
    string_array(6) = TRIM(sdf_bytes_compiler_flags)

    CALL sdf_write_namevalue(h, stitched_ids(n), &
        'SDF repository information', &
        (/'git_version         ', &
          'compile_date_string ', &
          'compile_date_seconds', &
          'compile_machine_info', &
          'compiler_info       ', &
          'compiler_flags      '/), string_array)

    CALL sdf_write_stitched(h, 'sdf_source', 'SDF source', &
        stitched_ids(1), c_stagger_cell_centre, stitched_ids, n)

  END SUBROUTINE sdf_write_source_info

END MODULE sdf_output_source
