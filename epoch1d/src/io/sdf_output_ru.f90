MODULE sdf_output_ru

  USE mpi
  USE sdf_common

  IMPLICIT NONE

CONTAINS

  SUBROUTINE sdf_write_header(h, code_name, code_io_version, step, time, &
      restart, jobid)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: code_name
    INTEGER, INTENT(IN) :: code_io_version, step
    REAL(num), INTENT(IN) :: time
    LOGICAL, INTENT(IN) :: restart
    TYPE(jobid_type), INTENT(IN), OPTIONAL :: jobid
    INTEGER(i4) :: int4
    INTEGER :: errcode
    INTEGER, PARAMETER :: sof8 = 8
    REAL(r8) :: real8
    CHARACTER(LEN=1) :: flag
    CHARACTER(LEN=6) :: padding

    IF (h%done_header) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'SDF header already written. Ignoring extra call.'
      ENDIF
      RETURN
    ENDIF

    IF (PRESENT(jobid)) THEN
      h%jobid = jobid
    ELSE
      h%jobid%start_seconds = 0
      h%jobid%start_milliseconds = 0
    ENDIF

    ! header length - must be updated if sdf_write_header changes
    h%first_block_location = c_header_length
    ! block header length - must be updated if sdf_write_block_header changes
    h%block_header_length = 3 * soi4 + 3 * soi8 + c_id_length &
        + h%string_length

    ! Currently no blocks written
    h%nblocks = 0
    h%summary_location = h%first_block_location
    h%summary_size = 0
    h%current_location = 0

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      ! Write the header
      CALL MPI_FILE_WRITE(h%filehandle, c_sdf_magic, 4, &
          MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, c_endianness, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, sdf_version, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, sdf_revision, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

      CALL sdf_safe_write_id(h, code_name)

      CALL MPI_FILE_WRITE(h%filehandle, h%first_block_location, 1, &
          MPI_INTEGER8, MPI_STATUS_IGNORE, errcode)

      ! Must be consistent with the c_summary_offset in sdf_common
      CALL MPI_FILE_WRITE(h%filehandle, h%summary_location, 1, &
          MPI_INTEGER8, MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, h%summary_size, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, h%nblocks, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, h%block_header_length, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

      int4 = INT(step,i4)
      CALL MPI_FILE_WRITE(h%filehandle, int4, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

      real8 = REAL(time,r8)
      CALL MPI_FILE_WRITE(h%filehandle, real8, 1, &
          MPI_REAL8, MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, h%jobid%start_seconds, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, h%jobid%start_milliseconds, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, h%string_length, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

      int4 = INT(code_io_version,i4)
      CALL MPI_FILE_WRITE(h%filehandle, int4, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

      IF (restart) THEN
        flag = ACHAR(1)
      ELSE
        flag = ACHAR(0)
      ENDIF

      CALL MPI_FILE_WRITE(h%filehandle, flag, 1, &
          MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)

      flag = ACHAR(0)
      CALL MPI_FILE_WRITE(h%filehandle, flag, 1, &
          MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)

      padding = ACHAR(0)
      CALL MPI_FILE_WRITE(h%filehandle, padding, 6, &
          MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_location = h%first_block_location
    h%done_header = .TRUE.

  END SUBROUTINE sdf_write_header



  SUBROUTINE write_block_header(h)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: b
    INTEGER :: errcode

    b => h%current_block
    IF (b%done_header) RETURN

    ! If this routine is changed then the value of h%block_header_length
    ! must be changed accordingly in sdf_write_header

    IF (h%rank .EQ. h%rank_master) THEN
      h%current_location = b%block_start

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      ! Write the block header
      CALL MPI_FILE_WRITE(h%filehandle, b%next_block_location, 1, &
          MPI_INTEGER8, MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%data_location, 1, &
          MPI_INTEGER8, MPI_STATUS_IGNORE, errcode)

      CALL sdf_safe_write_id(h, b%id)

      CALL MPI_FILE_WRITE(h%filehandle, b%data_length, 1, &
          MPI_INTEGER8, MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%blocktype, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%datatype, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%ndims, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

      CALL sdf_safe_write_string(h, b%name)
    ENDIF

    h%current_location = b%block_start + h%block_header_length
    b%done_header = .TRUE.

    !print*,'block header: b%id: ', TRIM(b%id)
    !print*,'   b%name: ', TRIM(b%name)
    !print*,'   b%blocktype: ', b%blocktype
    !print*,'   b%next_block_location: ', b%next_block_location
    !print*,'   b%data_location: ', b%data_location
    !print*,'   b%datatype: ', b%datatype
    !print*,'   b%ndims: ', b%ndims
    !print*,'   b%data_length: ', b%data_length

  END SUBROUTINE write_block_header



  SUBROUTINE sdf_write_block_header(h, id, name)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    CHARACTER(LEN=c_id_length) :: new_id
    TYPE(sdf_block_type), POINTER :: b

    b => h%current_block

    ! If this routine is changed then the value of h%block_header_length
    ! must be changed accordingly in sdf_open_clobber

    IF (.NOT. h%done_header) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'SDF header not yet written. Ignoring write call.'
      ENDIF
      RETURN
    ENDIF

    IF (b%done_header) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'SDF block header already written. Ignoring extra call.'
      ENDIF
      RETURN
    ENDIF

    b%data_location = b%block_start + b%info_length
    b%next_block_location = b%data_location + b%data_length

    CALL get_unique_id(h, id, new_id)
    CALL safe_copy_string(new_id, b%id)
    CALL safe_copy_string(name, b%name)

    CALL write_block_header(h)

    h%nblocks = h%nblocks + 1_4

  END SUBROUTINE sdf_write_block_header



  SUBROUTINE sdf_safe_write_string_len(h, string, length)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: length
    CHARACTER(LEN=length) :: output
    INTEGER :: len_s, errcode

    len_s = LEN(TRIM(string))

    IF (len_s .GT. length .AND. h%rank .EQ. h%rank_master) THEN
      PRINT*, '*** WARNING ***'
      PRINT*, 'Output string "' // TRIM(string) // '" has been truncated'
    ENDIF

    ! This subroutine expects that the record marker is in place and that
    ! the view is set correctly. Call it only on the node which is doing the
    ! writing. You still have to advance the file pointer yourself on all nodes

    output = ' '
    output(1:MIN(length, len_s)) = string(1:MIN(length, len_s))

    ! If this isn't the full string length then tag in a ACHAR(0) to help
    ! With C++ string handling
    IF (len_s + 1 .LT. length) output(len_s+1:length) = ACHAR(0)

    CALL MPI_FILE_WRITE(h%filehandle, output, length, &
        MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)

  END SUBROUTINE sdf_safe_write_string_len



  SUBROUTINE sdf_safe_write_string(h, string, length_in)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN), OPTIONAL :: length_in
    INTEGER :: length

    IF (PRESENT(length_in)) THEN
      length = length_in
    ELSE
      length = h%string_length
    ENDIF

    CALL sdf_safe_write_string_len(h, string, length)

  END SUBROUTINE sdf_safe_write_string



  SUBROUTINE sdf_safe_write_id(h, string)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: string

    CALL sdf_safe_write_string_len(h, string, c_id_length)

  END SUBROUTINE sdf_safe_write_id



  FUNCTION sdf_string_lowercase(string_in) RESULT(string_out)

    CHARACTER(LEN=*), PARAMETER :: lwr = 'abcdefghijklmnopqrstuvwxyz'
    CHARACTER(LEN=*), PARAMETER :: upr = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    CHARACTER(LEN=*), INTENT(IN) :: string_in
    CHARACTER(LEN=LEN(string_in)) :: string_out
    INTEGER :: i, idx

    string_out = string_in

    DO i = 1, LEN(string_out)
      idx = INDEX(upr, string_out(i:i))
      IF (idx .NE. 0) string_out(i:i) = lwr(idx:idx)
    ENDDO

  END FUNCTION sdf_string_lowercase



  SUBROUTINE safe_string_composite(string1, string2, output_string)

    CHARACTER(LEN=*), INTENT(IN) :: string1, string2
    CHARACTER(LEN=*), INTENT(OUT) :: output_string
    INTEGER :: len1, len2, olen

    len1 = LEN_TRIM(string1)
    len2 = LEN_TRIM(string2)
    olen = LEN(output_string)

    output_string = ''

    IF (olen < len1 + 1) THEN
      output_string(1:olen) = string1(1:olen)
    ELSE IF (olen < len1 + len2 + 1) THEN
      output_string(1:len1) = string1(1:len1)
      output_string(len1+1:len1+1) = '/'
      output_string(len1+2:olen) = string2(1:olen-len1-1)
    ELSE
      output_string(1:len1) = string1(1:len1)
      output_string(len1+1:len1+1) = '/'
      output_string(len1+2:len1+len2+1) = string2(1:len2)
    ENDIF

  END SUBROUTINE safe_string_composite



  SUBROUTINE write_run_info_meta(h, id, name)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    b => h%current_block

    b%datatype = c_datatype_other
    b%blocktype = c_blocktype_run_info
    b%ndims = 1

    ! Metadata is
    ! - version   INTEGER(i4)
    ! - revision  INTEGER(i4)
    ! - commit_id CHARACTER(string_length)
    ! - sha1sum   CHARACTER(string_length)
    ! - compmac   CHARACTER(string_length)
    ! - compflag  CHARACTER(string_length)
    ! - defines   INTEGER(i8)
    ! - compdate  INTEGER(i4)
    ! - rundate   INTEGER(i4)
    ! - iodate    INTEGER(i4)

    b%info_length = h%block_header_length + 5 * soi4 + soi8 &
        + 4 * h%string_length
    b%data_length = 0

    ! Write header
    IF (PRESENT(id)) THEN
      CALL sdf_write_block_header(h, id, name)
    ELSE
      CALL write_block_header(h)
    ENDIF

    IF (h%rank .EQ. h%rank_master) THEN
      ! Write metadata
      CALL MPI_FILE_WRITE(h%filehandle, b%run%version, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%run%revision, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      CALL sdf_safe_write_string(h, b%run%commit_id)

      CALL sdf_safe_write_string(h, b%run%sha1sum)

      CALL sdf_safe_write_string(h, b%run%compile_machine)

      CALL sdf_safe_write_string(h, b%run%compile_flags)

      CALL MPI_FILE_WRITE(h%filehandle, b%run%defines, 1, MPI_INTEGER8, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%run%compile_date, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%run%run_date, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%run%io_date, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_location = b%block_start + b%info_length
    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE write_run_info_meta



  SUBROUTINE sdf_write_run_info(h, version, revision, commit_id, &
      sha1sum, compile_machine, compile_flags, defines, compile_date, &
      run_date, rank_write)

    TYPE(sdf_file_handle) :: h
    INTEGER(i4), INTENT(IN) :: version, revision
    CHARACTER(LEN=*), INTENT(IN) :: commit_id, sha1sum
    CHARACTER(LEN=*), INTENT(IN) :: compile_machine, compile_flags
    INTEGER(i8), INTENT(IN) :: defines
    INTEGER(i4), INTENT(IN) :: compile_date, run_date
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    IF (.NOT. ASSOCIATED(b%run)) ALLOCATE(b%run)
    b%run%version = version
    b%run%revision = revision
    CALL safe_copy_string(commit_id, b%run%commit_id)
    CALL safe_copy_string(sha1sum, b%run%sha1sum)
    CALL safe_copy_string(compile_machine, b%run%compile_machine)
    CALL safe_copy_string(compile_flags, b%run%compile_flags)
    b%run%defines = defines
    b%run%compile_date = compile_date
    b%run%run_date = run_date
    b%run%io_date = get_unix_time()

    CALL write_run_info_meta(h, 'run_info', 'Run_info')

    h%rank_master = h%default_rank

  END SUBROUTINE sdf_write_run_info



  SUBROUTINE sdf_write_stitched_tensor(h, id, name, mesh_id, stagger, &
      variable_ids, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name, mesh_id
    INTEGER(i4), INTENT(IN), OPTIONAL :: stagger
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: variable_ids(:)
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (PRESENT(id)) THEN
      CALL sdf_get_next_block(h)
      b => h%current_block
      b%ndims = INT(SIZE(variable_ids),i4)
    ENDIF

    b => h%current_block

    b%datatype = c_datatype_other
    b%blocktype = c_blocktype_stitched_tensor

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Metadata is
    ! - stagger   INTEGER(i4)
    ! - meshid    CHARACTER(id_length)
    ! - varids    ndims*CHARACTER(id_length)

    b%info_length = h%block_header_length + soi4 + (b%ndims + 1) * c_id_length
    b%data_length = 0

    ! Write header
    IF (PRESENT(id)) THEN
      b%stagger = stagger
      CALL safe_copy_string(mesh_id, b%mesh_id)
      CALL sdf_write_block_header(h, id, name)
      ALLOCATE(b%variable_ids(b%ndims))
      DO i = 1, b%ndims
        CALL safe_copy_string(variable_ids(i), b%variable_ids(i))
      ENDDO
    ELSE
      CALL write_block_header(h)
    ENDIF

    IF (h%rank .EQ. h%rank_master) THEN
      ! Write metadata
      CALL MPI_FILE_WRITE(h%filehandle, b%stagger, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      CALL sdf_safe_write_id(h, b%mesh_id)

      DO i = 1, b%ndims
        CALL sdf_safe_write_id(h, b%variable_ids(i))
      ENDDO
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%block_start + b%info_length
    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_stitched_tensor



  SUBROUTINE sdf_write_stitched_material(h, id, name, mesh_id, stagger, &
      material_names, variable_ids, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name, mesh_id
    INTEGER(i4), INTENT(IN), OPTIONAL :: stagger
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: material_names(:), variable_ids(:)
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (PRESENT(id)) THEN
      CALL sdf_get_next_block(h)
      b => h%current_block
      b%ndims = INT(SIZE(variable_ids),i4)
    ENDIF

    b => h%current_block

    b%datatype = c_datatype_other
    b%blocktype = c_blocktype_stitched_material

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Metadata is
    ! - stagger   INTEGER(i4)
    ! - meshid    CHARACTER(id_length)
    ! - material_names ndims*CHARACTER(string_length)
    ! - varids    ndims*CHARACTER(id_length)

    b%info_length = h%block_header_length + soi4 + (b%ndims + 1) * c_id_length &
        + b%ndims * h%string_length
    b%data_length = 0

    ! Write header
    IF (PRESENT(id)) THEN
      b%stagger = stagger
      CALL safe_copy_string(mesh_id, b%mesh_id)
      CALL sdf_write_block_header(h, id, name)
      ALLOCATE(b%material_names(b%ndims))
      ALLOCATE(b%variable_ids(b%ndims))
      DO i = 1, b%ndims
        CALL safe_copy_string(material_names(i), b%material_names(i))
        CALL safe_copy_string(variable_ids(i), b%variable_ids(i))
      ENDDO
    ELSE
      CALL write_block_header(h)
    ENDIF

    IF (h%rank .EQ. h%rank_master) THEN
      ! Write metadata
      CALL MPI_FILE_WRITE(h%filehandle, b%stagger, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      CALL sdf_safe_write_id(h, b%mesh_id)

      DO i = 1, b%ndims
        CALL sdf_safe_write_string(h, b%material_names(i))
      ENDDO

      DO i = 1, b%ndims
        CALL sdf_safe_write_id(h, b%variable_ids(i))
      ENDDO
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%block_start + b%info_length
    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_stitched_material



  SUBROUTINE sdf_write_stitched_matvar(h, id, name, mesh_id, stagger, &
      material_id, variable_ids, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name, mesh_id
    INTEGER(i4), INTENT(IN), OPTIONAL :: stagger
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: material_id, variable_ids(:)
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (PRESENT(id)) THEN
      CALL sdf_get_next_block(h)
      b => h%current_block
      b%ndims = INT(SIZE(variable_ids),i4)
    ENDIF

    b => h%current_block

    b%datatype = c_datatype_other
    b%blocktype = c_blocktype_stitched_matvar

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Metadata is
    ! - stagger   INTEGER(i4)
    ! - meshid    CHARACTER(id_length)
    ! - matid     CHARACTER(id_length)
    ! - varids    ndims*CHARACTER(id_length)

    b%info_length = h%block_header_length + soi4 + (b%ndims + 2) * c_id_length
    b%data_length = 0

    ! Write header
    IF (PRESENT(id)) THEN
      b%stagger = stagger
      CALL safe_copy_string(mesh_id, b%mesh_id)
      CALL safe_copy_string(material_id, b%material_id)
      CALL sdf_write_block_header(h, id, name)
      ALLOCATE(b%variable_ids(b%ndims))
      DO i = 1, b%ndims
        CALL safe_copy_string(variable_ids(i), b%variable_ids(i))
      ENDDO
    ELSE
      CALL write_block_header(h)
    ENDIF

    IF (h%rank .EQ. h%rank_master) THEN
      ! Write metadata
      CALL MPI_FILE_WRITE(h%filehandle, b%stagger, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      CALL sdf_safe_write_id(h, b%mesh_id)

      CALL sdf_safe_write_id(h, b%material_id)

      DO i = 1, b%ndims
        CALL sdf_safe_write_id(h, b%variable_ids(i))
      ENDDO
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%block_start + b%info_length
    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_stitched_matvar



  SUBROUTINE sdf_write_stitched_species(h, id, name, mesh_id, stagger, &
      material_id, material_name, specnames, variable_ids, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name, mesh_id
    INTEGER(i4), INTENT(IN), OPTIONAL :: stagger
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: material_id, material_name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: specnames(:), variable_ids(:)
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (PRESENT(id)) THEN
      CALL sdf_get_next_block(h)
      b => h%current_block
      b%ndims = INT(SIZE(variable_ids),i4)
    ENDIF

    b => h%current_block

    b%datatype = c_datatype_other
    b%blocktype = c_blocktype_stitched_species

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Metadata is
    ! - stagger   INTEGER(i4)
    ! - meshid    CHARACTER(id_length)
    ! - matid     CHARACTER(id_length)
    ! - matname   CHARACTER(string_length)
    ! - specnames ndims*CHARACTER(string_length)
    ! - varids    ndims*CHARACTER(id_length)

    b%info_length = h%block_header_length + soi4 + (b%ndims + 2) * c_id_length &
        + (b%ndims + 1) * h%string_length
    b%data_length = 0

    ! Write header
    IF (PRESENT(id)) THEN
      b%stagger = stagger
      CALL safe_copy_string(mesh_id, b%mesh_id)
      CALL safe_copy_string(material_id, b%material_id)
      CALL safe_copy_string(material_name, b%material_name)
      ALLOCATE(b%material_names(b%ndims))
      ALLOCATE(b%variable_ids(b%ndims))
      DO i = 1, b%ndims
        CALL safe_copy_string(specnames(i), b%material_names(i))
        CALL safe_copy_string(variable_ids(i), b%variable_ids(i))
      ENDDO
      CALL sdf_write_block_header(h, id, name)
    ELSE
      CALL write_block_header(h)
    ENDIF

    IF (h%rank .EQ. h%rank_master) THEN
      ! Write metadata
      CALL MPI_FILE_WRITE(h%filehandle, b%stagger, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      CALL sdf_safe_write_id(h, b%mesh_id)

      CALL sdf_safe_write_id(h, b%material_id)

      CALL sdf_safe_write_string(h, b%material_name)

      DO i = 1, b%ndims
        CALL sdf_safe_write_string(h, b%material_names(i))
      ENDDO

      DO i = 1, b%ndims
        CALL sdf_safe_write_id(h, b%variable_ids(i))
      ENDDO
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%block_start + b%info_length
    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_stitched_species



  SUBROUTINE write_constant_meta(h, id, name)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name
    INTEGER :: errcode, var_len
    TYPE(sdf_block_type), POINTER :: b

    b => h%current_block

    b%blocktype = c_blocktype_constant
    b%ndims = 1

    var_len = 1
    b%nelements = var_len
    b%info_length = h%block_header_length + b%type_size
    b%data_length = 0

    ! Write header
    IF (PRESENT(id)) THEN
      CALL sdf_write_block_header(h, id, name)
    ELSE
      CALL write_block_header(h)
    ENDIF

    IF (h%rank .EQ. h%rank_master) THEN
      ! Write data (in metadata section)
      CALL MPI_FILE_WRITE(h%filehandle, b%const_value, var_len, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_location = b%block_start + b%info_length
    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE write_constant_meta



  SUBROUTINE sdf_write_constant_integer(h, id, name, value, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER, INTENT(IN) :: value
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%soi
    b%datatype = h%datatype_integer
    b%mpitype = h%mpitype_integer

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    b%const_value = TRANSFER(value, b%const_value)

    CALL write_constant_meta(h, id, name)

    h%rank_master = h%default_rank

  END SUBROUTINE sdf_write_constant_integer



  SUBROUTINE sdf_write_constant_logical(h, id, name, value, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    LOGICAL, INTENT(IN) :: value
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    CHARACTER(LEN=1) :: cvalue
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 1
    b%datatype = c_datatype_logical
    b%mpitype = MPI_CHARACTER

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    IF (value) THEN
      cvalue = ACHAR(1)
    ELSE
      cvalue = ACHAR(0)
    ENDIF

    b%const_value = TRANSFER(cvalue, b%const_value)

    CALL write_constant_meta(h, id, name)

    h%rank_master = h%default_rank

  END SUBROUTINE sdf_write_constant_logical



  SUBROUTINE sdf_write_source_code(h, id, name, array, last, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: array
    CHARACTER(LEN=*), INTENT(IN) :: last
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER(i8) :: sz
    INTEGER :: i, errcode, len1, len2
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 1
    b%datatype = c_datatype_character
    b%mpitype = MPI_CHARACTER
    b%blocktype = c_blocktype_source
    b%ndims = 0

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    IF (h%rank .EQ. h%rank_master) THEN
      sz   = SIZE(array)
      len1 = LEN(array)
      len2 = LEN(last)
      b%info_length = h%block_header_length
      b%data_length = sz*len1 + len2
    ENDIF

    CALL MPI_BCAST(b%data_length, 1, MPI_INTEGER8, 0, h%comm, errcode)

    ! Write header
    CALL sdf_write_block_header(h, id, name)

    h%current_location = b%data_location

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      ! Write data
      DO i = 1, sz
        CALL MPI_FILE_WRITE(h%filehandle, array(i), len1, &
            b%mpitype, MPI_STATUS_IGNORE, errcode)
      ENDDO
      CALL MPI_FILE_WRITE(h%filehandle, last, len2, &
          b%mpitype, MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_source_code



  SUBROUTINE write_array_meta(h, id, name)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name
    INTEGER :: errcode, ndims, i
    TYPE(sdf_block_type), POINTER :: b

    b => h%current_block

    b%blocktype = c_blocktype_array
    ndims = b%ndims

    b%nelements = 1
    DO i = 1,ndims
      b%nelements = b%nelements * b%dims(i)
    ENDDO

    ! Metadata is
    ! - dims      ndims*INTEGER(i4)

    b%info_length = h%block_header_length + b%ndims * soi4
    b%data_length = b%nelements * b%type_size

    ! Write header

    IF (PRESENT(id)) THEN
      CALL sdf_write_block_header(h, id, name)
    ELSE
      CALL write_block_header(h)
    ENDIF

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_WRITE(h%filehandle, b%dims, ndims, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_location = b%block_start + b%info_length
    b%done_info = .TRUE.

  END SUBROUTINE write_array_meta



  SUBROUTINE sdf_write_1d_array_integer(h, id, name, array, rank_write)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER, DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: errcode, n1
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%soi
    b%datatype = h%datatype_integer
    b%mpitype = h%mpitype_integer
    b%ndims = ndims

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    n1 = SIZE(array,1)
    b%dims(1) = n1

    ! Write header

    CALL write_array_meta(h, id, name)

    h%current_location = b%data_location

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      ! Actual array
      CALL MPI_FILE_WRITE(h%filehandle, array, n1, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_1d_array_integer



  SUBROUTINE sdf_write_2d_array_integer_spec(h, id, name, n1, n2, array, &
      rank_write)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER, INTENT(IN) :: n1, n2
    INTEGER, DIMENSION(:,:), INTENT(IN) :: array
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: errcode, var_len, i
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%soi
    b%datatype = h%datatype_integer
    b%mpitype = h%mpitype_integer
    b%ndims = ndims

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    b%dims(1) = n1
    b%dims(2) = n2

    ! Write header

    CALL write_array_meta(h, id, name)

    h%current_location = b%data_location

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      ! Actual array
      IF (n1 .EQ. SIZE(array,1)) THEN
        var_len = b%nelements
        CALL MPI_FILE_WRITE(h%filehandle, array, var_len, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
      ELSE
        DO i = 1,n2
          CALL MPI_FILE_WRITE(h%filehandle, array(1,i), n1, b%mpitype, &
              MPI_STATUS_IGNORE, errcode)
        ENDDO
      ENDIF
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_2d_array_integer_spec



  SUBROUTINE sdf_write_2d_array_integer(h, id, name, array, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER, DIMENSION(:,:), INTENT(IN) :: array
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: n1, n2

    n1 = SIZE(array,1)
    n2 = SIZE(array,2)
    CALL sdf_write_2d_array_integer_spec(h, id, name, n1, n2, &
        array, rank_write)

  END SUBROUTINE sdf_write_2d_array_integer



  SUBROUTINE sdf_write_1d_array_logical(h, id, name, array, rank_write)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    LOGICAL, DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    CHARACTER(LEN=1), DIMENSION(:), ALLOCATABLE :: carray
    INTEGER :: errcode, n1, i
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 1
    b%datatype = c_datatype_logical
    b%mpitype = MPI_CHARACTER
    b%ndims = ndims

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    n1 = SIZE(array,1)
    b%dims(1) = n1

    ! Write header

    CALL write_array_meta(h, id, name)

    h%current_location = b%data_location

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      ALLOCATE(carray(n1))
      DO i = 1,n1
        IF (array(i)) THEN
          carray(i) = ACHAR(1)
        ELSE
          carray(i) = ACHAR(0)
        ENDIF
      ENDDO

      ! Actual array
      CALL MPI_FILE_WRITE(h%filehandle, carray, n1, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)

      DEALLOCATE(carray)
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_1d_array_logical



  SUBROUTINE sdf_write_2d_array_character(h, id, name, array, rank_write)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(IN), OPTIONAL :: rank_write
    INTEGER :: errcode, var_len, n1, n2
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 1
    b%datatype = c_datatype_character
    b%mpitype = MPI_CHARACTER
    b%ndims = ndims

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    n1 = LEN(array)
    n2 = SIZE(array)
    b%dims(1) = n1
    b%dims(2) = n2

    ! Write header

    CALL write_array_meta(h, id, name)

    h%current_location = b%data_location

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      ! Actual array
      var_len = b%nelements
      CALL MPI_FILE_WRITE(h%filehandle, array, var_len, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_2d_array_character

END MODULE sdf_output_ru
