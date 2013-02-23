MODULE sdf_common

  USE mpi
  USE sdf_job_info

  IMPLICIT NONE

  INTEGER, PARAMETER :: i4  = SELECTED_INT_KIND(9)  ! 4-byte 2^31 ~ 10^9
  INTEGER, PARAMETER :: i8  = SELECTED_INT_KIND(18) ! 8-byte 2^63 ~ 10^18

  INTEGER, PARAMETER :: r4  = SELECTED_REAL_KIND(r=30)
  INTEGER, PARAMETER :: r8  = SELECTED_REAL_KIND(r=300)
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(r=3000)

  INTEGER, PARAMETER :: c_maxdims = 4
  INTEGER(i4), PARAMETER :: c_id_length = 32
  INTEGER(i4), PARAMETER :: c_long_id_length = 256
  INTEGER(i4), PARAMETER :: c_max_string_length = 64
  INTEGER(i8) :: npoint_per_iteration = 10000
  CHARACTER(LEN=4), PARAMETER :: c_sdf_magic = 'SDF1'

  TYPE sdf_run_type
    INTEGER(i4) :: version, revision, compile_date, run_date, io_date
    INTEGER(i8) :: defines
    CHARACTER(LEN=c_max_string_length) :: commit_id, sha1sum
    CHARACTER(LEN=c_max_string_length) :: compile_machine, compile_flags
  END TYPE sdf_run_type

  TYPE sdf_block_type
    REAL(r8), DIMENSION(2*c_maxdims) :: extents
    REAL(r8) :: mult
    REAL(r8), DIMENSION(:), POINTER :: dim_mults
    REAL(r8) :: const_value
    INTEGER(KIND=MPI_OFFSET_KIND) :: block_start
    INTEGER(i8) :: next_block_location, data_location
    INTEGER(i8) :: nelements, npoints, data_length, info_length
    INTEGER(i4) :: ndims, geometry, datatype, blocktype
    INTEGER(i4) :: mpitype, type_size, stagger
    INTEGER(i4), DIMENSION(c_maxdims) :: dims
    CHARACTER(LEN=c_id_length) :: id, units, mesh_id, material_id
    CHARACTER(LEN=c_id_length) :: vfm_id, obstacle_id
    CHARACTER(LEN=c_long_id_length) :: long_id
    CHARACTER(LEN=c_id_length), POINTER :: variable_ids(:)
    CHARACTER(LEN=c_id_length), POINTER :: dim_labels(:), dim_units(:)
    CHARACTER(LEN=c_max_string_length) :: name, material_name
    CHARACTER(LEN=c_max_string_length), POINTER :: material_names(:)
    LOGICAL :: done_header, done_info, done_data, truncated_id
    TYPE(sdf_run_type), POINTER :: run
    TYPE(sdf_block_type), POINTER :: next_block
  END TYPE sdf_block_type

  TYPE sdf_file_handle
    INTEGER(KIND=MPI_OFFSET_KIND) :: current_location
    REAL(r8) :: time
    INTEGER(i8) :: first_block_location, summary_location, start_location
    INTEGER(i8) :: soi ! large integer to prevent overflow in calculations
    INTEGER(i8) :: data_location
    INTEGER(i4) :: endianness, summary_size
    INTEGER(i4) :: block_header_length, string_length, nblocks, error_code
    INTEGER(i4) :: file_version, file_revision, code_io_version, step
    INTEGER(i4) :: datatype_integer, mpitype_integer
    INTEGER(i4) :: blocktype
    INTEGER :: filehandle, comm, rank, rank_master, default_rank, mode
    INTEGER :: errhandler
    LOGICAL :: done_header, restart_flag, other_domains, writing, handled_error
    CHARACTER(LEN=1), POINTER :: buffer(:)
    CHARACTER(LEN=c_id_length) :: code_name
    TYPE(jobid_type) :: jobid
    TYPE(sdf_block_type), POINTER :: blocklist, current_block
  END TYPE sdf_file_handle

  TYPE sdf_handle_type
    INTEGER :: filehandle
    TYPE(sdf_file_handle), POINTER :: handle
  END TYPE sdf_handle_type
  INTEGER, PARAMETER :: max_handles = 64
  TYPE(sdf_handle_type) :: sdf_handles(max_handles)

  INTEGER, PARAMETER :: c_sdf_read = 0
  INTEGER, PARAMETER :: c_sdf_write = 1

  INTEGER(i4), PARAMETER :: c_blocktype_scrubbed = -1
  INTEGER(i4), PARAMETER :: c_blocktype_null = 0
  INTEGER(i4), PARAMETER :: c_blocktype_plain_mesh = 1
  INTEGER(i4), PARAMETER :: c_blocktype_point_mesh = 2
  INTEGER(i4), PARAMETER :: c_blocktype_plain_variable = 3
  INTEGER(i4), PARAMETER :: c_blocktype_point_variable = 4
  INTEGER(i4), PARAMETER :: c_blocktype_constant = 5
  INTEGER(i4), PARAMETER :: c_blocktype_array = 6
  INTEGER(i4), PARAMETER :: c_blocktype_run_info = 7
  INTEGER(i4), PARAMETER :: c_blocktype_source = 8
  INTEGER(i4), PARAMETER :: c_blocktype_stitched_tensor = 9
  INTEGER(i4), PARAMETER :: c_blocktype_stitched_material = 10
  INTEGER(i4), PARAMETER :: c_blocktype_stitched_matvar = 11
  INTEGER(i4), PARAMETER :: c_blocktype_stitched_species = 12
  INTEGER(i4), PARAMETER :: c_blocktype_species = 13
  INTEGER(i4), PARAMETER :: c_blocktype_plain_derived = 14
  INTEGER(i4), PARAMETER :: c_blocktype_point_derived = 15
  INTEGER(i4), PARAMETER :: c_blocktype_contiguous_tensor = 16
  INTEGER(i4), PARAMETER :: c_blocktype_contiguous_material = 17
  INTEGER(i4), PARAMETER :: c_blocktype_contiguous_matvar = 18
  INTEGER(i4), PARAMETER :: c_blocktype_contiguous_species = 19
  INTEGER(i4), PARAMETER :: c_blocktype_cpu_split = 20
  INTEGER(i4), PARAMETER :: c_blocktype_stitched_obstacle_group = 21
  INTEGER(i4), PARAMETER :: c_blocktype_unstructured_mesh = 22
  INTEGER(i4), PARAMETER :: c_blocktype_stitched = 23
  INTEGER(i4), PARAMETER :: c_blocktype_contiguous = 24
  INTEGER(i4), PARAMETER :: c_blocktype_lagrangian_mesh = 25
  INTEGER(i4), PARAMETER :: c_blocktype_max = 25

  INTEGER(i4), PARAMETER :: c_datatype_null = 0
  INTEGER(i4), PARAMETER :: c_datatype_integer4 = 1
  INTEGER(i4), PARAMETER :: c_datatype_integer8 = 2
  INTEGER(i4), PARAMETER :: c_datatype_real4 = 3
  INTEGER(i4), PARAMETER :: c_datatype_real8 = 4
  INTEGER(i4), PARAMETER :: c_datatype_real16 = 5
  INTEGER(i4), PARAMETER :: c_datatype_character = 6
  INTEGER(i4), PARAMETER :: c_datatype_logical = 7
  INTEGER(i4), PARAMETER :: c_datatype_other = 8
  INTEGER(i4), PARAMETER :: c_datatype_max = 8

  INTEGER(i4), PARAMETER :: c_geometry_null = 0
  INTEGER(i4), PARAMETER :: c_geometry_cartesian = 1
  INTEGER(i4), PARAMETER :: c_geometry_cylindrical = 2
  INTEGER(i4), PARAMETER :: c_geometry_spherical = 3

  ! c_dimension_irrelevant is used where the dimensionality isn't needed, as
  ! with point variables still keep dimensionality as a common quantity
  ! because other than this, they really are very alike
  INTEGER(i4), PARAMETER :: c_dimension_irrelevant = 0
  INTEGER(i4), PARAMETER :: c_dimension_1d = 1
  INTEGER(i4), PARAMETER :: c_dimension_2d = 2
  INTEGER(i4), PARAMETER :: c_dimension_3d = 3

  INTEGER(i4), PARAMETER :: c_stagger_cell_centre = 0
  INTEGER(i4), PARAMETER :: c_stagger_face_x = 1
  INTEGER(i4), PARAMETER :: c_stagger_face_y = 2
  INTEGER(i4), PARAMETER :: c_stagger_face_z = 4
  INTEGER(i4), PARAMETER :: c_stagger_edge_x = &
      c_stagger_face_y + c_stagger_face_z
  INTEGER(i4), PARAMETER :: c_stagger_edge_y = &
      c_stagger_face_x + c_stagger_face_z
  INTEGER(i4), PARAMETER :: c_stagger_edge_z = &
      c_stagger_face_x + c_stagger_face_y
  INTEGER(i4), PARAMETER :: c_stagger_vertex = &
      c_stagger_face_x + c_stagger_face_y + c_stagger_face_z

  INTEGER(i4), PARAMETER :: sdf_version = 1, sdf_revision = 1

  INTEGER(i4), PARAMETER :: soi4 = 4 ! Size of 4-byte integer
  INTEGER(i4), PARAMETER :: soi8 = 8 ! Size of 8-byte integer
  INTEGER(i4), PARAMETER :: sof4 = 4 ! Size of 4-byte real
  INTEGER(i4), PARAMETER :: sof8 = 8 ! Size of 8-byte real

  ! header length (including padding) - must be updated if sdf_write_header
  ! changes
  INTEGER, PARAMETER :: c_header_length = 11 * soi4 + 2 * soi8 + sof8 + 12 &
      + c_id_length

  ! summary offset - must be updated if sdf_write_header changes
  INTEGER(i4), PARAMETER :: c_summary_offset = 4 + 3 * soi4 + c_id_length + soi8

  INTEGER(i4), PARAMETER :: c_endianness = 16911887

  INTEGER(KIND=MPI_OFFSET_KIND), PARAMETER :: c_off0 = 0
  INTEGER, PARAMETER :: max_mpi_error_codes = 20
  INTEGER, PARAMETER :: mpi_error_codes(max_mpi_error_codes) = (/ &
      MPI_ERR_ACCESS, MPI_ERR_AMODE, MPI_ERR_BAD_FILE, MPI_ERR_CONVERSION, &
      MPI_ERR_DUP_DATAREP, MPI_ERR_FILE, MPI_ERR_FILE_EXISTS, &
      MPI_ERR_FILE_IN_USE, MPI_ERR_INFO, MPI_ERR_INFO_KEY, MPI_ERR_INFO_NOKEY, &
      MPI_ERR_INFO_VALUE, MPI_ERR_IO, MPI_ERR_NOT_SAME, MPI_ERR_NO_SPACE, &
      MPI_ERR_NO_SUCH_FILE, MPI_ERR_QUOTA, MPI_ERR_READ_ONLY, &
      MPI_ERR_UNSUPPORTED_DATAREP, MPI_ERR_UNSUPPORTED_OPERATION /)

  INTEGER, PARAMETER :: c_err_success = 0
  INTEGER, PARAMETER :: c_err_access = 1
  INTEGER, PARAMETER :: c_err_amode = 2
  INTEGER, PARAMETER :: c_err_bad_file = 3
  INTEGER, PARAMETER :: c_err_conversion = 4
  INTEGER, PARAMETER :: c_err_dup_datarep = 5
  INTEGER, PARAMETER :: c_err_file = 6
  INTEGER, PARAMETER :: c_err_file_exists = 7
  INTEGER, PARAMETER :: c_err_file_in_use = 8
  INTEGER, PARAMETER :: c_err_info = 9
  INTEGER, PARAMETER :: c_err_info_key = 10
  INTEGER, PARAMETER :: c_err_info_nokey = 11
  INTEGER, PARAMETER :: c_err_info_value = 12
  INTEGER, PARAMETER :: c_err_io = 13
  INTEGER, PARAMETER :: c_err_not_same = 14
  INTEGER, PARAMETER :: c_err_no_space = 15
  INTEGER, PARAMETER :: c_err_no_such_file = 16
  INTEGER, PARAMETER :: c_err_quota = 17
  INTEGER, PARAMETER :: c_err_read_only = 18
  INTEGER, PARAMETER :: c_err_unsupported_datarep = 19
  INTEGER, PARAMETER :: c_err_unsupported_operation = 20
  INTEGER, PARAMETER :: c_err_unknown = 21
  INTEGER, PARAMETER :: c_err_unsupported_file = 22

  CHARACTER(LEN=*), PARAMETER :: c_blocktypes_char(-1:c_blocktype_max) = (/ &
      'SDF_BLOCKTYPE_SCRUBBED               ', &
      'SDF_BLOCKTYPE_NULL                   ', &
      'SDF_BLOCKTYPE_PLAIN_MESH             ', &
      'SDF_BLOCKTYPE_POINT_MESH             ', &
      'SDF_BLOCKTYPE_PLAIN_VARIABLE         ', &
      'SDF_BLOCKTYPE_POINT_VARIABLE         ', &
      'SDF_BLOCKTYPE_CONSTANT               ', &
      'SDF_BLOCKTYPE_ARRAY                  ', &
      'SDF_BLOCKTYPE_RUN_INFO               ', &
      'SDF_BLOCKTYPE_SOURCE                 ', &
      'SDF_BLOCKTYPE_STITCHED_TENSOR        ', &
      'SDF_BLOCKTYPE_STITCHED_MATERIAL      ', &
      'SDF_BLOCKTYPE_STITCHED_MATVAR        ', &
      'SDF_BLOCKTYPE_STITCHED_SPECIES       ', &
      'SDF_BLOCKTYPE_SPECIES                ', &
      'SDF_BLOCKTYPE_PLAIN_DERIVED          ', &
      'SDF_BLOCKTYPE_POINT_DERIVED          ', &
      'SDF_BLOCKTYPE_MULTI_TENSOR           ', &
      'SDF_BLOCKTYPE_MULTI_MATERIAL         ', &
      'SDF_BLOCKTYPE_MULTI_MATVAR           ', &
      'SDF_BLOCKTYPE_MULTI_SPECIES          ', &
      'SDF_BLOCKTYPE_CPU_SPLIT              ', &
      'SDF_BLOCKTYPE_STITCHED_OBSTACLE_GROUP', &
      'SDF_BLOCKTYPE_UNSTRUCTURED_MESH      ', &
      'SDF_BLOCKTYPE_STITCHED               ', &
      'SDF_BLOCKTYPE_CONTIGUOUS             ', &
      'SDF_BLOCKTYPE_LAGRANGIAN_MESH        ' /)

  CHARACTER(LEN=*), PARAMETER :: c_datatypes_char(0:c_datatype_max) = (/ &
      'SDF_DATATYPE_NULL     ', &
      'SDF_DATATYPE_INTEGER4 ', &
      'SDF_DATATYPE_INTEGER8 ', &
      'SDF_DATATYPE_REAL4    ', &
      'SDF_DATATYPE_REAL8    ', &
      'SDF_DATATYPE_REAL16   ', &
      'SDF_DATATYPE_CHARACTER', &
      'SDF_DATATYPE_LOGICAL  ', &
      'SDF_DATATYPE_OTHER    ' /)

CONTAINS

  SUBROUTINE sdf_get_next_block(h)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: next

    IF (ASSOCIATED(h%blocklist)) THEN
      IF (.NOT. ASSOCIATED(h%current_block)) THEN
        h%current_block => h%blocklist
        RETURN
      ELSE IF (ASSOCIATED(h%current_block%next_block)) THEN
        h%current_block => h%current_block%next_block
        RETURN
      ELSE
        ALLOCATE(h%current_block%next_block)
        CALL initialise_block_type(h%current_block%next_block)
        next => h%current_block%next_block
        next%block_start = h%current_block%next_block_location
      ENDIF
    ELSE
      ALLOCATE(h%blocklist)
      CALL initialise_block_type(h%blocklist)
      next => h%blocklist
      next%block_start = h%summary_location
    ENDIF

    next%done_header = .FALSE.
    next%done_info = .FALSE.
    next%done_data = .FALSE.
    NULLIFY(next%run)
    NULLIFY(next%next_block)
    h%current_block => next

  END SUBROUTINE sdf_get_next_block



  SUBROUTINE sdf_seek_start(h)

    TYPE(sdf_file_handle) :: h

    h%current_block => h%blocklist

  END SUBROUTINE sdf_seek_start



  FUNCTION sdf_find_block(h, b, block_id) RESULT(found)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: b
    CHARACTER(LEN=*), INTENT(IN) :: block_id
    INTEGER :: i
    LOGICAL :: found, use_truncated

    use_truncated = (LEN_TRIM(block_id) .GT. c_id_length)

    found = .TRUE.
    b => h%blocklist
    DO i = 1,h%nblocks
      IF (sdf_string_equal(block_id, b%id)) RETURN
      IF (use_truncated .AND. b%truncated_id) THEN
        IF (sdf_string_equal(block_id, b%long_id)) RETURN
      ENDIF
      b => b%next_block
    ENDDO

    found = .FALSE.
    NULLIFY(b)

  END FUNCTION sdf_find_block



  FUNCTION sdf_find_block_by_id(h, block_id) RESULT(found)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: b
    CHARACTER(LEN=*), INTENT(IN) :: block_id
    LOGICAL :: found

    found = sdf_find_block(h, b, block_id)
    IF (found) h%current_block => b

  END FUNCTION sdf_find_block_by_id



  FUNCTION sdf_get_block_id(h, long_id, block_id) RESULT(found)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: long_id
    CHARACTER(LEN=*), INTENT(OUT) :: block_id
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: found

    found = sdf_find_block(h, b, long_id)
    IF (found) THEN
      CALL safe_copy_string(b%id, block_id)
    ELSE
      CALL safe_copy_string(long_id, block_id)
    ENDIF

  END FUNCTION sdf_get_block_id



  FUNCTION sdf_seek_block(h, block_id) RESULT(found)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: block_id
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: found

    found = sdf_find_block(h, b, block_id)
    IF (found) h%current_block => b

  END FUNCTION sdf_seek_block



  FUNCTION sdf_get_data_location(h) RESULT(data_location)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8) :: data_location

    data_location = h%current_block%data_location

  END FUNCTION sdf_get_data_location



  SUBROUTINE sdf_set_data_location(h, data_location)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(IN) :: data_location

    h%data_location = data_location

  END SUBROUTINE sdf_set_data_location



  FUNCTION sdf_string_equal(str1, str2) RESULT(equal)

    CHARACTER(LEN=*), INTENT(IN) :: str1, str2
    INTEGER :: len1, len2
    LOGICAL :: equal

    len1 = LEN_TRIM(str1)
    len2 = LEN_TRIM(str2)

    IF (len1 .GT. 0) THEN
      IF (IACHAR(str1(len1:len1)) .EQ. 0) len1 = len1 - 1
    ENDIF
    IF (len2 .GT. 0) THEN
      IF (IACHAR(str2(len2:len2)) .EQ. 0) len2 = len2 - 1
    ENDIF

    IF (len1 .NE. len2) THEN
      equal = .FALSE.
      RETURN
    ENDIF

    equal = (str1(1:len1) .EQ. str2(1:len1))

  END FUNCTION sdf_string_equal



  SUBROUTINE safe_copy_string(s1, s2)

    CHARACTER(LEN=*), INTENT(IN) :: s1
    CHARACTER(LEN=*), INTENT(OUT) :: s2
    INTEGER :: len1, len2, olen, i

    len1 = LEN_TRIM(s1)
    len2 = LEN(s2)
    olen = MIN(len1,len2)
    IF (olen .GT. 0) THEN
      s2(1:olen) = s1(1:olen)
      DO i = olen+1,len2
        s2(i:i) = ACHAR(0)
      ENDDO
    ELSE
      DO i = 1,len2
        s2(i:i) = ACHAR(0)
      ENDDO
    ENDIF

  END SUBROUTINE safe_copy_string



  SUBROUTINE safe_copy_id(h, id, new_id)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id
    CHARACTER(LEN=c_id_length), INTENT(OUT) :: new_id

    IF (LEN_TRIM(id) .GT. c_id_length) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*, '*** WARNING ***'
        PRINT*, 'SDF ID string "' // TRIM(id) // '" was truncated.'
      ENDIF
    ENDIF

    CALL safe_copy_string(id, new_id)

  END SUBROUTINE safe_copy_id



  SUBROUTINE safe_copy_unique_id(h, b, id)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: b
    CHARACTER(LEN=*), INTENT(IN) :: id
    LOGICAL :: found
    CHARACTER(LEN=*), PARAMETER :: numbers = '0123456789'
    INTEGER :: i, n, num, pos
    TYPE(sdf_block_type), POINTER :: tmp

    IF (LEN_TRIM(id) .GT. c_id_length) THEN
      b%truncated_id = .TRUE.
      CALL safe_copy_string(id, b%long_id)
      IF (LEN_TRIM(id) .GT. c_long_id_length) THEN
        IF (h%rank .EQ. h%rank_master) THEN
          PRINT*, '*** WARNING ***'
          PRINT*, 'SDF ID string "' // TRIM(id) // '" was truncated.'
        ENDIF
      ENDIF
    ENDIF

    CALL safe_copy_string(id, b%id)
    found = sdf_find_block(h, tmp, b%id)
    i = 1
    DO WHILE(found)
      num = i

      ! Generate a new ID by adding ASCII digits to the end.
      pos = c_id_length
      n = MOD(num,10)
      b%id(pos:pos) = numbers(n+1:n+1)
      num = num / 10
      DO WHILE(num .GT. 0)
        pos = pos - 1
        n = MOD(num,10)
        b%id(pos:pos) = numbers(n+1:n+1)
        num = num / 10
      ENDDO

      found = sdf_find_block(h, tmp, b%id)
      i = i + 1
    ENDDO

  END SUBROUTINE safe_copy_unique_id



  SUBROUTINE initialise_block_type(var)

    TYPE(sdf_block_type) :: var

    NULLIFY(var%dim_mults)
    NULLIFY(var%variable_ids)
    NULLIFY(var%dim_labels)
    NULLIFY(var%dim_units)
    NULLIFY(var%material_names)
    NULLIFY(var%run)
    NULLIFY(var%next_block)
    var%done_header = .FALSE.
    var%done_info = .FALSE.
    var%done_data = .FALSE.
    var%truncated_id = .FALSE.
    var%data_location = 0
    var%blocktype = c_blocktype_null

  END SUBROUTINE initialise_block_type



  SUBROUTINE deallocate_block_type(var)

    TYPE(sdf_block_type) :: var

    IF (ASSOCIATED(var%dim_mults)) DEALLOCATE(var%dim_mults)
    IF (ASSOCIATED(var%variable_ids)) DEALLOCATE(var%variable_ids)
    IF (ASSOCIATED(var%dim_labels)) DEALLOCATE(var%dim_labels)
    IF (ASSOCIATED(var%dim_units)) DEALLOCATE(var%dim_units)
    IF (ASSOCIATED(var%material_names)) DEALLOCATE(var%material_names)
    IF (ASSOCIATED(var%run)) DEALLOCATE(var%run)

    CALL initialise_block_type(var)

  END SUBROUTINE deallocate_block_type



  SUBROUTINE initialise_file_handle(var)

    TYPE(sdf_file_handle) :: var

    NULLIFY(var%buffer)
    NULLIFY(var%blocklist)
    NULLIFY(var%current_block)
    ! Set filehandle to -1 to show that the file is closed
    var%filehandle = -1
    var%string_length = c_max_string_length
    var%done_header = .FALSE.
    var%restart_flag = .FALSE.
    var%other_domains = .FALSE.
    var%writing = .FALSE.
    var%handled_error = .FALSE.
    var%nblocks = 0
    var%error_code = 0
    var%errhandler = 0

  END SUBROUTINE initialise_file_handle



  SUBROUTINE deallocate_file_handle(var)

    TYPE(sdf_file_handle) :: var
    INTEGER :: errcode, i

    IF (ASSOCIATED(var%buffer)) DEALLOCATE(var%buffer)

    IF (var%errhandler .NE. 0) THEN
      CALL MPI_ERRHANDLER_FREE(var%errhandler, errcode)
    ENDIF

    DO i = 1, max_handles
      IF (sdf_handles(i)%filehandle .EQ. var%filehandle) THEN
        sdf_handles(i)%filehandle = 0
        EXIT
      ENDIF
    ENDDO

    CALL initialise_file_handle(var)

  END SUBROUTINE deallocate_file_handle



  FUNCTION map_error_code(error_code) RESULT(errcode)

    INTEGER, INTENT(IN) :: error_code
    INTEGER :: errcode, i

    errcode = c_err_unknown
    DO i = 1, max_mpi_error_codes
      IF (error_code .EQ. mpi_error_codes(i)) THEN
        errcode = i
        RETURN
      ENDIF
    ENDDO

  END FUNCTION map_error_code



  SUBROUTINE error_handler(filehandle, error_code)

    INTEGER :: filehandle, error_code
    TYPE(sdf_file_handle), POINTER :: h
    INTEGER :: i

    DO i = 1, max_handles
      IF (sdf_handles(i)%filehandle .EQ. filehandle) THEN
        h => sdf_handles(i)%handle
        IF (.NOT.h%handled_error) THEN
          h%error_code = map_error_code(error_code) + 64 * h%nblocks
          h%handled_error = .TRUE.
        ENDIF
        EXIT
      ENDIF
    ENDDO

  END SUBROUTINE error_handler

END MODULE sdf_common
