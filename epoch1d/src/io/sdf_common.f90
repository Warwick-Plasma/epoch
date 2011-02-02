MODULE sdf_common

  USE sdf_job_info
  USE mpi

  IMPLICIT NONE

  INTEGER, PARAMETER :: i4  = SELECTED_INT_KIND(9)  ! 4-byte 2^31 ~ 10^9
  INTEGER, PARAMETER :: i8  = SELECTED_INT_KIND(18) ! 8-byte 2^63 ~ 10^18

  INTEGER, PARAMETER :: r4  = SELECTED_REAL_KIND(r=30)
  INTEGER, PARAMETER :: r8  = SELECTED_REAL_KIND(r=300)
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(r=3000)

  INTEGER, PARAMETER :: c_maxdims = 4
  INTEGER(i4), PARAMETER :: c_id_length = 32
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
    INTEGER(KIND=MPI_OFFSET_KIND) :: block_start
    INTEGER(i8) :: next_block_location, data_location
    INTEGER(i8) :: nelements, npoints, data_length, info_length
    INTEGER(i4) :: ndims, geometry, datatype, blocktype
    INTEGER(i4) :: mpitype, type_size, stagger
    INTEGER(i4), DIMENSION(c_maxdims) :: dims
    CHARACTER(LEN=1) :: const_value(16)
    CHARACTER(LEN=c_id_length) :: id, units, mesh_id, material_id
    CHARACTER(LEN=c_id_length), POINTER :: variable_ids(:)
    CHARACTER(LEN=c_id_length), POINTER :: dim_labels(:), dim_units(:)
    CHARACTER(LEN=c_max_string_length) :: name, material_name
    CHARACTER(LEN=c_max_string_length), POINTER :: material_names(:)
    LOGICAL :: done_header, done_info, done_data
    TYPE(sdf_run_type), POINTER :: run
    TYPE(sdf_block_type), POINTER :: next_block
  END TYPE sdf_block_type

  TYPE sdf_file_handle
    INTEGER(KIND=MPI_OFFSET_KIND) :: current_location
    REAL(r8) :: time
    INTEGER(i8) :: first_block_location, summary_location, start_location
    INTEGER(i8) :: soi, sof ! large integer to prevent overflow in calculations
    INTEGER(i4) :: endianness, summary_size
    INTEGER(i4) :: block_header_length, string_length, nblocks
    INTEGER(i4) :: file_version, file_revision, code_io_version, step
    INTEGER(i4) :: datatype_real, datatype_integer
    INTEGER(i4) :: mpitype_real, mpitype_integer
    INTEGER :: filehandle, comm, rank, rank_master, default_rank, mode
    LOGICAL :: done_header, restart_flag, other_domains, writing
    CHARACTER(LEN=1), POINTER :: buffer(:)
    CHARACTER(LEN=c_id_length) :: code_name
    TYPE(jobid_type) :: jobid
    TYPE(sdf_block_type), POINTER :: blocklist, current_block
  END TYPE sdf_file_handle

  INTEGER, PARAMETER :: num = KIND(1.d0)

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
  INTEGER(i4), PARAMETER :: c_blocktype_family = 13

  INTEGER(i4), PARAMETER :: c_datatype_null = 0
  INTEGER(i4), PARAMETER :: c_datatype_integer4 = 1
  INTEGER(i4), PARAMETER :: c_datatype_integer8 = 2
  INTEGER(i4), PARAMETER :: c_datatype_real4 = 3
  INTEGER(i4), PARAMETER :: c_datatype_real8 = 4
  INTEGER(i4), PARAMETER :: c_datatype_real16 = 5
  INTEGER(i4), PARAMETER :: c_datatype_character = 6
  INTEGER(i4), PARAMETER :: c_datatype_logical = 7
  INTEGER(i4), PARAMETER :: c_datatype_other = 8

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

  INTEGER(i4), PARAMETER :: sdf_version = 1, sdf_revision = 0

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
    LOGICAL :: found

    found = .TRUE.
    b => h%blocklist
    DO i = 1,h%nblocks
      IF (sdf_string_equal(block_id, b%id)) RETURN
      b => b%next_block
    ENDDO

    found = .FALSE.
    NULLIFY(b)

  END FUNCTION sdf_find_block



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
    INTEGER :: len1, len2, olen

    len1 = LEN_TRIM(s1)
    len2 = LEN(s2)
    olen = MIN(len1,len2)
    IF (olen .GT. 0) THEN
      s2(olen:len2) = ACHAR(0)
      s2(1:olen) = s1(1:olen)
    ELSE
      s2(1:1) = ACHAR(0)
    ENDIF

  END SUBROUTINE safe_copy_string



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

  END SUBROUTINE initialise_block_type



  SUBROUTINE deallocate_block_type(var)

    TYPE(sdf_block_type) :: var

    IF (ASSOCIATED(var%dim_mults)) DEALLOCATE(var%dim_mults)
    IF (ASSOCIATED(var%variable_ids)) DEALLOCATE(var%variable_ids)
    IF (ASSOCIATED(var%dim_labels)) DEALLOCATE(var%dim_labels)
    IF (ASSOCIATED(var%dim_units)) DEALLOCATE(var%dim_units)
    IF (ASSOCIATED(var%material_names)) DEALLOCATE(var%material_names)

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

  END SUBROUTINE initialise_file_handle



  SUBROUTINE deallocate_file_handle(var)

    TYPE(sdf_file_handle) :: var

    IF (ASSOCIATED(var%buffer)) DEALLOCATE(var%buffer)

    CALL initialise_file_handle(var)

  END SUBROUTINE deallocate_file_handle

END MODULE sdf_common
