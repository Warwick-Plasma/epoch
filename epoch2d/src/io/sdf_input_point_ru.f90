MODULE sdf_input_point_ru

  USE sdf_input_ru

  IMPLICIT NONE

CONTAINS

  ! Mesh loading functions

  SUBROUTINE read_point_mesh_info_ru(h, npoints, geometry, species_id)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(OUT), OPTIONAL :: npoints
    INTEGER, INTENT(OUT), OPTIONAL :: geometry
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: species_id
    INTEGER :: i, clen
    TYPE(sdf_block_type), POINTER :: b

    ! Metadata is
    ! - mults     REAL(r8), DIMENSION(ndims)
    ! - labels    CHARACTER(id_length), DIMENSION(ndims)
    ! - units     CHARACTER(id_length), DIMENSION(ndims)
    ! - geometry  INTEGER(i4)
    ! - minval    REAL(r8), DIMENSION(ndims)
    ! - maxval    REAL(r8), DIMENSION(ndims)
    ! - npoints   INTEGER(i8)
    ! - speciesid CHARACTER(id_length)

    IF (sdf_info_init(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) THEN
      ALLOCATE(b%dim_labels(b%ndims))
      ALLOCATE(b%dim_units(b%ndims))
      ALLOCATE(b%dim_mults(b%ndims))

      CALL read_entry_array_real8(h, b%dim_mults, INT(b%ndims))

      DO i = 1,b%ndims
        CALL read_entry_id(h, b%dim_labels(i))
      ENDDO

      DO i = 1,b%ndims
        CALL read_entry_id(h, b%dim_units(i))
      ENDDO

      CALL read_entry_int4(h, b%geometry)

      CALL read_entry_array_real8(h, b%extents, INT(2*b%ndims))

      CALL read_entry_int8(h, b%npoints)

      ! species_id field added in version 1.3
      IF (1000 * h%file_version + h%file_revision .GT. 1002) THEN
        CALL read_entry_id(h, b%species_id)
      ELSE
        CALL safe_copy_id(h, '__unknown__', b%species_id)
      ENDIF
    ENDIF

    IF (PRESENT(geometry)) geometry = b%geometry
    IF (PRESENT(npoints)) npoints = b%npoints
    IF (PRESENT(species_id)) THEN
      clen = MIN(LEN(species_id),INT(c_id_length))
      species_id(1:clen) = b%species_id(1:clen)
    ENDIF

    h%current_location = b%block_start + h%block_header_length
    b%done_info = .TRUE.

  END SUBROUTINE read_point_mesh_info_ru



  SUBROUTINE read_point_mesh_info_i4_ru(h, npoints, geometry, species_id)

    TYPE(sdf_file_handle) :: h
    INTEGER(i4), INTENT(OUT) :: npoints
    INTEGER, INTENT(OUT), OPTIONAL :: geometry
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: species_id
    INTEGER(i8) :: npoints8

    CALL read_point_mesh_info_ru(h, npoints8, geometry, species_id)
    npoints = INT(npoints8,i4)

  END SUBROUTINE read_point_mesh_info_i4_ru



  ! Variable loading functions

  SUBROUTINE read_point_variable_info_ru(h, npoints, mesh_id, units, species_id)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(OUT), OPTIONAL :: npoints
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: mesh_id, units, species_id
    INTEGER :: clen
    TYPE(sdf_block_type), POINTER :: b

    ! Metadata is
    ! - mult      REAL(r8)
    ! - units     CHARACTER(id_length)
    ! - meshid    CHARACTER(id_length)
    ! - npoints   INTEGER(i8)
    ! - speciesid CHARACTER(id_length)

    IF (sdf_info_init(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) THEN
      CALL read_entry_real8(h, b%mult)

      CALL read_entry_id(h, b%units)

      CALL read_entry_id(h, b%mesh_id)

      CALL read_entry_int8(h, b%npoints)

      ! species_id field added in version 1.3
      IF (1000 * h%file_version + h%file_revision .GT. 1002) THEN
        CALL read_entry_id(h, b%species_id)
      ELSE
        CALL safe_copy_id(h, '__unknown__', b%species_id)
      ENDIF
    ENDIF

    IF (PRESENT(npoints)) npoints = b%npoints
    IF (PRESENT(units)) THEN
      clen = MIN(LEN(units),INT(c_id_length))
      units(1:clen) = b%units(1:clen)
    ENDIF
    IF (PRESENT(mesh_id)) THEN
      clen = MIN(LEN(mesh_id),INT(c_id_length))
      mesh_id(1:clen) = b%mesh_id(1:clen)
    ENDIF
    IF (PRESENT(species_id)) THEN
      clen = MIN(LEN(species_id),INT(c_id_length))
      species_id(1:clen) = b%species_id(1:clen)
    ENDIF

    h%current_location = b%block_start + h%block_header_length
    b%done_info = .TRUE.

  END SUBROUTINE read_point_variable_info_ru



  SUBROUTINE read_srl_pt_var_int_array(h, array)

    TYPE(sdf_file_handle) :: h
    INTEGER, DIMENSION(:), INTENT(OUT) :: array
    INTEGER :: errcode, npoints
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL read_point_variable_info_ru(h)

    h%current_location = b%data_location

    ! Read the real data

    npoints = INT(b%npoints)
    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_READ_AT(h%filehandle, h%current_location, array, npoints, &
          b%mpitype, MPI_STATUS_IGNORE, errcode)
    ENDIF

    CALL MPI_BCAST(array, npoints, b%mpitype, h%rank_master, h%comm, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_srl_pt_var_int_array



  SUBROUTINE read_srl_pt_var_logical_array(h, array)

    TYPE(sdf_file_handle) :: h
    LOGICAL, DIMENSION(:), INTENT(OUT) :: array
    INTEGER :: errcode, npoints, i
    TYPE(sdf_block_type), POINTER :: b
    CHARACTER(LEN=1), DIMENSION(:), ALLOCATABLE :: cvalues

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL read_point_variable_info_ru(h)

    h%current_location = b%data_location

    ! Read the real data

    npoints = INT(b%npoints)
    ALLOCATE(cvalues(npoints))

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_READ_AT(h%filehandle, h%current_location, cvalues, &
          npoints, b%mpitype, MPI_STATUS_IGNORE, errcode)
    ENDIF

    CALL MPI_BCAST(cvalues, npoints, b%mpitype, h%rank_master, h%comm, errcode)

    DO i = 1,npoints
      IF (cvalues(i) .EQ. ACHAR(0)) THEN
        array(i) = .FALSE.
      ELSE
        array(i) = .TRUE.
      ENDIF
    ENDDO

    DEALLOCATE(cvalues)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_srl_pt_var_logical_array



  SUBROUTINE read_point_variable_i4(h, npoint_local, distribution, &
      iterator, param)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(IN) :: npoint_local
    INTEGER, INTENT(IN) :: distribution
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER(i8) :: npoint_remain, npoint_per_it8, npoint_this_it8
    INTEGER :: errcode, npoint_per_it, npoint_this_it
    LOGICAL :: start
    INTEGER(i4), DIMENSION(:), ALLOCATABLE :: array
    INTEGER(i4) :: ret
    TYPE(sdf_block_type), POINTER :: b

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, param)
        USE sdf_common
        INTEGER(i4) :: iterator
        INTEGER(i4), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL read_point_variable_info_ru(h)

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)

    start = .TRUE.
    npoint_per_it8 = MIN(npoint_local, npoint_per_iteration)
    npoint_per_it  = INT(npoint_per_it8)
    ALLOCATE(array(1:npoint_per_it))
    npoint_remain = npoint_local
    npoint_this_it8 = MIN(npoint_remain, npoint_per_it8)
    npoint_this_it  = INT(npoint_this_it8)

    DO WHILE (npoint_this_it .GT. 0)
      npoint_this_it8 = MIN(npoint_remain, npoint_per_it8)
      npoint_this_it  = INT(npoint_this_it8)
      CALL MPI_FILE_READ(h%filehandle, array, npoint_this_it, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)

      npoint_remain = npoint_remain - npoint_this_it8
      ret = iterator(array, npoint_this_it, start, param)
      start = .FALSE.
    ENDDO

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    DEALLOCATE(array)

    h%current_location = b%next_block_location
    b%block_start = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_point_variable_i4



  SUBROUTINE read_point_variable_i8(h, npoint_local, distribution, &
      iterator, param)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(IN) :: npoint_local
    INTEGER, INTENT(IN) :: distribution
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER(i8) :: npoint_remain, npoint_per_it8, npoint_this_it8
    INTEGER :: errcode, npoint_per_it, npoint_this_it
    LOGICAL :: start
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: array
    INTEGER(i8) :: ret
    TYPE(sdf_block_type), POINTER :: b

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, param)
        USE sdf_common
        INTEGER(i8) :: iterator
        INTEGER(i8), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL read_point_variable_info_ru(h)

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)

    start = .TRUE.
    npoint_per_it8 = MIN(npoint_local, npoint_per_iteration)
    npoint_per_it  = INT(npoint_per_it8)
    ALLOCATE(array(1:npoint_per_it))
    npoint_remain = npoint_local
    npoint_this_it8 = MIN(npoint_remain, npoint_per_it8)
    npoint_this_it  = INT(npoint_this_it8)

    DO WHILE (npoint_this_it .GT. 0)
      npoint_this_it8 = MIN(npoint_remain, npoint_per_it8)
      npoint_this_it  = INT(npoint_this_it8)
      CALL MPI_FILE_READ(h%filehandle, array, npoint_this_it, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)

      npoint_remain = npoint_remain - npoint_this_it8
      ret = iterator(array, npoint_this_it, start, param)
      start = .FALSE.
    ENDDO

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    DEALLOCATE(array)

    h%current_location = b%next_block_location
    b%block_start = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_point_variable_i8

END MODULE sdf_input_point_ru
