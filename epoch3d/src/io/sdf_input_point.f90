MODULE sdf_input_point

  USE sdf_common
  USE sdf_input
  USE mpi

  IMPLICIT NONE

  INTERFACE sdf_read_srl_point_mesh
    MODULE PROCEDURE &
        sdf_read_srl_1d_pt_mesh_array, &
        sdf_read_srl_2d_pt_mesh_array, &
        sdf_read_srl_3d_pt_mesh_array
  END INTERFACE sdf_read_srl_point_mesh

  INTERFACE sdf_read_srl_point_variable
    MODULE PROCEDURE &
        sdf_read_srl_pt_var_int_array, &
        sdf_read_srl_pt_var_flt_array
  END INTERFACE sdf_read_srl_point_variable

CONTAINS

  ! Mesh loading functions

  SUBROUTINE sdf_info_init(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: errcode, ierr
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    h%current_location = b%block_start + h%block_header_length

    IF (.NOT. ALLOCATED(h%buffer)) THEN
      CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
          MPI_BYTE, "native", MPI_INFO_NULL, errcode)
    ENDIF

  END SUBROUTINE sdf_info_init



  SUBROUTINE sdf_read_point_mesh_info(h, npoints, geometry, extents, &
      dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(OUT), OPTIONAL :: npoints
    INTEGER, INTENT(OUT), OPTIONAL :: geometry
    REAL(num), DIMENSION(:), INTENT(OUT), OPTIONAL :: extents
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(OUT), OPTIONAL :: dim_mults
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

    CALL sdf_info_init(h)

    b => h%current_block
    IF (.NOT. b%done_info) THEN
      ALLOCATE(b%dim_labels(b%ndims))
      ALLOCATE(b%dim_units(b%ndims))
      ALLOCATE(b%dim_mults(b%ndims))

      CALL read_entry_array_real8(h, b%dim_mults, b%ndims)

      DO i = 1,b%ndims
        CALL read_entry_id(h, b%dim_labels(i))
      ENDDO

      DO i = 1,b%ndims
        CALL read_entry_id(h, b%dim_units(i))
      ENDDO

      CALL read_entry_int4(h, b%geometry)

      CALL read_entry_array_real8(h, b%extents, 2*b%ndims)

      CALL read_entry_int8(h, b%npoints)
    ENDIF

    IF (PRESENT(geometry)) geometry = b%geometry
    IF (PRESENT(npoints)) npoints = b%npoints
    IF (PRESENT(extents)) extents(1:2*b%ndims) = b%extents(1:2*b%ndims)
    IF (PRESENT(dim_mults)) dim_mults = b%dim_mults
    IF (PRESENT(dim_labels)) THEN
      DO i = 1,b%ndims
        clen = MIN(LEN(dim_labels(i)),c_id_length)
        dim_labels(i)(1:clen) = b%dim_labels(i)(1:clen)
      ENDDO
    ENDIF
    IF (PRESENT(dim_units)) THEN
      DO i = 1,b%ndims
        clen = MIN(LEN(dim_units(i)),c_id_length)
        dim_units(i)(1:clen) = b%dim_units(i)(1:clen)
      ENDDO
    ENDIF

    h%current_location = b%block_start + h%block_header_length
    b%done_info = .TRUE.

  END SUBROUTINE sdf_read_point_mesh_info



  SUBROUTINE sdf_read_srl_1d_pt_mesh_array(h, x)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:), INTENT(OUT) :: x
    INTEGER :: errcode, npoints
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_point_mesh_info(h)

    h%current_location = b%data_location

    ! Read the real data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        MPI_BYTE, "native", MPI_INFO_NULL, errcode)

    npoints = b%npoints
    CALL MPI_FILE_READ_ALL(h%filehandle, x, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_srl_1d_pt_mesh_array



  SUBROUTINE sdf_read_srl_2d_pt_mesh_array(h, x, y)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:), INTENT(OUT) :: x, y
    INTEGER :: errcode, npoints
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_point_mesh_info(h)

    h%current_location = b%data_location

    ! Read the real data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        MPI_BYTE, "native", MPI_INFO_NULL, errcode)

    npoints = b%npoints
    CALL MPI_FILE_READ_ALL(h%filehandle, x, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, y, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_srl_2d_pt_mesh_array



  SUBROUTINE sdf_read_srl_3d_pt_mesh_array(h, x, y, z)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:), INTENT(OUT) :: x, y, z
    INTEGER :: errcode, npoints
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_point_mesh_info(h)

    h%current_location = b%data_location

    ! Read the real data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        MPI_BYTE, "native", MPI_INFO_NULL, errcode)

    npoints = b%npoints
    CALL MPI_FILE_READ_ALL(h%filehandle, x, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, y, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, z, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_srl_3d_pt_mesh_array



  SUBROUTINE sdf_read_point_mesh(h, npoint_local, distribution, iterator)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(IN) :: npoint_local
    INTEGER, INTENT(IN) :: distribution
    INTEGER(i8) :: npoint_remain
    INTEGER :: direction, errcode, npoint_this_it
    LOGICAL :: start
    REAL(num), DIMENSION(:), ALLOCATABLE :: array
    REAL(num) :: ret
    TYPE(sdf_block_type), POINTER :: b

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, direction)
        USE sdf_common
        REAL(num) :: iterator
        REAL(num), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN) :: direction
      END FUNCTION iterator
    END INTERFACE

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_point_mesh_info(h)

    h%current_location = b%data_location

    ALLOCATE(array(1:npoint_per_iteration))

    DO direction = 1, b%ndims
      start = .TRUE.
      npoint_remain = npoint_local
      npoint_this_it = MIN(npoint_remain, npoint_per_iteration)

      CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, &
          b%mpitype, distribution, "native", MPI_INFO_NULL, errcode)

      DO WHILE (npoint_this_it .GT. 0)
        CALL MPI_FILE_READ(h%filehandle, array, npoint_this_it, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        npoint_remain = npoint_remain - npoint_this_it
        ret = iterator(array, npoint_this_it, start, direction)
        start = .FALSE.
        npoint_this_it = MIN(npoint_remain, npoint_per_iteration)
      ENDDO

      h%current_location = h%current_location + b%npoints * h%sof
    ENDDO

    DEALLOCATE(array)

    h%current_location = b%next_block_location
    b%block_start = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_point_mesh



  ! Variable loading functions

  SUBROUTINE sdf_read_point_variable_info(h, npoints, mesh_id, units, mult)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(OUT), OPTIONAL :: npoints
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: mesh_id, units
    REAL(num), INTENT(OUT), OPTIONAL :: mult
    INTEGER :: clen
    TYPE(sdf_block_type), POINTER :: b

    ! Metadata is
    ! - mult      REAL(r8)
    ! - units     CHARACTER(id_length)
    ! - meshid    CHARACTER(id_length)
    ! - npoints   INTEGER(i8)

    CALL sdf_info_init(h)

    b => h%current_block
    IF (.NOT. b%done_info) THEN
      CALL read_entry_real8(h, b%mult)

      CALL read_entry_id(h, b%units)

      CALL read_entry_id(h, b%mesh_id)

      CALL read_entry_int8(h, b%npoints)
    ENDIF

    IF (PRESENT(npoints)) npoints = b%npoints
    IF (PRESENT(mult)) mult = b%mult
    IF (PRESENT(units)) THEN
      clen = MIN(LEN(units),c_id_length)
      units(1:clen) = b%units(1:clen)
    ENDIF
    IF (PRESENT(mesh_id)) THEN
      clen = MIN(LEN(mesh_id),c_id_length)
      mesh_id(1:clen) = b%mesh_id(1:clen)
    ENDIF

    h%current_location = b%block_start + h%block_header_length
    b%done_info = .TRUE.

  END SUBROUTINE sdf_read_point_variable_info



  SUBROUTINE sdf_read_point_variable(h, npoint_local, distribution, iterator)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(IN) :: npoint_local
    INTEGER, INTENT(IN) :: distribution
    INTEGER(i8) :: npoint_remain
    INTEGER :: errcode, npoint_this_it
    LOGICAL :: start
    REAL(num), DIMENSION(:), ALLOCATABLE :: array
    REAL(num) :: ret
    TYPE(sdf_block_type), POINTER :: b

    INTERFACE
      FUNCTION iterator(array, npoint_it, start)
        USE sdf_common
        REAL(num) :: iterator
        REAL(num), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
      END FUNCTION iterator
    END INTERFACE

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_point_variable_info(h)

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, &
        b%mpitype, distribution, "native", MPI_INFO_NULL, errcode)

    start = .TRUE.
    ALLOCATE(array(1:npoint_per_iteration))
    npoint_remain = npoint_local
    npoint_this_it = MIN(npoint_remain, npoint_per_iteration)

    DO WHILE (npoint_this_it .GT. 0)
      npoint_this_it = MIN(npoint_remain, npoint_per_iteration)
      CALL MPI_FILE_READ(h%filehandle, array, npoint_this_it, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)

      npoint_remain = npoint_remain - npoint_this_it
      ret = iterator(array, npoint_this_it, start)
      start = .FALSE.
    ENDDO

    DEALLOCATE(array)

    h%current_location = b%next_block_location
    b%block_start = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_point_variable



  SUBROUTINE sdf_read_srl_pt_var_flt_array(h, array)

    TYPE(sdf_file_handle) :: h
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER :: errcode, npoints
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_point_variable_info(h)

    h%current_location = b%data_location

    ! Read the real data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        MPI_BYTE, "native", MPI_INFO_NULL, errcode)

    npoints = b%npoints
    CALL MPI_FILE_READ_ALL(h%filehandle, array, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_srl_pt_var_flt_array



  SUBROUTINE sdf_read_srl_pt_var_int_array(h, array)

    TYPE(sdf_file_handle) :: h
    INTEGER, DIMENSION(:), INTENT(OUT) :: array
    INTEGER :: errcode, npoints
    TYPE(sdf_block_type), POINTER :: b

    IF (.NOT.ASSOCIATED(h%current_block)) THEN
      IF (h%rank .EQ. h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      ENDIF
      RETURN
    ENDIF

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_point_variable_info(h)

    h%current_location = b%data_location

    ! Read the real data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        MPI_BYTE, "native", MPI_INFO_NULL, errcode)

    npoints = b%npoints
    CALL MPI_FILE_READ_ALL(h%filehandle, array, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_srl_pt_var_int_array

END MODULE sdf_input_point
