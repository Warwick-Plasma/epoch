MODULE sdf_input_point_r8

  USE mpi
  USE sdf_common
  USE sdf_input
  USE sdf_input_point_ru

  IMPLICIT NONE

CONTAINS

  ! Mesh loading functions

  SUBROUTINE read_point_mesh_info_r8(h, npoints, geometry, extents, &
      dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(OUT) :: npoints
    INTEGER, INTENT(OUT) :: geometry
    REAL(r8), DIMENSION(:), INTENT(OUT) :: extents
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r8), DIMENSION(:), INTENT(OUT), OPTIONAL :: dim_mults
    INTEGER :: i, clen
    TYPE(sdf_block_type), POINTER :: b

    CALL read_point_mesh_info_ru(h, npoints, geometry)
    b => h%current_block

    extents(1:2*b%ndims) = REAL(b%extents(1:2*b%ndims),r8)
    IF (PRESENT(dim_mults)) dim_mults = REAL(b%dim_mults,r8)
    IF (PRESENT(dim_labels)) THEN
      DO i = 1,b%ndims
        clen = MIN(LEN(dim_labels(i)),INT(c_id_length))
        dim_labels(i)(1:clen) = b%dim_labels(i)(1:clen)
      ENDDO
    ENDIF
    IF (PRESENT(dim_units)) THEN
      DO i = 1,b%ndims
        clen = MIN(LEN(dim_units(i)),INT(c_id_length))
        dim_units(i)(1:clen) = b%dim_units(i)(1:clen)
      ENDDO
    ENDIF

  END SUBROUTINE read_point_mesh_info_r8



  SUBROUTINE read_srl_1d_pt_mesh_array_r8(h, x)

    TYPE(sdf_file_handle) :: h
    REAL(r8), DIMENSION(:), INTENT(OUT) :: x
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
    IF (.NOT. b%done_info) CALL read_point_mesh_info_ru(h)

    h%current_location = b%data_location

    ! Read the real data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        MPI_BYTE, 'native', MPI_INFO_NULL, errcode)

    npoints = INT(b%npoints)
    CALL MPI_FILE_READ_ALL(h%filehandle, x, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_srl_1d_pt_mesh_array_r8



  SUBROUTINE read_srl_2d_pt_mesh_array_r8(h, x, y)

    TYPE(sdf_file_handle) :: h
    REAL(r8), DIMENSION(:), INTENT(OUT) :: x, y
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
    IF (.NOT. b%done_info) CALL read_point_mesh_info_ru(h)

    h%current_location = b%data_location

    ! Read the real data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        MPI_BYTE, 'native', MPI_INFO_NULL, errcode)

    npoints = INT(b%npoints)
    CALL MPI_FILE_READ_ALL(h%filehandle, x, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, y, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_srl_2d_pt_mesh_array_r8



  SUBROUTINE read_srl_3d_pt_mesh_array_r8(h, x, y, z)

    TYPE(sdf_file_handle) :: h
    REAL(r8), DIMENSION(:), INTENT(OUT) :: x, y, z
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
    IF (.NOT. b%done_info) CALL read_point_mesh_info_ru(h)

    h%current_location = b%data_location

    ! Read the real data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        MPI_BYTE, 'native', MPI_INFO_NULL, errcode)

    npoints = INT(b%npoints)
    CALL MPI_FILE_READ_ALL(h%filehandle, x, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, y, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, z, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_srl_3d_pt_mesh_array_r8



  SUBROUTINE read_point_mesh_r8(h, npoint_local, distribution, iterator)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(IN) :: npoint_local
    INTEGER, INTENT(IN) :: distribution
    INTEGER(i8) :: npoint_remain, npoint_per_it8, npoint_this_it8
    INTEGER :: direction, errcode, npoint_per_it, npoint_this_it
    LOGICAL :: start
    REAL(r8), DIMENSION(:), ALLOCATABLE :: array
    REAL(r8) :: ret
    TYPE(sdf_block_type), POINTER :: b

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, direction)
        USE sdf_common
        REAL(r8) :: iterator
        REAL(r8), DIMENSION(:), INTENT(INOUT) :: array
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
    IF (.NOT. b%done_info) CALL read_point_mesh_info_ru(h)

    h%current_location = b%data_location

    npoint_per_it8 = MIN(npoint_local, npoint_per_iteration)
    npoint_per_it  = INT(npoint_per_it8)
    ALLOCATE(array(1:npoint_per_it))

    DO direction = 1, b%ndims
      start = .TRUE.
      npoint_remain = npoint_local
      npoint_this_it8 = MIN(npoint_remain, npoint_per_it8)
      npoint_this_it  = INT(npoint_this_it8)

      CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
          distribution, 'native', MPI_INFO_NULL, errcode)

      DO WHILE (npoint_this_it .GT. 0)
        CALL MPI_FILE_READ(h%filehandle, array, npoint_this_it, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        npoint_remain = npoint_remain - npoint_this_it8
        ret = iterator(array, npoint_this_it, start, direction)
        start = .FALSE.
        npoint_this_it8 = MIN(npoint_remain, npoint_per_it8)
        npoint_this_it  = INT(npoint_this_it8)
      ENDDO

      h%current_location = h%current_location + b%npoints * b%type_size
    ENDDO

    DEALLOCATE(array)

    h%current_location = b%next_block_location
    b%block_start = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_point_mesh_r8



  ! Variable loading functions

  SUBROUTINE read_point_variable_info_r8(h, npoints, mesh_id, units, mult)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(OUT) :: npoints
    CHARACTER(LEN=*), INTENT(OUT) :: mesh_id, units
    REAL(r8), INTENT(OUT) :: mult
    TYPE(sdf_block_type), POINTER :: b

    CALL read_point_variable_info_ru(h, npoints, mesh_id, units)
    b => h%current_block
    mult = REAL(b%mult,r8)

  END SUBROUTINE read_point_variable_info_r8



  SUBROUTINE read_point_variable_r8(h, npoint_local, distribution, iterator)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(IN) :: npoint_local
    INTEGER, INTENT(IN) :: distribution
    INTEGER(i8) :: npoint_remain, npoint_per_it8, npoint_this_it8
    INTEGER :: errcode, npoint_per_it, npoint_this_it
    LOGICAL :: start
    REAL(r8), DIMENSION(:), ALLOCATABLE :: array
    REAL(r8) :: ret
    TYPE(sdf_block_type), POINTER :: b

    INTERFACE
      FUNCTION iterator(array, npoint_it, start)
        USE sdf_common
        REAL(r8) :: iterator
        REAL(r8), DIMENSION(:), INTENT(INOUT) :: array
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
      ret = iterator(array, npoint_this_it, start)
      start = .FALSE.
    ENDDO

    DEALLOCATE(array)

    h%current_location = b%next_block_location
    b%block_start = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_point_variable_r8



  SUBROUTINE read_srl_pt_var_flt_array_r8(h, array)

    TYPE(sdf_file_handle) :: h
    REAL(r8), DIMENSION(:), INTENT(OUT) :: array
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
    IF (.NOT. b%done_info) CALL read_point_variable_info_ru(h)

    h%current_location = b%data_location

    ! Read the real data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        MPI_BYTE, 'native', MPI_INFO_NULL, errcode)

    npoints = INT(b%npoints)
    CALL MPI_FILE_READ_ALL(h%filehandle, array, npoints, b%mpitype, &
        MPI_STATUS_IGNORE, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_srl_pt_var_flt_array_r8

END MODULE sdf_input_point_r8
