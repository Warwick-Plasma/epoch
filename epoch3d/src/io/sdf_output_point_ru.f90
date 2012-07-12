MODULE sdf_output_point_ru

  USE mpi
  USE sdf_common
  USE sdf_output

  IMPLICIT NONE

CONTAINS

  SUBROUTINE write_point_mesh_meta_r8(h, id, name, dim_labels, dim_units, &
      dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r8), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER :: ndims
    TYPE(sdf_block_type), POINTER :: b
    INTEGER :: i, errcode

    b => h%current_block

    b%blocktype = c_blocktype_point_mesh
    ndims = b%ndims
    b%nelements = b%ndims * b%npoints

    ! Metadata is
    ! - mults     REAL(r8), DIMENSION(ndims)
    ! - labels    CHARACTER(id_length), DIMENSION(ndims)
    ! - units     CHARACTER(id_length), DIMENSION(ndims)
    ! - geometry  INTEGER(i4)
    ! - minval    REAL(r8), DIMENSION(ndims)
    ! - maxval    REAL(r8), DIMENSION(ndims)
    ! - npoints   INTEGER(i8)

    b%info_length = h%block_header_length + soi4 + soi8 &
        + (3 * ndims) * sof8 + 2 * ndims * c_id_length
    b%data_length = b%nelements * b%type_size

    ! Write header

    IF (PRESENT(id)) THEN
      ALLOCATE(b%dim_labels(ndims), b%dim_units(ndims), b%dim_mults(ndims))

      IF (PRESENT(dim_labels)) THEN
        DO i = 1,ndims
          CALL safe_copy_string(dim_labels(i), b%dim_labels(i))
        ENDDO
      ELSE
        IF (ndims .GE. 1) CALL safe_copy_string('X', b%dim_labels(1))
        IF (ndims .GE. 2) CALL safe_copy_string('Y', b%dim_labels(2))
        IF (ndims .GE. 3) CALL safe_copy_string('Z', b%dim_labels(3))
      ENDIF

      IF (PRESENT(dim_units)) THEN
        DO i = 1,ndims
          CALL safe_copy_string(dim_units(i), b%dim_units(i))
        ENDDO
      ELSE
        DO i = 1,ndims
          CALL safe_copy_string('m', b%dim_units(i))
        ENDDO
      ENDIF

      IF (PRESENT(dim_mults)) THEN
        DO i = 1,ndims
          b%dim_mults(i) = REAL(dim_mults(i),r8)
        ENDDO
      ELSE
        DO i = 1,ndims
          b%dim_mults(i) = 1.d0
        ENDDO
      ENDIF

      CALL sdf_write_block_header(h, id, name)
    ELSE
      CALL write_block_header(h)
    ENDIF

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_WRITE(h%filehandle, b%dim_mults, ndims, MPI_REAL8, &
          MPI_STATUS_IGNORE, errcode)

      DO i = 1,ndims
        CALL sdf_safe_write_id(h, b%dim_labels(i))
      ENDDO

      DO i = 1,ndims
        CALL sdf_safe_write_id(h, b%dim_units(i))
      ENDDO

      CALL MPI_FILE_WRITE(h%filehandle, b%geometry, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%extents, 2 * ndims, MPI_REAL8, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%npoints, 1, MPI_INTEGER8, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_location = b%block_start + b%info_length
    b%done_info = .TRUE.

  END SUBROUTINE write_point_mesh_meta_r8



  SUBROUTINE write_point_mesh_meta_r4(h, id, name, dim_labels, dim_units, &
      dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    REAL(r8), DIMENSION(c_maxdims) :: dim_mults8
    TYPE(sdf_block_type), POINTER :: b
    INTEGER :: i

    IF (PRESENT(dim_mults)) THEN
      b => h%current_block
      DO i = 1,b%ndims
        dim_mults8(i) = REAL(dim_mults(i),r8)
      ENDDO
      CALL write_point_mesh_meta_r8(h, id, name, dim_labels, dim_units, &
          dim_mults8)
    ELSE
      CALL write_point_mesh_meta_r8(h, id, name, dim_labels, dim_units)
    ENDIF

  END SUBROUTINE write_point_mesh_meta_r4



  SUBROUTINE write_point_variable_meta_r8(h, id, name, units, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name, units
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: mesh_id
    REAL(r8), INTENT(IN), OPTIONAL :: mult
    INTEGER :: ndims
    TYPE(sdf_block_type), POINTER :: b
    INTEGER :: errcode

    b => h%current_block

    b%blocktype = c_blocktype_point_variable
    ndims = b%ndims
    b%nelements = b%ndims * b%npoints

    ! Metadata is
    ! - mult      REAL(r8)
    ! - units     CHARACTER(id_length)
    ! - meshid    CHARACTER(id_length)
    ! - npoints   INTEGER(i8)

    b%info_length = h%block_header_length + soi8 + sof8 + 2 * c_id_length
    b%data_length = b%nelements * b%type_size

    ! Write header

    IF (PRESENT(id)) THEN
      CALL safe_copy_string(units, b%units)
      CALL safe_copy_string(mesh_id, b%mesh_id)

      IF (PRESENT(mult)) THEN
        b%mult = REAL(mult,r8)
      ELSE
        b%mult = 1.d0
      ENDIF

      CALL sdf_write_block_header(h, id, name)
    ELSE
      CALL write_block_header(h)
    ENDIF

    IF (h%rank .EQ. h%rank_master) THEN
      CALL MPI_FILE_WRITE(h%filehandle, b%mult, 1, MPI_REAL8, &
          MPI_STATUS_IGNORE, errcode)

      CALL sdf_safe_write_id(h, b%units)

      CALL sdf_safe_write_id(h, b%mesh_id)

      CALL MPI_FILE_WRITE(h%filehandle, b%npoints, 1, MPI_INTEGER8, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_location = b%block_start + b%info_length
    b%done_info = .TRUE.

  END SUBROUTINE write_point_variable_meta_r8



  SUBROUTINE write_point_variable_meta_r4(h, id, name, units, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name, units
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: mesh_id
    REAL(r4), INTENT(IN), OPTIONAL :: mult

    IF (PRESENT(mult)) THEN
      CALL write_point_variable_meta_r8(h, id, name, units, mesh_id, &
          REAL(mult,r8))
    ELSE
      CALL write_point_variable_meta_r8(h, id, name, units, mesh_id)
    ENDIF

  END SUBROUTINE write_point_variable_meta_r4



  SUBROUTINE write_srl_pt_var_int_i8_r8(h, id, name, units, array, &
      npoint_global, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: array
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r8), INTENT(IN), OPTIONAL :: mult
    INTEGER(i8) :: i, idx, npoint_max, npoint_rem
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (npoint_global .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = INT(h%soi,r4)
    b%datatype = h%datatype_integer
    b%mpitype = h%mpitype_integer
    b%ndims = 1
    b%npoints = npoint_global

    ! Write header

    CALL write_point_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the real data

    IF (h%rank .EQ. h%rank_master) THEN
      h%current_location = b%data_location
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      ! This is all a bit messy, but it is necessary because MPI_FILE_WRITE
      ! accepts an INTEGER count of elements to write, which may not be
      ! big enough for npoint_global which is an INTEGER*8

      npoint_max = HUGE(npoint_max)
      npoint_rem = MOD(npoint_global, npoint_max)

      idx = 1
      DO i = 1, npoint_global / npoint_max
        CALL MPI_FILE_WRITE(h%filehandle, array(idx), npoint_max, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
        idx = idx + npoint_max
      ENDDO

      CALL MPI_FILE_WRITE(h%filehandle, array(idx), npoint_rem, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_pt_var_int_i8_r8



  SUBROUTINE write_srl_pt_var_int_i4_r8(h, id, name, units, array, &
      npoint_global, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r8), INTENT(IN), OPTIONAL :: mult

    CALL write_srl_pt_var_int_i8_r8(h, id, name, units, array, &
        INT(npoint_global,i8), mesh_id, mult)

  END SUBROUTINE write_srl_pt_var_int_i4_r8



  SUBROUTINE write_srl_pt_var_int_i8_r4(h, id, name, units, array, &
      npoint_global, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: array
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: mult

    CALL write_srl_pt_var_int_i8_r8(h, id, name, units, array, &
        npoint_global, mesh_id, REAL(mult,r8))

  END SUBROUTINE write_srl_pt_var_int_i8_r4



  SUBROUTINE write_srl_pt_var_int_i4_r4(h, id, name, units, array, &
      npoint_global, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: mult

    CALL write_srl_pt_var_int_i8_r8(h, id, name, units, array, &
        INT(npoint_global,i8), mesh_id, REAL(mult,r8))

  END SUBROUTINE write_srl_pt_var_int_i4_r4

END MODULE sdf_output_point_ru
