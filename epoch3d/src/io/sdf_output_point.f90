MODULE sdf_output_point

  USE sdf_common
  USE sdf_output
  USE mpi

  IMPLICIT NONE

  INTERFACE sdf_write_srl_point_mesh
    MODULE PROCEDURE &
        sdf_write_srl_1d_pt_mesh_array, &
        sdf_write_srl_2d_pt_mesh_array, &
        sdf_write_srl_3d_pt_mesh_array, &
        sdf_write_srl_1d_pt_mesh_array4, &
        sdf_write_srl_2d_pt_mesh_array4, &
        sdf_write_srl_3d_pt_mesh_array4
  END INTERFACE sdf_write_srl_point_mesh

  INTERFACE sdf_write_srl_point_variable
    MODULE PROCEDURE &
        sdf_write_srl_pt_var_int_array, &
        sdf_write_srl_pt_var_int_array4, &
        sdf_write_srl_pt_var_flt_array, &
        sdf_write_srl_pt_var_flt_array4
  END INTERFACE sdf_write_srl_point_variable

CONTAINS

  SUBROUTINE write_point_mesh_meta(h, id, name, dim_labels, dim_units, &
      dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
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

  END SUBROUTINE write_point_mesh_meta



  SUBROUTINE sdf_write_srl_1d_pt_mesh_array(h, id, name, x, &
      npoint_global, dim_labels, dim_units, dim_mults)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(num), DIMENSION(:), INTENT(IN) :: x
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER(i8) :: idx, npoint_max
    INTEGER :: errcode, npoint_rem, i
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    TYPE(sdf_block_type), POINTER :: b

    IF (npoint_global .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%geometry = c_geometry_cartesian
    b%ndims = ndims
    b%npoints = npoint_global

    gmn(1) = REAL(MINVAL(x),r8)
    gmx(1) = REAL(MAXVAL(x),r8)

    b%extents(1:ndims) = gmn(1:ndims)
    b%extents(ndims+1:2*ndims) = gmx(1:ndims)

    CALL MPI_ALLREDUCE(gmn, b%extents, ndims, MPI_REAL8, MPI_MIN, &
        h%comm, errcode)
    CALL MPI_ALLREDUCE(gmx, b%extents(ndims+1), ndims, MPI_REAL8, MPI_MAX, &
        h%comm, errcode)

    ! Write header

    CALL write_point_mesh_meta(h, id, name, dim_labels, dim_units, dim_mults)

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

      ! X-coordinates
      idx = 1
      DO i = 1, npoint_global / npoint_max
        CALL MPI_FILE_WRITE(h%filehandle, x(idx), npoint_max, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
        idx = idx + npoint_max
      ENDDO

      CALL MPI_FILE_WRITE(h%filehandle, x(idx), npoint_rem, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_srl_1d_pt_mesh_array



  SUBROUTINE sdf_write_srl_2d_pt_mesh_array(h, id, name, x, y, &
      npoint_global, dim_labels, dim_units, dim_mults)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(num), DIMENSION(:), INTENT(IN) :: x, y
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER(i8) :: idx, npoint_max, npoint_rem
    INTEGER :: errcode, i
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    TYPE(sdf_block_type), POINTER :: b

    IF (npoint_global .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%geometry = c_geometry_cartesian
    b%ndims = ndims
    b%npoints = npoint_global

    gmn(1) = REAL(MINVAL(x),r8)
    gmx(1) = REAL(MAXVAL(x),r8)
    gmn(2) = REAL(MINVAL(y),r8)
    gmx(2) = REAL(MAXVAL(y),r8)

    b%extents(1:ndims) = gmn(1:ndims)
    b%extents(ndims+1:2*ndims) = gmx(1:ndims)

    CALL MPI_ALLREDUCE(gmn, b%extents, ndims, MPI_REAL8, MPI_MIN, &
        h%comm, errcode)
    CALL MPI_ALLREDUCE(gmx, b%extents(ndims+1), ndims, MPI_REAL8, MPI_MAX, &
        h%comm, errcode)

    ! Write header

    CALL write_point_mesh_meta(h, id, name, dim_labels, dim_units, dim_mults)

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

      ! X-coordinates
      idx = 1
      DO i = 1, npoint_global / npoint_max
        CALL MPI_FILE_WRITE(h%filehandle, x(idx), npoint_max, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
        idx = idx + npoint_max
      ENDDO

      CALL MPI_FILE_WRITE(h%filehandle, x(idx), npoint_rem, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)

      ! Y-coordinates
      idx = 1
      DO i = 1, npoint_global / npoint_max
        CALL MPI_FILE_WRITE(h%filehandle, y(idx), npoint_max, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
        idx = idx + npoint_max
      ENDDO

      CALL MPI_FILE_WRITE(h%filehandle, y(idx), npoint_rem, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_srl_2d_pt_mesh_array



  SUBROUTINE sdf_write_srl_3d_pt_mesh_array(h, id, name, x, y, z, &
      npoint_global, dim_labels, dim_units, dim_mults)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(num), DIMENSION(:), INTENT(IN) :: x, y, z
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER(i8) :: idx, npoint_max, npoint_rem
    INTEGER :: errcode, i
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    TYPE(sdf_block_type), POINTER :: b

    IF (npoint_global .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%geometry = c_geometry_cartesian
    b%ndims = ndims
    b%npoints = npoint_global

    gmn(1) = REAL(MINVAL(x),r8)
    gmx(1) = REAL(MAXVAL(x),r8)
    gmn(2) = REAL(MINVAL(y),r8)
    gmx(2) = REAL(MAXVAL(y),r8)
    gmn(3) = REAL(MINVAL(z),r8)
    gmx(3) = REAL(MAXVAL(z),r8)

    b%extents(1:ndims) = gmn(1:ndims)
    b%extents(ndims+1:2*ndims) = gmx(1:ndims)

    CALL MPI_ALLREDUCE(gmn, b%extents, ndims, MPI_REAL8, MPI_MIN, &
        h%comm, errcode)
    CALL MPI_ALLREDUCE(gmx, b%extents(ndims+1), ndims, MPI_REAL8, MPI_MAX, &
        h%comm, errcode)

    ! Write header

    CALL write_point_mesh_meta(h, id, name, dim_labels, dim_units, dim_mults)

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

      ! X-coordinates
      idx = 1
      DO i = 1, npoint_global / npoint_max
        CALL MPI_FILE_WRITE(h%filehandle, x(idx), npoint_max, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
        idx = idx + npoint_max
      ENDDO

      CALL MPI_FILE_WRITE(h%filehandle, x(idx), npoint_rem, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)

      ! Y-coordinates
      idx = 1
      DO i = 1, npoint_global / npoint_max
        CALL MPI_FILE_WRITE(h%filehandle, y(idx), npoint_max, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
        idx = idx + npoint_max
      ENDDO

      CALL MPI_FILE_WRITE(h%filehandle, y(idx), npoint_rem, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)

      ! Z-coordinates
      idx = 1
      DO i = 1, npoint_global / npoint_max
        CALL MPI_FILE_WRITE(h%filehandle, z(idx), npoint_max, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
        idx = idx + npoint_max
      ENDDO

      CALL MPI_FILE_WRITE(h%filehandle, z(idx), npoint_rem, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_srl_3d_pt_mesh_array



  SUBROUTINE sdf_write_srl_1d_pt_mesh_array4(h, id, name, x, &
      npoint_global, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(num), DIMENSION(:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    CALL sdf_write_srl_1d_pt_mesh_array(h, id, name, x, &
        INT(npoint_global,i8), dim_labels, dim_units, dim_mults)

  END SUBROUTINE sdf_write_srl_1d_pt_mesh_array4



  SUBROUTINE sdf_write_srl_2d_pt_mesh_array4(h, id, name, x, y, &
      npoint_global, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(num), DIMENSION(:), INTENT(IN) :: x, y
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    CALL sdf_write_srl_2d_pt_mesh_array(h, id, name, x, y, &
        INT(npoint_global,i8), dim_labels, dim_units, dim_mults)

  END SUBROUTINE sdf_write_srl_2d_pt_mesh_array4



  SUBROUTINE sdf_write_srl_3d_pt_mesh_array4(h, id, name, x, y, z, &
      npoint_global, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(num), DIMENSION(:), INTENT(IN) :: x, y, z
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    CALL sdf_write_srl_3d_pt_mesh_array(h, id, name, x, y, z, &
        INT(npoint_global,i8), dim_labels, dim_units, dim_mults)

  END SUBROUTINE sdf_write_srl_3d_pt_mesh_array4



  !----------------------------------------------------------------------------
  ! Code to write a nD point mesh in parallel using an iterator
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_point_mesh(h, id, name, npoint_global, ndims, &
      iterator, offset, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER(i8), INTENT(IN) :: npoint_global
    INTEGER, INTENT(IN) :: ndims
    INTEGER(i8), INTENT(IN) :: offset
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(num), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, direction)
        USE sdf_common
        REAL(num) :: iterator
        REAL(num), DIMENSION(:), INTENT(OUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN) :: direction
      END FUNCTION iterator
    END INTERFACE

    INTEGER(MPI_OFFSET_KIND) :: offset_for_min_max
    INTEGER(i8) :: file_offset, nmax
    INTEGER :: errcode, idim, npoint_this_cycle
    LOGICAL :: start
    REAL(r8), DIMENSION(c_maxdims) :: gmn, gmx
    REAL(num), ALLOCATABLE, DIMENSION(:) :: array
    TYPE(sdf_block_type), POINTER :: b
    REAL(num) :: ret

    IF (npoint_global .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%geometry = c_geometry_cartesian
    b%ndims = ndims
    b%npoints = npoint_global

    gmn = 0.d0
    gmx = 0.d0

    ! Write header

    CALL write_point_mesh_meta(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the real data

    ALLOCATE(array(1:npoint_per_iteration))

    DO idim = 1, ndims
      npoint_this_cycle = npoint_per_iteration
      start = .TRUE.
      file_offset = h%current_location + offset * b%type_size

      DO
        ret = iterator(array, npoint_this_cycle, start, idim)
        nmax = npoint_this_cycle
        CALL MPI_ALLREDUCE(npoint_this_cycle, nmax, 1, h%mpitype_integer, &
            MPI_MAX, h%comm, errcode)
        IF (nmax .LE. 0) EXIT

        IF (start) THEN
          gmn(idim) = REAL(MINVAL(array(1:npoint_this_cycle)),r8)
          gmx(idim) = REAL(MAXVAL(array(1:npoint_this_cycle)),r8)
          start = .FALSE.
        ELSE
          gmn(idim) = REAL(MIN(gmn(idim),MINVAL(array(1:npoint_this_cycle))),r8)
          gmx(idim) = REAL(MAX(gmx(idim),MAXVAL(array(1:npoint_this_cycle))),r8)
        ENDIF

        CALL MPI_FILE_SET_VIEW(h%filehandle, file_offset, b%mpitype, &
            b%mpitype, "native", MPI_INFO_NULL, errcode)
        CALL MPI_FILE_WRITE_ALL(h%filehandle, array, npoint_this_cycle, &
            b%mpitype, MPI_STATUS_IGNORE, errcode)

        file_offset = file_offset + npoint_this_cycle * b%type_size
      ENDDO

      h%current_location = h%current_location &
          + npoint_global * b%type_size
    ENDDO

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    DEALLOCATE(array)

    ! Write updated values for the mesh extents into the metadata section

    b%extents(1:ndims) = gmn(1:ndims)
    b%extents(ndims+1:2*ndims) = gmx(1:ndims)

    CALL MPI_ALLREDUCE(gmn, b%extents, ndims, MPI_REAL8, MPI_MIN, &
        h%comm, errcode)
    CALL MPI_ALLREDUCE(gmx, b%extents(ndims+1), ndims, MPI_REAL8, MPI_MAX, &
        h%comm, errcode)

    IF (h%rank .EQ. h%rank_master) THEN
      offset_for_min_max = b%block_start + h%block_header_length &
          + (sof8 + 2 * c_id_length) * ndims + sof4

      CALL MPI_FILE_SEEK(h%filehandle, offset_for_min_max, MPI_SEEK_SET, &
          errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%extents, 2 * ndims, MPI_REAL8, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_point_mesh



  SUBROUTINE write_point_variable_meta(h, id, name, units, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name, units
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: mesh_id
    REAL(num), INTENT(IN), OPTIONAL :: mult
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

  END SUBROUTINE write_point_variable_meta



  SUBROUTINE sdf_write_srl_pt_var_flt_array(h, id, name, units, array, &
      npoint_global, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), INTENT(IN), OPTIONAL :: mult
    INTEGER(i8) :: idx, npoint_max, npoint_rem
    INTEGER :: errcode, i
    TYPE(sdf_block_type), POINTER :: b

    IF (npoint_global .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%ndims = 1
    b%npoints = npoint_global

    ! Write header

    CALL write_point_variable_meta(h, id, name, units, mesh_id, mult)

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

  END SUBROUTINE sdf_write_srl_pt_var_flt_array



  SUBROUTINE sdf_write_srl_pt_var_flt_array4(h, id, name, units, array, &
      npoint_global, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), INTENT(IN), OPTIONAL :: mult

    CALL sdf_write_srl_pt_var_flt_array(h, id, name, units, array, &
        INT(npoint_global,i8), mesh_id, mult)

  END SUBROUTINE sdf_write_srl_pt_var_flt_array4



  SUBROUTINE sdf_write_srl_pt_var_int_array(h, id, name, units, array, &
      npoint_global, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: array
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), INTENT(IN), OPTIONAL :: mult
    INTEGER(i8) :: idx, npoint_max, npoint_rem
    INTEGER :: errcode, i
    TYPE(sdf_block_type), POINTER :: b

    IF (npoint_global .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%soi
    b%datatype = h%datatype_integer
    b%mpitype = h%mpitype_integer
    b%ndims = 1
    b%npoints = npoint_global

    ! Write header

    CALL write_point_variable_meta(h, id, name, units, mesh_id, mult)

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

  END SUBROUTINE sdf_write_srl_pt_var_int_array



  SUBROUTINE sdf_write_srl_pt_var_int_array4(h, id, name, units, array, &
      npoint_global, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(num), INTENT(IN), OPTIONAL :: mult

    CALL sdf_write_srl_pt_var_int_array(h, id, name, units, array, &
        INT(npoint_global,i8), mesh_id, mult)

  END SUBROUTINE sdf_write_srl_pt_var_int_array4



  !----------------------------------------------------------------------------
  ! Code to write a point variable in parallel using an iterator
  !----------------------------------------------------------------------------

  SUBROUTINE sdf_write_point_variable(h, id, name, units, npoint_global, &
      mesh_id, iterator, offset, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), INTENT(IN) :: offset
    REAL(num), INTENT(IN), OPTIONAL :: mult

    INTERFACE
      FUNCTION iterator(array, npoint_it, start)
        USE sdf_common
        REAL(num) :: iterator
        REAL(num), DIMENSION(:), INTENT(OUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
      END FUNCTION iterator
    END INTERFACE

    INTEGER(i8) :: file_offset, nmax
    INTEGER :: errcode, npoint_this_cycle
    LOGICAL :: start
    REAL(num), ALLOCATABLE, DIMENSION(:) :: array
    TYPE(sdf_block_type), POINTER :: b
    REAL(num) :: ret

    IF (npoint_global .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = h%sof
    b%datatype = h%datatype_real
    b%mpitype = h%mpitype_real
    b%blocktype = c_blocktype_point_variable
    b%ndims = 1
    b%npoints = npoint_global

    ! Write header

    CALL write_point_variable_meta(h, id, name, units, mesh_id, mult)

    ! Write the real data

    ALLOCATE(array(1:npoint_per_iteration))

    npoint_this_cycle = npoint_per_iteration
    start = .TRUE.
    file_offset = h%current_location + offset * b%type_size

    DO
      ret = iterator(array, npoint_this_cycle, start)
      nmax = npoint_this_cycle
      CALL MPI_ALLREDUCE(npoint_this_cycle, nmax, 1, h%mpitype_integer, &
          MPI_MAX, h%comm, errcode)
      IF (nmax .LE. 0) EXIT

      IF (start) start = .FALSE.

      CALL MPI_FILE_SET_VIEW(h%filehandle, file_offset, b%mpitype, &
          b%mpitype, "native", MPI_INFO_NULL, errcode)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, array, npoint_this_cycle, &
          b%mpitype, MPI_STATUS_IGNORE, errcode)

      file_offset = file_offset + npoint_this_cycle * b%type_size
    ENDDO

    DEALLOCATE(array)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, "native", &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE sdf_write_point_variable

END MODULE sdf_output_point
