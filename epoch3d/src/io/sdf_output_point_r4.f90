MODULE sdf_output_point_r4

  USE sdf_output_point_ru

  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: sof = 4
  INTEGER, PARAMETER, PRIVATE :: datatype_real = c_datatype_real4
  INTEGER, PARAMETER, PRIVATE :: mpitype_real = MPI_REAL4

CONTAINS

  SUBROUTINE write_srl_1d_pt_mesh_i8_r4(h, id, name, species_id, x, &
      npoint_global, dim_labels, dim_units, dim_mults)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, species_id
    REAL(r4), DIMENSION(:), INTENT(IN) :: x
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER(i8) :: i, idx
    INTEGER :: npoint_max, npoint_rem, errcode
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    TYPE(sdf_block_type), POINTER :: b

    IF (npoint_global .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = sof
    b%datatype = datatype_real
    b%mpitype = mpitype_real
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

    CALL write_point_mesh_meta_r4(h, id, name, species_id, dim_labels, &
        dim_units, dim_mults)

    ! Write the real data

    IF (h%rank .EQ. h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      ! This is all a bit messy, but it is necessary because MPI_FILE_WRITE
      ! accepts an INTEGER count of elements to write, which may not be
      ! big enough for npoint_global which is an INTEGER*8

      npoint_max = HUGE(npoint_max)
      npoint_rem = INT(MOD(npoint_global, INT(npoint_max,i8)),i4)

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

  END SUBROUTINE write_srl_1d_pt_mesh_i8_r4



  SUBROUTINE write_srl_2d_pt_mesh_i8_r4(h, id, name, species_id, x, y, &
      npoint_global, dim_labels, dim_units, dim_mults)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, species_id
    REAL(r4), DIMENSION(:), INTENT(IN) :: x, y
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER(i8) :: i, idx
    INTEGER :: npoint_max, npoint_rem, errcode
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    TYPE(sdf_block_type), POINTER :: b

    IF (npoint_global .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = sof
    b%datatype = datatype_real
    b%mpitype = mpitype_real
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

    CALL write_point_mesh_meta_r4(h, id, name, species_id, dim_labels, &
        dim_units, dim_mults)

    ! Write the real data

    IF (h%rank .EQ. h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      ! This is all a bit messy, but it is necessary because MPI_FILE_WRITE
      ! accepts an INTEGER count of elements to write, which may not be
      ! big enough for npoint_global which is an INTEGER*8

      npoint_max = HUGE(npoint_max)
      npoint_rem = INT(MOD(npoint_global, INT(npoint_max,i8)),i4)

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

  END SUBROUTINE write_srl_2d_pt_mesh_i8_r4



  SUBROUTINE write_srl_3d_pt_mesh_i8_r4(h, id, name, species_id, x, y, z, &
      npoint_global, dim_labels, dim_units, dim_mults)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, species_id
    REAL(r4), DIMENSION(:), INTENT(IN) :: x, y, z
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER(i8) :: i, idx
    INTEGER :: npoint_max, npoint_rem, errcode
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    TYPE(sdf_block_type), POINTER :: b

    IF (npoint_global .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = sof
    b%datatype = datatype_real
    b%mpitype = mpitype_real
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

    CALL write_point_mesh_meta_r4(h, id, name, species_id, dim_labels, &
        dim_units, dim_mults)

    ! Write the real data

    IF (h%rank .EQ. h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      ! This is all a bit messy, but it is necessary because MPI_FILE_WRITE
      ! accepts an INTEGER count of elements to write, which may not be
      ! big enough for npoint_global which is an INTEGER*8

      npoint_max = HUGE(npoint_max)
      npoint_rem = INT(MOD(npoint_global, INT(npoint_max,i8)),i4)

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

  END SUBROUTINE write_srl_3d_pt_mesh_i8_r4



  SUBROUTINE write_srl_1d_pt_mesh_i4_r4(h, id, name, species_id, x, &
      npoint_global, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, species_id
    REAL(r4), DIMENSION(:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    CALL write_srl_1d_pt_mesh_i8_r4(h, id, name, species_id, x, &
        INT(npoint_global,i8), dim_labels, dim_units, dim_mults)

  END SUBROUTINE write_srl_1d_pt_mesh_i4_r4



  SUBROUTINE write_srl_2d_pt_mesh_i4_r4(h, id, name, species_id, x, y, &
      npoint_global, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, species_id
    REAL(r4), DIMENSION(:), INTENT(IN) :: x, y
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    CALL write_srl_2d_pt_mesh_i8_r4(h, id, name, species_id, x, y, &
        INT(npoint_global,i8), dim_labels, dim_units, dim_mults)

  END SUBROUTINE write_srl_2d_pt_mesh_i4_r4



  SUBROUTINE write_srl_3d_pt_mesh_i4_r4(h, id, name, species_id, x, y, z, &
      npoint_global, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, species_id
    REAL(r4), DIMENSION(:), INTENT(IN) :: x, y, z
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    CALL write_srl_3d_pt_mesh_i8_r4(h, id, name, species_id, x, y, z, &
        INT(npoint_global,i8), dim_labels, dim_units, dim_mults)

  END SUBROUTINE write_srl_3d_pt_mesh_i4_r4



  !----------------------------------------------------------------------------
  ! Code to write a nD point mesh in parallel using an iterator
  !----------------------------------------------------------------------------

  SUBROUTINE write_point_mesh_r4(h, id, name, species_id, npoint_global, &
      ndims, iterator, offset, convert_in, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, species_id
    INTEGER(i8), INTENT(IN) :: npoint_global
    INTEGER(i4), INTENT(IN) :: ndims
    INTEGER(i8), INTENT(IN) :: offset
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, direction, param)
        USE sdf_common
        REAL(r4) :: iterator
        REAL(r4), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN) :: direction
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    CALL write_point_mesh_gen_r4(h, id, name, species_id, npoint_global, &
        ndims, iterator, 0, offset, convert_in, dim_labels, dim_units, &
        dim_mults)

  END SUBROUTINE write_point_mesh_r4



  !----------------------------------------------------------------------------
  ! Code to write a nD point mesh in parallel using an iterator with parameter
  !----------------------------------------------------------------------------

  SUBROUTINE write_point_mesh_gen_r4(h, id, name, species_id, npoint_global, &
      ndims, iterator, param, offset, convert_in, dim_labels, dim_units, &
      dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, species_id
    INTEGER(i8), INTENT(IN) :: npoint_global
    INTEGER(i4), INTENT(IN) :: ndims, param
    INTEGER(i8), INTENT(IN) :: offset
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, direction, param)
        USE sdf_common
        REAL(r4) :: iterator
        REAL(r4), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN) :: direction
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    INTEGER(MPI_OFFSET_KIND) :: file_offset, offset_for_min_max
    INTEGER :: errcode, idim, npoint_this_cycle, nmax
    LOGICAL :: start, convert
    REAL(r8), DIMENSION(c_maxdims) :: gmn, gmx
    REAL(r4), ALLOCATABLE, DIMENSION(:) :: array
    REAL(r4), ALLOCATABLE, DIMENSION(:) :: r4array
    TYPE(sdf_block_type), POINTER :: b
    REAL(r4) :: ret

    IF (npoint_global .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    ENDIF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    ENDIF
    b%geometry = c_geometry_cartesian
    b%ndims = ndims
    b%npoints = npoint_global

    gmn =  HUGE(1.d0)
    gmx = -HUGE(1.d0)

    ! Write header

    CALL write_point_mesh_meta_r4(h, id, name, species_id, dim_labels, &
        dim_units, dim_mults)

    ! Write the real data

    ALLOCATE(array(1:npoint_per_iteration))
    IF (convert) ALLOCATE(r4array(1:npoint_per_iteration))

    DO idim = 1, ndims
      npoint_this_cycle = INT(npoint_per_iteration)
      start = .TRUE.
      file_offset = h%current_location + offset * b%type_size

      DO
        ret = iterator(array, npoint_this_cycle, start, idim, param)
        nmax = npoint_this_cycle
        CALL MPI_ALLREDUCE(npoint_this_cycle, nmax, 1, MPI_INTEGER, &
            MPI_MAX, h%comm, errcode)
        IF (nmax .LE. 0) EXIT

        start = .FALSE.
        gmn(idim) = MIN(gmn(idim),REAL(MINVAL(array(1:npoint_this_cycle)),r8))
        gmx(idim) = MAX(gmx(idim),REAL(MAXVAL(array(1:npoint_this_cycle)),r8))

        CALL MPI_FILE_SET_VIEW(h%filehandle, file_offset, MPI_BYTE, &
            b%mpitype, 'native', MPI_INFO_NULL, errcode)
        IF (convert) THEN
          r4array(1:npoint_this_cycle) = REAL(array(1:npoint_this_cycle),r4)
          CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, npoint_this_cycle, &
              b%mpitype, MPI_STATUS_IGNORE, errcode)
        ELSE
          CALL MPI_FILE_WRITE_ALL(h%filehandle, array, npoint_this_cycle, &
              b%mpitype, MPI_STATUS_IGNORE, errcode)
        ENDIF

        file_offset = file_offset + npoint_this_cycle * b%type_size
      ENDDO

      h%current_location = h%current_location &
          + npoint_global * b%type_size
    ENDDO

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    DEALLOCATE(array)
    IF (convert) DEALLOCATE(r4array)

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

  END SUBROUTINE write_point_mesh_gen_r4



  SUBROUTINE write_srl_pt_var_flt_i8_r4(h, id, name, species_id, units, array, &
      npoint_global, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, species_id, units
    REAL(r4), DIMENSION(:), INTENT(IN) :: array
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN), OPTIONAL :: mult
    INTEGER(i8) :: i, idx
    INTEGER :: npoint_max, npoint_rem, errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (npoint_global .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = sof
    b%datatype = datatype_real
    b%mpitype = mpitype_real
    b%ndims = 1
    b%npoints = npoint_global

    ! Write header

    CALL write_point_variable_meta_r4(h, id, name, species_id, units, &
        mesh_id, mult)

    ! Write the real data

    IF (h%rank .EQ. h%rank_master) THEN
      h%current_location = b%data_location
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      ! This is all a bit messy, but it is necessary because MPI_FILE_WRITE
      ! accepts an INTEGER count of elements to write, which may not be
      ! big enough for npoint_global which is an INTEGER*8

      npoint_max = HUGE(npoint_max)
      npoint_rem = INT(MOD(npoint_global, INT(npoint_max,i8)),i4)

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

  END SUBROUTINE write_srl_pt_var_flt_i8_r4



  SUBROUTINE write_srl_pt_var_flt_i4_r4(h, id, name, species_id, units, array, &
      npoint_global, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, species_id, units
    REAL(r4), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN), OPTIONAL :: mult

    CALL write_srl_pt_var_flt_i8_r4(h, id, name, species_id, units, array, &
        INT(npoint_global,i8), mesh_id, mult)

  END SUBROUTINE write_srl_pt_var_flt_i4_r4



  !----------------------------------------------------------------------------
  ! Code to write a point variable in parallel using an iterator
  !----------------------------------------------------------------------------

  SUBROUTINE write_point_variable_r4(h, id, name, species_id, units, &
      npoint_global, mesh_id, iterator, offset, convert_in, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, species_id, units
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), INTENT(IN) :: offset
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    REAL(r4), INTENT(IN), OPTIONAL :: mult

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, param)
        USE sdf_common
        REAL(r4) :: iterator
        REAL(r4), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    CALL write_point_variable_gen_r4(h, id, name, species_id, units, &
        npoint_global, mesh_id, iterator, 0, offset, convert_in, mult)

  END SUBROUTINE write_point_variable_r4



  !----------------------------------------------------------------------------
  ! Code to write a point variable in parallel using an iterator and parameter
  !----------------------------------------------------------------------------

  SUBROUTINE write_point_variable_gen_r4(h, id, name, species_id, units, &
      npoint_global, mesh_id, iterator, param, offset, convert_in, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, species_id, units
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER, INTENT(IN) :: param
    INTEGER(i8), INTENT(IN) :: offset
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    REAL(r4), INTENT(IN), OPTIONAL :: mult

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, param)
        USE sdf_common
        REAL(r4) :: iterator
        REAL(r4), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    INTEGER(MPI_OFFSET_KIND) :: file_offset
    INTEGER :: errcode, npoint_this_cycle, nmax
    LOGICAL :: start, convert
    REAL(r4), ALLOCATABLE, DIMENSION(:) :: array
    REAL(r4), ALLOCATABLE, DIMENSION(:) :: r4array
    TYPE(sdf_block_type), POINTER :: b
    REAL(r4) :: ret

    IF (npoint_global .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    ENDIF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    ENDIF
    b%blocktype = c_blocktype_point_variable
    b%ndims = 1
    b%npoints = npoint_global

    ! Write header

    CALL write_point_variable_meta_r4(h, id, name, species_id, units, &
        mesh_id, mult)

    ! Write the real data

    ALLOCATE(array(1:npoint_per_iteration))
    IF (convert) ALLOCATE(r4array(1:npoint_per_iteration))

    npoint_this_cycle = INT(npoint_per_iteration)
    start = .TRUE.
    file_offset = h%current_location + offset * b%type_size

    DO
      ret = iterator(array, npoint_this_cycle, start, param)
      nmax = npoint_this_cycle
      CALL MPI_ALLREDUCE(npoint_this_cycle, nmax, 1, MPI_INTEGER, &
          MPI_MAX, h%comm, errcode)
      IF (nmax .LE. 0) EXIT

      IF (start) start = .FALSE.

      CALL MPI_FILE_SET_VIEW(h%filehandle, file_offset, MPI_BYTE, &
          b%mpitype, 'native', MPI_INFO_NULL, errcode)
      IF (convert) THEN
        r4array(1:npoint_this_cycle) = REAL(array(1:npoint_this_cycle),r4)
        CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, npoint_this_cycle, &
            b%mpitype, MPI_STATUS_IGNORE, errcode)
      ELSE
        CALL MPI_FILE_WRITE_ALL(h%filehandle, array, npoint_this_cycle, &
            b%mpitype, MPI_STATUS_IGNORE, errcode)
      ENDIF

      file_offset = file_offset + npoint_this_cycle * b%type_size
    ENDDO

    DEALLOCATE(array)
    IF (convert) DEALLOCATE(r4array)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_point_variable_gen_r4



  ! Calls without the species_id argument for backwards compatibility
  SUBROUTINE write_nospec_srl_1d_pt_mesh_i8_r4(h, id, name, x, &
      npoint_global, dim_labels, dim_units, dim_mults)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    CALL write_srl_1d_pt_mesh_i8_r4(h, id, name, 'nospec', x, &
        npoint_global, dim_labels, dim_units, dim_mults)

  END SUBROUTINE write_nospec_srl_1d_pt_mesh_i8_r4



  SUBROUTINE write_nospec_srl_2d_pt_mesh_i8_r4(h, id, name, x, y, &
      npoint_global, dim_labels, dim_units, dim_mults)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x, y
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    CALL write_srl_2d_pt_mesh_i8_r4(h, id, name, 'nospec', x, y, &
        npoint_global, dim_labels, dim_units, dim_mults)

  END SUBROUTINE write_nospec_srl_2d_pt_mesh_i8_r4



  SUBROUTINE write_nospec_srl_3d_pt_mesh_i8_r4(h, id, name, x, y, z, &
      npoint_global, dim_labels, dim_units, dim_mults)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x, y, z
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    CALL write_srl_3d_pt_mesh_i8_r4(h, id, name, 'nospec', x, y, z, &
        npoint_global, dim_labels, dim_units, dim_mults)

  END SUBROUTINE write_nospec_srl_3d_pt_mesh_i8_r4



  SUBROUTINE write_nospec_srl_1d_pt_mesh_i4_r4(h, id, name, x, &
      npoint_global, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    CALL write_srl_1d_pt_mesh_i4_r4(h, id, name, 'nospec', x, &
        npoint_global, dim_labels, dim_units, dim_mults)

  END SUBROUTINE write_nospec_srl_1d_pt_mesh_i4_r4



  SUBROUTINE write_nospec_srl_2d_pt_mesh_i4_r4(h, id, name, x, y, &
      npoint_global, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x, y
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    CALL write_srl_2d_pt_mesh_i4_r4(h, id, name, 'nospec', x, y, &
        npoint_global, dim_labels, dim_units, dim_mults)

  END SUBROUTINE write_nospec_srl_2d_pt_mesh_i4_r4



  SUBROUTINE write_nospec_srl_3d_pt_mesh_i4_r4(h, id, name, x, y, z, &
      npoint_global, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x, y, z
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    CALL write_srl_3d_pt_mesh_i4_r4(h, id, name, 'nospec', x, y, z, &
        npoint_global, dim_labels, dim_units, dim_mults)

  END SUBROUTINE write_nospec_srl_3d_pt_mesh_i4_r4



  SUBROUTINE write_nospec_point_mesh_r4(h, id, name, npoint_global, &
      ndims, iterator, offset, convert_in, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER(i8), INTENT(IN) :: npoint_global
    INTEGER(i4), INTENT(IN) :: ndims
    INTEGER(i8), INTENT(IN) :: offset
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, direction, param)
        USE sdf_common
        REAL(r4) :: iterator
        REAL(r4), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN) :: direction
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    CALL write_point_mesh_r4(h, id, name, 'nospec', npoint_global, &
        ndims, iterator, offset, convert_in, dim_labels, dim_units, dim_mults)

  END SUBROUTINE write_nospec_point_mesh_r4



  SUBROUTINE write_nospec_point_mesh_gen_r4(h, id, name, npoint_global, &
      ndims, iterator, param, offset, convert_in, dim_labels, dim_units, &
      dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER(i8), INTENT(IN) :: npoint_global
    INTEGER(i4), INTENT(IN) :: ndims, param
    INTEGER(i8), INTENT(IN) :: offset
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, direction, param)
        USE sdf_common
        REAL(r4) :: iterator
        REAL(r4), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN) :: direction
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    CALL write_point_mesh_gen_r4(h, id, name, 'nospec', npoint_global, &
        ndims, iterator, param, offset, convert_in, dim_labels, dim_units, &
        dim_mults)

  END SUBROUTINE write_nospec_point_mesh_gen_r4



  SUBROUTINE write_nospec_srl_pt_var_flt_i8_r4(h, id, name, units, array, &
      npoint_global, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    REAL(r4), DIMENSION(:), INTENT(IN) :: array
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN), OPTIONAL :: mult

    CALL write_srl_pt_var_flt_i8_r4(h, id, name, 'nospec', units, array, &
        npoint_global, mesh_id, mult)

  END SUBROUTINE write_nospec_srl_pt_var_flt_i8_r4



  SUBROUTINE write_nospec_srl_pt_var_flt_i4_r4(h, id, name, units, array, &
      npoint_global, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    REAL(r4), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN), OPTIONAL :: mult

    CALL write_srl_pt_var_flt_i4_r4(h, id, name, 'nospec', units, array, &
        npoint_global, mesh_id, mult)

  END SUBROUTINE write_nospec_srl_pt_var_flt_i4_r4



  SUBROUTINE write_nospec_point_variable_r4(h, id, name, units, &
      npoint_global, mesh_id, iterator, offset, convert_in, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), INTENT(IN) :: offset
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    REAL(r4), INTENT(IN), OPTIONAL :: mult

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, param)
        USE sdf_common
        REAL(r4) :: iterator
        REAL(r4), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    CALL write_point_variable_r4(h, id, name, 'nospec', units, &
        npoint_global, mesh_id, iterator, offset, convert_in, mult)

  END SUBROUTINE write_nospec_point_variable_r4



  SUBROUTINE write_nospec_point_variable_gen_r4(h, id, name, units, &
      npoint_global, mesh_id, iterator, param, offset, convert_in, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER(i8), INTENT(IN) :: npoint_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER, INTENT(IN) :: param
    INTEGER(i8), INTENT(IN) :: offset
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    REAL(r4), INTENT(IN), OPTIONAL :: mult

    INTERFACE
      FUNCTION iterator(array, npoint_it, start, param)
        USE sdf_common
        REAL(r4) :: iterator
        REAL(r4), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER, INTENT(INOUT) :: npoint_it
        LOGICAL, INTENT(IN) :: start
        INTEGER, INTENT(IN), OPTIONAL :: param
      END FUNCTION iterator
    END INTERFACE

    CALL write_point_variable_gen_r4(h, id, name, 'nospec', units, &
        npoint_global, mesh_id, iterator, param, offset, convert_in, mult)

  END SUBROUTINE write_nospec_point_variable_gen_r4

END MODULE sdf_output_point_r4
