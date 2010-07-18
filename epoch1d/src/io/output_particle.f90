MODULE output_particle

  USE iocommon
  USE output

  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! Code to write a nD particle grid in parallel
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_nd_particle_grid_all(name, class, particles, &
      npart_global, particle_coord_type, particle_type)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), DIMENSION(:,:), INTENT(IN) :: particles
    INTEGER(8), INTENT(IN) :: npart_global
    INTEGER(4), INTENT(IN) :: particle_coord_type
    INTEGER, INTENT(IN) :: particle_type
    INTEGER(8) :: npart_local, block_length, md_length
    INTEGER(4) :: ndims, idim
    INTEGER :: sizes(2)
    REAL(num) :: mn, mx

    IF (npart_global .LE. 0) RETURN

    sizes = SHAPE(particles)
    npart_local = sizes(2)
    ndims = INT(sizes(1),4)

    ! Metadata is
    ! - meshtype (INTEGER(4)) All mesh blocks contain this
    ! - nd     INTEGER(4)
    ! - sof    INTEGER(4)
    ! Specific to particle mesh
    ! - ct     INTEGER(4)
    ! - npart  INTEGER(8)
    ! - d1min  REAL(num)
    ! - d1max  REAL(num)
    ! - d2min  REAL(num)
    ! - d2max  REAL(num)
    ! .
    ! .
    ! .

    md_length = meshtype_header_offset + soi + soi8 + 2 * ndims * sof
    block_length = md_length + npart_global * ndims * sof

    ! Write header
    CALL cfd_write_block_header(name, class, c_type_mesh, &
        block_length, md_length, default_rank)
    CALL cfd_write_meshtype_header(c_mesh_particle, ndims, &
        sof, default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER4, MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, particle_coord_type, 1, &
          MPI_INTEGER4, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER8, MPI_INTEGER8, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, npart_global, 1, &
          MPI_INTEGER8, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + soi8

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    DO idim = 1, ndims
      CALL MPI_ALLREDUCE(MINVAL(particles(:,idim)), mn, 1, mpireal, MPI_MIN, &
          cfd_comm, cfd_errcode)
      CALL MPI_ALLREDUCE(MAXVAL(particles(:,idim)), mx, 1, mpireal, MPI_MAX, &
          cfd_comm, cfd_errcode)

      IF (cfd_rank .EQ. default_rank) THEN
        CALL MPI_FILE_WRITE(cfd_filehandle, mn, 1, mpireal, &
            cfd_status, cfd_errcode)
        CALL MPI_FILE_WRITE(cfd_filehandle, mx, 1, mpireal, &
            cfd_status, cfd_errcode)
      ENDIF

      current_displacement = current_displacement + 2 * sof
    ENDDO

    ! Write the real data

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        particle_type, "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_WRITE_ALL(cfd_filehandle, particles, npart_local * ndims, &
        mpireal, cfd_status, cfd_errcode)

    current_displacement = current_displacement + npart_global * ndims * sof

  END SUBROUTINE cfd_write_nd_particle_grid_all



  !----------------------------------------------------------------------------
  ! Code to write a nD particle grid in parallel using an iterator
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_nd_particle_grid_with_iterator_all(name, class, &
      iterator, ndims, npart_local, npart_global, npart_per_iteration, &
      particle_coord_type, lengths, offsets)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    INTEGER(4), INTENT(IN) :: ndims
    INTEGER(8), INTENT(IN) :: npart_local
    INTEGER(8), INTENT(IN) :: npart_global
    INTEGER(8), INTENT(IN) :: npart_per_iteration
    INTEGER(4), INTENT(IN) :: particle_coord_type
    INTEGER(8), DIMENSION(:), INTENT(IN) :: lengths
    INTEGER(8), DIMENSION(:), INTENT(IN) :: offsets

    INTERFACE
      SUBROUTINE iterator(data, npart_it, direction, start)
        USE shared_data
        REAL(num), DIMENSION(:), INTENT(INOUT) :: data
        INTEGER(8), INTENT(INOUT) :: npart_it
        INTEGER, INTENT(IN) :: direction
        LOGICAL, INTENT(IN) :: start
      END SUBROUTINE iterator
    END INTERFACE

    INTEGER(MPI_OFFSET_KIND) :: offset_for_min_max
    INTEGER(8) :: block_length, md_length, npart_this_cycle
    INTEGER(8) :: file_offset, nmax, nsec_left, nwrite_left, nelements, off
    INTEGER :: idim, osec, sec
    LOGICAL :: start
    REAL(num) :: mn, mx
    REAL(num), ALLOCATABLE, DIMENSION(:) :: gmn, gmx
    REAL(num), ALLOCATABLE, DIMENSION(:) :: data

    IF (npart_global .LE. 0) RETURN

    ! Metadata is
    ! - meshtype (INTEGER(4)) All mesh blocks contain this
    ! - nd     INTEGER(4)
    ! - sof    INTEGER(4)
    ! Specific to particle mesh
    ! - ct     INTEGER(4)
    ! - npart  INTEGER(8)
    ! - d1min  REAL(num)
    ! - d1max  REAL(num)
    ! - d2min  REAL(num)
    ! - d2max  REAL(num)
    ! .
    ! .
    ! .

    md_length = meshtype_header_offset + soi + soi8 + 2 * ndims * sof
    block_length = md_length + npart_global * ndims * sof

    ! Write header
    CALL cfd_write_block_header(name, class, c_type_mesh, &
        block_length, md_length, default_rank)
    CALL cfd_write_meshtype_header(c_mesh_particle, ndims, &
        sof, default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER4, MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, particle_coord_type, 1, &
          MPI_INTEGER4, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER8, MPI_INTEGER8, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, npart_global, 1, &
          MPI_INTEGER8, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + soi8

    ! This is to skip past the location for the min/max values (Just write
    ! zeros). They will be filled in later

    ALLOCATE(gmn(ndims), gmx(ndims))
    gmn = 0.0_num
    gmx = 0.0_num

    offset_for_min_max = current_displacement

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, gmn, ndims, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, gmx, ndims, mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 2 * ndims * sof

    ! Write the real data

    ALLOCATE(data(1:npart_per_iteration))

    DO idim = 1, ndims
      npart_this_cycle = npart_per_iteration
      start = .TRUE.
      sec = 1
      osec = 1
      nsec_left = lengths(sec)
      file_offset = current_displacement + offsets(sec) * num

      DO
        CALL iterator(data, npart_this_cycle, idim, start)
        CALL MPI_ALLREDUCE(npart_this_cycle, nmax, 1, MPI_INTEGER8, MPI_MAX, &
            cfd_comm, cfd_errcode)
        IF (nmax .LE. 0) EXIT

        off = 1
        nwrite_left = npart_this_cycle

        IF (start) THEN
          gmn(idim) = MINVAL(data(1:npart_this_cycle))
          gmx(idim) = MAXVAL(data(1:npart_this_cycle))
          start = .FALSE.
        ELSE
          gmn(idim) = MIN(gmn(idim), MINVAL(data(1:npart_this_cycle)))
          gmx(idim) = MAX(gmx(idim), MAXVAL(data(1:npart_this_cycle)))
        ENDIF

        DO
          IF (nwrite_left .LE. nsec_left) THEN
            nelements = nwrite_left
            nsec_left = nsec_left - nelements
          ELSE
            nelements = nsec_left
            sec = sec + 1
            nsec_left = lengths(sec)
          ENDIF

          CALL MPI_FILE_SET_VIEW(cfd_filehandle, file_offset, mpireal, &
              mpireal, "native", MPI_INFO_NULL, cfd_errcode)
          CALL MPI_FILE_WRITE_ALL(cfd_filehandle, data(off), nelements, &
              mpireal, cfd_status, cfd_errcode)

          nwrite_left = nwrite_left - nelements
          CALL MPI_ALLREDUCE(nwrite_left, nmax, 1, MPI_INTEGER8, MPI_MAX, &
              cfd_comm, cfd_errcode)
          IF (nmax .LE. 0) EXIT

          IF (sec .NE. osec) THEN
            file_offset = current_displacement + offsets(sec) * num
          ELSE
            file_offset = file_offset + nelements * num
          ENDIF
          off = off + nelements
          osec = sec
        ENDDO
      ENDDO

      current_displacement = current_displacement + npart_global * sof
    ENDDO

    DEALLOCATE(data)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, offset_for_min_max, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    DO idim = 1, ndims
      CALL MPI_ALLREDUCE(gmn(idim), mn, 1, mpireal, MPI_MIN, &
          cfd_comm, cfd_errcode)
      CALL MPI_ALLREDUCE(gmx(idim), mx, 1, mpireal, MPI_MAX, &
          cfd_comm, cfd_errcode)

      IF (cfd_rank .EQ. default_rank) THEN
        CALL MPI_FILE_WRITE(cfd_filehandle, mn, 1, mpireal, &
            cfd_status, cfd_errcode)
        CALL MPI_FILE_WRITE(cfd_filehandle, mx, 1, mpireal, &
            cfd_status, cfd_errcode)
      ENDIF
    ENDDO

    DEALLOCATE(gmn, gmx)

  END SUBROUTINE cfd_write_nd_particle_grid_with_iterator_all



  !----------------------------------------------------------------------------
  ! Code to write a nD particle variable in parallel
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_nd_particle_variable_all(name, class, particles, &
      npart_global, mesh_name, mesh_class, particle_type)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), DIMENSION(:), INTENT(IN) :: particles
    INTEGER(8), INTENT(IN) :: npart_global
    CHARACTER(LEN=*), INTENT(IN) :: mesh_name, mesh_class
    INTEGER, INTENT(IN) :: particle_type
    INTEGER(8) :: npart_local, block_length, md_length
    REAL(num) :: mn, mx

    IF (npart_global .LE. 0) RETURN

    npart_local = SIZE(particles)

    ! Metadata is
    ! - meshtype (INTEGER(4)) All mesh blocks contain this
    ! - nd     INTEGER(4)
    ! - sof    INTEGER(4)
    ! Specific to particle variable
    ! - npart  INTEGER(8)
    ! - vmin   REAL(num)
    ! - vmax   REAL(num)
    ! - mesh   CHARACTER
    ! - mclass CHARACTER

    md_length = meshtype_header_offset + soi8 + 2 * sof + 2 * max_string_len
    block_length = md_length + npart_global * sof

    ! Write header
    CALL cfd_write_block_header(name, class, c_type_mesh_variable, &
        block_length, md_length, default_rank)
    CALL cfd_write_meshtype_header(c_var_particle, c_dimension_irrelevant, &
        sof, default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER8, MPI_INTEGER8, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, npart_global, 1, &
          MPI_INTEGER8, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + soi8

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_ALLREDUCE(MINVAL(particles), mn, 1, mpireal, MPI_MIN, &
        cfd_comm, cfd_errcode)
    CALL MPI_ALLREDUCE(MAXVAL(particles), mx, 1, mpireal, MPI_MAX, &
        cfd_comm, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, mn, 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx, 1, mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 2 * sof

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF

    current_displacement = current_displacement + 2 * max_string_len

    ! Write the real data

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        particle_type, "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_WRITE_ALL(cfd_filehandle, particles, npart_local, &
        mpireal, cfd_status, cfd_errcode)

    current_displacement = current_displacement + npart_global * sof

  END SUBROUTINE cfd_write_nd_particle_variable_all



  !----------------------------------------------------------------------------
  ! Code to write a nD particle variable in parallel using an iterator
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_nd_particle_variable_with_iterator_all(name, class, &
      iterator, npart_global, npart_per_iteration, mesh_name, mesh_class, &
      lengths, offsets)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    INTEGER(8), INTENT(IN) :: npart_global
    INTEGER(8), INTENT(IN) :: npart_per_iteration
    CHARACTER(LEN=*), INTENT(IN) :: mesh_name, mesh_class
    INTEGER(8), DIMENSION(:), INTENT(IN) :: lengths
    INTEGER(8), DIMENSION(:), INTENT(IN) :: offsets

    INTERFACE
      SUBROUTINE iterator(data, npart_it, start)
        USE shared_data
        REAL(num), DIMENSION(:), INTENT(INOUT) :: data
        INTEGER(8), INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
      END SUBROUTINE iterator
    END INTERFACE

    INTEGER(MPI_OFFSET_KIND) :: offset_for_min_max
    INTEGER(8) :: block_length, md_length, npart_this_cycle
    INTEGER(8) :: file_offset, nmax, nsec_left, nwrite_left, nelements, off
    INTEGER :: sec, osec
    LOGICAL :: start
    REAL(num) :: mn, mx, gmn, gmx
    REAL(num), ALLOCATABLE, DIMENSION(:) :: data

    IF (npart_global .LE. 0) RETURN

    ! Metadata is
    ! - meshtype (INTEGER(4)) All mesh blocks contain this
    ! - nd     INTEGER(4)
    ! - sof    INTEGER(4)
    ! Specific to particle variable
    ! - npart  INTEGER(8)
    ! - vmin   REAL(num)
    ! - vmax   REAL(num)
    ! - mesh   CHARACTER
    ! - mclass CHARACTER

    md_length = meshtype_header_offset + soi8 + 2 * sof + 2 * max_string_len
    block_length = md_length + npart_global * sof

    ! Write header
    CALL cfd_write_block_header(name, class, c_type_mesh_variable, &
        block_length, md_length, default_rank)
    CALL cfd_write_meshtype_header(c_var_particle, c_dimension_irrelevant, &
        sof, default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER8, MPI_INTEGER8, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, npart_global, 1, &
          MPI_INTEGER8, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + soi8

    ! This is to skip past the location for the min/max values (Just write
    ! zeros). They will be filled in later

    offset_for_min_max = current_displacement

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, 0.0_num, 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, 0.0_num, 1, mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 2 * sof

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF

    current_displacement = current_displacement + 2 * max_string_len

    ! Write the real data

    ALLOCATE(data(1:npart_per_iteration))

    npart_this_cycle = npart_per_iteration
    start = .TRUE.
    sec = 1
    osec = 1
    nsec_left = lengths(sec)
    file_offset = current_displacement + offsets(sec) * num

    DO
      CALL iterator(data, npart_this_cycle, start)
      CALL MPI_ALLREDUCE(npart_this_cycle, nmax, 1, MPI_INTEGER8, MPI_MAX, &
          cfd_comm, cfd_errcode)
      IF (nmax .LE. 0) EXIT

      off = 1
      nwrite_left = npart_this_cycle

      IF (start) THEN
        gmn = MINVAL(data(1:npart_this_cycle))
        gmx = MAXVAL(data(1:npart_this_cycle))
        start = .FALSE.
      ELSE
        gmn = MIN(gmn, MINVAL(data(1:npart_this_cycle)))
        gmx = MAX(gmx, MAXVAL(data(1:npart_this_cycle)))
      ENDIF

      DO
        IF (nwrite_left .LE. nsec_left) THEN
          nelements = nwrite_left
          nsec_left = nsec_left - nelements
        ELSE
          nelements = nsec_left
          sec = sec + 1
          nsec_left = lengths(sec)
        ENDIF

        CALL MPI_FILE_SET_VIEW(cfd_filehandle, file_offset, mpireal, &
            mpireal, "native", MPI_INFO_NULL, cfd_errcode)
        CALL MPI_FILE_WRITE_ALL(cfd_filehandle, data(off), nelements, &
            mpireal, cfd_status, cfd_errcode)

        nwrite_left = nwrite_left - nelements
        CALL MPI_ALLREDUCE(nwrite_left, nmax, 1, MPI_INTEGER8, MPI_MAX, &
            cfd_comm, cfd_errcode)
        IF (nmax .LE. 0) EXIT

        IF (sec .NE. osec) THEN
          file_offset = current_displacement + offsets(sec) * num
        ELSE
          file_offset = file_offset + nelements * num
        ENDIF
        off = off + nelements
        osec = sec
      ENDDO
    ENDDO

    current_displacement = current_displacement + npart_global * sof

    DEALLOCATE(data)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, offset_for_min_max, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_ALLREDUCE(gmn, mn, 1, mpireal, MPI_MIN, cfd_comm, cfd_errcode)
    CALL MPI_ALLREDUCE(gmx, mx, 1, mpireal, MPI_MAX, cfd_comm, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, mn, 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx, 1, mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

  END SUBROUTINE cfd_write_nd_particle_variable_with_iterator_all

END MODULE output_particle
