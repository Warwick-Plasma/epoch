MODULE cfd_output_particle

  USE cfd_common
  USE cfd_output
  USE mpi

  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! Code to write a nD particle grid in parallel using an iterator
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_nd_particle_grid_with_iterator_all(h, name, class, &
      iterator, ndims, npart_local, npart_global, npart_per_iteration, &
      particle_coord_type, lengths, offsets)

    TYPE(cfd_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    INTEGER(4), INTENT(IN) :: ndims
    INTEGER(8), INTENT(IN) :: npart_local
    INTEGER(8), INTENT(IN) :: npart_global
    INTEGER(8), INTENT(IN) :: npart_per_iteration
    INTEGER(4), INTENT(IN) :: particle_coord_type
    INTEGER(8), DIMENSION(:), INTENT(IN) :: lengths
    INTEGER(8), DIMENSION(:), INTENT(IN) :: offsets

    INTERFACE
      SUBROUTINE iterator(array, npart_it, direction, start)
        USE cfd_common
        REAL(num), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER(8), INTENT(INOUT) :: npart_it
        INTEGER, INTENT(IN) :: direction
        LOGICAL, INTENT(IN) :: start
      END SUBROUTINE iterator
    END INTERFACE

    INTEGER(MPI_OFFSET_KIND) :: offset_for_min_max
    INTEGER(8) :: block_length, md_length, npart_this_cycle
    INTEGER(8) :: file_offset, nmax, nsec_left, nwrite_left, nelements, off
    INTEGER :: idim, sec, osec, errcode
    LOGICAL :: start
    REAL(num) :: mn, mx
    REAL(num), ALLOCATABLE, DIMENSION(:) :: gmn, gmx
    REAL(num), ALLOCATABLE, DIMENSION(:) :: array

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

    md_length = c_meshtype_header_offset + soi + soi8 + 2 * ndims * sof
    block_length = md_length + npart_global * ndims * sof

    ! Write header
    CALL cfd_write_block_header(h, name, class, c_type_mesh, &
        block_length, md_length, h%default_rank)
    CALL cfd_write_meshtype_header(h, c_mesh_particle, ndims, &
        sof, h%default_rank)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_INTEGER4, MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. h%default_rank) THEN
      CALL MPI_FILE_WRITE(h%filehandle, particle_coord_type, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_displacement = h%current_displacement + soi

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_INTEGER8, MPI_INTEGER8, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. h%default_rank) THEN
      CALL MPI_FILE_WRITE(h%filehandle, npart_global, 1, &
          MPI_INTEGER8, MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_displacement = h%current_displacement + soi8

    ! This is to skip past the location for the min/max values (Just write
    ! zeros). They will be filled in later

    ALLOCATE(gmn(ndims), gmx(ndims))
    gmn = 0.0_num
    gmx = 0.0_num

    offset_for_min_max = h%current_displacement

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        h%mpireal, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. h%default_rank) THEN
      CALL MPI_FILE_WRITE(h%filehandle, gmn, ndims, h%mpireal, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, gmx, ndims, h%mpireal, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_displacement = h%current_displacement + 2 * ndims * sof

    ! Write the real data

    ALLOCATE(array(1:npart_per_iteration))

    DO idim = 1, ndims
      npart_this_cycle = npart_per_iteration
      start = .TRUE.
      sec = 1
      osec = 1
      nsec_left = lengths(sec)
      file_offset = h%current_displacement + offsets(sec) * num

      DO
        CALL iterator(array, npart_this_cycle, idim, start)
        CALL MPI_ALLREDUCE(npart_this_cycle, nmax, 1, MPI_INTEGER8, MPI_MAX, &
            h%comm, errcode)
        IF (nmax .LE. 0) EXIT

        off = 1
        nwrite_left = npart_this_cycle

        IF (start) THEN
          gmn(idim) = MINVAL(array(1:npart_this_cycle))
          gmx(idim) = MAXVAL(array(1:npart_this_cycle))
          start = .FALSE.
        ELSE
          gmn(idim) = MIN(gmn(idim), MINVAL(array(1:npart_this_cycle)))
          gmx(idim) = MAX(gmx(idim), MAXVAL(array(1:npart_this_cycle)))
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

          CALL MPI_FILE_SET_VIEW(h%filehandle, file_offset, h%mpireal, &
              h%mpireal, "native", MPI_INFO_NULL, errcode)
          CALL MPI_FILE_WRITE_ALL(h%filehandle, array(off), nelements, &
              h%mpireal, MPI_STATUS_IGNORE, errcode)

          nwrite_left = nwrite_left - nelements
          CALL MPI_ALLREDUCE(nwrite_left, nmax, 1, MPI_INTEGER8, MPI_MAX, &
              h%comm, errcode)

          IF (sec .NE. osec) THEN
            file_offset = h%current_displacement + offsets(sec) * num
          ELSE
            file_offset = file_offset + nelements * num
          ENDIF

          IF (nmax .LE. 0) EXIT

          off = off + nelements
          osec = sec
        ENDDO
      ENDDO

      h%current_displacement = h%current_displacement + npart_global * sof
    ENDDO

    DEALLOCATE(array)

    CALL MPI_FILE_SET_VIEW(h%filehandle, offset_for_min_max, h%mpireal, &
        h%mpireal, "native", MPI_INFO_NULL, errcode)

    DO idim = 1, ndims
      CALL MPI_ALLREDUCE(gmn(idim), mn, 1, h%mpireal, MPI_MIN, &
          h%comm, errcode)
      CALL MPI_ALLREDUCE(gmx(idim), mx, 1, h%mpireal, MPI_MAX, &
          h%comm, errcode)

      IF (h%rank .EQ. h%default_rank) THEN
        CALL MPI_FILE_WRITE(h%filehandle, mn, 1, h%mpireal, &
            MPI_STATUS_IGNORE, errcode)
        CALL MPI_FILE_WRITE(h%filehandle, mx, 1, h%mpireal, &
            MPI_STATUS_IGNORE, errcode)
      ENDIF
    ENDDO

    DEALLOCATE(gmn, gmx)

  END SUBROUTINE cfd_write_nd_particle_grid_with_iterator_all



  !----------------------------------------------------------------------------
  ! Code to write a nD particle variable in parallel using an iterator
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_nd_particle_variable_with_iterator_all(h, name, class, &
      iterator, npart_global, npart_per_iteration, mesh_name, mesh_class, &
      lengths, offsets)

    TYPE(cfd_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    INTEGER(8), INTENT(IN) :: npart_global
    INTEGER(8), INTENT(IN) :: npart_per_iteration
    CHARACTER(LEN=*), INTENT(IN) :: mesh_name, mesh_class
    INTEGER(8), DIMENSION(:), INTENT(IN) :: lengths
    INTEGER(8), DIMENSION(:), INTENT(IN) :: offsets

    INTERFACE
      SUBROUTINE iterator(array, npart_it, start)
        USE cfd_common
        REAL(num), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER(8), INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
      END SUBROUTINE iterator
    END INTERFACE

    INTEGER(MPI_OFFSET_KIND) :: offset_for_min_max
    INTEGER(8) :: block_length, md_length, npart_this_cycle
    INTEGER(8) :: file_offset, nmax, nsec_left, nwrite_left, nelements, off
    INTEGER :: sec, osec, errcode
    LOGICAL :: start
    REAL(num) :: mn, mx, gmn, gmx
    REAL(num), ALLOCATABLE, DIMENSION(:) :: array

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

    md_length = c_meshtype_header_offset + soi8 + 2 * sof + 2 * h%max_string_len
    block_length = md_length + npart_global * sof

    ! Write header
    CALL cfd_write_block_header(h, name, class, c_type_mesh_variable, &
        block_length, md_length, h%default_rank)
    CALL cfd_write_meshtype_header(h, c_var_particle, c_dimension_irrelevant, &
        sof, h%default_rank)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_INTEGER8, MPI_INTEGER8, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. h%default_rank) THEN
      CALL MPI_FILE_WRITE(h%filehandle, npart_global, 1, &
          MPI_INTEGER8, MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_displacement = h%current_displacement + soi8

    ! This is to skip past the location for the min/max values (Just write
    ! zeros). They will be filled in later

    offset_for_min_max = h%current_displacement

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        h%mpireal, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. h%default_rank) THEN
      CALL MPI_FILE_WRITE(h%filehandle, 0.0_num, 1, h%mpireal, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, 0.0_num, 1, h%mpireal, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_displacement = h%current_displacement + 2 * sof

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. h%default_rank) THEN
      CALL cfd_safe_write_string(h, mesh_name)
      CALL cfd_safe_write_string(h, mesh_class)
    ENDIF

    h%current_displacement = h%current_displacement + 2 * h%max_string_len

    ! Write the real data

    ALLOCATE(array(1:npart_per_iteration))

    npart_this_cycle = npart_per_iteration
    start = .TRUE.
    sec = 1
    osec = 1
    nsec_left = lengths(sec)
    file_offset = h%current_displacement + offsets(sec) * num

    DO
      CALL iterator(array, npart_this_cycle, start)
      CALL MPI_ALLREDUCE(npart_this_cycle, nmax, 1, MPI_INTEGER8, MPI_MAX, &
          h%comm, errcode)
      IF (nmax .LE. 0) EXIT

      off = 1
      nwrite_left = npart_this_cycle

      IF (start) THEN
        gmn = MINVAL(array(1:npart_this_cycle))
        gmx = MAXVAL(array(1:npart_this_cycle))
        start = .FALSE.
      ELSE
        gmn = MIN(gmn, MINVAL(array(1:npart_this_cycle)))
        gmx = MAX(gmx, MAXVAL(array(1:npart_this_cycle)))
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

        CALL MPI_FILE_SET_VIEW(h%filehandle, file_offset, h%mpireal, &
            h%mpireal, "native", MPI_INFO_NULL, errcode)
        CALL MPI_FILE_WRITE_ALL(h%filehandle, array(off), nelements, &
            h%mpireal, MPI_STATUS_IGNORE, errcode)

        nwrite_left = nwrite_left - nelements
        CALL MPI_ALLREDUCE(nwrite_left, nmax, 1, MPI_INTEGER8, MPI_MAX, &
            h%comm, errcode)

        IF (sec .NE. osec) THEN
          file_offset = h%current_displacement + offsets(sec) * num
        ELSE
          file_offset = file_offset + nelements * num
        ENDIF

        IF (nmax .LE. 0) EXIT

        off = off + nelements
        osec = sec
      ENDDO
    ENDDO

    h%current_displacement = h%current_displacement + npart_global * sof

    DEALLOCATE(array)

    CALL MPI_FILE_SET_VIEW(h%filehandle, offset_for_min_max, h%mpireal, &
        h%mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_ALLREDUCE(gmn, mn, 1, h%mpireal, MPI_MIN, h%comm, errcode)
    CALL MPI_ALLREDUCE(gmx, mx, 1, h%mpireal, MPI_MAX, h%comm, errcode)

    IF (h%rank .EQ. h%default_rank) THEN
      CALL MPI_FILE_WRITE(h%filehandle, mn, 1, h%mpireal, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, mx, 1, h%mpireal, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

  END SUBROUTINE cfd_write_nd_particle_variable_with_iterator_all

END MODULE cfd_output_particle
