MODULE output_particle

  USE iocommon
  USE output

  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! Code to write a 2D Cartesian grid in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify nx, ny
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_nd_particle_grid_all(name, class, particles, &
      npart_global, particle_coord_type, particle_type)

    REAL(num), DIMENSION(:,:), INTENT(IN) :: particles
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    INTEGER(8), INTENT(IN) :: npart_global
    INTEGER(4), INTENT(IN) :: particle_coord_type
    INTEGER, INTENT(IN) :: particle_type
    INTEGER(8) :: npart_local
    INTEGER(8) :: block_length, md_length
    INTEGER(4) :: ndim, i, disp0
    INTEGER(4) :: sizes(2)
    REAL(num) :: mn, mx

    sizes = SHAPE(particles)
    npart_local = sizes(2)
    ndim = sizes(1)

    ! Metadata is
    ! * ) meshtype (INTEGER(4)) All mesh blocks contain this
    ! * ) nd    INTEGER(4)
    ! * ) sof   INTEGER(4)
    ! Specific to particle mesh
    ! 1 ) ct    INTEGER(4)
    ! 2 ) npart INTEGER(8)
    ! 3 ) d1min REAL(num)
    ! 4 ) d1max REAL(num)
    ! 5 ) d2min REAL(num)
    ! 6 ) d2max REAL(num)
    ! .
    ! .
    ! .
    ! n ) dnmin REAL(num)
    ! n+1) dnmax REAL(num)

    ! 1 INT, 1 INT8, 2REAL per Dim
    md_length = meshtype_header_offset + 1 * soi + 1 * soi8  + ndim * 2 * num
    block_length = md_length + num * ndim * npart_global

    ! Now written header, write metadata
    CALL cfd_write_block_header(name, class, c_type_mesh, block_length, &
        md_length, default_rank)

    disp0 = current_displacement
    CALL cfd_write_meshtype_header(c_mesh_particle, ndim, num, default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, &
        MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, particle_coord_type, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 1 * soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER8, MPI_INTEGER8, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, npart_global, 1, MPI_INTEGER8, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 1 * soi8

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    DO i = 1, ndim
      CALL MPI_ALLREDUCE(MINVAL(particles(:,i)), mn, 1, mpireal, MPI_MIN, &
          cfd_comm, cfd_errcode)
      CALL MPI_ALLREDUCE(MAXVAL(particles(:,i)), mx, 1, mpireal, MPI_MAX, &
          cfd_comm, cfd_errcode)

      IF (cfd_rank .EQ. default_rank) THEN
        CALL MPI_FILE_WRITE(cfd_filehandle, mn, 1, mpireal, cfd_status, &
            cfd_errcode)
        CALL MPI_FILE_WRITE(cfd_filehandle, mx, 1, mpireal, cfd_status, &
            cfd_errcode)
      ENDIF

      current_displacement = current_displacement + 2 * num
    ENDDO

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        particle_type, "native", MPI_INFO_NULL, cfd_errcode)

    ! Write the real data
    CALL MPI_FILE_WRITE_ALL(cfd_filehandle, particles, npart_local * ndim, &
        mpireal, cfd_status, cfd_errcode)

    current_displacement = current_displacement + ndim * npart_global * num

  END SUBROUTINE cfd_write_nd_particle_grid_all



  !----------------------------------------------------------------------------
  ! Code to write a 2D Cartesian grid in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify nx, ny
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_nd_particle_grid_with_iterator_all(name, class, &
      iterator, ndims, npart_local, npart_global, npart_per_iteration, &
      particle_coord_type, particle_type)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    INTEGER(8), INTENT(IN) :: npart_global
    INTEGER(8), INTENT(IN) :: npart_local
    INTEGER(8), INTENT(IN) :: npart_per_iteration
    INTEGER(4), INTENT(IN) :: ndims
    INTEGER(4), INTENT(IN) :: particle_coord_type
    INTEGER, INTENT(IN) :: particle_type
    REAL(num), ALLOCATABLE, DIMENSION(:) :: data

    INTERFACE
      SUBROUTINE iterator(data, npart_it, direction, start)
        USE shared_data
        REAL(num), DIMENSION(:), INTENT(INOUT) :: data
        INTEGER, INTENT(IN) :: direction
        INTEGER(8), INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
      END SUBROUTINE iterator
    END INTERFACE

    INTEGER(8) :: block_length, md_length, npart_this_cycle, npart_sent
    INTEGER(4) :: idim
    INTEGER(MPI_OFFSET_KIND) :: offset_for_min_max
    REAL(num) :: mn, mx
    REAL(num), ALLOCATABLE, DIMENSION(:,:) :: min_max
    LOGICAL :: start
    INTEGER :: mpi_request

    ! Metadata is
    ! * ) meshtype (INTEGER(4)) All mesh blocks contain this
    ! * ) nd    INTEGER(4)
    ! * ) sof   INTEGER(4)
    ! Specific to particle mesh
    ! 1 ) ct    INTEGER(4)
    ! 2 ) npart INTEGER(8)
    ! 3 ) d1min REAL(num)
    ! 4 ) d1max REAL(num)
    ! 5 ) d2min REAL(num)
    ! 6 ) d2max REAL(num)
    ! .
    ! .
    ! .
    ! n ) dnmin REAL(num)
    ! n+1) dnmax REAL(num)

    ! 1 INT, 1 INT8, 2REAL per Dim
    md_length = meshtype_header_offset + 1 * soi + 1 * soi8  + ndims * 2 * num
    block_length = md_length + num * ndims * npart_global

    ALLOCATE(min_max(1:ndims, 1:2))
    min_max = 0.0_num

    ! Now written header, write metadata
    CALL cfd_write_block_header(name, class, c_type_mesh, block_length, &
        md_length, default_rank)
    CALL cfd_write_meshtype_header(c_mesh_particle, ndims, num, default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, particle_coord_type, 1, &
          MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 1 * soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER8, MPI_INTEGER8, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, npart_global, 1, &
          MPI_INTEGER8, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 1 * soi8

    ! This is to skip past the location for the min/max values (Just write
    ! zeros). They will be filled in later
    offset_for_min_max = current_displacement
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, offset_for_min_max, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, min_max, ndims * 2, &
          mpireal, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 2 * ndims * num

    ! Write the real data

    start = .TRUE.
    ALLOCATE(data(1:npart_per_iteration))
    npart_sent = 0

    DO idim = 1, ndims
      CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
          particle_type, "native", MPI_INFO_NULL, cfd_errcode)
      npart_this_cycle = npart_per_iteration
      start = .TRUE.

      DO
        CALL iterator(data, npart_this_cycle, idim, start)
        IF (npart_this_cycle .LE. 0) EXIT

        IF (start) THEN
          min_max(idim, 1) = MINVAL(data(1:npart_this_cycle))
          min_max(idim, 2) = MAXVAL(data(1:npart_this_cycle))
        ELSE
          min_max(idim, 1) = MIN(min_max(idim, 1), &
              MINVAL(data(1:npart_this_cycle)))
          min_max(idim, 2) = MAX(min_max(idim, 2), &
              MAXVAL(data(1:npart_this_cycle)))
        ENDIF

        start = .FALSE.
        npart_sent = npart_sent + npart_this_cycle
        CALL MPI_FILE_IWRITE(cfd_filehandle, data, npart_this_cycle, &
            mpireal, mpi_request, cfd_errcode)
      ENDDO

      current_displacement = current_displacement +  npart_global * num
    ENDDO

    CALL MPI_WAIT(mpi_request, cfd_status, cfd_errcode)
    DEALLOCATE(data)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, offset_for_min_max, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    DO idim = 1, ndims
      CALL MPI_ALLREDUCE(min_max(idim, 1), mn, 1, mpireal, MPI_MIN, &
          cfd_comm, cfd_errcode)
      CALL MPI_ALLREDUCE(min_max(idim, 2), mx, 1, mpireal, MPI_MAX, &
          cfd_comm, cfd_errcode)

      IF (cfd_rank .EQ. default_rank) THEN
        CALL MPI_FILE_WRITE(cfd_filehandle, mn, 1, mpireal, &
            cfd_status, cfd_errcode)
        CALL MPI_FILE_WRITE(cfd_filehandle, mx, 1, mpireal, &
            cfd_status, cfd_errcode)
      ENDIF
    ENDDO

    DEALLOCATE(min_max)

    CALL MPI_BARRIER(comm, errcode)

  END SUBROUTINE cfd_write_nd_particle_grid_with_iterator_all



  !----------------------------------------------------------------------------
  ! Code to write a 2D Cartesian grid in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify nx, ny
  !----------------------------------------------------------------------------

  SUBROUTINE cfd_write_nd_particle_variable_all(name, class, particles, &
      npart_global, mesh_name, mesh_class, particle_type)

    REAL(num), DIMENSION(:), INTENT(IN) :: particles
    CHARACTER(LEN=*), INTENT(IN) :: name, class, mesh_name, mesh_class
    INTEGER, INTENT(IN) :: particle_type
    INTEGER(8), INTENT(IN) :: npart_global
    INTEGER(8) :: npart_local
    INTEGER(8) :: block_length, md_length
    REAL(num) :: mn, mx

    npart_local = SIZE(particles)

    ! Metadata is
    ! * ) meshtype (INTEGER(4)) All mesh blocks contain this
    ! * ) nd     INTEGER(4)
    ! * ) sof    INTEGER(4)
    ! Specific to particle variable
    ! 1 ) npart  INTEGER(8)
    ! 2 ) vmin   REAL(num)
    ! 3 ) vmax   REAL(num)
    ! 4 ) mesh   CHARACTER
    ! 5 ) mclass CHARACTER

    md_length = meshtype_header_offset + 1 * soi8 + 2 * num + 2 * max_string_len
    block_length = md_length + num * npart_global

    ! Now written header, write metadata
    CALL cfd_write_block_header(name, class, c_type_mesh_variable, &
        block_length, md_length, default_rank)
    CALL cfd_write_meshtype_header(c_var_particle, c_dimension_irrelevant, &
        num, default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER8, MPI_INTEGER8, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, npart_global, 1, MPI_INTEGER8, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 1 * soi8

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_ALLREDUCE(MINVAL(particles), mn, 1, mpireal, MPI_MIN, cfd_comm, &
        cfd_errcode)
    CALL MPI_ALLREDUCE(MAXVAL(particles), mx, 1, mpireal, MPI_MAX, cfd_comm, &
        cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, mn, 1, mpireal, cfd_status, &
          cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx, 1, mpireal, cfd_status, &
          cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 2 * num

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF

    current_displacement = current_displacement + 2 * max_string_len

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        particle_type, "native", MPI_INFO_NULL, cfd_errcode)

    ! Write the real data
    CALL MPI_FILE_WRITE_ALL(cfd_filehandle, particles, npart_local, mpireal, &
        cfd_status, cfd_errcode)

    current_displacement = current_displacement + npart_global * num

  END SUBROUTINE cfd_write_nd_particle_variable_all



  SUBROUTINE cfd_write_nd_particle_variable_with_iterator_all(name, class, &
      iterator, npart_global, npart_per_iteration, mesh_name, mesh_class, &
      particle_type)

    CHARACTER(LEN=*), INTENT(IN) :: name, class, mesh_name, mesh_class
    INTEGER, INTENT(IN) :: particle_type
    INTEGER(8), INTENT(IN) :: npart_global, npart_per_iteration
    INTEGER(8) :: npart_this_cycle
    REAL(num), ALLOCATABLE, DIMENSION(:) :: data

    INTERFACE
      SUBROUTINE iterator(data, npart_it, start)
        USE shared_data
        REAL(num), DIMENSION(:), INTENT(INOUT) :: data
        INTEGER(8), INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
      END SUBROUTINE iterator
    END INTERFACE

    INTEGER(8) :: block_length, md_length
    REAL(num) :: mn, mx, mn_g, mx_g
    INTEGER(MPI_OFFSET_KIND) :: offset_for_min_max
    LOGICAL :: start
    INTEGER :: mpi_request

    ! Metadata is
    ! * ) meshtype (INTEGER(4)) All mesh blocks contain this
    ! * ) nd     INTEGER(4)
    ! * ) sof    INTEGER(4)
    ! Specific to particle variable
    ! 1 ) npart  INTEGER(8)
    ! 2 ) vmin   REAL(num)
    ! 3 ) vmax   REAL(num)
    ! 4 ) mesh   CHARACTER
    ! 5 ) mclass CHARACTER

    md_length = meshtype_header_offset + 1 * soi8 + 2 * num + 2 * max_string_len
    block_length = md_length + num * npart_global

    ! Now written header, write metadata
    CALL cfd_write_block_header(name, class, c_type_mesh_variable, &
        block_length, md_length, default_rank)
    CALL cfd_write_meshtype_header(c_var_particle, c_dimension_irrelevant, &
        num, default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER8, MPI_INTEGER8, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, npart_global, 1, MPI_INTEGER8, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 1 * soi8

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        mpireal, mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, 0.0_num, 1, mpireal, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, 0.0_num, 1, mpireal, &
          cfd_status, cfd_errcode)
    ENDIF

    offset_for_min_max = current_displacement
    current_displacement = current_displacement + 2 * num

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF

    current_displacement = current_displacement + 2 * max_string_len

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        mpireal, particle_type, "native", MPI_INFO_NULL, cfd_errcode)

    start = .TRUE.
    npart_this_cycle = npart_per_iteration
    ALLOCATE(data(1:npart_per_iteration))

    DO
      data = 27.224_num
      CALL iterator(data, npart_this_cycle, start)
      IF (npart_this_cycle .LE. 0) EXIT

      IF (start) THEN
        mn = MINVAL(data(1:npart_this_cycle))
        mx = MAXVAL(data(1:npart_this_cycle))
      ELSE
        mn = MIN(mn, MINVAL(data(1:npart_this_cycle)))
        mx = MAX(mx, MAXVAL(data(1:npart_this_cycle)))
      ENDIF
      start = .FALSE.

      CALL MPI_FILE_IWRITE(cfd_filehandle, data, npart_this_cycle, mpireal, &
          mpi_request, cfd_errcode)
    ENDDO

    CALL MPI_WAIT(mpi_request, cfd_status, cfd_errcode)
    DEALLOCATE(data)

    current_displacement = current_displacement + npart_global * num

    CALL MPI_ALLREDUCE(mn, mn_g, 1, mpireal, MPI_MIN, cfd_comm, cfd_errcode)
    CALL MPI_ALLREDUCE(mx, mx_g, 1, mpireal, MPI_MAX, cfd_comm, cfd_errcode)
    mn = mn_g
    mx = mx_g

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, offset_for_min_max, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, mn, 1, mpireal, cfd_status, &
          cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, mx, 1, mpireal, cfd_status, &
          cfd_errcode)
    ENDIF

    CALL MPI_BARRIER(comm, errcode)

  END SUBROUTINE cfd_write_nd_particle_variable_with_iterator_all

END MODULE output_particle
