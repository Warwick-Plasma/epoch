MODULE boundary

  USE partlist
  USE particle_temperature
  USE deck_io_block

  IMPLICIT NONE

CONTAINS

  SUBROUTINE setup_particle_boundaries

    ! For some types of boundary, fields and particles are treated in
    ! different ways, deal with that here

    IF (bc_x_min_particle .EQ. c_bc_other) bc_x_min_particle = c_bc_reflect
    IF (bc_x_max_particle .EQ. c_bc_other) bc_x_max_particle = c_bc_reflect
    IF (bc_y_min_particle .EQ. c_bc_other) bc_y_min_particle = c_bc_reflect
    IF (bc_y_max_particle .EQ. c_bc_other) bc_y_max_particle = c_bc_reflect
    IF (bc_z_min_particle .EQ. c_bc_other) bc_z_min_particle = c_bc_reflect
    IF (bc_z_max_particle .EQ. c_bc_other) bc_z_max_particle = c_bc_reflect

    IF (bc_x_min_field .EQ. c_bc_other) bc_x_min_field = c_bc_clamp
    IF (bc_x_max_field .EQ. c_bc_other) bc_x_max_field = c_bc_clamp
    IF (bc_y_min_field .EQ. c_bc_other) bc_y_min_field = c_bc_clamp
    IF (bc_y_max_field .EQ. c_bc_other) bc_y_max_field = c_bc_clamp
    IF (bc_z_min_field .EQ. c_bc_other) bc_z_min_field = c_bc_clamp
    IF (bc_z_max_field .EQ. c_bc_other) bc_z_max_field = c_bc_clamp

    ! Note, for laser bcs to work, the main bcs must be set IN THE CODE to
    ! simple_laser (or outflow) and the field bcs to c_bc_clamp. Particles
    ! can then be set separately. IN THE DECK, laser bcs are chosen either
    ! by seting the main bcs OR by setting the field bcs to simple_laser
    ! (or outflow).

    ! Laser boundaries assume open particles unless otherwise specified.
    IF (bc_x_min_particle .EQ. c_bc_simple_laser &
        .OR. bc_x_min_particle .EQ. c_bc_simple_outflow) &
            bc_x_min_particle = c_bc_open
    IF (bc_x_max_particle .EQ. c_bc_simple_laser &
        .OR. bc_x_max_particle .EQ. c_bc_simple_outflow) &
            bc_x_max_particle = c_bc_open
    IF (bc_y_min_particle .EQ. c_bc_simple_laser &
        .OR. bc_y_min_particle .EQ. c_bc_simple_outflow) &
            bc_y_min_particle = c_bc_open
    IF (bc_y_max_particle .EQ. c_bc_simple_laser &
        .OR. bc_y_max_particle .EQ. c_bc_simple_outflow) &
            bc_y_max_particle = c_bc_open
    IF (bc_z_min_particle .EQ. c_bc_simple_laser &
        .OR. bc_z_min_particle .EQ. c_bc_simple_outflow) &
            bc_z_min_particle = c_bc_open
    IF (bc_z_max_particle .EQ. c_bc_simple_laser &
        .OR. bc_z_max_particle .EQ. c_bc_simple_outflow) &
            bc_z_max_particle = c_bc_open

    ! Note: reflecting EM boundaries not yet implemented.
    IF (bc_x_min_field .EQ. c_bc_reflect) bc_x_min_field = c_bc_clamp
    IF (bc_x_max_field .EQ. c_bc_reflect) bc_x_max_field = c_bc_clamp
    IF (bc_y_min_field .EQ. c_bc_reflect) bc_y_min_field = c_bc_clamp
    IF (bc_y_max_field .EQ. c_bc_reflect) bc_y_max_field = c_bc_clamp
    IF (bc_z_min_field .EQ. c_bc_reflect) bc_z_min_field = c_bc_clamp
    IF (bc_z_max_field .EQ. c_bc_reflect) bc_z_max_field = c_bc_clamp

    IF (bc_x_min_field .EQ. c_bc_open) bc_x_min_field = c_bc_simple_outflow
    IF (bc_x_max_field .EQ. c_bc_open) bc_x_min_field = c_bc_simple_outflow
    IF (bc_y_min_field .EQ. c_bc_open) bc_y_min_field = c_bc_simple_outflow
    IF (bc_y_max_field .EQ. c_bc_open) bc_y_max_field = c_bc_simple_outflow
    IF (bc_z_min_field .EQ. c_bc_open) bc_z_min_field = c_bc_simple_outflow
    IF (bc_z_max_field .EQ. c_bc_open) bc_z_max_field = c_bc_simple_outflow

    IF (bc_x_min_particle .EQ. c_bc_open &
        .OR. bc_x_max_particle .EQ. c_bc_open &
        .OR. bc_y_min_particle .EQ. c_bc_open &
        .OR. bc_y_max_particle .EQ. c_bc_open &
        .OR. bc_z_min_particle .EQ. c_bc_open &
        .OR. bc_z_max_particle .EQ. c_bc_open) &
            CALL create_empty_partlist(ejected_particles)

  END SUBROUTINE setup_particle_boundaries



  ! Exchanges field values at processor boundaries and applies field
  ! boundary conditions
  SUBROUTINE field_bc(field)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: field

    CALL do_field_mpi_with_lengths(field, nx, ny, nz)

  END SUBROUTINE field_bc



  SUBROUTINE do_field_mpi_with_lengths(field, nx_local, ny_local, nz_local)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: nx_local, ny_local, nz_local
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER :: subarray

    sizes(1) = nx_local + 6
    sizes(2) = ny_local + 6
    sizes(3) = nz_local + 6
    subsizes(1) = 3
    subsizes(2) = ny_local + 6
    subsizes(3) = nz_local + 6
    starts = 0

    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, subarray, errcode)
    CALL MPI_TYPE_COMMIT(subarray, errcode)

    CALL MPI_SENDRECV(field(1,-2,-2), 1, subarray, proc_x_min, tag, &
        field(nx_local+1,-2,-2), 1, subarray, proc_x_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV(field(nx_local-2,-2,-2), 1, subarray, proc_x_max, tag, &
        field(-2,-2,-2), 1, subarray, proc_x_min, tag, &
        comm, status, errcode)

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = nx_local + 6
    subsizes(2) = 3
    subsizes(3) = nz_local + 6

    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, subarray, errcode)
    CALL MPI_TYPE_COMMIT(subarray, errcode)

    CALL MPI_SENDRECV(field(-2,1,-2), 1, subarray, proc_y_min, tag, &
        field(-2,ny_local+1,-2), 1, subarray, proc_y_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV(field(-2,ny_local-2,-2), 1, subarray, proc_y_max, tag, &
        field(-2,-2,-2), 1, subarray, proc_y_min, tag, &
        comm, status, errcode)

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = nx_local + 6
    subsizes(2) = ny_local + 6
    subsizes(3) = 3

    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, subarray, errcode)
    CALL MPI_TYPE_COMMIT(subarray, errcode)

    CALL MPI_SENDRECV(field(-2,-2,1), 1, subarray, proc_z_min, tag, &
        field(-2,-2,nz_local+1), 1, subarray, proc_z_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV(field(-2,-2,nz_local-2), 1, subarray, proc_z_max, tag, &
        field(-2,-2,-2), 1, subarray, proc_z_min, tag, &
        comm, status, errcode)

    CALL MPI_TYPE_FREE(subarray, errcode)

  END SUBROUTINE do_field_mpi_with_lengths



  SUBROUTINE field_zero_gradient(field, force)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: field
    LOGICAL, INTENT(IN) :: force

    IF ((bc_x_min_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. proc_x_min .EQ. MPI_PROC_NULL) THEN
      field(-1,:,:) = field(2,:,:)
      field( 0,:,:) = field(1,:,:)
    ENDIF

    IF ((bc_x_max_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. proc_x_max .EQ. MPI_PROC_NULL) THEN
      field(nx+1,:,:) = field(nx,:,:)
      field(nx+2,:,:) = field(nx-1,:,:)
    ENDIF

    IF ((bc_y_min_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. proc_y_min .EQ. MPI_PROC_NULL) THEN
      field(:,-1,:) = field(:,2,:)
      field(:, 0,:) = field(:,1,:)
    ENDIF

    IF ((bc_y_max_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. proc_y_max .EQ. MPI_PROC_NULL) THEN
      field(:,ny+1,:) = field(:,ny,:)
      field(:,ny+2,:) = field(:,ny-1,:)
    ENDIF

    IF ((bc_z_min_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. proc_z_min .EQ. MPI_PROC_NULL) THEN
      field(:,:,-1) = field(:,:,2)
      field(:,:, 0) = field(:,:,1)
    ENDIF

    IF ((bc_z_max_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. proc_z_max .EQ. MPI_PROC_NULL) THEN
      field(:,:,nz+1) = field(:,:,nz)
      field(:,:,nz+2) = field(:,:,nz-1)
    ENDIF

  END SUBROUTINE field_zero_gradient



  SUBROUTINE field_clamp_zero(field, s1, s2, s3)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: s1, s2, s3

    ! Use clamp when the laser is on.

    IF ((bc_x_min_field .EQ. c_bc_clamp &
        .OR. bc_x_min_field .EQ. c_bc_simple_laser &
        .OR. bc_x_min_field .EQ. c_bc_simple_outflow) &
        .AND. proc_x_min .EQ. MPI_PROC_NULL) THEN
      IF (s1 .EQ. 1) THEN
        field(-1,:,:) = -field(1,:,:)
        field( 0,:,:) = 0.0_num
      ELSE
        field(-1,:,:) = -field(2,:,:)
        field( 0,:,:) = -field(1,:,:)
      ENDIF
    ENDIF

    IF ((bc_x_max_field .EQ. c_bc_clamp &
        .OR. bc_x_max_field .EQ. c_bc_simple_laser &
        .OR. bc_x_max_field .EQ. c_bc_simple_outflow) &
        .AND. proc_x_max .EQ. MPI_PROC_NULL) THEN
      IF (s1 .EQ. 1) THEN
        field(nx,  :,:) = 0.0_num
        field(nx+1,:,:) = -field(nx-1,:,:)
      ELSE
        field(nx+1,:,:) = -field(nx,  :,:)
        field(nx+2,:,:) = -field(nx-1,:,:)
      ENDIF
    ENDIF

    IF ((bc_y_min_field .EQ. c_bc_clamp &
        .OR. bc_y_min_field .EQ. c_bc_simple_laser &
        .OR. bc_y_min_field .EQ. c_bc_simple_outflow) &
        .AND. proc_y_min .EQ. MPI_PROC_NULL) THEN
      IF (s2 .EQ. 1) THEN
        field(:,-1,:) = -field(:,1,:)
        field(:, 0,:) = 0.0_num
      ELSE
        field(:,-1,:) = -field(:,2,:)
        field(:, 0,:) = -field(:,1,:)
      ENDIF
    ENDIF

    IF ((bc_y_max_field .EQ. c_bc_clamp &
        .OR. bc_y_max_field .EQ. c_bc_simple_laser &
        .OR. bc_y_max_field .EQ. c_bc_simple_outflow) &
        .AND. proc_y_max .EQ. MPI_PROC_NULL) THEN
      IF (s2 .EQ. 1) THEN
        field(:,ny,  :) = 0.0_num
        field(:,ny+1,:) = -field(:,ny-1,:)
      ELSE
        field(:,ny+1,:) = -field(:,ny,  :)
        field(:,ny+2,:) = -field(:,ny-1,:)
      ENDIF
    ENDIF

    IF ((bc_z_min_field .EQ. c_bc_clamp &
        .OR. bc_z_min_field .EQ. c_bc_simple_laser &
        .OR. bc_z_min_field .EQ. c_bc_simple_outflow) &
        .AND. proc_z_min .EQ. MPI_PROC_NULL) THEN
      IF (s3 .EQ. 1) THEN
        field(:,:,-1) = -field(:,:,1)
        field(:,:, 0) = 0.0_num
      ELSE
        field(:,:,-1) = -field(:,:,2)
        field(:,:, 0) = -field(:,:,1)
      ENDIF
    ENDIF

    IF ((bc_z_max_field .EQ. c_bc_clamp &
        .OR. bc_z_max_field .EQ. c_bc_simple_laser &
        .OR. bc_z_max_field .EQ. c_bc_simple_outflow) &
        .AND. proc_z_max .EQ. MPI_PROC_NULL) THEN
      IF (s3 .EQ. 1) THEN
        field(:,:,nz  ) = 0.0_num
        field(:,:,nz+1) = -field(:,:,nz-1)
      ELSE
        field(:,:,nz+1) = -field(:,:,nz  )
        field(:,:,nz+2) = -field(:,:,nz-1)
      ENDIF
    ENDIF

  END SUBROUTINE field_clamp_zero



  SUBROUTINE processor_summation_bcs(array)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER :: subarray

    sizes(1) = nx + 6
    sizes(2) = ny + 6
    sizes(3) = nz + 6
    subsizes(1) = 3
    subsizes(2) = ny + 6
    subsizes(3) = nz + 6
    starts = 0

    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, subarray, errcode)
    CALL MPI_TYPE_COMMIT(subarray, errcode)

    ALLOCATE(temp(3, ny+6, nz+6))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-2,-2,-2), 1, subarray, &
        neighbour(-1,0,0), tag, temp, 3*(ny+6)*(nz+6), mpireal, &
        neighbour( 1,0,0), tag, comm, status, errcode)
    array(-2+nx:nx,:,:) = array(-2+nx:nx,:,:) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(nx+1,-2,-2), 1, subarray, &
        neighbour( 1,0,0), tag, temp, 3*(ny+6)*(nz+6), mpireal, &
        neighbour(-1,0,0), tag, comm, status, errcode)
    array(1:3,:,:) = array(1:3,:,:) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = nx + 6
    subsizes(2) = 3
    subsizes(3) = nz + 6

    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, subarray, errcode)
    CALL MPI_TYPE_COMMIT(subarray, errcode)

    ALLOCATE(temp(nx+6, 3, nz+6))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-2,-2,-2), 1, subarray, &
        neighbour(0,-1,0), tag, temp, 3*(nx+6)*(nz+6), mpireal, &
        neighbour(0, 1,0), tag, comm, status, errcode)
    array(:,-2+ny:ny,:) = array(:,-2+ny:ny,:) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-2,ny+1,-2), 1, subarray, &
        neighbour(0, 1,0), tag, temp, 3*(nx+6)*(nz+6), mpireal, &
        neighbour(0,-1,0), tag, comm, status, errcode)
    array(:,1:3,:) = array(:,1:3,:) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = nx + 6
    subsizes(2) = ny + 6
    subsizes(3) = 3

    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, subarray, errcode)
    CALL MPI_TYPE_COMMIT(subarray, errcode)
    ALLOCATE(temp(nx+6, ny+6, 3))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-2,-2,-2), 1, subarray, &
        neighbour(0,0,-1), tag, temp, 3*(nx+6)*(ny+6), mpireal, &
        neighbour(0,0, 1), tag, comm, status, errcode)
    array(:,:,-2+nz:nz) = array(:,:,-2+nz:nz) + temp

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-2,-2,nz+1), 1, subarray, &
        neighbour(0,0, 1), tag, temp, 3*(nx+6)*(ny+6), mpireal, &
        neighbour(0,0,-1), tag, comm, status, errcode)
    array(:,:,1:3) = array(:,:,1:3) + temp

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    CALL field_bc(array)

  END SUBROUTINE processor_summation_bcs



  SUBROUTINE efield_bcs

    ! These are the MPI boundaries
    CALL field_bc(ex)
    CALL field_bc(ey)
    CALL field_bc(ez)

    ! These apply zero field boundary conditions on the edges
    CALL field_clamp_zero(ex, 1, 0, 0)
    CALL field_clamp_zero(ey, 0, 1, 0)
    CALL field_clamp_zero(ez, 0, 0, 1)

    ! These apply zero field gradient boundary conditions on the edges
    CALL field_zero_gradient(ex, .FALSE.)
    CALL field_zero_gradient(ey, .FALSE.)
    CALL field_zero_gradient(ez, .FALSE.)

  END SUBROUTINE efield_bcs



  SUBROUTINE bfield_bcs(mpi_only)

    LOGICAL, INTENT(IN) :: mpi_only

    ! These are the MPI boundaries
    CALL field_bc(bx)
    CALL field_bc(by)
    CALL field_bc(bz)

    IF (.NOT. mpi_only) THEN
      ! These apply zero field boundary conditions on the edges
      CALL field_clamp_zero(bx, 0, 1, 1)
      CALL field_clamp_zero(by, 1, 0, 1)
      CALL field_clamp_zero(bz, 1, 1, 0)
      ! These apply zero field boundary conditions on the edges
      CALL field_zero_gradient(bx, .FALSE.)
      CALL field_zero_gradient(by, .FALSE.)
      CALL field_zero_gradient(bz, .FALSE.)
    ENDIF

  END SUBROUTINE bfield_bcs



  SUBROUTINE particle_bcs

    TYPE(particle), POINTER :: cur, next
    TYPE(particle_list), DIMENSION(-1:1, -1:1, -1:1) :: send, recv
    INTEGER :: xbd, ybd, zbd
    INTEGER(KIND=8) :: ixp, iyp, izp
    LOGICAL :: out_of_bounds
    INTEGER :: ispecies

    DO ispecies = 1, n_species
      cur=>particle_species(ispecies)%attached_list%head

      DO iz = -1, 1
        DO iy = -1, 1
          DO ix = -1, 1
            IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
            CALL create_empty_partlist(send(ix, iy, iz))
            CALL create_empty_partlist(recv(ix, iy, iz))
          ENDDO
        ENDDO
      ENDDO

      DO WHILE (ASSOCIATED(cur))
        next=>cur%next

        xbd = 0
        ybd = 0
        zbd = 0
        out_of_bounds = .FALSE.

        ! Particle has left this processor
        IF (cur%part_pos(1) .LT. x_min_local - dx / 2.0_num) THEN
          xbd = -1
          ! Particle has left the system
          IF (coordinates(c_ndims) .EQ. 0) THEN
            IF (bc_x_min_particle .EQ. c_bc_open) THEN
              out_of_bounds = .TRUE.
            ELSE IF (bc_x_min_particle .EQ. c_bc_reflect) THEN
              cur%part_pos(1) = 2.0_num * x_min - dx - cur%part_pos(1)
              cur%part_p(1) = -cur%part_p(1)
            ELSE IF (bc_x_min_particle .EQ. c_bc_periodic) THEN
              cur%part_pos(1) = cur%part_pos(1) + (length_x + dx)
            ENDIF
          ENDIF

        ! Particle has left this processor
        ELSE IF (cur%part_pos(1) .GE. x_max_local + dx / 2.0_num) THEN
          xbd = 1
          ! Particle has left the system
          IF (coordinates(c_ndims) .EQ. nprocx - 1) THEN
            IF (bc_x_max_particle .EQ. c_bc_open) THEN
              out_of_bounds = .TRUE.
            ELSE IF (bc_x_max_particle .EQ. c_bc_reflect) THEN
              cur%part_pos(1) = 2.0_num * x_max + dx - cur%part_pos(1)
              cur%part_p(1) = -cur%part_p(1)
            ELSE IF (bc_x_max_particle .EQ. c_bc_periodic) THEN
              cur%part_pos(1) = cur%part_pos(1) - (length_x + dx)
            ENDIF
          ENDIF
        ENDIF

        ! Particle has left this processor
        IF (cur%part_pos(2) .LT. y_min_local - dy / 2.0_num) THEN
          ybd = -1
          ! Particle has left the system
          IF (coordinates(c_ndims-1) .EQ. 0) THEN
            IF (bc_y_min_particle .EQ. c_bc_open) THEN
              out_of_bounds = .TRUE.
            ELSE IF (bc_y_min_particle .EQ. c_bc_reflect) THEN
              cur%part_pos(2) = 2.0_num * y_min - dy - cur%part_pos(2)
              cur%part_p(2) = -cur%part_p(2)
            ELSE IF (bc_y_min_particle .EQ. c_bc_periodic) THEN
              cur%part_pos(2) = cur%part_pos(2) + (length_y + dy)
            ENDIF
          ENDIF

        ! Particle has left this processor
        ELSE IF (cur%part_pos(2) .GE. y_max_local + dy / 2.0_num) THEN
          ybd = 1
          ! Particle has left the system
          IF (coordinates(c_ndims-1) .EQ. nprocy - 1) THEN
            IF (bc_y_max_particle .EQ. c_bc_open) THEN
              out_of_bounds = .TRUE.
            ELSE IF (bc_y_max_particle .EQ. c_bc_reflect) THEN
              cur%part_pos(2) = 2.0_num * y_max + dy - cur%part_pos(2)
              cur%part_p(2) = -cur%part_p(2)
            ELSE IF (bc_y_max_particle .EQ. c_bc_periodic) THEN
              cur%part_pos(2) = cur%part_pos(2) - (length_y + dy)
            ENDIF
          ENDIF
        ENDIF

        ! Particle has left this processor
        IF (cur%part_pos(3) .LT. z_min_local - dz / 2.0_num) THEN
          zbd = -1
          ! Particle has left the system
          IF (coordinates(c_ndims-2) .EQ. 0) THEN
            IF (bc_z_min_particle .EQ. c_bc_open) THEN
              out_of_bounds = .TRUE.
            ELSE IF (bc_z_min_particle .EQ. c_bc_reflect) THEN
              cur%part_pos(3) = 2.0_num * z_min - dz - cur%part_pos(3)
              cur%part_p(3) = -cur%part_p(3)
            ELSE IF (bc_z_min_particle .EQ. c_bc_periodic) THEN
              cur%part_pos(3) = cur%part_pos(3) + (length_z + dz)
            ENDIF
          ENDIF

        ! Particle has left this processor
        ELSE IF (cur%part_pos(3) .GE. z_max_local + dz / 2.0_num) THEN
          zbd = 1
          ! Particle has left the system
          IF (coordinates(c_ndims-2) .EQ. nprocz - 1) THEN
            IF (bc_z_max_particle .EQ. c_bc_open) THEN
              out_of_bounds = .TRUE.
            ELSE IF (bc_z_max_particle .EQ. c_bc_reflect) THEN
              cur%part_pos(3) = 2.0_num * z_max + dz - cur%part_pos(3)
              cur%part_p(3) = -cur%part_p(3)
            ELSE IF (bc_z_max_particle .EQ. c_bc_periodic) THEN
              cur%part_pos(3) = cur%part_pos(3) - (length_z + dz)
            ENDIF
          ENDIF
        ENDIF

        IF (out_of_bounds) THEN
          ! Particle has gone forever
          CALL remove_particle_from_partlist(&
              particle_species(ispecies)%attached_list, cur)
          IF (dumpmask(c_dump_ejected_particles) .NE. c_io_never) THEN
            CALL add_particle_to_partlist(ejected_particles, cur)
          ELSE
            DEALLOCATE(cur)
          ENDIF
        ELSE IF (ABS(xbd) + ABS(ybd) + ABS(zbd) .GT. 0) THEN
          ! Particle has left processor, send it to its neighbour
          CALL remove_particle_from_partlist(&
              particle_species(ispecies)%attached_list, cur)
          CALL add_particle_to_partlist(send(xbd, ybd, zbd), cur)
        ENDIF

        ! Move to next particle
        cur=>next
      ENDDO

      ! swap Particles
      DO iz = -1, 1
        DO iy = -1, 1
          DO ix = -1, 1
            IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
            ixp = -ix
            iyp = -iy
            izp = -iz
            CALL partlist_sendrecv(send(ix, iy, iz), recv(ixp, iyp, izp), &
                neighbour(ix, iy, iz), neighbour(ixp, iyp, izp))
            CALL append_partlist(particle_species(ispecies)%attached_list, &
                recv(ixp, iyp, izp))
          ENDDO
        ENDDO
      ENDDO

      DO iz = -1, 1
        DO iy = -1, 1
          DO ix = -1, 1
            IF (ABS(ix) + ABS(iy) + ABS(iz) .EQ. 0) CYCLE
            CALL destroy_partlist(send(ix, iy, iz))
            CALL destroy_partlist(recv(ix, iy, iz))
          ENDDO
        ENDDO
      ENDDO

    ENDDO

  END SUBROUTINE particle_bcs

END MODULE boundary
