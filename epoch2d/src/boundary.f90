MODULE boundary

  USE mpi
  USE partlist
  USE particle_temperature
  USE deck_io_block
  USE laser

  IMPLICIT NONE

CONTAINS

  SUBROUTINE setup_particle_boundaries

    INTEGER :: i
    LOGICAL :: error
    CHARACTER(LEN=5), DIMENSION(2*c_ndims) :: &
        boundary = (/ "x_min", "x_max", "y_min", "y_max" /)

    ! For some types of boundary, fields and particles are treated in
    ! different ways, deal with that here

    DO i = 1, 2*c_ndims
      IF (bc_particle(i) .EQ. c_bc_other) bc_particle(i) = c_bc_reflect
      IF (bc_field(i) .EQ. c_bc_other) bc_field(i) = c_bc_clamp
      IF (bc_field(i) .EQ. c_bc_simple_laser) add_laser(i) = .TRUE.
    ENDDO

    ! Note, for laser bcs to work, the main bcs must be set IN THE CODE to
    ! simple_laser (or outflow) and the field bcs to c_bc_clamp. Particles
    ! can then be set separately. IN THE DECK, laser bcs are chosen either
    ! by seting the main bcs OR by setting the field bcs to simple_laser
    ! (or outflow).

    ! Laser boundaries assume open particles unless otherwise specified.
    DO i = 1, 2*c_ndims
      IF (bc_particle(i) .EQ. c_bc_simple_laser &
          .OR. bc_particle(i) .EQ. c_bc_simple_outflow) &
              bc_particle(i) = c_bc_open
    ENDDO

    ! Note: reflecting EM boundaries not yet implemented.
    DO i = 1, 2*c_ndims
      IF (bc_field(i) .EQ. c_bc_reflect) bc_field(i) = c_bc_clamp
      IF (bc_field(i) .EQ. c_bc_open) bc_field(i) = c_bc_simple_outflow
    ENDDO

    ! Sanity check on particle boundaries
    error = .FALSE.
    DO i = 1, 2*c_ndims
      IF (bc_particle(i) .EQ. c_bc_periodic &
          .OR. bc_particle(i) .EQ. c_bc_reflect &
          .OR. bc_particle(i) .EQ. c_bc_thermal &
          .OR. bc_particle(i) .EQ. c_bc_open) CYCLE
      IF (rank .EQ. 0) THEN
        WRITE(*,*)
        WRITE(*,*) '*** ERROR ***'
        WRITE(*,*) 'Unrecognised particle boundary condition on "', &
            boundary(i), '" boundary.'
      ENDIF
      error = .TRUE.
    ENDDO

    IF (error) CALL MPI_ABORT(MPI_COMM_WORLD, errcode, errcode)

  END SUBROUTINE setup_particle_boundaries



  ! Exchanges field values at processor boundaries and applies field
  ! boundary conditions
  SUBROUTINE field_bc(field)

    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: field

    CALL do_field_mpi_with_lengths(field, nx, ny)

  END SUBROUTINE field_bc



  SUBROUTINE do_field_mpi_with_lengths(field, nx_local, ny_local)

    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: nx_local, ny_local
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER :: subarray

    sizes(1) = nx_local + 6
    sizes(2) = ny_local + 6
    subsizes(1) = 3
    subsizes(2) = ny_local + 6
    starts = 0

    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, subarray, errcode)
    CALL MPI_TYPE_COMMIT(subarray, errcode)

    CALL MPI_SENDRECV(field(1,-2), 1, subarray, proc_x_min, tag, &
        field(nx_local+1,-2), 1, subarray, proc_x_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV(field(nx_local-2,-2), 1, subarray, proc_x_max, tag, &
        field(-2,-2), 1, subarray, proc_x_min, tag, &
        comm, status, errcode)

    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = nx_local + 6
    subsizes(2) = 3

    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, subarray, errcode)
    CALL MPI_TYPE_COMMIT(subarray, errcode)

    CALL MPI_SENDRECV(field(-2,1), 1, subarray, proc_y_min, tag, &
        field(-2,ny_local+1), 1, subarray, proc_y_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV(field(-2,ny_local-2), 1, subarray, proc_y_max, tag, &
        field(-2,-2), 1, subarray, proc_y_min, tag, &
        comm, status, errcode)

    CALL MPI_TYPE_FREE(subarray, errcode)

  END SUBROUTINE do_field_mpi_with_lengths



  SUBROUTINE field_zero_gradient(field, stagger_type, boundary)

    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: stagger_type, boundary

    IF (bc_field(boundary) .EQ. c_bc_periodic) RETURN

    IF (boundary .EQ. c_bd_x_min .AND. x_min_boundary) THEN
      IF (stagger(c_dir_x,stagger_type)) THEN
        field(-1,:) = field(1,:)
      ELSE
        field(-1,:) = field(2,:)
        field( 0,:) = field(1,:)
      ENDIF
    ELSE IF (boundary .EQ. c_bd_x_max .AND. x_max_boundary) THEN
      IF (stagger(c_dir_x,stagger_type)) THEN
        field(nx+1,:) = field(nx-1,:)
      ELSE
        field(nx+1,:) = field(nx  ,:)
        field(nx+2,:) = field(nx-1,:)
      ENDIF
    ELSE IF (boundary .EQ. c_bd_y_min .AND. y_min_boundary) THEN
      IF (stagger(c_dir_y,stagger_type)) THEN
        field(:,-1) = field(:,1)
      ELSE
        field(:,-1) = field(:,2)
        field(:, 0) = field(:,1)
      ENDIF
    ELSE IF (boundary .EQ. c_bd_y_max .AND. y_max_boundary) THEN
      IF (stagger(c_dir_y,stagger_type)) THEN
        field(:,ny+1) = field(:,ny-1)
      ELSE
        field(:,ny+1) = field(:,ny  )
        field(:,ny+2) = field(:,ny-1)
      ENDIF
    ENDIF

  END SUBROUTINE field_zero_gradient



  SUBROUTINE field_clamp_zero(field, stagger_type, boundary)

    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: stagger_type, boundary

    IF (bc_field(boundary) .EQ. c_bc_periodic) RETURN

    IF (boundary .EQ. c_bd_x_min .AND. x_min_boundary) THEN
      IF (stagger(c_dir_x,stagger_type)) THEN
        field(-1,:) = -field(1,:)
        field( 0,:) = 0.0_num
      ELSE
        field(-1,:) = -field(2,:)
        field( 0,:) = -field(1,:)
      ENDIF
    ELSE IF (boundary .EQ. c_bd_x_max .AND. x_max_boundary) THEN
      IF (stagger(c_dir_x,stagger_type)) THEN
        field(nx  ,:) = 0.0_num
        field(nx+1,:) = -field(nx-1,:)
      ELSE
        field(nx+1,:) = -field(nx  ,:)
        field(nx+2,:) = -field(nx-1,:)
      ENDIF
    ELSE IF (boundary .EQ. c_bd_y_min .AND. y_min_boundary) THEN
      IF (stagger(c_dir_y,stagger_type)) THEN
        field(:,-1) = -field(:,1)
        field(:, 0) = 0.0_num
      ELSE
        field(:,-1) = -field(:,2)
        field(:, 0) = -field(:,1)
      ENDIF
    ELSE IF (boundary .EQ. c_bd_y_max .AND. y_max_boundary) THEN
      IF (stagger(c_dir_y,stagger_type)) THEN
        field(:,ny  ) = 0.0_num
        field(:,ny+1) = -field(:,ny-1)
      ELSE
        field(:,ny+1) = -field(:,ny  )
        field(:,ny+2) = -field(:,ny-1)
      ENDIF
    ENDIF

  END SUBROUTINE field_clamp_zero



  SUBROUTINE processor_summation_bcs(array, flip_direction)

    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: array
    INTEGER, INTENT(IN), OPTIONAL :: flip_direction
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temp
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes, starts
    INTEGER :: subarray, sgn

    sizes(1) = nx + 6
    sizes(2) = ny + 6
    subsizes(1) = 3
    subsizes(2) = ny + 6
    starts = 0

    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, subarray, errcode)
    CALL MPI_TYPE_COMMIT(subarray, errcode)

    ALLOCATE(temp(3, ny+6))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(nx+1,-2), 1, subarray, &
        neighbour( 1,0), tag, temp, 3*(ny+6), mpireal, &
        neighbour(-1,0), tag, comm, status, errcode)

    ! Deal with reflecting boundaries differently
    IF ((bc_particle(c_bd_x_min) .EQ. c_bc_reflect &
        .OR. bc_particle(c_bd_x_min) .EQ. c_bc_thermal) &
        .AND. x_min_boundary) THEN
      sgn = 1
      IF (PRESENT(flip_direction)) THEN
        ! Currents get reversed in the direction of the boundary
        IF (flip_direction .EQ. c_dir_x) sgn = -1
      ENDIF
      array(1,:) = array(1,:) + sgn * array( 0,:)
      array(2,:) = array(2,:) + sgn * array(-1,:)
      array(3,:) = array(3,:) + sgn * array(-2,:)
    ELSE
      array(1:3,:) = array(1:3,:) + temp
    ENDIF

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-2,-2), 1, subarray, &
        neighbour(-1,0), tag, temp, 3*(ny+6), mpireal, &
        neighbour( 1,0), tag, comm, status, errcode)

    ! Deal with reflecting boundaries differently
    IF ((bc_particle(c_bd_x_max) .EQ. c_bc_reflect &
        .OR. bc_particle(c_bd_x_max) .EQ. c_bc_thermal) &
        .AND. x_max_boundary) THEN
      sgn = 1
      IF (PRESENT(flip_direction)) THEN
        ! Currents get reversed in the direction of the boundary
        IF (flip_direction .EQ. c_dir_x) sgn = -1
      ENDIF
      array(nx-2,:) = array(nx-2,:) + sgn * array(nx+3,:)
      array(nx-1,:) = array(nx-1,:) + sgn * array(nx+2,:)
      array(nx  ,:) = array(nx  ,:) + sgn * array(nx+1,:)
    ELSE
      array(nx-2:nx,:) = array(nx-2:nx,:) + temp
    ENDIF

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    subsizes(1) = nx + 6
    subsizes(2) = 3

    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, subarray, errcode)
    CALL MPI_TYPE_COMMIT(subarray, errcode)

    ALLOCATE(temp(nx+6, 3))

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-2,ny+1), 1, subarray, &
        neighbour(0, 1), tag, temp, 3*(nx+6), mpireal, &
        neighbour(0,-1), tag, comm, status, errcode)

    ! Deal with reflecting boundaries differently
    IF ((bc_particle(c_bd_y_min) .EQ. c_bc_reflect &
        .OR. bc_particle(c_bd_y_min) .EQ. c_bc_thermal) &
        .AND. y_min_boundary) THEN
      sgn = 1
      IF (PRESENT(flip_direction)) THEN
        ! Currents get reversed in the direction of the boundary
        IF (flip_direction .EQ. c_dir_y) sgn = -1
      ENDIF
      array(:,1) = array(:,1) + sgn * array(:, 0)
      array(:,2) = array(:,2) + sgn * array(:,-1)
      array(:,3) = array(:,3) + sgn * array(:,-2)
    ELSE
      array(:,1:3) = array(:,1:3) + temp
    ENDIF

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-2,-2), 1, subarray, &
        neighbour(0,-1), tag, temp, 3*(nx+6), mpireal, &
        neighbour(0, 1), tag, comm, status, errcode)

    ! Deal with reflecting boundaries differently
    IF ((bc_particle(c_bd_y_max) .EQ. c_bc_reflect &
        .OR. bc_particle(c_bd_y_max) .EQ. c_bc_thermal) &
        .AND. y_max_boundary) THEN
      sgn = 1
      IF (PRESENT(flip_direction)) THEN
        ! Currents get reversed in the direction of the boundary
        IF (flip_direction .EQ. c_dir_y) sgn = -1
      ENDIF
      array(:,ny-2) = array(:,ny-2) + sgn * array(:,ny+3)
      array(:,ny-1) = array(:,ny-1) + sgn * array(:,ny+2)
      array(:,ny  ) = array(:,ny  ) + sgn * array(:,ny+1)
    ELSE
      array(:,ny-2:ny) = array(:,ny-2:ny) + temp
    ENDIF

    DEALLOCATE(temp)
    CALL MPI_TYPE_FREE(subarray, errcode)

    CALL field_bc(array)

  END SUBROUTINE processor_summation_bcs



  SUBROUTINE efield_bcs

    INTEGER :: i

    ! These are the MPI boundaries
    CALL field_bc(ex)
    CALL field_bc(ey)
    CALL field_bc(ez)

    ! Perfectly conducting boundaries
    DO i = c_bd_x_min, c_bd_x_max, c_bd_x_max - c_bd_x_min
      IF (bc_field(i) .EQ. c_bc_conduct) THEN
        CALL field_clamp_zero(ey, c_stagger_ey, i)
        CALL field_clamp_zero(ez, c_stagger_ez, i)
      ENDIF
    ENDDO

    DO i = c_bd_y_min, c_bd_y_max, c_bd_y_max - c_bd_y_min
      IF (bc_field(i) .EQ. c_bc_conduct) THEN
        CALL field_clamp_zero(ex, c_stagger_ex, i)
        CALL field_clamp_zero(ez, c_stagger_ez, i)
      ENDIF
    ENDDO

    DO i = 1, 2*c_ndims
      IF (bc_field(i) .EQ. c_bc_clamp &
          .OR. bc_field(i) .EQ. c_bc_simple_laser &
          .OR. bc_field(i) .EQ. c_bc_simple_outflow) THEN
        ! These apply zero field boundary conditions on the edges
        CALL field_clamp_zero(ex, c_stagger_ex, i)
        CALL field_clamp_zero(ey, c_stagger_ey, i)
        CALL field_clamp_zero(ez, c_stagger_ez, i)
      ENDIF

      ! These apply zero field gradient boundary conditions on the edges
      IF (bc_field(i) .EQ. c_bc_zero_gradient) THEN
        CALL field_zero_gradient(ex, c_stagger_ex, i)
        CALL field_zero_gradient(ey, c_stagger_ey, i)
        CALL field_zero_gradient(ez, c_stagger_ez, i)
      ENDIF
    ENDDO

  END SUBROUTINE efield_bcs



  SUBROUTINE bfield_bcs(mpi_only)

    LOGICAL, INTENT(IN) :: mpi_only
    INTEGER :: i

    ! These are the MPI boundaries
    CALL field_bc(bx)
    CALL field_bc(by)
    CALL field_bc(bz)

    IF (mpi_only) RETURN

    ! Perfectly conducting boundaries
    DO i = c_bd_x_min, c_bd_x_max, c_bd_x_max - c_bd_x_min
      IF (bc_field(i) .EQ. c_bc_conduct) THEN
        CALL field_clamp_zero(bx, c_stagger_bx, i)
        CALL field_zero_gradient(by, c_stagger_by, i)
        CALL field_zero_gradient(bz, c_stagger_bz, i)
      ENDIF
    ENDDO

    DO i = c_bd_y_min, c_bd_y_max, c_bd_y_max - c_bd_y_min
      IF (bc_field(i) .EQ. c_bc_conduct) THEN
        CALL field_clamp_zero(by, c_stagger_by, i)
        CALL field_zero_gradient(bx, c_stagger_bx, i)
        CALL field_zero_gradient(bz, c_stagger_bz, i)
      ENDIF
    ENDDO

    DO i = 1, 2*c_ndims
      ! These apply zero field boundary conditions on the edges
      IF (bc_field(i) .EQ. c_bc_clamp &
          .OR. bc_field(i) .EQ. c_bc_simple_laser &
          .OR. bc_field(i) .EQ. c_bc_simple_outflow) THEN
        CALL field_clamp_zero(bx, c_stagger_bx, i)
        CALL field_clamp_zero(by, c_stagger_by, i)
        CALL field_clamp_zero(bz, c_stagger_bz, i)
      ENDIF

      ! These apply zero field boundary conditions on the edges
      IF (bc_field(i) .EQ. c_bc_zero_gradient) THEN
        CALL field_zero_gradient(bx, c_stagger_bx, i)
        CALL field_zero_gradient(by, c_stagger_by, i)
        CALL field_zero_gradient(bz, c_stagger_bz, i)
      ENDIF
    ENDDO

  END SUBROUTINE bfield_bcs



  SUBROUTINE bfield_final_bcs

    INTEGER :: i

    CALL bfield_bcs(.FALSE.)

    IF (x_min_boundary) THEN
      i = c_bd_x_min
      IF (add_laser(i) .OR. bc_field(i) .EQ. c_bc_simple_outflow) &
          CALL outflow_bcs_x_min
    ENDIF

    IF (x_max_boundary) THEN
      i = c_bd_x_max
      IF (add_laser(i) .OR. bc_field(i) .EQ. c_bc_simple_outflow) &
          CALL outflow_bcs_x_max
    ENDIF

    IF (y_min_boundary) THEN
      i = c_bd_y_min
      IF (add_laser(i) .OR. bc_field(i) .EQ. c_bc_simple_outflow) &
          CALL outflow_bcs_y_min
    ENDIF

    IF (y_max_boundary) THEN
      i = c_bd_y_max
      IF (add_laser(i) .OR. bc_field(i) .EQ. c_bc_simple_outflow) &
          CALL outflow_bcs_y_max
    ENDIF

    CALL bfield_bcs(.TRUE.)

  END SUBROUTINE bfield_final_bcs



  SUBROUTINE particle_bcs

    TYPE(particle), POINTER :: cur, next
    TYPE(particle_list), DIMENSION(-1:1,-1:1) :: send, recv
    INTEGER :: xbd, ybd
    INTEGER(KIND=8) :: ixp, iyp
    LOGICAL :: out_of_bounds
    INTEGER :: ispecies, i, ix, iy
    INTEGER :: cell_x, cell_y
    REAL(num), DIMENSION(-1:1) :: gx, gy
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    REAL(num) :: cf2, temp(3)
    REAL(num) :: part_pos

    DO ispecies = 1, n_species
      cur=>species_list(ispecies)%attached_list%head

      DO iy = -1, 1
        DO ix = -1, 1
          IF (ABS(ix) + ABS(iy) .EQ. 0) CYCLE
          CALL create_empty_partlist(send(ix, iy))
          CALL create_empty_partlist(recv(ix, iy))
        ENDDO
      ENDDO

      DO WHILE (ASSOCIATED(cur))
        next=>cur%next

        xbd = 0
        ybd = 0
        out_of_bounds = .FALSE.

        part_pos = cur%part_pos(1)
        ! Particle has left this processor
        IF (part_pos .LT. x_min_local - dx / 2.0_num) THEN
          xbd = -1
          ! Particle has left the system
          IF (x_min_boundary) THEN
            xbd = 0
            IF (bc_particle(c_bd_x_min) .EQ. c_bc_reflect) THEN
              cur%part_pos(1) = 2.0_num * x_min - dx - part_pos
              cur%part_p(1) = -cur%part_p(1)
            ELSE IF (bc_particle(c_bd_x_min) .EQ. c_bc_thermal) THEN
              ! Always use the triangle particle weighting for simplicity
              cell_y_r = (cur%part_pos(2) - y_min_local) / dy
              cell_y = FLOOR(cell_y_r + 0.5_num)
              cell_frac_y = REAL(cell_y, num) - cell_y_r
              cell_y = cell_y + 1

              cf2 = cell_frac_y**2
              gy(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_y)
              gy( 0) = 0.75_num - cf2
              gy( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_y)

              DO i = 1, 3
                temp(i) = 0.0_num
                DO iy = -1, 1
                  temp(i) = temp(i) + gy(iy) &
                      * species_list(ispecies)%ext_temp_x_min(cell_y+iy, i)
                ENDDO
              ENDDO

              ! x-direction
              i = 1
              cur%part_p(i) = ABS(momentum_from_temperature(&
                  species_list(ispecies)%mass, temp(i), 0.0_num))

              ! y-direction
              i = 2
              cur%part_p(i) = momentum_from_temperature(&
                  species_list(ispecies)%mass, temp(i), 0.0_num)

              ! z-direction
              i = 3
              cur%part_p(i) = momentum_from_temperature(&
                  species_list(ispecies)%mass, temp(i), 0.0_num)

              cur%part_pos(1) = 2.0_num * x_min - dx - part_pos

            ELSE IF (bc_particle(c_bd_x_min) .EQ. c_bc_periodic) THEN
              xbd = -1
              cur%part_pos(1) = part_pos + length_x
            ELSE
              ! Default to open boundary conditions - remove particle
              out_of_bounds = .TRUE.
            ENDIF
          ENDIF
        ENDIF

        ! Particle has left this processor
        IF (part_pos .GE. x_max_local + dx / 2.0_num) THEN
          xbd = 1
          ! Particle has left the system
          IF (x_max_boundary) THEN
            xbd = 0
            IF (bc_particle(c_bd_x_max) .EQ. c_bc_reflect) THEN
              cur%part_pos(1) = 2.0_num * x_max + dx - part_pos
              cur%part_p(1) = -cur%part_p(1)
            ELSE IF (bc_particle(c_bd_x_max) .EQ. c_bc_thermal) THEN
              ! Always use the triangle particle weighting for simplicity
              cell_y_r = (cur%part_pos(2) - y_min_local) / dy
              cell_y = FLOOR(cell_y_r + 0.5_num)
              cell_frac_y = REAL(cell_y, num) - cell_y_r
              cell_y = cell_y + 1

              cf2 = cell_frac_y**2
              gy(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_y)
              gy( 0) = 0.75_num - cf2
              gy( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_y)

              DO i = 1, 3
                temp(i) = 0.0_num
                DO iy = -1, 1
                  temp(i) = temp(i) + gy(iy) &
                      * species_list(ispecies)%ext_temp_x_max(cell_y+iy, i)
                ENDDO
              ENDDO

              ! x-direction
              i = 1
              cur%part_p(i) = -ABS(momentum_from_temperature(&
                  species_list(ispecies)%mass, temp(i), 0.0_num))

              ! y-direction
              i = 2
              cur%part_p(i) = momentum_from_temperature(&
                  species_list(ispecies)%mass, temp(i), 0.0_num)

              ! z-direction
              i = 3
              cur%part_p(i) = momentum_from_temperature(&
                  species_list(ispecies)%mass, temp(i), 0.0_num)

              cur%part_pos(1) = 2.0_num * x_max + dx - part_pos

            ELSE IF (bc_particle(c_bd_x_max) .EQ. c_bc_periodic) THEN
              xbd = 1
              cur%part_pos(1) = part_pos - length_x
            ELSE
              ! Default to open boundary conditions - remove particle
              out_of_bounds = .TRUE.
            ENDIF
          ENDIF
        ENDIF

        part_pos = cur%part_pos(2)
        ! Particle has left this processor
        IF (part_pos .LT. y_min_local - dy / 2.0_num) THEN
          ybd = -1
          ! Particle has left the system
          IF (y_min_boundary) THEN
            ybd = 0
            IF (bc_particle(c_bd_y_min) .EQ. c_bc_reflect) THEN
              cur%part_pos(2) = 2.0_num * y_min - dy - part_pos
              cur%part_p(2) = -cur%part_p(2)
            ELSE IF (bc_particle(c_bd_y_min) .EQ. c_bc_thermal) THEN
              ! Always use the triangle particle weighting for simplicity
              cell_x_r = (cur%part_pos(1) - x_min_local) / dx
              cell_x = FLOOR(cell_x_r + 0.5_num)
              cell_frac_x = REAL(cell_x, num) - cell_x_r
              cell_x = cell_x + 1

              cf2 = cell_frac_x**2
              gx(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_x)
              gx( 0) = 0.75_num - cf2
              gx( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_x)

              DO i = 1, 3
                temp(i) = 0.0_num
                DO ix = -1, 1
                  temp(i) = temp(i) + gx(ix) &
                      * species_list(ispecies)%ext_temp_y_min(cell_x+ix, i)
                ENDDO
              ENDDO

              ! x-direction
              i = 1
              cur%part_p(i) = momentum_from_temperature(&
                  species_list(ispecies)%mass, temp(i), 0.0_num)

              ! y-direction
              i = 2
              cur%part_p(i) = ABS(momentum_from_temperature(&
                  species_list(ispecies)%mass, temp(i), 0.0_num))

              ! z-direction
              i = 3
              cur%part_p(i) = momentum_from_temperature(&
                  species_list(ispecies)%mass, temp(i), 0.0_num)

              cur%part_pos(2) = 2.0_num * y_min - dy - part_pos

            ELSE IF (bc_particle(c_bd_y_min) .EQ. c_bc_periodic) THEN
              ybd = -1
              cur%part_pos(2) = part_pos + length_y
            ELSE
              ! Default to open boundary conditions - remove particle
              out_of_bounds = .TRUE.
            ENDIF
          ENDIF
        ENDIF

        ! Particle has left this processor
        IF (part_pos .GE. y_max_local + dy / 2.0_num) THEN
          ybd = 1
          ! Particle has left the system
          IF (y_max_boundary) THEN
            ybd = 0
            IF (bc_particle(c_bd_y_max) .EQ. c_bc_reflect) THEN
              cur%part_pos(2) = 2.0_num * y_max + dy - part_pos
              cur%part_p(2) = -cur%part_p(2)
            ELSE IF (bc_particle(c_bd_y_max) .EQ. c_bc_thermal) THEN
              ! Always use the triangle particle weighting for simplicity
              cell_x_r = (cur%part_pos(1) - x_min_local) / dx
              cell_x = FLOOR(cell_x_r + 0.5_num)
              cell_frac_x = REAL(cell_x, num) - cell_x_r
              cell_x = cell_x + 1

              cf2 = cell_frac_x**2
              gx(-1) = 0.5_num * (0.25_num + cf2 + cell_frac_x)
              gx( 0) = 0.75_num - cf2
              gx( 1) = 0.5_num * (0.25_num + cf2 - cell_frac_x)

              DO i = 1, 3
                temp(i) = 0.0_num
                DO ix = -1, 1
                  temp(i) = temp(i) + gx(ix) &
                      * species_list(ispecies)%ext_temp_y_max(cell_x+ix, i)
                ENDDO
              ENDDO

              ! x-direction
              i = 1
              cur%part_p(i) = momentum_from_temperature(&
                  species_list(ispecies)%mass, temp(i), 0.0_num)

              ! y-direction
              i = 2
              cur%part_p(i) = -ABS(momentum_from_temperature(&
                  species_list(ispecies)%mass, temp(i), 0.0_num))

              ! z-direction
              i = 3
              cur%part_p(i) = momentum_from_temperature(&
                  species_list(ispecies)%mass, temp(i), 0.0_num)

              cur%part_pos(2) = 2.0_num * y_max + dy - part_pos

            ELSE IF (bc_particle(c_bd_y_max) .EQ. c_bc_periodic) THEN
              ybd = 1
              cur%part_pos(2) = part_pos - length_y
            ELSE
              ! Default to open boundary conditions - remove particle
              out_of_bounds = .TRUE.
            ENDIF
          ENDIF
        ENDIF

        IF (out_of_bounds) THEN
          ! Particle has gone forever
          CALL remove_particle_from_partlist(&
              species_list(ispecies)%attached_list, cur)
          IF (dumpmask(c_dump_ejected_particles) .NE. c_io_never) THEN
            CALL add_particle_to_partlist(&
                ejected_list(ispecies)%attached_list, cur)
          ELSE
            DEALLOCATE(cur)
          ENDIF
        ELSE IF (ABS(xbd) + ABS(ybd) .GT. 0) THEN
          ! Particle has left processor, send it to its neighbour
          CALL remove_particle_from_partlist(&
              species_list(ispecies)%attached_list, cur)
          CALL add_particle_to_partlist(send(xbd, ybd), cur)
        ENDIF

        ! Move to next particle
        cur=>next
      ENDDO

      ! swap Particles
      DO iy = -1, 1
        DO ix = -1, 1
          IF (ABS(ix) + ABS(iy) .EQ. 0) CYCLE
          ixp = -ix
          iyp = -iy
          CALL partlist_sendrecv(send(ix, iy), recv(ixp, iyp), &
              neighbour(ix, iy), neighbour(ixp, iyp))
          CALL append_partlist(species_list(ispecies)%attached_list, &
              recv(ixp, iyp))
        ENDDO
      ENDDO

      DO iy = -1, 1
        DO ix = -1, 1
          IF (ABS(ix) + ABS(iy) .EQ. 0) CYCLE
          CALL destroy_partlist(send(ix, iy))
          CALL destroy_partlist(recv(ix, iy))
        ENDDO
      ENDDO

    ENDDO

  END SUBROUTINE particle_bcs



  SUBROUTINE current_bcs

    INTEGER :: i

    ! domain is decomposed. Just add currents at edges
    CALL processor_summation_bcs(jx, c_dir_x)
    CALL processor_summation_bcs(jy, c_dir_y)
    CALL processor_summation_bcs(jz, c_dir_z)

    DO i = 1, 2*c_ndims
      IF (bc_particle(i) .EQ. c_bc_reflect) THEN
        CALL field_clamp_zero(jx, c_stagger_jx, i)
        CALL field_clamp_zero(jy, c_stagger_jy, i)
        CALL field_clamp_zero(jz, c_stagger_jz, i)
      ENDIF
    ENDDO

  END SUBROUTINE current_bcs

END MODULE boundary
