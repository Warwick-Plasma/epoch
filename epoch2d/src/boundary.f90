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

    cpml_boundaries = .FALSE.
    DO i = 1, 2*c_ndims
      IF (bc_particle(i) .EQ. c_bc_other) bc_particle(i) = c_bc_reflect
      IF (bc_field(i) .EQ. c_bc_other) bc_field(i) = c_bc_clamp
      IF (bc_field(i) .EQ. c_bc_cpml_laser &
          .OR. bc_field(i) .EQ. c_bc_cpml_outflow) cpml_boundaries = .TRUE.
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
          .OR. bc_particle(i) .EQ. c_bc_simple_outflow &
          .OR. bc_particle(i) .EQ. c_bc_cpml_laser &
          .OR. bc_particle(i) .EQ. c_bc_cpml_outflow) &
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
        IF (bc_field(c_bd_x_min) .EQ. c_bc_cpml_laser &
            .OR. bc_field(c_bd_x_min) .EQ. c_bc_cpml_outflow) THEN
          IF (x_min_boundary) THEN
            ! Particle has left the system
            IF (part_pos .LT. x_min + dx * (cpml_thickness - 0.5_num)) THEN
              xbd = 0
              out_of_bounds = .TRUE.
            ENDIF
          ELSE
            ! Particle has left this processor
            IF (part_pos .LT. x_min_local - dx / 2.0_num) xbd = -1
          ENDIF
        ELSE
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
        ENDIF

        IF (bc_field(c_bd_x_max) .EQ. c_bc_cpml_laser &
            .OR. bc_field(c_bd_x_max) .EQ. c_bc_cpml_outflow) THEN
          IF (x_max_boundary) THEN
            ! Particle has left the system
            IF (part_pos .GE. x_max - dx * (cpml_thickness - 0.5_num)) THEN
              xbd = 0
              out_of_bounds = .TRUE.
            ENDIF
          ELSE
            ! Particle has left this processor
            IF (part_pos .GE. x_max_local + dx / 2.0_num) xbd =  1
          ENDIF
        ELSE
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
        ENDIF

        part_pos = cur%part_pos(2)
        IF (bc_field(c_bd_y_min) .EQ. c_bc_cpml_laser &
            .OR. bc_field(c_bd_y_min) .EQ. c_bc_cpml_outflow) THEN
          IF (y_min_boundary) THEN
            ! Particle has left the system
            IF (part_pos .LT. y_min + dy * (cpml_thickness - 0.5_num)) THEN
              ybd = 0
              out_of_bounds = .TRUE.
            ENDIF
          ELSE
            ! Particle has left this processor
            IF (part_pos .LT. y_min_local - dy / 2.0_num) ybd = -1
          ENDIF
        ELSE
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
        ENDIF

        IF (bc_field(c_bd_y_max) .EQ. c_bc_cpml_laser &
            .OR. bc_field(c_bd_y_max) .EQ. c_bc_cpml_outflow) THEN
          IF (y_max_boundary) THEN
            ! Particle has left the system
            IF (part_pos .GE. y_max - dy * (cpml_thickness - 0.5_num)) THEN
              ybd = 0
              out_of_bounds = .TRUE.
            ENDIF
          ELSE
            ! Particle has left this processor
            IF (part_pos .GE. y_max_local + dy / 2.0_num) ybd =  1
          ENDIF
        ELSE
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



  SUBROUTINE set_cpml_helpers(nx, nx_global_min, nx_global_max, &
      ny, ny_global_min, ny_global_max)

    INTEGER, INTENT(IN) :: nx, nx_global_min, nx_global_max
    INTEGER, INTENT(IN) :: ny, ny_global_min, ny_global_max
    INTEGER :: i
    INTEGER :: ix, ix_glob
    INTEGER :: iy, iy_glob
    INTEGER, PARAMETER :: cpml_m = 3
    INTEGER, PARAMETER :: cpml_ma = 1
    REAL(num) :: x_pos, x_pos_m, x_pos_ma
    REAL(num) :: y_pos, y_pos_m, y_pos_ma
    REAL(num) :: sigma, kappa, acoeff, bcoeff, ccoeff

    ALLOCATE(cpml_kappa_ex(-2:nx+3), cpml_kappa_bx(-2:nx+3))
    ALLOCATE(cpml_a_ex(-2:nx+3), cpml_a_bx(-2:nx+3))
    ALLOCATE(cpml_sigma_ex(-2:nx+3), cpml_sigma_bx(-2:nx+3))

    ALLOCATE(cpml_kappa_ey(-2:ny+3), cpml_kappa_by(-2:ny+3))
    ALLOCATE(cpml_a_ey(-2:ny+3), cpml_a_by(-2:ny+3))
    ALLOCATE(cpml_sigma_ey(-2:ny+3), cpml_sigma_by(-2:ny+3))

    cpml_kappa_ex = 1.0_num
    cpml_kappa_bx = 1.0_num

    cpml_a_ex = 0.0_num
    cpml_sigma_ex = 0.0_num
    cpml_a_bx = 0.0_num
    cpml_sigma_bx = 0.0_num

    cpml_kappa_ey = 1.0_num
    cpml_kappa_by = 1.0_num

    cpml_a_ey = 0.0_num
    cpml_sigma_ey = 0.0_num
    cpml_a_by = 0.0_num
    cpml_sigma_by = 0.0_num

    cpml_sigma_max = cpml_sigma_max * c * 0.8_num * (cpml_m + 1.0_num) / dx

    ! ============= x_min boundary =============

    i = c_bd_x_min
    IF (bc_field(i) .EQ. c_bc_cpml_laser &
        .OR. bc_field(i) .EQ. c_bc_cpml_outflow) THEN
      cpml_x_min_start = nx+1
      cpml_x_min_end = 0

      IF (nx_global_min .LE. cpml_thickness) THEN
        cpml_x_min = .TRUE.
        cpml_x_min_start = 1 ! in local grid coordinates

        ! The following distinction is necessary because, in principle, it is
        ! possible for the local domain to lie completely within the boundary
        ! layer.
        IF (nx_global_max .GE. cpml_thickness) THEN
          ! in local grid coordinates
          ! global -> local: ixl = ixg - nx_global_min + 1
          cpml_x_min_end = cpml_thickness - nx_global_min + 1
        ELSE
          cpml_x_min_end = nx ! in local grid coordinates
        ENDIF

        DO ix = cpml_x_min_start,cpml_x_min_end
          ! runs from 1 to cpml_thickness in global coordinates
          ! local -> global: ixg = ixl + nx_global_min - 1
          ix_glob = ix + nx_global_min - 1

          ! runs from 1.0 to nearly 0.0 (actually 0.0 at cpml_thickness+1)
          x_pos = 1.0_num - REAL(ix_glob-1,num) / REAL(cpml_thickness,num)
          x_pos_m = x_pos**cpml_m
          x_pos_ma = (1.0_num - x_pos)**cpml_ma

          cpml_kappa_ex(ix) = 1.0_num + (cpml_kappa_max - 1.0_num) * x_pos_m
          cpml_sigma_ex(ix) = cpml_sigma_max * x_pos_m
          cpml_a_ex(ix) = cpml_a_max * x_pos_ma

          ! runs from nearly 1.0 to nearly 0.0 on the half intervals
          ! 1.0 at ix_glob=1-1/2 and 0.0 at ix_glob=cpml_thickness+1/2
          x_pos = 1.0_num - (REAL(ix_glob,num) - 0.5_num) &
              / REAL(cpml_thickness,num)
          x_pos_m = x_pos**cpml_m
          x_pos_ma = (1.0_num - x_pos)**cpml_ma

          cpml_kappa_bx(ix) = 1.0_num + (cpml_kappa_max - 1.0_num) * x_pos_m
          cpml_sigma_bx(ix) = cpml_sigma_max * x_pos_m
          cpml_a_bx(ix) = cpml_a_max * x_pos_ma
        ENDDO
      ENDIF

      ! Ghost cells start at the edge of the CPML boundary
      IF (nx_global_min .LE. cpml_thickness + ng + 1 &
          .AND. nx_global_max .GE. cpml_thickness + ng + 1) THEN
        add_laser(i) = .TRUE.
        cpml_x_min_laser_idx = cpml_thickness + ng - nx_global_min + 1
      ENDIF
    ENDIF

    ! ============= x_max boundary =============

    i = c_bd_x_max
    ! Same as x_min using the transformation ix -> nx_global - ix + 1
    IF (bc_field(i) .EQ. c_bc_cpml_laser &
        .OR. bc_field(i) .EQ. c_bc_cpml_outflow) THEN
      cpml_x_max_start = nx+1
      cpml_x_max_end = 0

      IF (nx_global_max .GE. nx_global - cpml_thickness + 1) THEN
        cpml_x_max = .TRUE.
        cpml_x_max_end = nx ! in local grid coordinates

        ! The following distinction is necessary because, in principle, it is
        ! possible for the local domain to lie completely within the boundary
        ! layer.
        IF (nx_global_min .LE. nx_global - cpml_thickness + 1) THEN
          ! in local grid coordinates
          ! global -> local: ixl = ixg - nx_global_min + 1
          cpml_x_max_start = nx_global - cpml_thickness + 1 - nx_global_min + 1
        ELSE
          cpml_x_max_start = 1 ! in local grid coordinates
        ENDIF

        DO ix = cpml_x_max_start,cpml_x_max_end
          ! runs from cpml_thickness to 1 in global coordinates
          ! local -> global: ixg = ixl + nx_global_min - 1
          ix_glob = nx_global - (ix + nx_global_min - 1) + 1

          ! runs from nearly 0.0 (actually 0.0 at cpml_thickness+1) to 1.0
          x_pos = 1.0_num - REAL(ix_glob-1,num) / REAL(cpml_thickness,num)
          x_pos_m = x_pos**cpml_m
          x_pos_ma = (1.0_num - x_pos)**cpml_ma

          cpml_kappa_ex(ix) = 1.0_num + (cpml_kappa_max - 1.0_num) * x_pos_m
          cpml_sigma_ex(ix) = cpml_sigma_max * x_pos_m
          cpml_a_ex(ix) = cpml_a_max * x_pos_ma

          ! runs from nearly 0.0 to nearly 1.0 on the half intervals
          ! 0.0 at ix_glob=cpml_thickness+1/2 and 1.0 at ix_glob=1-1/2
          x_pos = 1.0_num - (REAL(ix_glob,num) - 0.5_num) &
              / REAL(cpml_thickness,num)
          x_pos_m = x_pos**cpml_m
          x_pos_ma = (1.0_num - x_pos)**cpml_ma

          cpml_kappa_bx(ix-1) = 1.0_num + (cpml_kappa_max - 1.0_num) * x_pos_m
          cpml_sigma_bx(ix-1) = cpml_sigma_max * x_pos_m
          cpml_a_bx(ix-1) = cpml_a_max * x_pos_ma
        ENDDO
      ENDIF

      ! Ghost cells start at the edge of the CPML boundary
      IF (nx_global_min .LE. nx_global - cpml_thickness - ng + 1 &
          .AND. nx_global_max .GE. nx_global - cpml_thickness - ng + 1) THEN
        add_laser(i) = .TRUE.
        cpml_x_max_laser_idx = &
            nx_global - cpml_thickness - ng + 1 - nx_global_min + 1
      ENDIF
    ENDIF

    ! ============= y_min boundary =============

    i = c_bd_y_min
    IF (bc_field(i) .EQ. c_bc_cpml_laser &
        .OR. bc_field(i) .EQ. c_bc_cpml_outflow) THEN
      cpml_y_min_start = ny+1
      cpml_y_min_end = 0

      IF (ny_global_min .LE. cpml_thickness) THEN
        cpml_y_min = .TRUE.
        cpml_y_min_start = 1 ! in local grid coordinates

        ! The following distinction is necessary because, in principle, it is
        ! possible for the local domain to lie completely within the boundary
        ! layer.
        IF (ny_global_max .GE. cpml_thickness) THEN
          ! in local grid coordinates
          ! global -> local: iyl = iyg - ny_global_min + 1
          cpml_y_min_end = cpml_thickness - ny_global_min + 1
        ELSE
          cpml_y_min_end = ny ! in local grid coordinates
        ENDIF

        DO iy = cpml_y_min_start,cpml_y_min_end
          ! runs from 1 to cpml_thickness in global coordinates
          ! local -> global: iyg = iyl + ny_global_min - 1
          iy_glob = iy + ny_global_min - 1

          ! runs from 1.0 to nearly 0.0 (actually 0.0 at cpml_thickness+1)
          y_pos = 1.0_num - REAL(iy_glob-1,num) / REAL(cpml_thickness,num)
          y_pos_m = y_pos**cpml_m
          y_pos_ma = (1.0_num - y_pos)**cpml_ma

          cpml_kappa_ey(iy) = 1.0_num + (cpml_kappa_max - 1.0_num) * y_pos_m
          cpml_sigma_ey(iy) = cpml_sigma_max * y_pos_m
          cpml_a_ey(iy) = cpml_a_max * y_pos_ma

          ! runs from nearly 1.0 to nearly 0.0 on the half intervals
          ! 1.0 at iy_glob=1-1/2 and 0.0 at iy_glob=cpml_thickness+1/2
          y_pos = 1.0_num - (REAL(iy_glob,num) - 0.5_num) &
              / REAL(cpml_thickness,num)
          y_pos_m = y_pos**cpml_m
          y_pos_ma = (1.0_num - y_pos)**cpml_ma

          cpml_kappa_by(iy) = 1.0_num + (cpml_kappa_max - 1.0_num) * y_pos_m
          cpml_sigma_by(iy) = cpml_sigma_max * y_pos_m
          cpml_a_by(iy) = cpml_a_max * y_pos_ma
        ENDDO
      ENDIF

      ! Ghost cells start at the edge of the CPML boundary
      IF (ny_global_min .LE. cpml_thickness + ng + 1 &
          .AND. ny_global_max .GE. cpml_thickness + ng + 1) THEN
        add_laser(i) = .TRUE.
        cpml_y_min_laser_idx = cpml_thickness + ng - ny_global_min + 1
      ENDIF
    ENDIF

    ! ============= y_max boundary =============

    i = c_bd_y_max
    ! Same as y_min using the transformation iy -> ny_global - iy + 1
    IF (bc_field(i) .EQ. c_bc_cpml_laser &
        .OR. bc_field(i) .EQ. c_bc_cpml_outflow) THEN
      cpml_y_max_start = ny+1
      cpml_y_max_end = 0

      IF (ny_global_max .GE. ny_global - cpml_thickness + 1) THEN
        cpml_y_max = .TRUE.
        cpml_y_max_end = ny ! in local grid coordinates

        ! The following distinction is necessary because, in principle, it is
        ! possible for the local domain to lie completely within the boundary
        ! layer.
        IF (ny_global_min .LE. ny_global - cpml_thickness + 1) THEN
          ! in local grid coordinates
          ! global -> local: iyl = iyg - ny_global_min + 1
          cpml_y_max_start = ny_global - cpml_thickness + 1 - ny_global_min + 1
        ELSE
          cpml_y_max_start = 1 ! in local grid coordinates
        ENDIF

        DO iy = cpml_y_max_start,cpml_y_max_end
          ! runs from cpml_thickness to 1 in global coordinates
          ! local -> global: iyg = iyl + ny_global_min - 1
          iy_glob = ny_global - (iy + ny_global_min - 1) + 1

          ! runs from nearly 0.0 (actually 0.0 at cpml_thickness+1) to 1.0
          y_pos = 1.0_num - REAL(iy_glob-1,num) / REAL(cpml_thickness,num)
          y_pos_m = y_pos**cpml_m
          y_pos_ma = (1.0_num - y_pos)**cpml_ma

          cpml_kappa_ey(iy) = 1.0_num + (cpml_kappa_max - 1.0_num) * y_pos_m
          cpml_sigma_ey(iy) = cpml_sigma_max * y_pos_m
          cpml_a_ey(iy) = cpml_a_max * y_pos_ma

          ! runs from nearly 0.0 to nearly 1.0 on the half intervals
          ! 0.0 at iy_glob=cpml_thickness+1/2 and 1.0 at iy_glob=1-1/2
          y_pos = 1.0_num - (REAL(iy_glob,num) - 0.5_num) &
              / REAL(cpml_thickness,num)
          y_pos_m = y_pos**cpml_m
          y_pos_ma = (1.0_num - y_pos)**cpml_ma

          cpml_kappa_by(iy-1) = 1.0_num + (cpml_kappa_max - 1.0_num) * y_pos_m
          cpml_sigma_by(iy-1) = cpml_sigma_max * y_pos_m
          cpml_a_by(iy-1) = cpml_a_max * y_pos_ma
        ENDDO
      ENDIF

      ! Ghost cells start at the edge of the CPML boundary
      IF (ny_global_min .LE. ny_global - cpml_thickness - ng + 1 &
          .AND. ny_global_max .GE. ny_global - cpml_thickness - ng + 1) THEN
        add_laser(i) = .TRUE.
        cpml_y_max_laser_idx = &
            ny_global - cpml_thickness - ng + 1 - ny_global_min + 1
      ENDIF
    ENDIF

  END SUBROUTINE set_cpml_helpers



  SUBROUTINE allocate_cpml_fields

    ! I will ignore memory consumption issues and, for simplicity,
    ! allocate the boundary fields throughout the whole simulation box.
    ALLOCATE(cpml_psi_eyx(-2:nx+3,-2:ny+3), cpml_psi_ezx(-2:nx+3,-2:ny+3))
    ALLOCATE(cpml_psi_byx(-2:nx+3,-2:ny+3), cpml_psi_bzx(-2:nx+3,-2:ny+3))
    ALLOCATE(cpml_psi_exy(-2:nx+3,-2:ny+3), cpml_psi_ezy(-2:nx+3,-2:ny+3))
    ALLOCATE(cpml_psi_bxy(-2:nx+3,-2:ny+3), cpml_psi_bzy(-2:nx+3,-2:ny+3))

    cpml_psi_eyx = 0.0_num
    cpml_psi_ezx = 0.0_num
    cpml_psi_byx = 0.0_num
    cpml_psi_bzx = 0.0_num

    cpml_psi_exy = 0.0_num
    cpml_psi_ezy = 0.0_num
    cpml_psi_bxy = 0.0_num
    cpml_psi_bzy = 0.0_num

  END SUBROUTINE allocate_cpml_fields



  SUBROUTINE deallocate_cpml_helpers

    DEALLOCATE(cpml_kappa_ex, cpml_kappa_bx)
    DEALLOCATE(cpml_a_ex, cpml_a_bx)
    DEALLOCATE(cpml_sigma_ex, cpml_sigma_bx)

    DEALLOCATE(cpml_kappa_ey, cpml_kappa_by)
    DEALLOCATE(cpml_a_ey, cpml_a_by)
    DEALLOCATE(cpml_sigma_ey, cpml_sigma_by)

  END SUBROUTINE deallocate_cpml_helpers



  SUBROUTINE cpml_advance_e_currents(tstep)

    REAL(num), INTENT(IN) :: tstep
    INTEGER :: ipos, ix, iy
    REAL(num) :: acoeff, bcoeff, ccoeff_d
    REAL(num) :: kappa, sigma

    ! ============= x_min boundary =============

    IF (bc_field(c_bd_x_min) .EQ. c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_min) .EQ. c_bc_cpml_outflow) THEN
      DO iy = 1,ny
        DO ipos = cpml_x_min_start,cpml_x_min_end
          kappa = cpml_kappa_ex(ipos)
          sigma = cpml_sigma_ex(ipos)
          acoeff = cpml_a_ex(ipos)
          bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
          ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
              / (sigma + kappa * acoeff) / dx

          cpml_psi_eyx(ipos,iy) = bcoeff * cpml_psi_eyx(ipos,iy) &
              + ccoeff_d * (bz(ipos,iy) - bz(ipos-1,iy))
          cpml_psi_ezx(ipos,iy) = bcoeff * cpml_psi_ezx(ipos,iy) &
              + ccoeff_d * (by(ipos,iy) - by(ipos-1,iy))
        ENDDO
      ENDDO
    ENDIF

    ! ============= x_max boundary =============

    IF (bc_field(c_bd_x_max) .EQ. c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_max) .EQ. c_bc_cpml_outflow) THEN
      DO iy = 1,ny
        DO ipos = cpml_x_max_start,cpml_x_max_end
          kappa = cpml_kappa_ex(ipos)
          sigma = cpml_sigma_ex(ipos)
          acoeff = cpml_a_ex(ipos)
          bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
          ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
              / (sigma + kappa * acoeff) / dx

          cpml_psi_eyx(ipos,iy) = bcoeff * cpml_psi_eyx(ipos,iy) &
              + ccoeff_d * (bz(ipos,iy) - bz(ipos-1,iy))
          cpml_psi_ezx(ipos,iy) = bcoeff * cpml_psi_ezx(ipos,iy) &
              + ccoeff_d * (by(ipos,iy) - by(ipos-1,iy))
        ENDDO
      ENDDO
    ENDIF

    ! ============= y_min boundary =============

    IF (bc_field(c_bd_y_min) .EQ. c_bc_cpml_laser &
        .OR. bc_field(c_bd_y_min) .EQ. c_bc_cpml_outflow) THEN
      DO ipos = cpml_y_min_start,cpml_y_min_end
        kappa = cpml_kappa_ey(ipos)
        sigma = cpml_sigma_ey(ipos)
        acoeff = cpml_a_ey(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dy

        DO ix = 1,nx
          cpml_psi_exy(ix,ipos) = bcoeff * cpml_psi_exy(ix,ipos) &
              + ccoeff_d * (bz(ix,ipos) - bz(ix,ipos-1))
          cpml_psi_ezy(ix,ipos) = bcoeff * cpml_psi_ezy(ix,ipos) &
              + ccoeff_d * (bx(ix,ipos) - bx(ix,ipos-1))
        ENDDO
      ENDDO
    ENDIF

    ! ============= y_max boundary =============

    IF (bc_field(c_bd_y_max) .EQ. c_bc_cpml_laser &
        .OR. bc_field(c_bd_y_max) .EQ. c_bc_cpml_outflow) THEN
      DO ipos = cpml_y_max_start,cpml_y_max_end
        kappa = cpml_kappa_ey(ipos)
        sigma = cpml_sigma_ey(ipos)
        acoeff = cpml_a_ey(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dy

        DO ix = 1,nx
          cpml_psi_exy(ix,ipos) = bcoeff * cpml_psi_exy(ix,ipos) &
              + ccoeff_d * (bz(ix,ipos) - bz(ix,ipos-1))
          cpml_psi_ezy(ix,ipos) = bcoeff * cpml_psi_ezy(ix,ipos) &
              + ccoeff_d * (bx(ix,ipos) - bx(ix,ipos-1))
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE cpml_advance_e_currents



  SUBROUTINE cpml_advance_b_currents(tstep)

    REAL(num), INTENT(IN) :: tstep
    INTEGER :: ipos, ix, iy
    REAL(num) :: acoeff, bcoeff, ccoeff_d
    REAL(num) :: kappa, sigma

    ! ============= x_min boundary =============

    IF (bc_field(c_bd_x_min) .EQ. c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_min) .EQ. c_bc_cpml_outflow) THEN
      DO iy = 1,ny
        DO ipos = cpml_x_min_start,cpml_x_min_end
          kappa = cpml_kappa_bx(ipos)
          sigma = cpml_sigma_bx(ipos)
          acoeff = cpml_a_bx(ipos)
          bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
          ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
              / (sigma + kappa * acoeff) / dx

          cpml_psi_byx(ipos,iy) = bcoeff * cpml_psi_byx(ipos,iy) &
              + ccoeff_d * (ez(ipos+1,iy) - ez(ipos,iy))
          cpml_psi_bzx(ipos,iy) = bcoeff * cpml_psi_bzx(ipos,iy) &
              + ccoeff_d * (ey(ipos+1,iy) - ey(ipos,iy))
        ENDDO
      ENDDO
    ENDIF

    ! ============= x_max boundary =============

    IF (bc_field(c_bd_x_max) .EQ. c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_max) .EQ. c_bc_cpml_outflow) THEN
      DO iy = 1,ny
        DO ipos = cpml_x_max_start-1,cpml_x_max_end-1
          kappa = cpml_kappa_bx(ipos)
          sigma = cpml_sigma_bx(ipos)
          acoeff = cpml_a_bx(ipos)
          bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
          ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
              / (sigma + kappa * acoeff) / dx

          cpml_psi_byx(ipos,iy) = bcoeff * cpml_psi_byx(ipos,iy) &
              + ccoeff_d * (ez(ipos+1,iy) - ez(ipos,iy))
          cpml_psi_bzx(ipos,iy) = bcoeff * cpml_psi_bzx(ipos,iy) &
              + ccoeff_d * (ey(ipos+1,iy) - ey(ipos,iy))
        ENDDO
      ENDDO
    ENDIF

    ! ============= y_min boundary =============

    IF (bc_field(c_bd_y_min) .EQ. c_bc_cpml_laser &
        .OR. bc_field(c_bd_y_min) .EQ. c_bc_cpml_outflow) THEN
      DO ipos = cpml_y_min_start,cpml_y_min_end
        kappa = cpml_kappa_by(ipos)
        sigma = cpml_sigma_by(ipos)
        acoeff = cpml_a_by(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dy

        DO ix = 1,nx
          cpml_psi_bxy(ix,ipos) = bcoeff * cpml_psi_bxy(ix,ipos) &
              + ccoeff_d * (ez(ix,ipos+1) - ez(ix,ipos))
          cpml_psi_bzy(ix,ipos) = bcoeff * cpml_psi_bzy(ix,ipos) &
              + ccoeff_d * (ex(ix,ipos+1) - ex(ix,ipos))
        ENDDO
      ENDDO
    ENDIF

    ! ============= y_max boundary =============

    IF (bc_field(c_bd_y_max) .EQ. c_bc_cpml_laser &
        .OR. bc_field(c_bd_y_max) .EQ. c_bc_cpml_outflow) THEN
      DO ipos = cpml_y_max_start-1,cpml_y_max_end-1
        kappa = cpml_kappa_by(ipos)
        sigma = cpml_sigma_by(ipos)
        acoeff = cpml_a_by(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dy

        DO ix = 1,nx
          cpml_psi_bxy(ix,ipos) = bcoeff * cpml_psi_bxy(ix,ipos) &
              + ccoeff_d * (ez(ix,ipos+1) - ez(ix,ipos))
          cpml_psi_bzy(ix,ipos) = bcoeff * cpml_psi_bzy(ix,ipos) &
              + ccoeff_d * (ex(ix,ipos+1) - ex(ix,ipos))
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE cpml_advance_b_currents

END MODULE boundary
