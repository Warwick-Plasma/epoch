MODULE boundary

  USE partlist

  IMPLICIT NONE

CONTAINS

  SUBROUTINE setup_particle_boundaries

    ! For some types of boundary, fields and particles are treated in
    ! different ways, deal with that here

    IF (bc_x_min_particle .EQ. c_bc_other) bc_x_min_particle = c_bc_reflect
    IF (bc_x_max_particle .EQ. c_bc_other) bc_x_max_particle = c_bc_reflect
    IF (bc_y_min_particle .EQ. c_bc_other) bc_y_min_particle = c_bc_reflect
    IF (bc_y_max_particle .EQ. c_bc_other) bc_y_max_particle = c_bc_reflect

    IF (bc_x_min_field .EQ. c_bc_other) bc_x_min_field = c_bc_clamp
    IF (bc_x_max_field .EQ. c_bc_other) bc_x_max_field = c_bc_clamp
    IF (bc_y_min_field .EQ. c_bc_other) bc_y_min_field = c_bc_clamp
    IF (bc_y_max_field .EQ. c_bc_other) bc_y_max_field = c_bc_clamp

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

    ! Note: reflecting EM boundaries not yet implemented.
    IF (bc_x_min_field .EQ. c_bc_reflect) bc_x_min_field = c_bc_clamp
    IF (bc_x_max_field .EQ. c_bc_reflect) bc_x_max_field = c_bc_clamp
    IF (bc_y_min_field .EQ. c_bc_reflect) bc_y_min_field = c_bc_clamp
    IF (bc_y_max_field .EQ. c_bc_reflect) bc_y_max_field = c_bc_clamp

    IF (bc_x_min_field .EQ. c_bc_open) bc_x_min_field = c_bc_simple_outflow
    IF (bc_x_max_field .EQ. c_bc_open) bc_x_min_field = c_bc_simple_outflow
    IF (bc_y_min_field .EQ. c_bc_open) bc_y_min_field = c_bc_simple_outflow
    IF (bc_y_max_field .EQ. c_bc_open) bc_y_max_field = c_bc_simple_outflow

    IF (bc_x_min_particle .EQ. c_bc_open &
        .OR. bc_x_max_particle .EQ. c_bc_open &
        .OR. bc_y_min_particle .EQ. c_bc_open &
        .OR. bc_y_max_particle .EQ. c_bc_open) &
            CALL create_empty_partlist(ejected_particles)

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

    CALL MPI_SENDRECV(field(1:3,:), &
        3*(ny_local+6), mpireal, proc_x_min, tag, &
        field(nx_local+1:nx_local+3,:), &
        3*(ny_local+6), mpireal, proc_x_max, tag, comm, status, errcode)
    CALL MPI_SENDRECV(field(nx_local-2:nx_local,:), &
        3*(ny_local+6), mpireal, proc_x_max, tag, &
        field(-2:0,:), &
        3*(ny_local+6), mpireal, proc_x_min, tag, comm, status, errcode)

    CALL MPI_SENDRECV(field(:,1:3), &
        3*(nx_local+6), mpireal, proc_y_min, tag, &
        field(:,ny_local+1:ny_local+3), &
        3*(nx_local+6), mpireal, proc_y_max, tag, comm, status, errcode)
    CALL MPI_SENDRECV(field(:,ny_local-2:ny_local), &
        3*(nx_local+6), mpireal, proc_y_max, tag, &
        field(:,-2:0), &
        3*(nx_local+6), mpireal, proc_y_min, tag, comm, status, errcode)

  END SUBROUTINE do_field_mpi_with_lengths



  SUBROUTINE field_zero_gradient(field, force)

    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: field
    LOGICAL, INTENT(IN) :: force

    IF ((bc_x_min_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. proc_x_min .EQ. MPI_PROC_NULL) THEN
      field(-1,:) = field(2,:)
      field( 0,:) = field(1,:)
    ENDIF

    IF ((bc_x_max_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. proc_x_max .EQ. MPI_PROC_NULL) THEN
      field(nx+1,:) = field(nx,  :)
      field(nx+2,:) = field(nx-1,:)
    ENDIF

    IF ((bc_y_min_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. proc_y_min .EQ. MPI_PROC_NULL) THEN
      field(:,-1) = field(:,2)
      field(:, 0) = field(:,1)
    ENDIF

    IF ((bc_y_max_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. proc_y_max .EQ. MPI_PROC_NULL) THEN
      field(:,ny+1) = field(:,ny)
      field(:,ny+2) = field(:,ny-1)
    ENDIF

  END SUBROUTINE field_zero_gradient



  SUBROUTINE field_clamp_zero(field, stagger)

    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: field
    INTEGER, DIMENSION(2), INTENT(IN) :: stagger

    ! Use clamp when the laser is on.

    IF ((bc_x_min_field .EQ. c_bc_clamp &
        .OR. bc_x_min_field .EQ. c_bc_simple_laser &
        .OR. bc_x_min_field .EQ. c_bc_simple_outflow) &
        .AND. proc_x_min .EQ. MPI_PROC_NULL) THEN
      IF (stagger(1) .EQ. 1) THEN
        field(-1,:) = -field(1,:)
        field( 0,:) = 0.0_num
      ELSE
        field(-1,:) = -field(2,:)
        field( 0,:) = -field(1,:)
      ENDIF
    ENDIF

    IF ((bc_x_max_field .EQ. c_bc_clamp &
        .OR. bc_x_max_field .EQ. c_bc_simple_laser &
        .OR. bc_x_max_field .EQ. c_bc_simple_outflow) &
        .AND. proc_x_max .EQ. MPI_PROC_NULL) THEN
      IF (stagger(1) .EQ. 1) THEN
        field(nx,  :) = 0.0_num
        field(nx+1,:) = -field(nx-1,:)
      ELSE
        field(nx+1,:) = -field(nx,  :)
        field(nx+2,:) = -field(nx-1,:)
      ENDIF
    ENDIF

    IF ((bc_y_min_field .EQ. c_bc_clamp &
        .OR. bc_y_min_field .EQ. c_bc_simple_laser &
        .OR. bc_y_min_field .EQ. c_bc_simple_outflow) &
        .AND. proc_y_min .EQ. MPI_PROC_NULL) THEN
      IF (stagger(2) .EQ. 1) THEN
        field(:,-1) = -field(:,1)
        field(:, 0) = 0.0_num
      ELSE
        field(:,-1) = -field(:,2)
        field(:, 0) = -field(:,1)
      ENDIF
    ENDIF

    IF ((bc_y_max_field .EQ. c_bc_clamp &
        .OR. bc_y_max_field .EQ. c_bc_simple_laser &
        .OR. bc_y_max_field .EQ. c_bc_simple_outflow) &
        .AND. proc_y_max .EQ. MPI_PROC_NULL) THEN
      IF (stagger(2) .EQ. 1) THEN
        field(:,ny  ) = 0.0_num
        field(:,ny+1) = -field(:,ny-1)
      ELSE
        field(:,ny+1) = -field(:,ny)
        field(:,ny+2) = -field(:,ny-1)
      ENDIF
    ENDIF

  END SUBROUTINE field_clamp_zero



  SUBROUTINE processor_summation_bcs(array)

    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temp

    INTEGER, DIMENSION(-1:1, -1:1) :: sizes, x_min, x_max, x_shift
    INTEGER, DIMENSION(-1:1, -1:1) :: y_min, y_max, y_shift
    INTEGER :: xs, xe, xf, ys, ye, yf

    sizes = 0
    x_min = 0
    y_min = 0
    x_max = 0
    y_max = 0
    x_shift = 0
    y_shift = 0

    DO iy = -1, 1
      DO ix = -1, 1
        sizes(ix, iy) = 1
        IF (ix .EQ. 0) THEN
          sizes(ix, iy) = sizes(ix, iy) * (nx+6)
          x_min(ix, iy) = -2
          x_max(ix, iy) = nx+3
        ELSE IF (ix .EQ. 1) THEN
          sizes(ix, iy) = sizes(ix, iy) * 3
          x_min(ix, iy) = nx+1
          x_max(ix, iy) = nx+3
          x_shift(ix, iy) = -nx
        ELSE
          sizes(ix, iy) = sizes(ix, iy) * 3
          x_min(ix, iy) = -2
          x_max(ix, iy) = 0
          x_shift(ix, iy) = nx
        ENDIF
        IF (iy .EQ. 0) THEN
          sizes(ix, iy) = sizes(ix, iy) * (ny+6)
          y_min(ix, iy) = -2
          y_max(ix, iy) = ny+3
        ELSE IF (iy .EQ. 1) THEN
          sizes(ix, iy) = sizes(ix, iy) * 3
          y_min(ix, iy) = ny+1
          y_max(ix, iy) = ny+3
          y_shift(ix, iy) = -ny
        ELSE
          sizes(ix, iy) = sizes(ix, iy) * 3
          y_min(ix, iy) = -2
          y_max(ix, iy) = 0
          y_shift(ix, iy) = ny
        ENDIF
      ENDDO
    ENDDO

    DO iy = -1, 1
      DO ix = -1, 1
        IF (ix .EQ. 0 .AND. iy .EQ. 0) CYCLE
        IF (ABS(ix)+ABS(iy) .NE. 1) CYCLE

        xs = x_min(ix, iy)
        ys = y_min(ix, iy)
        xe = x_max(ix, iy)
        ye = y_max(ix, iy)
        xf = x_shift(ix, iy)
        yf = y_shift(ix, iy)

        ALLOCATE(temp(xs:xe, ys:ye))
        temp = 0.0_num
        CALL MPI_SENDRECV(array(xs:xe, ys:ye), sizes(ix, iy), mpireal, &
            neighbour(ix, iy), tag, temp, sizes(-ix, -iy), mpireal, &
            neighbour(-ix, -iy), tag, comm, status, errcode)
        array(xs+xf:xe+xf, ys+yf:ye+yf) = array(xs+xf:xe+xf, ys+yf:ye+yf) + temp
        DEALLOCATE(temp)
      ENDDO
    ENDDO

    CALL field_bc(array)

  END SUBROUTINE processor_summation_bcs



  SUBROUTINE efield_bcs

    ! These are the MPI boundaries
    CALL field_bc(ex)
    CALL field_bc(ey)
    CALL field_bc(ez)

    ! These apply zero field boundary conditions on the edges
    CALL field_clamp_zero(ex, (/1, 0/))
    CALL field_clamp_zero(ey, (/0, 1/))
    CALL field_clamp_zero(ez, (/0, 0/))

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
      CALL field_clamp_zero(bx, (/0, 1/))
      CALL field_clamp_zero(by, (/1, 0/))
      CALL field_clamp_zero(bz, (/1, 1/))
      ! These apply zero field boundary conditions on the edges
      CALL field_zero_gradient(bx, .FALSE.)
      CALL field_zero_gradient(by, .FALSE.)
      CALL field_zero_gradient(bz, .FALSE.)
    ENDIF

  END SUBROUTINE bfield_bcs



  SUBROUTINE particle_bcs

    TYPE(particle), POINTER :: cur, next
    TYPE(particle_list), DIMENSION(-1:1, -1:1) :: send, recv
    INTEGER :: xbd, ybd
    INTEGER(KIND=8) :: ixp, iyp
    LOGICAL :: out_of_bounds
    INTEGER :: ispecies

    DO ispecies = 1, n_species
      cur=>particle_species(ispecies)%attached_list%head

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

        ! Particle has left this processor
        IF (cur%part_pos(1) .LT. x_min_local - dx / 2.0_num) THEN
          xbd = -1
          ! Particle has left the system
          IF (cur%part_pos(1) .LT. x_min - dx / 2.0_num) THEN
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
          IF (cur%part_pos(1) .GE. x_max + dx / 2.0_num) THEN
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
          IF (cur%part_pos(2) .LT. y_min - dy / 2.0_num) THEN
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
          IF (cur%part_pos(2) .GE. y_max + dy / 2.0_num) THEN
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

        IF (out_of_bounds) THEN
          ! Particle has gone forever
          CALL remove_particle_from_partlist(&
              particle_species(ispecies)%attached_list, cur)
          IF (dumpmask(29) .NE. c_io_never) THEN
            CALL add_particle_to_partlist(ejected_particles, cur)
          ELSE
            DEALLOCATE(cur)
          ENDIF
        ELSE IF (ABS(xbd) + ABS(ybd) .GT. 0) THEN
          ! Particle has left processor, send it to its neighbour
          CALL remove_particle_from_partlist(&
              particle_species(ispecies)%attached_list, cur)
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
          CALL append_partlist(particle_species(ispecies)%attached_list, &
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

END MODULE boundary
