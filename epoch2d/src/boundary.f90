MODULE boundary

  USE partlist

  IMPLICIT NONE

CONTAINS

  SUBROUTINE setup_particle_boundaries

    any_open = .FALSE.

    ! For some types of boundary, fields and particles are treated in
    ! different ways, deal with that here
    IF (xbc_right .EQ. c_bc_periodic) THEN
      xbc_right_particle = c_bc_periodic
      xbc_right_field = c_bc_periodic
    ENDIF
    IF (xbc_left .EQ. c_bc_periodic) THEN
      xbc_left_particle = c_bc_periodic
      xbc_left_field = c_bc_periodic
    ENDIF

    IF (ybc_up .EQ. c_bc_periodic) THEN
      ybc_up_particle = c_bc_periodic
      ybc_up_field = c_bc_periodic
    ENDIF
    IF (ybc_down .EQ. c_bc_periodic) THEN
      ybc_down_particle = c_bc_periodic
      ybc_down_field = c_bc_periodic
    ENDIF

    ! For some types of boundary, fields and particles are treated in
    ! different ways, deal with that here
    IF (xbc_right .EQ. c_bc_other) THEN
      xbc_right_particle = c_bc_reflect
      xbc_right_field = c_bc_clamp
    ENDIF
    IF (xbc_left .EQ. c_bc_other) THEN
      xbc_left_particle = c_bc_reflect
      xbc_left_field = c_bc_clamp
    ENDIF

    IF (ybc_up .EQ. c_bc_other) THEN
      ybc_up_particle = c_bc_reflect
      ybc_up_field = c_bc_clamp
    ENDIF
    IF (ybc_down .EQ. c_bc_other) THEN
      ybc_down_particle = c_bc_reflect
      ybc_down_field = c_bc_clamp
    ENDIF

    ! laser boundaries reflect particles off a hard wall
    IF (xbc_left .EQ. c_bc_simple_laser &
        .OR. xbc_left .EQ. c_bc_simple_outflow) THEN
      xbc_left_particle = c_bc_open
      xbc_left_field = c_bc_clamp
      any_open = .TRUE.
    ENDIF
    IF (xbc_right .EQ. c_bc_simple_laser &
        .OR. xbc_right .EQ. c_bc_simple_outflow) THEN
      xbc_right_particle = c_bc_open
      xbc_right_field = c_bc_clamp
      any_open = .TRUE.
    ENDIF

    IF (ybc_up .EQ. c_bc_simple_laser &
        .OR. ybc_up .EQ. c_bc_simple_outflow) THEN
      ybc_up_particle = c_bc_open
      ybc_up_field = c_bc_clamp
      any_open = .TRUE.
    ENDIF
    IF (ybc_down .EQ. c_bc_simple_laser &
        .OR. ybc_down .EQ. c_bc_simple_outflow) THEN
      ybc_down_particle = c_bc_open
      ybc_down_field = c_bc_clamp
      any_open = .TRUE.
    ENDIF

    IF (any_open) CALL create_empty_partlist(ejected_particles)

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

    CALL MPI_SENDRECV(field(1:3,:), 3*(ny_local+6), mpireal, left, tag, &
        field(nx_local+1:nx_local+3,:), 3*(ny_local+6), mpireal, right, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV(field(nx_local-2:nx_local,:), 3*(ny_local+6), mpireal, &
        right, tag, field(-2:0,:), 3*(ny_local+6), mpireal, left, tag, &
        comm, status, errcode)

    CALL MPI_SENDRECV(field(:,1:3), 3*(nx_local+6), mpireal, down, tag, &
        field(:,ny_local+1:ny_local+3), 3*(nx_local+6), mpireal, up, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV(field(:,ny_local-2:ny_local), 3*(nx_local+6), mpireal, &
        up, tag, field(:,-2:0), 3*(nx_local+6), mpireal, down, tag, &
        comm, status, errcode)

  END SUBROUTINE do_field_mpi_with_lengths



  SUBROUTINE field_zero_gradient(field, force)

    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: field
    LOGICAL, INTENT(IN) :: force

    IF ((xbc_left_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. left .EQ. MPI_PROC_NULL) THEN
      field(0,:) = field(1,:)
      field(-1,:) = field(2,:)
    ENDIF

    IF ((xbc_right_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. right .EQ. MPI_PROC_NULL) THEN
      field(nx+1,:) = field(nx,:)
      field(nx+2,:) = field(nx-1,:)
    ENDIF

    IF ((ybc_down_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. down .EQ. MPI_PROC_NULL) THEN
      field(:,0) = field(:,1)
      field(:,-1) = field(:,2)
    ENDIF

    IF ((ybc_up_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. up .EQ. MPI_PROC_NULL) THEN
      field(:,ny+1) = field(:,ny)
      field(:,ny+2) = field(:,ny-1)
    ENDIF

  END SUBROUTINE field_zero_gradient



  SUBROUTINE field_clamp_zero(field, stagger)

    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: field
    INTEGER, DIMENSION(2), INTENT(IN) :: stagger

    IF (xbc_left_field .EQ. c_bc_clamp .AND. left .EQ. MPI_PROC_NULL) THEN
      IF (stagger(1) .EQ. 1) THEN
        field(0,:) = 0.0_num
        field(-1,:) = -field(1,:)
      ELSE
        field(0,:) = -field(1,:)
        field(-1,:) = -field(2,:)
      ENDIF
    ENDIF

    IF (xbc_right_field .EQ. c_bc_clamp .AND. right .EQ. MPI_PROC_NULL) THEN
      IF (stagger(1) .EQ. 1) THEN
        field(nx,:) = 0.0_num
      ELSE
        field(nx+1,:) = -field(nx,:)
        field(nx+2,:) = -field(nx-1,:)
      ENDIF
    ENDIF

    IF (ybc_down_field .EQ. c_bc_clamp .AND. down .EQ. MPI_PROC_NULL) THEN
      IF (stagger(2) .EQ. 1) THEN
        field(:,-1:0) = 0.0_num
      ELSE
        field(:,0) = -field(:,1)
        field(:,-1) = -field(:,2)
      ENDIF
    ENDIF

    IF (ybc_up_field .EQ. c_bc_clamp .AND. up .EQ. MPI_PROC_NULL) THEN
      IF (stagger(2) .EQ. 1) THEN
        field(:,ny:ny+2) = 0.0_num
      ELSE
        field(:,ny+1) = -field(:,ny)
        field(:,ny+2) = -field(:,ny-1)
      ENDIF
    ENDIF

  END SUBROUTINE field_clamp_zero



  SUBROUTINE processor_summation_bcs(array)

    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temp

    INTEGER, DIMENSION(-1:1, -1:1) :: sizes, x_start, x_end, x_shift
    INTEGER, DIMENSION(-1:1, -1:1) :: y_start, y_end, y_shift
    INTEGER :: xs, xe, xf, ys, ye, yf

    sizes = 0
    x_start = 0
    y_start = 0
    x_end = 0
    y_end = 0
    x_shift = 0
    y_shift = 0

    DO iy = -1, 1
      DO ix = -1, 1
        sizes(ix, iy) = 1
        IF (ix .EQ. 0) THEN
          sizes(ix, iy) = sizes(ix, iy) * (nx+6)
          x_start(ix, iy) = -2
          x_end(ix, iy) = nx+3
        ELSE IF (ix .EQ. 1) THEN
          sizes(ix, iy) = sizes(ix, iy) *3
          x_start(ix, iy) = nx+1
          x_end(ix, iy) = nx+3
          x_shift(ix, iy) = -nx
        ELSE
          sizes(ix, iy) = sizes(ix, iy) *3
          x_start(ix, iy) = -2
          x_end(ix, iy) = 0
          x_shift(ix, iy) = nx
        ENDIF
        IF (iy .EQ. 0) THEN
          sizes(ix, iy) = sizes(ix, iy) * (ny+6)
          y_start(ix, iy) = -2
          y_end(ix, iy) = ny+3
        ELSE IF (iy .EQ. 1) THEN
          sizes(ix, iy) = sizes(ix, iy) *3
          y_start(ix, iy) = ny+1
          y_end(ix, iy) = ny+3
          y_shift(ix, iy) = -ny
        ELSE
          sizes(ix, iy) = sizes(ix, iy) *3
          y_start(ix, iy) = -2
          y_end(ix, iy) = 0
          y_shift(ix, iy) = ny
        ENDIF
      ENDDO
    ENDDO

    DO iy = -1, 1
      DO ix = -1, 1
        IF (ix .EQ. 0 .AND. iy .EQ. 0) CYCLE
        IF (ABS(ix)+ABS(iy) .NE. 1) CYCLE
        ! Copy the starts into variables with shorter names, or this is
        ! HORRIFIC to read
        xs = x_start(ix, iy)
        ys = y_start(ix, iy)
        xe = x_end(ix, iy)
        ye = y_end(ix, iy)
        xf = x_shift(ix, iy)
        yf = y_shift(ix, iy)
        ALLOCATE(temp(xs:xe, ys:ye))
        temp = 0.0_num
        ! IF (neighbour(ix, iy) .EQ. MPI_PROC_NULL) &
        !     PRINT *, "BAD NEIGHBOUR", ix, iy
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
          CALL create_empty_partlist(send(ix, iy))
          CALL create_empty_partlist(recv(ix, iy))
        ENDDO
      ENDDO

      DO WHILE (ASSOCIATED(cur))
        out_of_bounds = .FALSE.
        next=>cur%next

        xbd = 0
        ybd = 0

        ! These conditions apply if a particle has passed a physical boundary
        ! Not a processor boundary or a periodic boundary
        IF (cur%part_pos(1) .LE. x_start-dx/2.0_num &
            .AND. left .EQ. MPI_PROC_NULL &
            .AND. xbc_left_particle .EQ. c_bc_reflect) THEN
          ! particle has crossed left boundary
          cur%part_pos(1) =  2.0_num * (x_start-dx/2.0_num) - cur%part_pos(1)
          cur%part_p(1) = - cur%part_p(1)
        ENDIF

        IF (cur%part_pos(1) .GE. x_end+dx/2.0_num &
            .AND. right .EQ. MPI_PROC_NULL &
            .AND. xbc_right_particle .EQ. c_bc_reflect) THEN
          ! particle has crossed right boundary
          cur%part_pos(1) =  2.0_num *(x_end+dx/2.0_num) - cur%part_pos(1)
          cur%part_p(1) = - cur%part_p(1)
        ENDIF

        IF (cur%part_pos(2) .LE. y_start-dy/2.0_num &
            .AND. down .EQ. MPI_PROC_NULL &
            .AND. ybc_down_particle .EQ. c_bc_reflect) THEN
          ! particle has crossed bottom boundary
          cur%part_pos(2) =  2.0_num * (y_start-dy/2.0_num) - cur%part_pos(2)
          cur%part_p(2) = - cur%part_p(2)
        ENDIF

        IF (cur%part_pos(2) .GE. y_end+dy/2.0_num &
            .AND. up .EQ. MPI_PROC_NULL &
            .AND. ybc_up_particle .EQ. c_bc_reflect) THEN
          ! particle has crossed top boundary
          cur%part_pos(2) =  2.0_num * (y_end + dy/2.0_num) - cur%part_pos(2)
          cur%part_p(2) = - cur%part_p(2)
        ENDIF

        IF (cur%part_pos(1) .LT. x_start_local - dx/2.0_num) xbd = -1
        IF (cur%part_pos(1) .GT. x_end_local + dx/2.0_num )  xbd = 1
        IF (cur%part_pos(2) .LT. y_start_local - dy/2.0_num) ybd = -1
        IF (cur%part_pos(2) .GT. y_end_local + dy/2.0_num)   ybd = 1

        IF ((cur%part_pos(1) .LT. x_start - dx/2.0_num) &
            .AND. (xbc_left_particle .EQ. c_bc_open)) out_of_bounds = .TRUE.
        IF ((cur%part_pos(1) .GT. x_end + dx/2.0_num) &
            .AND. (xbc_right_particle .EQ. c_bc_open)) out_of_bounds = .TRUE.
        IF ((cur%part_pos(2) .LT. y_start - dy/2.0_num) &
            .AND. (ybc_down_particle .EQ. c_bc_open)) out_of_bounds = .TRUE.
        IF ((cur%part_pos(2) .GT. y_end +dy/2.0_num) &
            .AND. (ybc_up_particle .EQ. c_bc_open)) out_of_bounds = .TRUE.

        IF (ABS(xbd) + ABS(ybd) .GT. 0) THEN
          ! particle has left box
          CALL remove_particle_from_partlist(&
              particle_species(ispecies)%attached_list, cur)
          IF (.NOT. out_of_bounds) THEN
            CALL add_particle_to_partlist(send(xbd, ybd), cur)
          ELSE
            IF (dumpmask(29) .NE. c_io_never) THEN
              CALL add_particle_to_partlist(ejected_particles, cur)
            ELSE
              DEALLOCATE(cur)
            ENDIF
          ENDIF
        ENDIF

        ! Move to next particle
        cur=>next
      ENDDO

      ! swap Particles
      DO iy = -1, 1
        DO ix = -1, 1
          IF (ABS(ix)+ABS(iy) .EQ. 0) CYCLE
          ixp = -ix
          iyp = -iy
          CALL partlist_sendrecv(send(ix, iy), recv(ixp, iyp), &
              neighbour(ix, iy), neighbour(ixp, iyp))
          CALL append_partlist(particle_species(ispecies)%attached_list, &
              recv(ixp, iyp))
        ENDDO
      ENDDO

      ! Particles should only lie outside boundaries if the periodic boundaries
      ! are turned on. This now moves them to within the boundaries
      cur=>particle_species(ispecies)%attached_list%head
      ct = 0
      DO WHILE(ASSOCIATED(cur))
        IF (cur%part_pos(1) .GT. x_end+dx/2.0_num &
            .AND. xbc_left_particle .EQ. c_bc_periodic) &
                cur%part_pos(1) = cur%part_pos(1)-length_x - dx
        IF (cur%part_pos(1) .LT. x_start-dx/2.0_num &
            .AND. xbc_right_particle .EQ. c_bc_periodic) &
                cur%part_pos(1) = cur%part_pos(1)+length_x + dx
        IF (cur%part_pos(2) .GT. y_end+dy/2.0_num &
            .AND. ybc_up_particle .EQ. c_bc_periodic) &
                cur%part_pos(2) = cur%part_pos(2)-length_y - dy
        IF (cur%part_pos(2) .LT. y_start-dy/2.0_num &
            .AND. ybc_down_particle .EQ. c_bc_periodic) &
                cur%part_pos(2) = cur%part_pos(2)+length_y + dy
        cur=>cur%next
      ENDDO
    ENDDO

  END SUBROUTINE particle_bcs

END MODULE boundary
