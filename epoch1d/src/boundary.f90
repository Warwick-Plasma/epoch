MODULE boundary

  USE partlist

  IMPLICIT NONE

CONTAINS

  SUBROUTINE setup_particle_boundaries

    any_open = .FALSE.

    ! For some types of boundary, fields and particles are treated in
    ! different ways, deal with that here
    IF (bc_x_min .EQ. c_bc_periodic) THEN
      bc_x_min_particle = c_bc_periodic
      bc_x_min_field = c_bc_periodic
    ENDIF
    IF (bc_x_max .EQ. c_bc_periodic) THEN
      bc_x_max_particle = c_bc_periodic
      bc_x_max_field = c_bc_periodic
    ENDIF

    ! For some types of boundary, fields and particles are treated in
    ! different ways, deal with that here
    IF (bc_x_min .EQ. c_bc_other) THEN
      bc_x_min_particle = c_bc_reflect
      bc_x_min_field = c_bc_clamp
    ENDIF
    IF (bc_x_max .EQ. c_bc_other) THEN
      bc_x_max_particle = c_bc_reflect
      bc_x_max_field = c_bc_clamp
    ENDIF

    ! laser boundaries reflect particles off a hard wall
    IF (bc_x_min .EQ. c_bc_simple_laser &
        .OR. bc_x_min .EQ. c_bc_simple_outflow) THEN
      bc_x_min_particle = c_bc_open
      bc_x_min_field = c_bc_clamp
      any_open = .TRUE.
    ENDIF
    IF (bc_x_max .EQ. c_bc_simple_laser &
        .OR. bc_x_max .EQ. c_bc_simple_outflow) THEN
      bc_x_max_particle = c_bc_open
      bc_x_max_field = c_bc_clamp
      any_open = .TRUE.
    ENDIF

    IF (any_open) CALL create_empty_partlist(ejected_particles)

  END SUBROUTINE setup_particle_boundaries



  ! Exchanges field values at processor boundaries and applies field
  ! boundary conditions
  SUBROUTINE field_bc(field)

    REAL(num), DIMENSION(-2:), INTENT(INOUT) :: field

    CALL do_field_mpi_with_lengths(field, nx)

  END SUBROUTINE field_bc



  SUBROUTINE do_field_mpi_with_lengths(field, nx_local)

    REAL(num), DIMENSION(-2:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: nx_local

    CALL MPI_SENDRECV(field(1:3), &
        3, mpireal, proc_x_min, tag, field(nx_local+1:nx_local+3), &
        3, mpireal, proc_x_max, tag, comm, status, errcode)
    CALL MPI_SENDRECV(field(nx_local-2:nx_local), &
        3, mpireal, proc_x_max, tag, field(-2:0), &
        3, mpireal, proc_x_min, tag, comm, status, errcode)

  END SUBROUTINE do_field_mpi_with_lengths



  SUBROUTINE field_zero_gradient(field, force)

    REAL(num), DIMENSION(-2:), INTENT(INOUT) :: field
    LOGICAL, INTENT(IN) :: force

    IF ((bc_x_min_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. proc_x_min .EQ. MPI_PROC_NULL) THEN
      field(-1) = field(2)
      field( 0) = field(1)
    ENDIF

    IF ((bc_x_max_field .EQ. c_bc_zero_gradient .OR. force) &
        .AND. proc_x_max .EQ. MPI_PROC_NULL) THEN
      field(nx+1) = field(nx)
      field(nx+2) = field(nx-1)
    ENDIF

  END SUBROUTINE field_zero_gradient



  SUBROUTINE field_clamp_zero(field, stagger)

    REAL(num), DIMENSION(-2:), INTENT(INOUT) :: field
    INTEGER, DIMENSION(1), INTENT(IN) :: stagger

    IF (bc_x_min_field .EQ. c_bc_clamp .AND. proc_x_min .EQ. MPI_PROC_NULL) THEN
      IF (stagger(1) .EQ. 1) THEN
        field(-1) = -field(1)
        field( 0) = 0.0_num
      ELSE
        field(-1) = -field(2)
        field( 0) = -field(1)
      ENDIF
    ENDIF

    IF (bc_x_max_field .EQ. c_bc_clamp .AND. proc_x_max .EQ. MPI_PROC_NULL) THEN
      IF (stagger(1) .EQ. 1) THEN
        field(nx  ) = 0.0_num
        field(nx+1) = -field(nx-1)
      ELSE
        field(nx+1) = -field(nx)
        field(nx+2) = -field(nx-1)
      ENDIF
    ENDIF

  END SUBROUTINE field_clamp_zero



  SUBROUTINE processor_summation_bcs(array)

    REAL(num), DIMENSION(-2:), INTENT(INOUT) :: array
    REAL(num), DIMENSION(:), ALLOCATABLE :: temp

    INTEGER, DIMENSION(-1:1) :: sizes, x_min, x_max, x_shift
    INTEGER :: xs, xe, xf

    sizes = 0
    x_min = 0
    x_max = 0
    x_shift = 0

    DO ix = -1, 1
      sizes(ix) = 1
      IF (ix .EQ. 0) THEN
        sizes(ix) = sizes(ix) * (nx+6)
        x_min(ix) = -2
        x_max(ix) = nx+3
      ELSE IF (ix .EQ. 1) THEN
        sizes(ix) = sizes(ix) * 3
        x_min(ix) = nx+1
        x_max(ix) = nx+3
        x_shift(ix) = -nx
      ELSE
        sizes(ix) = sizes(ix) * 3
        x_min(ix) = -2
        x_max(ix) = 0
        x_shift(ix) = nx
      ENDIF
    ENDDO

    DO ix = -1, 1, 2
      xs = x_min(ix)
      xe = x_max(ix)
      xf = x_shift(ix)

      ALLOCATE(temp(xs:xe))
      temp = 0.0_num
      CALL MPI_SENDRECV(array(xs:xe), sizes(ix), mpireal, neighbour(ix), tag, &
          temp, sizes(-ix), mpireal, neighbour(-ix), tag, comm, status, errcode)
      array(xs+xf:xe+xf) = array(xs+xf:xe+xf) + temp
      DEALLOCATE(temp)
    ENDDO

    CALL field_bc(array)

  END SUBROUTINE processor_summation_bcs



  SUBROUTINE efield_bcs

    ! These are the MPI boundaries
    CALL field_bc(ex)
    CALL field_bc(ey)
    CALL field_bc(ez)

    ! These apply zero field boundary conditions on the edges
    CALL field_clamp_zero(ex, (/1/))
    CALL field_clamp_zero(ey, (/0/))
    CALL field_clamp_zero(ez, (/0/))

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
      CALL field_clamp_zero(bx, (/0/))
      CALL field_clamp_zero(by, (/1/))
      CALL field_clamp_zero(bz, (/1/))
      ! These apply zero field boundary conditions on the edges
      CALL field_zero_gradient(bx, .FALSE.)
      CALL field_zero_gradient(by, .FALSE.)
      CALL field_zero_gradient(bz, .FALSE.)
    ENDIF

  END SUBROUTINE bfield_bcs



  SUBROUTINE particle_bcs

    TYPE(particle), POINTER :: cur, next
    TYPE(particle_list), DIMENSION(-1:1) :: send, recv
    INTEGER :: xbd
    INTEGER(KIND=8) :: ixp
    LOGICAL :: out_of_bounds
    INTEGER :: ispecies

    DO ispecies = 1, n_species
      cur=>particle_species(ispecies)%attached_list%head

      DO ix = -1, 1
        CALL create_empty_partlist(send(ix))
        CALL create_empty_partlist(recv(ix))
      ENDDO

      DO WHILE (ASSOCIATED(cur))
        next=>cur%next

        out_of_bounds = .FALSE.

        xbd = 0

        ! These conditions apply if a particle has passed a physical boundary
        ! Not a processor boundary or a periodic boundary
        IF (cur%part_pos .LE. x_min-dx/2.0_num &
            .AND. proc_x_min .EQ. MPI_PROC_NULL &
            .AND. bc_x_min_particle .EQ. c_bc_reflect) THEN
          ! particle has crossed left boundary
          cur%part_pos =  2.0_num * (x_min-dx/2.0_num) - cur%part_pos
          cur%part_p(1) = - cur%part_p(1)
        ENDIF

        IF (cur%part_pos .GE. x_max+dx/2.0_num &
            .AND. proc_x_max .EQ. MPI_PROC_NULL &
            .AND. bc_x_max_particle .EQ. c_bc_reflect) THEN
          ! particle has crossed right boundary
          cur%part_pos =  2.0_num *(x_max+dx/2.0_num) - cur%part_pos
          cur%part_p(1) = - cur%part_p(1)
        ENDIF

        IF (cur%part_pos .LT. x_min_local - dx/2.0_num) xbd = -1
        IF (cur%part_pos .GT. x_max_local + dx/2.0_num) xbd =  1

        IF ((cur%part_pos .LT. x_min - dx/2.0_num) &
            .AND. (bc_x_min_particle .EQ. c_bc_open)) out_of_bounds = .TRUE.
        IF ((cur%part_pos .GT. x_max + dx/2.0_num) &
            .AND. (bc_x_max_particle .EQ. c_bc_open)) out_of_bounds = .TRUE.

        IF (ABS(xbd) .GT. 0) THEN
          ! particle has left box
          CALL remove_particle_from_partlist(&
              particle_species(ispecies)%attached_list, cur)
          IF (.NOT. out_of_bounds) THEN
            CALL add_particle_to_partlist(send(xbd), cur)
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
      DO ix = -1, 1
        IF (ABS(ix) .EQ. 0) CYCLE
        ixp = -ix
        CALL partlist_sendrecv(send(ix), recv(ixp), neighbour(ix), &
            neighbour(ixp))
        CALL append_partlist(particle_species(ispecies)%attached_list, &
            recv(ixp))
      ENDDO

      ! Particles should only lie outside boundaries if the periodic boundaries
      ! are turned on. This now moves them to within the boundaries
      cur=>particle_species(ispecies)%attached_list%head
      ct = 0
      DO WHILE(ASSOCIATED(cur))
        IF (cur%part_pos .LT. x_min-dx/2.0_num) THEN
          cur%part_pos = cur%part_pos + length_x + dx
        ELSE IF (cur%part_pos .GT. x_max+dx/2.0_num) THEN
          cur%part_pos = cur%part_pos - length_x - dx
        ENDIF
        cur=>cur%next
      ENDDO
    ENDDO

  END SUBROUTINE particle_bcs

END MODULE boundary
