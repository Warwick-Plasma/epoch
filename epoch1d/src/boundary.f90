MODULE boundary

  USE partlist
  USE particle_temperature
  USE deck_io_block

  IMPLICIT NONE

CONTAINS

  SUBROUTINE setup_particle_boundaries

    INTEGER :: i

    ! For some types of boundary, fields and particles are treated in
    ! different ways, deal with that here

    DO i = 1, 2*c_ndims
      IF (bc_particle(i) .EQ. c_bc_other) bc_particle(i) = c_bc_reflect
      IF (bc_field(i) .EQ. c_bc_other) bc_field(i) = c_bc_clamp
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

    CALL MPI_SENDRECV(field(1), 3, mpireal, proc_x_min, tag, &
        field(nx_local+1), 3, mpireal, proc_x_max, tag, &
        comm, status, errcode)
    CALL MPI_SENDRECV(field(nx_local-2), 3, mpireal, proc_x_max, tag, &
        field(-2), 3, mpireal, proc_x_min, tag, &
        comm, status, errcode)

  END SUBROUTINE do_field_mpi_with_lengths



  SUBROUTINE field_zero_gradient(field, stagger_type, boundary)

    REAL(num), DIMENSION(-2:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: stagger_type, boundary

    IF (bc_field(boundary) .EQ. c_bc_periodic) RETURN

    IF (boundary .EQ. c_bd_x_min .AND. x_min_boundary) THEN
      field(-1) = field(2)
      field( 0) = field(1)
    ELSE IF (boundary .EQ. c_bd_x_max .AND. x_max_boundary) THEN
      field(nx+1) = field(nx  )
      field(nx+2) = field(nx-1)
    ENDIF

  END SUBROUTINE field_zero_gradient



  SUBROUTINE field_clamp_zero(field, stagger_type, boundary)

    REAL(num), DIMENSION(-2:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: stagger_type, boundary

    IF (bc_field(boundary) .EQ. c_bc_periodic) RETURN

    IF (boundary .EQ. c_bd_x_min .AND. x_min_boundary) THEN
      IF (stagger(1,stagger_type) .EQ. 1) THEN
        field(-1) = -field(1)
        field( 0) = 0.0_num
      ELSE
        field(-1) = -field(2)
        field( 0) = -field(1)
      ENDIF
    ELSE IF (boundary .EQ. c_bd_x_max .AND. x_max_boundary) THEN
      IF (stagger(1,stagger_type) .EQ. 1) THEN
        field(nx  ) = 0.0_num
        field(nx+1) = -field(nx-1)
      ELSE
        field(nx+1) = -field(nx  )
        field(nx+2) = -field(nx-1)
      ENDIF
    ENDIF

  END SUBROUTINE field_clamp_zero



  SUBROUTINE processor_summation_bcs(array, flip_direction)

    REAL(num), DIMENSION(-2:), INTENT(INOUT) :: array
    INTEGER, INTENT(IN), OPTIONAL :: flip_direction
    REAL(num), DIMENSION(3) :: temp
    INTEGER :: sgn

    temp = 0.0_num
    CALL MPI_SENDRECV(array(nx+1), 3, mpireal, &
        neighbour( 1), tag, temp, 3, mpireal, &
        neighbour(-1), tag, comm, status, errcode)

    ! Deal with reflecting boundaries differently
    IF (bc_particle(c_bd_x_min) .EQ. c_bc_reflect .AND. x_min_boundary) THEN
      sgn = 1
      IF (PRESENT(flip_direction)) THEN
        ! Currents get reversed in the direction of the boundary
        IF (flip_direction .EQ. c_dir_x) sgn = -1
      ENDIF
      array(1) = array(1) + sgn * array( 0)
      array(2) = array(2) + sgn * array(-1)
      array(3) = array(3) + sgn * array(-2)
    ELSE
      array(1:3) = array(1:3) + temp
    ENDIF

    temp = 0.0_num
    CALL MPI_SENDRECV(array(-2), 3, mpireal, &
        neighbour(-1), tag, temp, 3, mpireal, &
        neighbour( 1), tag, comm, status, errcode)

    ! Deal with reflecting boundaries differently
    IF (bc_particle(c_bd_x_max) .EQ. c_bc_reflect .AND. x_max_boundary) THEN
      sgn = 1
      IF (PRESENT(flip_direction)) THEN
        ! Currents get reversed in the direction of the boundary
        IF (flip_direction .EQ. c_dir_x) sgn = -1
      ENDIF
      array(nx-2) = array(nx-2) + sgn * array(nx+3)
      array(nx-1) = array(nx-1) + sgn * array(nx+2)
      array(nx  ) = array(nx  ) + sgn * array(nx+1)
    ELSE
      array(nx-2:nx) = array(nx-2:nx) + temp
    ENDIF

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



  SUBROUTINE particle_bcs

    TYPE(particle), POINTER :: cur, next
    TYPE(particle_list), DIMENSION(-1:1) :: send, recv
    INTEGER :: xbd
    INTEGER(KIND=8) :: ixp
    LOGICAL :: out_of_bounds
    INTEGER :: ispecies, ix
    REAL(num) :: part_pos

    DO ispecies = 1, n_species
      cur=>species_list(ispecies)%attached_list%head

      DO ix = -1, 1, 2
        CALL create_empty_partlist(send(ix))
        CALL create_empty_partlist(recv(ix))
      ENDDO

      DO WHILE (ASSOCIATED(cur))
        next=>cur%next

        xbd = 0
        out_of_bounds = .FALSE.

        part_pos = cur%part_pos
        ! Particle has left this processor
        IF (part_pos .LT. x_min_local - dx / 2.0_num) THEN
          xbd = -1
          ! Particle has left the system
          IF (x_min_boundary) THEN
            xbd = 0
            IF (bc_particle(c_bd_x_min) .EQ. c_bc_open) THEN
              out_of_bounds = .TRUE.
            ELSE IF (bc_particle(c_bd_x_min) .EQ. c_bc_reflect) THEN
              cur%part_pos = 2.0_num * x_min - dx - part_pos
              cur%part_p(1) = -cur%part_p(1)
            ELSE IF (bc_particle(c_bd_x_min) .EQ. c_bc_periodic) THEN
              xbd = -1
              cur%part_pos = part_pos + length_x
            ELSE
              IF (rank .EQ. 0) THEN
                WRITE(*,*) '*** ERROR ***'
                WRITE(*,*) 'Unrecognised particle boundary condition.'
              ENDIF
              CALL MPI_ABORT(MPI_COMM_WORLD, errcode, errcode)
            ENDIF
          ENDIF
        ENDIF

        ! Particle has left this processor
        IF (part_pos .GE. x_max_local + dx / 2.0_num) THEN
          xbd = 1
          ! Particle has left the system
          IF (x_max_boundary) THEN
            xbd = 0
            IF (bc_particle(c_bd_x_max) .EQ. c_bc_open) THEN
              out_of_bounds = .TRUE.
            ELSE IF (bc_particle(c_bd_x_max) .EQ. c_bc_reflect) THEN
              cur%part_pos = 2.0_num * x_max + dx - part_pos
              cur%part_p(1) = -cur%part_p(1)
            ELSE IF (bc_particle(c_bd_x_max) .EQ. c_bc_periodic) THEN
              xbd = 1
              cur%part_pos = part_pos - length_x
            ELSE
              IF (rank .EQ. 0) THEN
                WRITE(*,*) '*** ERROR ***'
                WRITE(*,*) 'Unrecognised particle boundary condition.'
              ENDIF
              CALL MPI_ABORT(MPI_COMM_WORLD, errcode, errcode)
            ENDIF
          ENDIF
        ENDIF

        IF (out_of_bounds) THEN
          ! Particle has gone forever
          CALL remove_particle_from_partlist(&
              species_list(ispecies)%attached_list, cur)
          IF (dumpmask(c_dump_ejected_particles) .NE. c_io_never) THEN
            CALL add_particle_to_partlist(ejected_particles, cur)
          ELSE
            DEALLOCATE(cur)
          ENDIF
        ELSE IF (ABS(xbd) .GT. 0) THEN
          ! Particle has left processor, send it to its neighbour
          CALL remove_particle_from_partlist(&
              species_list(ispecies)%attached_list, cur)
          CALL add_particle_to_partlist(send(xbd), cur)
        ENDIF

        ! Move to next particle
        cur=>next
      ENDDO

      ! swap Particles
      DO ix = -1, 1, 2
        ixp = -ix
        CALL partlist_sendrecv(send(ix), recv(ixp), &
            neighbour(ix), neighbour(ixp))
        CALL append_partlist(species_list(ispecies)%attached_list, &
            recv(ixp))
      ENDDO

      DO ix = -1, 1, 2
        CALL destroy_partlist(send(ix))
        CALL destroy_partlist(recv(ix))
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
