MODULE boundary

  USE mpi
  USE partlist
  USE particle_temperature
  USE deck_io_block
  USE laser

  IMPLICIT NONE

CONTAINS

  SUBROUTINE setup_particle_boundaries

    INTEGER :: i, ierr
    LOGICAL :: error
    CHARACTER(LEN=5), DIMENSION(2*c_ndims) :: &
        boundary = (/ 'x_min', 'x_max' /)

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

    IF (error) CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)

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
      IF (stagger(c_dir_x,stagger_type)) THEN
        field(-1) = field(1)
      ELSE
        field(-1) = field(2)
        field( 0) = field(1)
      ENDIF
    ELSE IF (boundary .EQ. c_bd_x_max .AND. x_max_boundary) THEN
      IF (stagger(c_dir_x,stagger_type)) THEN
        field(nx+1) = field(nx-1)
      ELSE
        field(nx+1) = field(nx  )
        field(nx+2) = field(nx-1)
      ENDIF
    ENDIF

  END SUBROUTINE field_zero_gradient



  SUBROUTINE field_clamp_zero(field, stagger_type, boundary)

    REAL(num), DIMENSION(-2:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: stagger_type, boundary

    IF (bc_field(boundary) .EQ. c_bc_periodic) RETURN

    IF (boundary .EQ. c_bd_x_min .AND. x_min_boundary) THEN
      IF (stagger(c_dir_x,stagger_type)) THEN
        field(-1) = -field(1)
        field( 0) = 0.0_num
      ELSE
        field(-1) = -field(2)
        field( 0) = -field(1)
      ENDIF
    ELSE IF (boundary .EQ. c_bd_x_max .AND. x_max_boundary) THEN
      IF (stagger(c_dir_x,stagger_type)) THEN
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
    IF ((bc_particle(c_bd_x_min) .EQ. c_bc_reflect &
        .OR. bc_particle(c_bd_x_min) .EQ. c_bc_thermal) &
        .AND. x_min_boundary) THEN
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
    IF ((bc_particle(c_bd_x_max) .EQ. c_bc_reflect &
        .OR. bc_particle(c_bd_x_max) .EQ. c_bc_thermal) &
        .AND. x_max_boundary) THEN
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
      ! These apply zero field boundary conditions on the edges
      IF (bc_field(i) .EQ. c_bc_clamp &
          .OR. bc_field(i) .EQ. c_bc_simple_laser &
          .OR. bc_field(i) .EQ. c_bc_simple_outflow) THEN
        CALL field_clamp_zero(ex, c_stagger_ex, i)
        CALL field_clamp_zero(ey, c_stagger_ey, i)
        CALL field_clamp_zero(ez, c_stagger_ez, i)
      ENDIF

      ! These apply zero gradient boundary conditions on the edges
      IF (bc_field(i) .EQ. c_bc_zero_gradient &
          .OR. bc_field(i) .EQ. c_bc_cpml_laser &
          .OR. bc_field(i) .EQ. c_bc_cpml_outflow) THEN
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

      ! These apply zero gradient boundary conditions on the edges
      IF (bc_field(i) .EQ. c_bc_zero_gradient &
          .OR. bc_field(i) .EQ. c_bc_cpml_laser &
          .OR. bc_field(i) .EQ. c_bc_cpml_outflow) THEN
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

    CALL bfield_bcs(.TRUE.)

  END SUBROUTINE bfield_final_bcs



  SUBROUTINE particle_bcs

    TYPE(particle), POINTER :: cur, next
    TYPE(particle_list), DIMENSION(-1:1) :: send, recv
    INTEGER :: xbd
    INTEGER(KIND=8) :: ixp
    LOGICAL :: out_of_bounds
    INTEGER :: ispecies, i, ix
    REAL(num) :: temp(3)
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
                cur%part_pos = 2.0_num * x_min - dx - part_pos
                cur%part_p(1) = -cur%part_p(1)
              ELSE IF (bc_particle(c_bd_x_min) .EQ. c_bc_thermal) THEN
                DO i = 1, 3
                  temp(i) = species_list(ispecies)%ext_temp_x_min(i)
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

                cur%part_pos = 2.0_num * x_min - dx - part_pos

              ELSE IF (bc_particle(c_bd_x_min) .EQ. c_bc_periodic) THEN
                xbd = -1
                cur%part_pos = part_pos + length_x
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
                cur%part_pos = 2.0_num * x_max + dx - part_pos
                cur%part_p(1) = -cur%part_p(1)
              ELSE IF (bc_particle(c_bd_x_max) .EQ. c_bc_thermal) THEN
                DO i = 1, 3
                  temp(i) = species_list(ispecies)%ext_temp_x_max(i)
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

                cur%part_pos = 2.0_num * x_max + dx - part_pos

              ELSE IF (bc_particle(c_bd_x_max) .EQ. c_bc_periodic) THEN
                xbd = 1
                cur%part_pos = part_pos - length_x
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



  SUBROUTINE set_cpml_helpers(nx, nx_global_min, nx_global_max)

    INTEGER, INTENT(IN) :: nx, nx_global_min, nx_global_max
    INTEGER :: i
    INTEGER :: ix, ix_glob
    INTEGER, PARAMETER :: cpml_m = 3
    INTEGER, PARAMETER :: cpml_ma = 1
    REAL(num) :: x_pos, x_pos_m, x_pos_ma

    ALLOCATE(cpml_kappa_ex(-2:nx+3), cpml_kappa_bx(-2:nx+3))
    ALLOCATE(cpml_a_ex(-2:nx+3), cpml_a_bx(-2:nx+3))
    ALLOCATE(cpml_sigma_ex(-2:nx+3), cpml_sigma_bx(-2:nx+3))

    cpml_kappa_ex = 1.0_num
    cpml_kappa_bx = 1.0_num

    cpml_a_ex = 0.0_num
    cpml_sigma_ex = 0.0_num
    cpml_a_bx = 0.0_num
    cpml_sigma_bx = 0.0_num

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
        cpml_x_min_laser_idx = cpml_thickness + ng + 1 - nx_global_min
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
      IF (nx_global_min .LE. nx_global - cpml_thickness - ng + 2 &
          .AND. nx_global_max .GE. nx_global - cpml_thickness - ng + 2) THEN
        add_laser(i) = .TRUE.
        cpml_x_max_laser_idx = &
            nx_global - cpml_thickness - ng + 2 - nx_global_min
      ENDIF
    ENDIF

  END SUBROUTINE set_cpml_helpers



  SUBROUTINE allocate_cpml_fields

    ! I will ignore memory consumption issues and, for simplicity,
    ! allocate the boundary fields throughout the whole simulation box.
    ALLOCATE(cpml_psi_eyx(-2:nx+3), cpml_psi_ezx(-2:nx+3))
    ALLOCATE(cpml_psi_byx(-2:nx+3), cpml_psi_bzx(-2:nx+3))

    cpml_psi_eyx = 0.0_num
    cpml_psi_ezx = 0.0_num
    cpml_psi_byx = 0.0_num
    cpml_psi_bzx = 0.0_num

  END SUBROUTINE allocate_cpml_fields



  SUBROUTINE deallocate_cpml_helpers

    DEALLOCATE(cpml_kappa_ex, cpml_kappa_bx)
    DEALLOCATE(cpml_a_ex, cpml_a_bx)
    DEALLOCATE(cpml_sigma_ex, cpml_sigma_bx)

  END SUBROUTINE deallocate_cpml_helpers



  SUBROUTINE cpml_advance_e_currents(tstep)

    REAL(num), INTENT(IN) :: tstep
    INTEGER :: ipos, ix
    REAL(num) :: acoeff, bcoeff, ccoeff_d
    REAL(num) :: kappa, sigma

    ! ============= x_min boundary =============

    IF (bc_field(c_bd_x_min) .EQ. c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_min) .EQ. c_bc_cpml_outflow) THEN
      DO ipos = cpml_x_min_start,cpml_x_min_end
        kappa = cpml_kappa_ex(ipos)
        sigma = cpml_sigma_ex(ipos)
        acoeff = cpml_a_ex(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dx

        cpml_psi_eyx(ipos) = bcoeff * cpml_psi_eyx(ipos) &
            + ccoeff_d * (bz(ipos) - bz(ipos-1))
        cpml_psi_ezx(ipos) = bcoeff * cpml_psi_ezx(ipos) &
            + ccoeff_d * (by(ipos) - by(ipos-1))
      ENDDO
    ENDIF

    ! ============= x_max boundary =============

    IF (bc_field(c_bd_x_max) .EQ. c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_max) .EQ. c_bc_cpml_outflow) THEN
      DO ipos = cpml_x_max_start,cpml_x_max_end
        kappa = cpml_kappa_ex(ipos)
        sigma = cpml_sigma_ex(ipos)
        acoeff = cpml_a_ex(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dx

        cpml_psi_eyx(ipos) = bcoeff * cpml_psi_eyx(ipos) &
            + ccoeff_d * (bz(ipos) - bz(ipos-1))
        cpml_psi_ezx(ipos) = bcoeff * cpml_psi_ezx(ipos) &
            + ccoeff_d * (by(ipos) - by(ipos-1))
      ENDDO
    ENDIF

  END SUBROUTINE cpml_advance_e_currents



  SUBROUTINE cpml_advance_b_currents(tstep)

    REAL(num), INTENT(IN) :: tstep
    INTEGER :: ipos, ix
    REAL(num) :: acoeff, bcoeff, ccoeff_d
    REAL(num) :: kappa, sigma

    ! ============= x_min boundary =============

    IF (bc_field(c_bd_x_min) .EQ. c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_min) .EQ. c_bc_cpml_outflow) THEN
      DO ipos = cpml_x_min_start,cpml_x_min_end
        kappa = cpml_kappa_bx(ipos)
        sigma = cpml_sigma_bx(ipos)
        acoeff = cpml_a_bx(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dx

        cpml_psi_byx(ipos) = bcoeff * cpml_psi_byx(ipos) &
            + ccoeff_d * (ez(ipos+1) - ez(ipos))
        cpml_psi_bzx(ipos) = bcoeff * cpml_psi_bzx(ipos) &
            + ccoeff_d * (ey(ipos+1) - ey(ipos))
      ENDDO
    ENDIF

    ! ============= x_max boundary =============

    IF (bc_field(c_bd_x_max) .EQ. c_bc_cpml_laser &
        .OR. bc_field(c_bd_x_max) .EQ. c_bc_cpml_outflow) THEN
      DO ipos = cpml_x_max_start-1,cpml_x_max_end-1
        kappa = cpml_kappa_bx(ipos)
        sigma = cpml_sigma_bx(ipos)
        acoeff = cpml_a_bx(ipos)
        bcoeff = EXP(-(sigma / kappa + acoeff) * tstep)
        ccoeff_d = (bcoeff - 1.0_num) * sigma / kappa &
            / (sigma + kappa * acoeff) / dx

        cpml_psi_byx(ipos) = bcoeff * cpml_psi_byx(ipos) &
            + ccoeff_d * (ez(ipos+1) - ez(ipos))
        cpml_psi_bzx(ipos) = bcoeff * cpml_psi_bzx(ipos) &
            + ccoeff_d * (ey(ipos+1) - ey(ipos))
      ENDDO
    ENDIF

  END SUBROUTINE cpml_advance_b_currents

END MODULE boundary
