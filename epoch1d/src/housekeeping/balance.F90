! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2012      Martin Ramsay <M.G.Ramsay@warwick.ac.uk>
! Copyright (C) 2009      Chris Brady <C.S.Brady@warwick.ac.uk>
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE balance

  USE boundary
  USE mpi_subtype_control
  USE redblack_module
  USE timer
  USE utilities

  IMPLICIT NONE

  INTEGER, DIMENSION(:), ALLOCATABLE :: new_cell_x_min, new_cell_x_max
  LOGICAL :: overriding
  REAL(num) :: load_av

CONTAINS

  SUBROUTINE balance_workload(over_ride)

    ! This subroutine determines whether or not the code needs rebalancing,
    ! calculates where to split the domain and calls other subroutines to
    ! actually rearrange the fields and particles onto the new processors

    LOGICAL, INTENT(IN) :: over_ride
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: load_x
    REAL(num) :: balance_frac, balance_frac_final, balance_improvement
    REAL(num) :: load_local, load_sum, load_max
    INTEGER(i8) :: npart_local
    INTEGER, SAVE :: balance_check_frequency = 1
    INTEGER, SAVE :: last_check = -HUGE(1) / 2
    INTEGER, SAVE :: last_full_check = -HUGE(1) / 2
    LOGICAL, SAVE :: first_flag = .TRUE.
    LOGICAL :: first_message, restarting, full_check
    LOGICAL :: use_redistribute_domain, use_redistribute_particles
#ifdef PARTICLE_DEBUG
    TYPE(particle), POINTER :: current
    INTEGER :: ispecies
#endif

    ! On one processor do nothing to save time
    IF (nproc == 1) RETURN

    full_check = over_ride
    IF (step - last_full_check < dlb_force_interval) THEN
      IF (step - last_check < balance_check_frequency) RETURN
    ELSE
      full_check = .TRUE.
    END IF

    restarting = .FALSE.
    use_redistribute_domain = .FALSE.

    IF (first_flag) THEN
      first_flag = .FALSE.
      first_message = .TRUE.
      use_redistribute_particles = .TRUE.
      IF (ic_from_restart) restarting = .TRUE.
    ELSE
      first_message = .FALSE.
      use_redistribute_particles = .FALSE.
    END IF
    last_check = step

    ! count particles
    npart_local = get_total_local_particles()
    load_local = REAL(push_per_field * npart_local + nx, num)

    CALL MPI_ALLREDUCE(load_local, load_max, 1, mpireal, MPI_MAX, comm, errcode)
    CALL MPI_ALLREDUCE(load_local, load_sum, 1, mpireal, MPI_SUM, comm, errcode)

    load_av = load_sum / nproc

    balance_frac = (load_av + SQRT(load_av)) / (load_max + SQRT(load_max))

    ! The over_ride flag allows the code to force a load balancing sweep
    ! at t = 0
    IF (.NOT. full_check .AND. balance_frac > dlb_threshold) THEN
      balance_check_frequency = &
          MIN(balance_check_frequency * 2, dlb_maximum_interval)
      IF (rank == 0) THEN
        PRINT'(''Skipping redistribution. Balance:'', F6.3, &
              &'', threshold:'', F6.3, '', next: '', i9)', &
              balance_frac, dlb_threshold, &
              MIN(step + balance_check_frequency, &
                  last_full_check + dlb_force_interval)
      END IF
      RETURN
    END IF

    IF (timer_collect) CALL timer_start(c_timer_balance)

    IF (.NOT.use_exact_restart) THEN
      overriding = full_check

      ALLOCATE(new_cell_x_min(nprocx), new_cell_x_max(nprocx))

      new_cell_x_min = cell_x_min
      new_cell_x_max = cell_x_max

      ! Sweep in X
      IF (nprocx > 1) THEN
        ! Rebalancing in X
        ALLOCATE(load_x(nx_global + 2 * ng))
        CALL get_load_in_x(load_x)
        CALL calculate_breaks(load_x, nprocx, new_cell_x_min, new_cell_x_max)
        DEALLOCATE(load_x)
      END IF

      IF (.NOT.restarting) THEN
        CALL calculate_new_load_imbalance(balance_frac, balance_frac_final)
        balance_improvement = (balance_frac_final - balance_frac) / balance_frac

        last_full_check = step

        ! Consider load balancing a success if the load imbalance improved by
        ! more than 5 percent
        IF (balance_improvement > 0.05_num) THEN
          use_redistribute_domain = .TRUE.
          use_redistribute_particles = .TRUE.
          first_message = .FALSE.
          balance_check_frequency = 1
          IF (rank == 0) THEN
            PRINT'(''Redistributing.          Balance:'', F6.3, &
                  &'',     after:'', F6.3, '', next: '', i9)', &
                  balance_frac, balance_frac_final, &
                  MIN(step + balance_check_frequency, &
                      last_full_check + dlb_force_interval)
          END IF
        ELSE
          IF (.NOT.first_message) THEN
            balance_check_frequency = &
                MIN(balance_check_frequency * 2, dlb_maximum_interval)
            IF (rank == 0) THEN
              PRINT'(''Skipping redistribution. Balance:'', F6.3, &
                    &'',     after:'', F6.3, '', next: '', i9)', &
                    balance_frac, balance_frac_final, &
                    MIN(step + balance_check_frequency, &
                        last_full_check + dlb_force_interval)
            END IF
          END IF
        END IF
      END IF
    END IF

    IF (use_redistribute_domain) THEN
      CALL redistribute_domain
    END IF

    IF (ALLOCATED(new_cell_x_min)) THEN
      DEALLOCATE(new_cell_x_min, new_cell_x_max)
    END IF

    ! Redistribute the particles onto their new processors
    IF (use_redistribute_particles) CALL distribute_particles

    ! If running with particle debugging then set the t = 0 processor if
    ! over_ride = true
#ifdef PARTICLE_DEBUG
    IF (full_check) THEN
      DO ispecies = 1, n_species
        current => species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          current%processor_at_t0 = rank
          current => current%next
        END DO
      END DO
    END IF
#endif

    IF (first_message) THEN
      npart_local = get_total_local_particles()
      load_local = REAL(push_per_field * npart_local + nx, num)

      CALL MPI_ALLREDUCE(load_local, load_max, 1, mpireal, MPI_MAX, &
          comm, errcode)

      balance_frac_final = (load_av + SQRT(load_av)) &
          / (load_max + SQRT(load_max))
      balance_check_frequency = 1

      IF (rank == 0) THEN
        PRINT'(''Redistributing.          Balance:'', F6.3, &
              &'',     next: '', i9)', &
              balance_frac_final, (step + balance_check_frequency)
      END IF
    END IF

    use_exact_restart = .FALSE.

    IF (timer_collect) CALL timer_stop(c_timer_balance)

  END SUBROUTINE balance_workload



  SUBROUTINE redistribute_domain

    INTEGER, DIMENSION(c_ndims,2) :: domain
    INTEGER :: iproc

    IF (.NOT.ALLOCATED(new_cell_x_min)) RETURN

    ! Now need to calculate the start and end points for the new domain on
    ! the current processor

    domain(1,:) = (/new_cell_x_min(x_coords+1), new_cell_x_max(x_coords+1)/)

    ! Redistribute the field variables
    CALL redistribute_fields(domain)

    ! Copy the new lengths into the permanent variables
    cell_x_min = new_cell_x_min
    cell_x_max = new_cell_x_max

    ! Set the new nx
    nx_global_min = cell_x_min(x_coords+1)
    nx_global_max = cell_x_max(x_coords+1)
    n_global_min(1) = nx_global_min
    n_global_max(1) = nx_global_max

    nx = nx_global_max - nx_global_min + 1

    ! Do X array separately because we already have global copies
    DEALLOCATE(x)
    ALLOCATE(x(1-ng:nx+ng))
    x(1-ng:nx+ng) = x_global(nx_global_min-ng:nx_global_max+ng)

    DEALLOCATE(xb)
    ALLOCATE(xb(1-ng:nx+ng))
    xb(1-ng:nx+ng) = xb_global(nx_global_min-ng:nx_global_max+ng)

    ! Recalculate x_grid_mins/maxs so that rebalancing works next time
    DO iproc = 0, nprocx - 1
      x_grid_mins(iproc) = x_global(cell_x_min(iproc+1))
      x_grid_maxs(iproc) = x_global(cell_x_max(iproc+1))
    END DO

    ! Set the lengths of the current domain so that the particle balancer
    ! works properly
    x_grid_min_local = x_grid_mins(x_coords)
    x_grid_max_local = x_grid_maxs(x_coords)

    x_min_local = x_grid_min_local + (cpml_x_min_offset - 0.5_num) * dx
    x_max_local = x_grid_max_local - (cpml_x_max_offset - 0.5_num) * dx

  END SUBROUTINE redistribute_domain



  SUBROUTINE redistribute_fields(new_domain)

    ! This subroutine redistributes the field variables over the new
    ! processor layout. If using a field of your own then set the
    ! redistribute_field subroutine to implement it.

    INTEGER :: nx_new
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: new_domain
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temp_sum
    REAL(r4), DIMENSION(:,:), ALLOCATABLE :: r4temp_sum
    REAL(num), DIMENSION(:), ALLOCATABLE :: temp, temp2
    INTEGER :: i, ispecies, io, id, nspec_local, mask

    nx_new = new_domain(1,2) - new_domain(1,1) + 1

    ! The following code is quite messy and repetitive. Unfortunately, the
    ! F90 standard does not allow the ALLOCATABLE attribute for subroutine
    ! arguments and POINTER arrays are not as fast.

    ! Full domain arrays

    ALLOCATE(temp(1-ng:nx_new+ng))

    ! Current will be recalculated during the particle push, so there
    ! is no need to copy the contents of the old arrays.
    ! If overriding, then we may not be doing a particle push next
    ! so we still have to balance the arrays.
    ! It is done slightly differently since the arrays may be
    ! a different size.

    IF (overriding) THEN
      ALLOCATE(temp2(1-ng:nx+ng))

      temp2(0:nx+1) = jx(0:nx+1)
      CALL remap_field(temp2, temp)
      DEALLOCATE(jx)
      ALLOCATE(jx(1-jng:nx_new+jng))
      jx(0:nx_new+1) = temp(0:nx_new+1)

      temp2(0:nx+1) = jy(0:nx+1)
      CALL remap_field(temp2, temp)
      DEALLOCATE(jy)
      ALLOCATE(jy(1-jng:nx_new+jng))
      jy(0:nx_new+1) = temp(0:nx_new+1)

      temp2(0:nx+1) = jz(0:nx+1)
      CALL remap_field(temp2, temp)
      DEALLOCATE(jz)
      ALLOCATE(jz(1-jng:nx_new+jng))
      jz(0:nx_new+1) = temp(0:nx_new+1)

      DEALLOCATE(temp2)
    ELSE
      DEALLOCATE(jx)
      DEALLOCATE(jy)
      DEALLOCATE(jz)
      ALLOCATE(jx(1-jng:nx_new+jng))
      ALLOCATE(jy(1-jng:nx_new+jng))
      ALLOCATE(jz(1-jng:nx_new+jng))
    END IF

    CALL remap_field(ex, temp)
    DEALLOCATE(ex)
    ALLOCATE(ex(1-ng:nx_new+ng))
    ex = temp

    CALL remap_field(ey, temp)
    DEALLOCATE(ey)
    ALLOCATE(ey(1-ng:nx_new+ng))
    ey = temp

    CALL remap_field(ez, temp)
    DEALLOCATE(ez)
    ALLOCATE(ez(1-ng:nx_new+ng))
    ez = temp

    CALL remap_field(bx, temp)
    DEALLOCATE(bx)
    ALLOCATE(bx(1-ng:nx_new+ng))
    bx = temp

    CALL remap_field(by, temp)
    DEALLOCATE(by)
    ALLOCATE(by(1-ng:nx_new+ng))
    by = temp

    CALL remap_field(bz, temp)
    DEALLOCATE(bz)
    ALLOCATE(bz(1-ng:nx_new+ng))
    bz = temp

    DO ispecies = 1, n_species
      IF (species_list(ispecies)%migrate%fluid) THEN
        CALL remap_field(species_list(ispecies)%migrate%fluid_energy, temp)
        DEALLOCATE(species_list(ispecies)%migrate%fluid_energy)
        ALLOCATE(species_list(ispecies)%migrate%fluid_energy(1-ng:nx_new+ng))
        species_list(ispecies)%migrate%fluid_energy = temp

        CALL remap_field(species_list(ispecies)%migrate%fluid_density, temp)
        DEALLOCATE(species_list(ispecies)%migrate%fluid_density)
        ALLOCATE(species_list(ispecies)%migrate%fluid_density(1-ng:nx_new+ng))
        species_list(ispecies)%migrate%fluid_density = temp
      END IF
    END DO

    IF (cpml_boundaries) THEN
      CALL remap_field(cpml_psi_eyx, temp)
      DEALLOCATE(cpml_psi_eyx)
      ALLOCATE(cpml_psi_eyx(1-ng:nx_new+ng))
      cpml_psi_eyx = temp

      CALL remap_field(cpml_psi_byx, temp)
      DEALLOCATE(cpml_psi_byx)
      ALLOCATE(cpml_psi_byx(1-ng:nx_new+ng))
      cpml_psi_byx = temp

      CALL remap_field(cpml_psi_ezx, temp)
      DEALLOCATE(cpml_psi_ezx)
      ALLOCATE(cpml_psi_ezx(1-ng:nx_new+ng))
      cpml_psi_ezx = temp

      CALL remap_field(cpml_psi_bzx, temp)
      DEALLOCATE(cpml_psi_bzx)
      ALLOCATE(cpml_psi_bzx(1-ng:nx_new+ng))
      cpml_psi_bzx = temp

      CALL deallocate_cpml_helpers
      CALL set_cpml_helpers(nx_new, new_domain(1,1), new_domain(1,2))
    END IF

    DEALLOCATE(temp)

    ! Full domain arrays with an additional index

    DO id = 1, num_vars_to_dump
      io = averaged_var_block(id)
      IF (io == 0) CYCLE

      mask = io_block_list(io)%dumpmask(id)
      nspec_local = 0
      IF (IAND(mask, c_io_no_sum) == 0) &
          nspec_local = 1
      IF (IAND(mask, c_io_species) /= 0) &
          nspec_local = nspec_local + n_species

      IF (nspec_local <= 0) CYCLE

      IF (io_block_list(io)%averaged_data(id)%dump_single) THEN
        IF (.NOT. ASSOCIATED(io_block_list(io)%averaged_data(id)%r4array)) CYCLE

        ALLOCATE(r4temp_sum(1-ng:nx_new+ng, nspec_local))

        DO i = 1, nspec_local
          CALL remap_field_r4(&
              io_block_list(io)%averaged_data(id)%r4array(:,i), &
              r4temp_sum(:,i))
        END DO

        DEALLOCATE(io_block_list(io)%averaged_data(id)%r4array)
        ALLOCATE(io_block_list(io)%averaged_data(id)&
            %r4array(1-ng:nx_new+ng, nspec_local))

        io_block_list(io)%averaged_data(id)%r4array = r4temp_sum

        DEALLOCATE(r4temp_sum)
      ELSE
        IF (.NOT. ASSOCIATED(io_block_list(io)%averaged_data(id)%array)) CYCLE

        ALLOCATE(temp_sum(1-ng:nx_new+ng, nspec_local))

        DO i = 1, nspec_local
          CALL remap_field(&
              io_block_list(io)%averaged_data(id)%array(:,i), &
              temp_sum(:,i))
        END DO

        DEALLOCATE(io_block_list(io)%averaged_data(id)%array)
        ALLOCATE(io_block_list(io)%averaged_data(id)&
            %array(1-ng:nx_new+ng, nspec_local))

        io_block_list(io)%averaged_data(id)%array = temp_sum

        DEALLOCATE(temp_sum)
      END IF
    END DO

  END SUBROUTINE redistribute_fields



  SUBROUTINE remap_field(field_in, field_out)

    ! This is a wrapper for the field redistribution routine
    REAL(num), DIMENSION(:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(c_ndims) :: n_new, cdim
    INTEGER :: i

    n_new = SHAPE(field_out) - 2 * ng

    DO i = 1, c_ndims
      cdim(i) = c_ndims + 1 - i
    END DO

    CALL redistribute_field_1d(field_in, field_out, cdim, &
        cell_x_min, cell_x_max, new_cell_x_min, new_cell_x_max)

    CALL do_field_mpi_with_lengths(field_out, ng, n_new(1))

  END SUBROUTINE remap_field



  SUBROUTINE remap_field_r4(field_in, field_out)

    ! This is a wrapper for the field redistribution routine
    REAL(r4), DIMENSION(:), INTENT(IN) :: field_in
    REAL(r4), DIMENSION(:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(c_ndims) :: n_new, cdim
    INTEGER :: i

    n_new = SHAPE(field_out) - 2 * ng

    DO i = 1, c_ndims
      cdim(i) = c_ndims + 1 - i
    END DO

    CALL redistribute_field_1d_r4(field_in, field_out, cdim, &
        cell_x_min, cell_x_max, new_cell_x_min, new_cell_x_max)

    CALL do_field_mpi_with_lengths_r4(field_out, ng, n_new(1))

  END SUBROUTINE remap_field_r4



  SUBROUTINE redistribute_field_1d(field_in, field_out, cdim, &
      old_cell_min1, old_cell_max1, new_cell_min1, new_cell_max1)

    ! This subroutine redistributes the fields over the new processor layout
    INTEGER, PARAMETER :: nd = 1
    REAL(num), DIMENSION(1-ng:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(1-ng:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(nd), INTENT(IN) :: cdim
    INTEGER, DIMENSION(:), INTENT(IN) :: old_cell_min1, old_cell_max1
    INTEGER, DIMENSION(:), INTENT(IN) :: new_cell_min1, new_cell_max1
    INTEGER :: irank, basetype, n, ng0, ng1
    INTEGER :: i, iproc, inew
    INTEGER, DIMENSION(nd) :: type_min, type_max, old_0, old_1, new_0
    INTEGER, DIMENSION(nd) :: n_global, n_local, start, nprocs
    INTEGER, DIMENSION(nd) :: old_min, old_max, new_min, new_max
    INTEGER, DIMENSION(c_ndims) :: coord
    INTEGER, DIMENSION(nd) :: our_coords, nmin, nmax
    INTEGER, DIMENSION(:), ALLOCATABLE :: sendtypes, recvtypes

    basetype = mpireal

    ALLOCATE(sendtypes(0:nproc-1))
    ALLOCATE(recvtypes(0:nproc-1))

    DO i = 1, nd
      our_coords(i) = coordinates(cdim(i))
    END DO

    nprocs(1) = SIZE(old_cell_min1)

    old_min(1) = old_cell_min1(our_coords(1)+1)
    old_max(1) = old_cell_max1(our_coords(1)+1)
    new_min(1) = new_cell_min1(our_coords(1)+1)
    new_max(1) = new_cell_max1(our_coords(1)+1)

    nmin(1) = new_cell_min1(1)
    nmax(1) = new_cell_max1(nprocs(1))

    tag = 0
    sendtypes = 0
    recvtypes = 0
    coord = coordinates

    ! Create array of sendtypes

    DO i = 1,nd
      n_global(i) = old_max(i) - old_min(i) + 2 * ng + 1
    END DO

    n = 1
    type_min(n) = old_min(n)
    type_max(n) = old_min(n)

    ! Find the new processor on which the old x_min resides
    ! This could be sped up by using bisection.
    DO iproc = 1, nprocs(n)-1
      IF (new_cell_min1(iproc) <= old_min(n) &
          .AND. new_cell_max1(iproc) >= old_min(n)) EXIT
    END DO

    DO WHILE(type_max(n) <= old_max(n))
      coord(cdim(n)) = iproc - 1
      type_max(n) = new_cell_max1(iproc)
      IF (type_max(n) > old_max(n)) type_max(n) = old_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) == nmin(n)) ng0 = ng
      IF (type_max(n) == nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - old_min(n) + ng - ng0 + 1

      CALL MPI_CART_RANK(comm, coord, irank, errcode)

      IF (rank /= irank) THEN
        sendtypes(irank) = create_1d_array_subtype(basetype, n_local, &
            n_global, start)
      ELSE
        ! New domain is on the same processor as the old domain.
        ! Just copy the region rather than using MPI.
        DO i = 1,nd
          old_0(i) = start(i) - ng
          old_1(i) = old_0(i) + n_local(i) - 1
        END DO
      END IF

      n = 1
      IF (type_max(n) == old_max(n)) EXIT
      iproc = iproc + 1
      type_min(n) = new_cell_min1(iproc)
    END DO

    ! Create array of recvtypes

    DO i = 1,nd
      n_global(i) = new_max(i) - new_min(i) + 2 * ng + 1
    END DO

    n = 1
    type_min(n) = new_min(n)
    type_max(n) = new_min(n)

    ! Find the old processor on which the new x_min resides
    ! This could be sped up by using bisection.
    DO iproc = 1, nprocs(n)-1
      IF (old_cell_min1(iproc) <= new_min(n) &
          .AND. old_cell_max1(iproc) >= new_min(n)) EXIT
    END DO

    DO WHILE(type_max(n) <= new_max(n))
      coord(cdim(n)) = iproc - 1
      type_max(n) = old_cell_max1(iproc)
      IF (type_max(n) > new_max(n)) type_max(n) = new_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) == nmin(n)) ng0 = ng
      IF (type_max(n) == nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - new_min(n) + ng - ng0 + 1

      CALL MPI_CART_RANK(comm, coord, irank, errcode)

      IF (rank /= irank) THEN
        recvtypes(irank) = create_1d_array_subtype(basetype, n_local, &
            n_global, start)
      ELSE
        ! New domain is on the same processor as the old domain.
        ! Just copy the region rather than using MPI.
        DO i = 1,nd
          new_0(i) = start(i) - ng
        END DO
        DO i = old_0(1),old_1(1)
          inew = new_0(1) + i - old_0(1)
          field_out(inew) = field_in(i)
        END DO
      END IF

      n = 1
      IF (type_max(n) == new_max(n)) EXIT
      iproc = iproc + 1
      type_min(n) = old_cell_min1(iproc)
    END DO

    CALL redblack(field_in, field_out, sendtypes, recvtypes)

    DO i = 0,nproc-1
      IF (sendtypes(i) /= 0) CALL MPI_TYPE_FREE(sendtypes(i), errcode)
      IF (recvtypes(i) /= 0) CALL MPI_TYPE_FREE(recvtypes(i), errcode)
    END DO

    DEALLOCATE(sendtypes)
    DEALLOCATE(recvtypes)

  END SUBROUTINE redistribute_field_1d



  SUBROUTINE redistribute_field_1d_r4(field_in, field_out, cdim, &
      old_cell_min1, old_cell_max1, new_cell_min1, new_cell_max1)

    ! This subroutine redistributes the fields over the new processor layout
    INTEGER, PARAMETER :: nd = 1
    REAL(r4), DIMENSION(1-ng:), INTENT(IN) :: field_in
    REAL(r4), DIMENSION(1-ng:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(nd), INTENT(IN) :: cdim
    INTEGER, DIMENSION(:), INTENT(IN) :: old_cell_min1, old_cell_max1
    INTEGER, DIMENSION(:), INTENT(IN) :: new_cell_min1, new_cell_max1
    INTEGER :: irank, basetype, n, ng0, ng1
    INTEGER :: i, iproc, inew
    INTEGER, DIMENSION(nd) :: type_min, type_max, old_0, old_1, new_0
    INTEGER, DIMENSION(nd) :: n_global, n_local, start, nprocs
    INTEGER, DIMENSION(nd) :: old_min, old_max, new_min, new_max
    INTEGER, DIMENSION(c_ndims) :: coord
    INTEGER, DIMENSION(nd) :: our_coords, nmin, nmax
    INTEGER, DIMENSION(:), ALLOCATABLE :: sendtypes, recvtypes

    basetype = MPI_REAL4

    ALLOCATE(sendtypes(0:nproc-1))
    ALLOCATE(recvtypes(0:nproc-1))

    DO i = 1, nd
      our_coords(i) = coordinates(cdim(i))
    END DO

    nprocs(1) = SIZE(old_cell_min1)

    old_min(1) = old_cell_min1(our_coords(1)+1)
    old_max(1) = old_cell_max1(our_coords(1)+1)
    new_min(1) = new_cell_min1(our_coords(1)+1)
    new_max(1) = new_cell_max1(our_coords(1)+1)

    nmin(1) = new_cell_min1(1)
    nmax(1) = new_cell_max1(nprocs(1))

    tag = 0
    sendtypes = 0
    recvtypes = 0
    coord = coordinates

    ! Create array of sendtypes

    DO i = 1,nd
      n_global(i) = old_max(i) - old_min(i) + 2 * ng + 1
    END DO

    n = 1
    type_min(n) = old_min(n)
    type_max(n) = old_min(n)

    ! Find the new processor on which the old x_min resides
    ! This could be sped up by using bisection.
    DO iproc = 1, nprocs(n)-1
      IF (new_cell_min1(iproc) <= old_min(n) &
          .AND. new_cell_max1(iproc) >= old_min(n)) EXIT
    END DO

    DO WHILE(type_max(n) <= old_max(n))
      coord(cdim(n)) = iproc - 1
      type_max(n) = new_cell_max1(iproc)
      IF (type_max(n) > old_max(n)) type_max(n) = old_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) == nmin(n)) ng0 = ng
      IF (type_max(n) == nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - old_min(n) + ng - ng0 + 1

      CALL MPI_CART_RANK(comm, coord, irank, errcode)

      IF (rank /= irank) THEN
        sendtypes(irank) = create_1d_array_subtype(basetype, n_local, &
            n_global, start)
      ELSE
        ! New domain is on the same processor as the old domain.
        ! Just copy the region rather than using MPI.
        DO i = 1,nd
          old_0(i) = start(i) - ng
          old_1(i) = old_0(i) + n_local(i) - 1
        END DO
      END IF

      n = 1
      IF (type_max(n) == old_max(n)) EXIT
      iproc = iproc + 1
      type_min(n) = new_cell_min1(iproc)
    END DO

    ! Create array of recvtypes

    DO i = 1,nd
      n_global(i) = new_max(i) - new_min(i) + 2 * ng + 1
    END DO

    n = 1
    type_min(n) = new_min(n)
    type_max(n) = new_min(n)

    ! Find the old processor on which the new x_min resides
    ! This could be sped up by using bisection.
    DO iproc = 1, nprocs(n)-1
      IF (old_cell_min1(iproc) <= new_min(n) &
          .AND. old_cell_max1(iproc) >= new_min(n)) EXIT
    END DO

    DO WHILE(type_max(n) <= new_max(n))
      coord(cdim(n)) = iproc - 1
      type_max(n) = old_cell_max1(iproc)
      IF (type_max(n) > new_max(n)) type_max(n) = new_max(n)

      ng0 = 0
      ng1 = 0
      IF (type_min(n) == nmin(n)) ng0 = ng
      IF (type_max(n) == nmax(n)) ng1 = ng

      n_local(n) = type_max(n) - type_min(n) + ng0 + ng1 + 1
      start(n) = type_min(n) - new_min(n) + ng - ng0 + 1

      CALL MPI_CART_RANK(comm, coord, irank, errcode)

      IF (rank /= irank) THEN
        recvtypes(irank) = create_1d_array_subtype(basetype, n_local, &
            n_global, start)
      ELSE
        ! New domain is on the same processor as the old domain.
        ! Just copy the region rather than using MPI.
        DO i = 1,nd
          new_0(i) = start(i) - ng
        END DO
        DO i = old_0(1),old_1(1)
          inew = new_0(1) + i - old_0(1)
          field_out(inew) = field_in(i)
        END DO
      END IF

      n = 1
      IF (type_max(n) == new_max(n)) EXIT
      iproc = iproc + 1
      type_min(n) = old_cell_min1(iproc)
    END DO

    CALL redblack(field_in, field_out, sendtypes, recvtypes)

    DO i = 0,nproc-1
      IF (sendtypes(i) /= 0) CALL MPI_TYPE_FREE(sendtypes(i), errcode)
      IF (recvtypes(i) /= 0) CALL MPI_TYPE_FREE(recvtypes(i), errcode)
    END DO

    DEALLOCATE(sendtypes)
    DEALLOCATE(recvtypes)

  END SUBROUTINE redistribute_field_1d_r4



  SUBROUTINE get_load_in_x(load)

    ! Calculate total load across the X direction

    INTEGER(i8), DIMENSION(:), INTENT(OUT) :: load
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: temp
    TYPE(particle), POINTER :: current
    INTEGER :: cell, ispecies, sz

    load = 0

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so x_grid_min, NOT x_grid_min_local
        cell = FLOOR((current%part_pos - x_grid_min) / dx + 1.5_num) + ng

        load(cell) = load(cell) + 1
        current => current%next
      END DO
    END DO

    ! Now have local densities, so add using MPI
    sz = SIZE(load)
    ALLOCATE(temp(sz))
    CALL MPI_ALLREDUCE(load, temp, sz, MPI_INTEGER8, MPI_SUM, comm, errcode)

    ! Adjust the load of pushing one particle relative to the load
    ! of updating one field cell, then add on the field load.
    ! The push_per_field factor will be updated automatically in future.
    load = push_per_field * temp
    load(ng+1:sz-ng) = load(ng+1:sz-ng) + 1

    DEALLOCATE(temp)

  END SUBROUTINE get_load_in_x



  SUBROUTINE calculate_breaks(load, nproc, mins, maxs)

    ! This subroutine calculates the places in a given load profile to split
    ! The domain to give the most even subdivision possible

    INTEGER(i8), INTENT(IN), DIMENSION(1-ng:) :: load
    INTEGER, INTENT(IN) :: nproc
    INTEGER, DIMENSION(:), INTENT(OUT) :: mins, maxs
    INTEGER :: sz, idim, proc, old, nextra
    INTEGER(i8) :: total, total_old, load_per_proc_ideal

    sz = SIZE(load) - 2 * ng
    maxs = sz

    load_per_proc_ideal = FLOOR(REAL(SUM(load(1:sz)), num) / nproc + 0.5d0, i8)

    proc = 1
    total = 0
    old = 1
    nextra = 0
    DO idim = 1, sz
      IF (nextra > 0) THEN
        nextra = nextra - 1
        CYCLE
      END IF
      total_old = total
      total = total + load(idim)
      IF (total >= load_per_proc_ideal) THEN
        ! Pick the split that most closely matches the load
        IF (load_per_proc_ideal - total_old &
            < total - load_per_proc_ideal) THEN
          maxs(proc) = idim - 1
        ELSE
          maxs(proc) = idim
        END IF
        ! To communicate ghost cell information correctly, each domain must
        ! contain at least ng cells.
        nextra = old - maxs(proc) + ng
        IF (nextra > 0) THEN
          maxs(proc) = maxs(proc) + nextra
        END IF
        proc = proc + 1
        IF (proc == nproc) EXIT
        total = 0
        old = maxs(proc-1)
        load_per_proc_ideal = FLOOR(REAL(SUM(load(MIN(idim+1,sz):)), num) &
            / (nproc - proc + 1) + 0.5d0, i8)
      END IF
    END DO
    maxs(nproc) = sz

    ! Sanity check. Must be one cell of separation between each endpoint.
    ! Backwards
    old = sz
    DO proc = nproc-1, 1, -1
      IF (old - maxs(proc) < ng) THEN
        maxs(proc) = old - ng
      END IF
      old = maxs(proc)
    END DO

    ! Forwards (unnecessary?)
    old = 0
    DO proc = 1, nproc-1
      IF (maxs(proc) - old < ng) THEN
        maxs(proc) = old + ng
      END IF
      old = maxs(proc)
    END DO

    ! Set mins
    mins(1) = 1
    DO proc = 2, nproc
      mins(proc) = maxs(proc-1) + 1
    END DO

  END SUBROUTINE calculate_breaks



  FUNCTION get_particle_processor(part)

    ! This subroutine calculates which processor a given particles resides on

    TYPE(particle), INTENT(IN) :: part
    INTEGER :: get_particle_processor
    INTEGER :: iproc, coords(c_ndims)
    REAL(num) :: minpos, maxpos

    get_particle_processor = -1
    coords = -1

    ! This could be replaced by a bisection method, but for the moment I
    ! just don't care

    DO iproc = 0, nprocx - 1
      IF (iproc == 0) THEN
        minpos = x_grid_mins(iproc) - dx * (0.5_num + png)
      ELSE
        minpos = x_grid_mins(iproc) - dx * 0.5_num
      END IF
      IF (iproc == nprocx - 1) THEN
        maxpos = x_grid_maxs(iproc) + dx * (0.5_num + png)
      ELSE
        maxpos = x_grid_maxs(iproc) + dx * 0.5_num
      END IF
      IF (part%part_pos >= minpos .AND. part%part_pos < maxpos) THEN
        coords(c_ndims) = iproc
        EXIT
      END IF
    END DO

    IF (MINVAL(coords) < 0) THEN
      WRITE(*,*) 'UNLOCATABLE PARTICLE', coords
      RETURN
    END IF
    CALL MPI_CART_RANK(comm, coords, get_particle_processor, errcode)
    ! IF (get_particle_processor /= rank) PRINT *,

  END FUNCTION get_particle_processor



  ! This subroutine is used to rearrange particles over processors
  SUBROUTINE distribute_particles

    ! This subroutine moves particles which are on the wrong processor
    ! to the correct processor.

    TYPE(particle_list), DIMENSION(:), ALLOCATABLE :: pointers_send
    TYPE(particle_list), DIMENSION(:), ALLOCATABLE :: pointers_recv
    TYPE(particle), POINTER :: current, next
    INTEGER :: part_proc, iproc, ispecies
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: sendcounts, recvcounts

    ALLOCATE(pointers_send(0:nproc-1), pointers_recv(0:nproc-1))
    ALLOCATE(sendcounts(0:nproc-1), recvcounts(0:nproc-1))

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO iproc = 0, nproc - 1
        CALL create_empty_partlist(pointers_send(iproc))
        CALL create_empty_partlist(pointers_recv(iproc))
      END DO

      DO WHILE(ASSOCIATED(current))
        next => current%next
        part_proc = get_particle_processor(current)
        IF (part_proc < 0) THEN
          PRINT *, 'Unlocatable particle on processor', rank, current%part_pos
          CALL abort_code(c_err_bad_value)
          STOP
        END IF
#ifdef PARTICLE_DEBUG
        current%processor = part_proc
#endif
        IF (part_proc /= rank) THEN
          CALL remove_particle_from_partlist(&
              species_list(ispecies)%attached_list, current)
          CALL add_particle_to_partlist(pointers_send(part_proc), current)
        END IF
        current => next
      END DO

      DO iproc = 0, nproc - 1
        sendcounts(iproc) = pointers_send(iproc)%count
      END DO

      CALL MPI_ALLTOALL(sendcounts, 1, MPI_INTEGER8, recvcounts, 1, &
          MPI_INTEGER8, comm, errcode)

      CALL redblack(pointers_send, pointers_recv, sendcounts, recvcounts)

      DO iproc = 0, nproc - 1
        CALL append_partlist(species_list(ispecies)%attached_list, &
            pointers_recv(iproc))
      END DO
    END DO

    DEALLOCATE(sendcounts, recvcounts)
    DEALLOCATE(pointers_send, pointers_recv)

  END SUBROUTINE distribute_particles



  SUBROUTINE create_npart_per_cell(npart_per_cell)

    INTEGER(i8), ALLOCATABLE, INTENT(OUT) :: npart_per_cell(:)
    INTEGER :: ispecies
    INTEGER :: cell_x
    TYPE(particle), POINTER :: current, next
    INTEGER :: i0, i1

    i0 = 1 - ng
    IF (use_field_ionisation) i0 = -ng
    i1 = 1 - i0

    ALLOCATE(npart_per_cell(i0:nx+i1))
    npart_per_cell(:) = 0

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        next => current%next
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x = FLOOR((current%part_pos - x_grid_min_local) / dx) + 1
#else
        cell_x = FLOOR((current%part_pos - x_grid_min_local) / dx + 1.5_num)
#endif
        npart_per_cell(cell_x) = npart_per_cell(cell_x) + 1

        current => next
      END DO
    END DO

  END SUBROUTINE create_npart_per_cell



  SUBROUTINE calculate_new_load_imbalance(balance_frac, balance_frac_final)

    REAL(num), INTENT(INOUT) :: balance_frac, balance_frac_final
    REAL(num), ALLOCATABLE :: load_per_cpu(:)
    INTEGER(i8), ALLOCATABLE :: npart_per_cell(:)
    INTEGER(i8) :: npart_local
    INTEGER :: i, i0, i1, ix
    INTEGER :: ierr
    REAL(num) :: load_local, load_sum, load_max

    CALL create_npart_per_cell(npart_per_cell)

    IF (use_injectors) THEN
      npart_local = SUM(npart_per_cell(1:nx))
      load_local = REAL(push_per_field * npart_local + nx, num)

      CALL MPI_ALLREDUCE(load_local, load_max, 1, mpireal, MPI_MAX, comm, ierr)
      CALL MPI_ALLREDUCE(load_local, load_sum, 1, mpireal, MPI_SUM, comm, ierr)

      load_av = load_sum / nproc

      balance_frac = (load_av + SQRT(load_av)) / (load_max + SQRT(load_max))
    END IF

    ALLOCATE(load_per_cpu(nprocx))
    load_per_cpu = 0.0_num

    DO i = 1, nprocx
      i0 = new_cell_x_min(i) - nx_global_min + 1
      i1 = new_cell_x_max(i) - nx_global_min + 1

      IF (i1 < 1 .OR. i0 > nx) CYCLE

      i0 = MAX(i0, 1)
      i1 = MIN(i1, nx)

      DO ix = i0, i1
        load_per_cpu(i) = load_per_cpu(i) &
            + REAL(push_per_field * npart_per_cell(ix) + 1, num)
      END DO
    END DO

    DEALLOCATE(npart_per_cell)

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, load_per_cpu, nprocx, &
                       mpireal, MPI_SUM, comm, ierr)

    load_max = MAXVAL(load_per_cpu)

    balance_frac_final = (load_av + SQRT(load_av)) / (load_max + SQRT(load_max))

    DEALLOCATE(load_per_cpu)

  END SUBROUTINE calculate_new_load_imbalance

END MODULE balance
