! Copyright (C) 2009-2019 University of Warwick
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

MODULE injectors

  USE partlist
  USE particle_temperature
  USE evaluator
  USE random_generator
  USE utilities

  IMPLICIT NONE

  REAL(num) :: flow_limit_val = 10.0_num

CONTAINS

  SUBROUTINE init_injector(boundary, injector)

    INTEGER, INTENT(IN) :: boundary
    TYPE(injector_block), INTENT(INOUT) :: injector

    injector%npart_per_cell = -1.0_num
    injector%species = -1
    injector%boundary = boundary
    injector%t_start = 0.0_num
    injector%t_end = t_end
    injector%has_t_end = .FALSE.
    injector%density_min = 0.0_num
    injector%density_max = HUGE(1.0_num)
    injector%use_flux_injector = .TRUE.
    NULLIFY(injector%depth)
    NULLIFY(injector%next)

    IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
      ALLOCATE(injector%depth(1-ng:ny+ng))
    END IF

    IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
      ALLOCATE(injector%depth(1-ng:nx+ng))
    END IF

    injector%depth = 1.0_num
    need_random_state = .TRUE.

    injector%inject_from_file = .FALSE.
    injector%file_finished = .FALSE.
    injector%x_data_given = .FALSE.
    injector%y_data_given = .FALSE.
    injector%px_data_given = .FALSE.
    injector%py_data_given = .FALSE.
    injector%pz_data_given = .FALSE.
    injector%t_data_given = .FALSE.
#ifndef PER_SPECIES_WEIGHT
    injector%w_data_given = .FALSE.
#endif
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
    injector%id_data_given = .FALSE.
#endif

  END SUBROUTINE init_injector



  SUBROUTINE attach_injector(injector)

    TYPE(injector_block), POINTER :: injector
    INTEGER :: boundary

    boundary = injector%boundary

    IF (boundary == c_bd_x_min) THEN
      CALL attach_injector_to_list(injector_x_min, injector)
    ELSE IF (boundary == c_bd_x_max) THEN
      CALL attach_injector_to_list(injector_x_max, injector)
    ELSE IF (boundary == c_bd_y_min) THEN
      CALL attach_injector_to_list(injector_y_min, injector)
    ELSE IF (boundary == c_bd_y_max) THEN
      CALL attach_injector_to_list(injector_y_max, injector)
    END IF

  END SUBROUTINE attach_injector



  ! Actually does the attaching of the injector to the correct list
  SUBROUTINE attach_injector_to_list(list, injector)

    TYPE(injector_block), POINTER :: list
    TYPE(injector_block), POINTER :: injector
    TYPE(injector_block), POINTER :: current

    NULLIFY(injector%next)

    IF (ASSOCIATED(list)) THEN
      current => list
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      END DO
      current%next => injector
    ELSE
      list => injector
    END IF

  END SUBROUTINE attach_injector_to_list



  SUBROUTINE deallocate_injectors

    CALL deallocate_injector_list(injector_x_min)
    CALL deallocate_injector_list(injector_x_max)
    CALL deallocate_injector_list(injector_y_min)
    CALL deallocate_injector_list(injector_y_max)

  END SUBROUTINE deallocate_injectors



  SUBROUTINE deallocate_injector_list(list)

    TYPE(injector_block), POINTER :: list
    TYPE(injector_block), POINTER :: current, next
    INTEGER :: i

    current => list
    DO WHILE(ASSOCIATED(current))
      next => current%next
      IF (current%density_function%init) &
          CALL deallocate_stack(current%density_function)
      IF (ASSOCIATED(current%depth)) DEALLOCATE(current%depth)
      DO i = 1, 3
        IF (current%temperature_function(i)%init) &
            CALL deallocate_stack(current%temperature_function(i))
        IF (current%drift_function(i)%init) &
            CALL deallocate_stack(current%drift_function(i))
      END DO
      DEALLOCATE(current)
      current => next
    END DO

  END SUBROUTINE deallocate_injector_list



  SUBROUTINE open_injector_files(injector)

    ! Called in deck_injector_block after we have read the injector variables,
    ! and only for injectors with the inject_from_file flag. The file units are
    ! chosen such that each variable for each injector has a unique file unit.

    TYPE(injector_block), POINTER :: injector
    INTEGER :: inj_base_unit, iostat
    LOGICAL :: file_exists

    inj_base_unit = custom_base_unit + (injector%custom_id-1)*custom_var_num

    ! Injection times (mandatory)
    INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%t_data),&
        EXIST=file_exists)
    IF (file_exists) THEN
      OPEN(UNIT=inj_base_unit + unit_t, &
          FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%t_data))
      IF (rank == 0) THEN
        PRINT*, ''
        PRINT*, 'Successfully opened time file: ', TRIM(data_dir) // '/' // &
            TRIM(injector_filenames%t_data)
      END IF
    ELSE
      IF (rank == 0) THEN
        PRINT*, ''
        PRINT*, '*** ERROR ***'
        PRINT*, 'Unable to locate time file: ', &
            TRIM(data_dir) // '/' // TRIM(injector_filenames%t_data)
      END IF
    END IF

#ifndef PER_SPECIES_WEIGHT
    ! Particle weights (mandatory)
    INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%w_data),&
        EXIST=file_exists)
    IF (file_exists) THEN
      OPEN(UNIT=inj_base_unit + unit_w, &
          FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%w_data))
      IF (rank == 0) THEN
        PRINT*, ''
        PRINT*, 'Successfully opened weight file: ', TRIM(data_dir) // '/' // &
            TRIM(injector_filenames%w_data)
      END IF
    ELSE
      IF (rank == 0) THEN
        PRINT*, ''
        PRINT*, '*** ERROR ***'
        PRINT*, 'Unable to locate weight file: ', &
            TRIM(data_dir) // '/' // TRIM(injector_filenames%w_data)
      END IF
    END IF
#endif

    ! Position (positions that aren't in boundary dimension are mandatory)
    IF (injector%x_data_given) THEN
      INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%x_data),&
          EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(UNIT=inj_base_unit + unit_x, &
            FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%x_data))
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, 'Successfully opened x file: ', TRIM(data_dir) // '/' // &
              TRIM(injector_filenames%x_data)
        END IF
      ELSE
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, '*** ERROR ***'
          PRINT*, 'Unable to locate x file: ', &
              TRIM(data_dir) // '/' // TRIM(injector_filenames%x_data)
        END IF
      END IF
    END IF

    IF (injector%y_data_given) THEN
      INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%y_data),&
          EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(UNIT=inj_base_unit + unit_y, &
            FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%y_data))
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, 'Successfully opened y file: ', TRIM(data_dir) // '/' // &
              TRIM(injector_filenames%y_data)
        END IF
      ELSE
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, '*** ERROR ***'
          PRINT*, 'Unable to locate y file: ', &
              TRIM(data_dir) // '/' // TRIM(injector_filenames%y_data)
        END IF
      END IF
    END IF

    ! Momentum, optional - if missing, will be set to 0
    IF (injector%px_data_given) THEN
      INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%px_data),&
          EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(UNIT=inj_base_unit + unit_px, &
            FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%px_data))
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, 'Successfully opened px file: ', TRIM(data_dir) // '/' // &
              TRIM(injector_filenames%px_data)
        END IF
      ELSE
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, '*** ERROR ***'
          PRINT*, 'Unable to locate px file: ', &
              TRIM(data_dir) // '/' // TRIM(injector_filenames%px_data)
        END IF
      END IF
    END IF

    IF (injector%py_data_given) THEN
      INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%py_data),&
          EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(UNIT=inj_base_unit + unit_py, &
            FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%py_data))
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, 'Successfully opened py file: ', TRIM(data_dir) // '/' // &
              TRIM(injector_filenames%py_data)
        END IF
      ELSE
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, '*** ERROR ***'
          PRINT*, 'Unable to locate py file: ', &
              TRIM(data_dir) // '/' // TRIM(injector_filenames%py_data)
        END IF
      END IF
    END IF

    IF (injector%pz_data_given) THEN
      INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%pz_data),&
          EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(UNIT=inj_base_unit + unit_pz, &
            FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%pz_data))
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, 'Successfully opened pz file: ', TRIM(data_dir) // '/' // &
              TRIM(injector_filenames%pz_data)
        END IF
      ELSE
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, '*** ERROR ***'
          PRINT*, 'Unable to locate pz file: ', &
              TRIM(data_dir) // '/' // TRIM(injector_filenames%pz_data)
        END IF
      END IF
    END IF

#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
    ! Particle ID, optional - if missing, particle ID will be left empty
    IF (injector%id_data_given) THEN
      INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%id_data),&
          EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(UNIT=inj_base_unit + unit_id, &
            FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%id_data))
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, 'Successfully opened id file: ', TRIM(data_dir) // '/' // &
              TRIM(injector_filenames%id_data)
        END IF
      ELSE
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, '*** ERROR ***'
          PRINT*, 'Unable to locate id file: ', &
              TRIM(data_dir) // '/' // TRIM(injector_filenames%id_data)
        END IF
      END IF
    END IF
#endif

  END SUBROUTINE open_injector_files



  SUBROUTINE run_injectors

    TYPE(injector_block), POINTER :: current

    IF (x_min_boundary) THEN
      current => injector_x_min
      DO WHILE(ASSOCIATED(current))
        IF (.NOT. current%inject_from_file) THEN
          CALL run_single_injector(current, c_bd_x_min)
        END IF
        current => current%next
      END DO
    END IF

    IF (x_max_boundary) THEN
      current => injector_x_max
      DO WHILE(ASSOCIATED(current))
        IF (.NOT. current%inject_from_file) THEN
          CALL run_single_injector(current, c_bd_x_max)
        END IF
        current => current%next
      END DO
    END IF

    IF (y_min_boundary) THEN
      current => injector_y_min
      DO WHILE(ASSOCIATED(current))
        IF (.NOT. current%inject_from_file) THEN
          CALL run_single_injector(current, c_bd_y_min)
        END IF
        current => current%next
      END DO
    END IF

    IF (y_max_boundary) THEN
      current => injector_y_max
      DO WHILE(ASSOCIATED(current))
        IF (.NOT. current%inject_from_file) THEN
          CALL run_single_injector(current, c_bd_y_max)
        END IF
        current => current%next
      END DO
    END IF

    ! Injectors exist on all ranks. Above, the IF (..._boundary) lines prevent
    ! non-boundary ranks from generating particles. When injecting from file,
    ! all ranks must know which particles have been added, so different ranks
    ! can pick up where the others left off if they get switched around in the
    ! load balancer. Hence, all processors must read the files (only one rank
    ! adds the particles though)
    current => injector_x_max
    DO WHILE(ASSOCIATED(current))
      IF (current%inject_from_file) THEN
        CALL run_file_injection(current)
      END IF
      current => current%next
    END DO

    current => injector_x_min
    DO WHILE(ASSOCIATED(current))
      IF (current%inject_from_file) THEN
        CALL run_file_injection(current)
      END IF
      current => current%next
    END DO

    current => injector_y_max
    DO WHILE(ASSOCIATED(current))
      IF (current%inject_from_file) THEN
        CALL run_file_injection(current)
      END IF
      current => current%next
    END DO

    current => injector_y_min
    DO WHILE(ASSOCIATED(current))
      IF (current%inject_from_file) THEN
        CALL run_file_injection(current)
      END IF
      current => current%next
    END DO

  END SUBROUTINE run_injectors



  SUBROUTINE run_single_injector(injector, direction)

    TYPE(injector_block), POINTER :: injector
    INTEGER, INTENT(IN) :: direction
    REAL(num) :: bdy_pos, cell_size
    TYPE(particle), POINTER :: new
    TYPE(particle_list) :: plist
    REAL(num) :: mass, typical_mc2, p_therm, p_inject_drift
    REAL(num) :: gamma_mass, v_inject, density, vol, p_drift, p_ratio
    REAL(num) :: npart_ideal, itemp, v_inject_s, density_correction, dir_mult
    REAL(num) :: v_inject_dt
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: weight_fac
#endif
    REAL(num), DIMENSION(3) :: temperature, drift
    INTEGER :: parts_this_time, ipart, idir, dir_index, flux_dir, flux_dir_cell
    INTEGER :: ii
    INTEGER :: perp_dir_index, nperp
    REAL(num) :: perp_cell_size, cur_cell
    TYPE(parameter_pack) :: parameters
    REAL(num), PARAMETER :: sqrt2 = SQRT(2.0_num)
    REAL(num), PARAMETER :: sqrt2_inv = 1.0_num / sqrt2
    REAL(num), PARAMETER :: sqrt2pi_inv = 1.0_num / SQRT(2.0_num * pi)

    IF (time < injector%t_start .OR. time > injector%t_end) RETURN

    ! If you have a moving window that has started moving then unless you
    ! EXPLICITLY give a t_end value to the injector stop the injector
    IF (move_window .AND. window_started .AND. .NOT. injector%has_t_end) &
        RETURN

    IF (direction == c_bd_x_min) THEN
      bdy_pos = x_min
      parameters%pack_ix = 0
      dir_mult = 1.0_num
      ! x-direction
      dir_index = 1
      cell_size = dx
      perp_dir_index = 2
      perp_cell_size = dy
      nperp = ny
    ELSE IF (direction == c_bd_x_max) THEN
      bdy_pos = x_max
      parameters%pack_ix = nx
      dir_mult = -1.0_num
      ! x-direction
      dir_index = 1
      cell_size = dx
      perp_dir_index = 2
      perp_cell_size = dy
      nperp = ny
    ELSE IF (direction == c_bd_y_min) THEN
      bdy_pos = y_min
      parameters%pack_iy = 0
      dir_mult = 1.0_num
      ! y-direction
      dir_index = 2
      cell_size = dy
      perp_dir_index = 1
      perp_cell_size = dx
      nperp = nx
    ELSE IF (direction == c_bd_y_max) THEN
      bdy_pos = y_max
      parameters%pack_iy = ny
      dir_mult = -1.0_num
      ! y-direction
      dir_index = 2
      cell_size = dy
      perp_dir_index = 1
      perp_cell_size = dx
      nperp = nx
    ELSE
      RETURN
    END IF

    IF (injector%use_flux_injector) THEN
      flux_dir = dir_index
    ELSE
      flux_dir = -1
    END IF

    vol = dx * dy
    bdy_pos = bdy_pos - 0.5_num * dir_mult * cell_size * png

    mass = species_list(injector%species)%mass
    typical_mc2 = (mass * c)**2
    cur_cell = 0.0_num
#ifndef PER_SPECIES_WEIGHT
    weight_fac = vol / injector%npart_per_cell
#endif

    CALL create_empty_partlist(plist)

    DO ii = 1, nperp
      IF (perp_dir_index == 1) THEN
        cur_cell = x(ii)
        parameters%pack_ix = ii
      ELSE
        cur_cell = y(ii)
        parameters%pack_iy = ii
      END IF

      parameters%use_grid_position = .TRUE.

      CALL populate_injector_properties(injector, parameters, density=density)

      IF (density < injector%density_min) CYCLE

      CALL populate_injector_properties(injector, parameters, &
          temperature=temperature, drift=drift)

      ! Assume agressive maximum thermal momentum, all components
      ! like hottest component
      p_therm = SQRT(mass * kb * MAXVAL(temperature))
      p_inject_drift = drift(dir_index)
      flux_dir_cell = flux_dir

      IF (flux_dir_cell /= -1) THEN
        ! Drift adjusted so that +ve is 'inwards' through boundary
        p_drift = p_inject_drift * dir_mult

        ! Average momentum of inflowing part
        ! For large inwards drift, is asymptotic to drift
        ! Otherwise it is a complicated expression
        ! Inwards drift - lhs terms are same sign -> +ve
        IF (p_drift > flow_limit_val * p_therm) THEN
          ! For sufficiently large drifts, net inflow -> p_drift
          gamma_mass = SQRT(p_inject_drift**2 + typical_mc2) / c
          v_inject_s = p_inject_drift / gamma_mass
          density_correction = 1.0_num
          ! Large drift flux Maxwellian can be approximated by a
          ! non-flux Maxwellian
          flux_dir_cell = -1
        ELSE IF (p_drift < -flow_limit_val * p_therm) THEN
          ! Net is outflow - inflow velocity is zero
          CYCLE
        ELSE IF (ABS(p_therm) < c_tiny) THEN
          CYCLE
        ELSE IF (ABS(p_drift) < p_therm * 1.0e-9_num) THEN
          v_inject_s = 2.0_num * sqrt2pi_inv * p_therm &
              + (1.0_num - 2.0_num * sqrt2 / pi) * p_drift
          gamma_mass = SQRT(v_inject_s**2 + typical_mc2) / c
          v_inject_s = v_inject_s / gamma_mass
          density_correction = 0.5_num
        ELSE
          p_ratio = sqrt2_inv * p_drift / p_therm

          ! Fraction of the drifting Maxwellian distribution inflowing
          density_correction = 0.5_num * (1.0_num + erf_func(p_ratio))
          IF (density_correction < c_tiny) CYCLE

          ! Below is actually MOMENTUM, will correct on next line
          v_inject_s = dir_mult * (p_drift &
              + sqrt2pi_inv * p_therm * EXP(-p_ratio**2) / density_correction)

          gamma_mass = SQRT(v_inject_s**2 + typical_mc2) / c
          v_inject_s = v_inject_s / gamma_mass
        END IF
      ELSE
        ! User asked for Maxwellian only - no correction to apply
        gamma_mass = SQRT(p_inject_drift**2 + typical_mc2) / c
        v_inject_s = p_inject_drift / gamma_mass
        density_correction = 1.0_num
      END IF

      v_inject = ABS(v_inject_s)
      v_inject_dt = dt * v_inject_s

      npart_ideal = injector%npart_per_cell * v_inject * density_correction &
          * dt / cell_size
      itemp = random_box_muller(0.5_num * SQRT(npart_ideal &
          * (1.0_num - npart_ideal / injector%npart_per_cell))) + npart_ideal
      injector%depth(ii) = injector%depth(ii) - itemp

      IF (injector%depth(ii) >= 0.0_num) CYCLE

      parts_this_time = FLOOR(ABS(injector%depth(ii) - 1.0_num))
      injector%depth(ii) = injector%depth(ii) + REAL(parts_this_time, num)

      DO ipart = 1, parts_this_time
        CALL create_particle(new)

        new%part_pos(perp_dir_index) = &
            (random() - 0.5_num) * perp_cell_size + cur_cell

        new%part_pos(dir_index) = bdy_pos - random() * v_inject_dt
        parameters%pack_pos = new%part_pos
        parameters%use_grid_position = .FALSE.

#ifdef PER_SPECIES_WEIGHT
        CALL populate_injector_properties(injector, parameters, &
            temperature=temperature, drift=drift)
#else
        CALL populate_injector_properties(injector, parameters, density, &
            temperature, drift)
#endif

        DO idir = 1, 3
          IF (idir == flux_dir_cell) THEN
            ! Drift is signed - dir mult is the direciton we want to get
            new%part_p(idir) = flux_momentum_from_temperature(&
                mass, temperature(idir), drift(idir), dir_mult)
          ELSE
            new%part_p(idir) = momentum_from_temperature(mass, &
                temperature(idir), drift(idir))
          END IF
        END DO
#ifdef PER_PARTICLE_CHARGE_MASS
        new%charge = species_list(injector%species)%charge
        new%mass = mass
#endif
#ifndef PER_SPECIES_WEIGHT
        density = MIN(density, injector%density_max)
        new%weight = weight_fac * density
#endif
        CALL add_particle_to_partlist(plist, new)
      END DO
    END DO

    CALL append_partlist(species_list(injector%species)%attached_list, plist)

  END SUBROUTINE run_single_injector



  SUBROUTINE populate_injector_properties(injector, parameters, density, &
      temperature, drift)

    TYPE(injector_block), POINTER :: injector
    TYPE(parameter_pack), INTENT(IN) :: parameters
    REAL(num), INTENT(OUT), OPTIONAL :: density
    REAL(num), DIMENSION(3), INTENT(OUT), OPTIONAL :: temperature, drift
    INTEGER :: errcode, i

    errcode = 0
    IF (PRESENT(density)) THEN
      density = 0.0_num
      IF (injector%density_function%init) THEN
        density = MAX(evaluate_with_parameters(injector%density_function, &
            parameters, errcode), 0.0_num)
      END IF
    END IF

    ! Stack can only be time varying if valid. Change if this isn't true
    IF (PRESENT(temperature)) THEN
      temperature(:) = 0.0_num
      DO i = 1, 3
        IF (injector%temperature_function(i)%init) THEN
          temperature(i) = &
              MAX(evaluate_with_parameters(injector%temperature_function(i), &
                  parameters, errcode), 0.0_num)
        END IF
      END DO
    END IF

    IF (PRESENT(drift)) THEN
      drift(:) = 0.0_num
      DO i = 1, 3
        IF (injector%drift_function(i)%init) THEN
          drift(i) = &
              evaluate_with_parameters(injector%drift_function(i), &
                                       parameters, errcode)
        END IF
      END DO
    END IF

    IF (errcode /= c_err_none) CALL abort_code(errcode)

  END SUBROUTINE populate_injector_properties



  SUBROUTINE finish_injector_setup

    TYPE(injector_block), POINTER :: current

    IF (x_min_boundary) THEN
      current => injector_x_min
      DO WHILE(ASSOCIATED(current))
        CALL finish_single_injector_setup(current, c_bd_x_min)
        current => current%next
      END DO
    END IF

    IF (x_max_boundary) THEN
      current => injector_x_max
      DO WHILE(ASSOCIATED(current))
        CALL finish_single_injector_setup(current, c_bd_x_max)
        current => current%next
      END DO
    END IF

    IF (y_min_boundary) THEN
      current => injector_y_min
      DO WHILE(ASSOCIATED(current))
        CALL finish_single_injector_setup(current, c_bd_y_min)
        current => current%next
      END DO
    END IF

    IF (y_max_boundary) THEN
      current => injector_y_max
      DO WHILE(ASSOCIATED(current))
        CALL finish_single_injector_setup(current, c_bd_y_max)
        current => current%next
      END DO
    END IF

  END SUBROUTINE finish_injector_setup



  SUBROUTINE finish_single_injector_setup(injector, boundary)

    TYPE(injector_block), POINTER :: injector
    INTEGER, INTENT(IN) :: boundary
    TYPE(particle_species), POINTER :: species
    INTEGER :: i

    species => species_list(injector%species)
    IF (injector%npart_per_cell < 0.0_num) THEN
      injector%npart_per_cell = species%npart_per_cell
    END IF

    IF (.NOT.injector%density_function%init) THEN
      CALL copy_stack(species%density_function, injector%density_function)
    END IF

    DO i = 1, 3
      IF (.NOT.injector%drift_function(i)%init) THEN
        CALL copy_stack(species%drift_function(i), injector%drift_function(i))
      END IF
      IF (.NOT.injector%temperature_function(i)%init) THEN
        CALL copy_stack(species%temperature_function(i), &
            injector%temperature_function(i))
      END IF
    END DO

  END SUBROUTINE finish_single_injector_setup



  SUBROUTINE create_boundary_injector(ispecies, bnd)

    INTEGER, INTENT(IN) :: ispecies, bnd
    TYPE(injector_block), POINTER :: working_injector

    species_list(ispecies)%bc_particle(bnd) = c_bc_open
    use_injectors = .TRUE.

    ALLOCATE(working_injector)

    CALL init_injector(bnd, working_injector)
    working_injector%species = ispecies

    CALL attach_injector(working_injector)

  END SUBROUTINE create_boundary_injector



  SUBROUTINE setup_injector_depths(inj_init, depths, inj_count)

    TYPE(injector_block), POINTER :: inj_init
    REAL(num), DIMENSION(:,:), INTENT(IN) :: depths
    INTEGER, INTENT(OUT) :: inj_count
    TYPE(injector_block), POINTER :: inj
    INTEGER :: iinj

    iinj = 1
    inj => inj_init

    DO WHILE(ASSOCIATED(inj))
      ! Exclude ghost cells
      IF (inj%boundary == c_bd_x_min .OR. inj%boundary == c_bd_x_max) THEN
        inj%depth(1:ny) = depths(:,iinj)
      ELSE
        inj%depth(1:nx) = depths(:,iinj)
      END IF
      iinj = iinj + 1
      inj => inj%next
    END DO

    inj_count = iinj - 1

  END SUBROUTINE setup_injector_depths



  SUBROUTINE run_file_injection(injector)

    ! This subroutine reads particles from files (opened after the full injector
    ! block has been read, in a call from deck_injector_block.f90). We read from
    ! the injection time file (t_data) until we hit a particle injected after
    ! the NEXT timestep (time + dt). We position injected particles such that
    ! they pass the injection boundary at the time specified in t_data (assuming
    ! constant velocity), only injecting with ranks at the appropriate position

    TYPE(injector_block), POINTER :: injector
    REAL(num) :: mass, inv_m2c2
    REAL(num) :: x_in, y_in, px_in, py_in, pz_in
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: w_in
#endif
#ifdef PARTICLE_ID4
    INTEGER :: id_in
#elif PARTICLE_ID
    INTEGER(i8) :: id_in
#endif
    INTEGER :: boundary
    REAL(num) :: next_time, time_to_bdy
    REAL(num) :: vx, vy, gamma, inv_gamma_mass, x_start, y_start, z_start
    TYPE(particle), POINTER :: new
    TYPE(particle_list) :: plist
    LOGICAL :: no_particles_added, skip_processor

    ! Set to true if any of the associated files have reached the end of file
    IF (injector%file_finished) RETURN

    mass = species_list(injector%species)%mass
    inv_m2c2 = 1.0_num/(mass*c)**2
    no_particles_added = .TRUE.

    ! Add particles until we reach an injection time greater than the next
    ! timestep (particles must pass the injection boundary in the following
    ! particle push), or until there are no more particles to add
    DO
      ! We always start with injector%next_time known. Global time is a half
      ! timestep ahead of particle time when this is called
      IF (.NOT. injector%next_time < time + 0.5_num*dt ) EXIT

      ! If on x boundary read the y position, and vice versa
      boundary = injector%boundary
      IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
        CALL read_injector_real(unit_y, y_start, injector)
      ELSE IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
        CALL read_injector_real(unit_x, x_start, injector)
      END IF

#ifndef PER_SPECIES_WEIGHT
      ! Read weight data
      CALL read_injector_real(unit_w, w_in, injector)
#endif

      ! Read momentum data
      IF (injector%px_data_given) THEN
        CALL read_injector_real(unit_px, px_in, injector)
      ELSE
        px_in = 0.0_num
      END IF

      IF (injector%py_data_given) THEN
        CALL read_injector_real(unit_py, py_in, injector)
      ELSE
        py_in = 0.0_num
      END IF

      IF (injector%pz_data_given) THEN
        CALL read_injector_real(unit_pz, pz_in, injector)
      ELSE
        pz_in = 0.0_num
      END IF

#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
      IF (injector%id_data_given) THEN
#ifdef PARTICLE_ID4
        CALL read_injector_int4(unit_id, id_in, injector)
#elif PARTICLE_ID
        CALL read_injector_int(unit_id, id_in, injector)
#endif
      ELSE
        id_in = 0
      END IF
#endif

      ! Ensure we still have values for each variable here
      IF (injector%file_finished) EXIT

      ! Identify processors which aren't adding this particle to the simulation
      boundary = injector%boundary
      skip_processor = .FALSE.
      IF (boundary == c_bd_x_min) THEN
        IF (.NOT. x_min_boundary) skip_processor = .TRUE.
      ELSE IF (boundary == c_bd_x_max) THEN
        IF (.NOT. x_max_boundary) skip_processor = .TRUE.
      ELSE IF (boundary == c_bd_y_min) THEN
        IF (.NOT. y_min_boundary) skip_processor = .TRUE.
      ELSE IF (boundary == c_bd_y_max) THEN
        IF (.NOT. y_max_boundary) skip_processor = .TRUE.
      END IF

      ! Skip the following phases if this rank isn't on the right boundary
      IF (skip_processor) THEN
        ! Calculate time for the next particle to be injected
        CALL read_injector_real(unit_t, next_time, injector)
        ! If there are no more particles to add, exit the loop, otherwise save
        ! the next injection time
        IF (injector%file_finished) EXIT
        injector%next_time = next_time
        CYCLE
      END IF

      ! If code has been restarted from an output dump, the position in our file
      ! will be lost. All particles with injection time < simulation time were
      ! injected in the previous run, and shouldn't be added again
      !
      ! Note that particles injected in the step before the restart will have
      ! been deleted, and must be re-added
      IF (injector%next_time < time - 0.5_num*dt) THEN
        CALL read_injector_real(unit_t, next_time, injector)
        IF (.NOT. injector%file_finished) injector%next_time = next_time
        CYCLE
      END IF

      ! Only ranks on the same boundary as the particle can reach here
      ! Calculate particle velocity
      gamma = SQRT(1.0_num + (px_in**2 + py_in**2 + pz_in**2)*inv_m2c2)
      inv_gamma_mass = 1.0_num/(gamma*mass)
      vx = px_in*inv_gamma_mass
      vy = py_in*inv_gamma_mass

      ! Calculate position of injection such that paritlces reach the boundary
      ! at next_time. Note that global time is a half timestep ahead of the time
      ! our particles are at
      time_to_bdy = (injector%next_time - (time-0.5_num*dt))
      IF (boundary == c_bd_x_min) THEN
        x_in = x_min - time_to_bdy * vx
        y_in = y_start - time_to_bdy * vy
      ELSE IF (boundary == c_bd_x_max) THEN
        x_in = x_max - time_to_bdy * vx
        y_in = y_start - time_to_bdy * vy
      ELSE IF (boundary == c_bd_y_min) THEN
        x_in = x_start - time_to_bdy * vx
        y_in = y_min - time_to_bdy * vy
      ELSE IF (boundary == c_bd_y_max) THEN
        x_in = x_start - time_to_bdy * vx
        y_in = y_max - time_to_bdy * vy
      END IF

      ! Now we have start position, find the processor which will be adding the
      ! particle
      IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
        ! Skip all processors which are at the wrong y position
        IF (y_in <= y_min_local .OR. y_in > y_max_local) THEN
          skip_processor = .TRUE.
        END IF

      ELSE IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
        ! Skip all processors which are at the wrong x position
        IF (x_in <= x_min_local .OR. x_in > x_max_local) THEN
          skip_processor = .TRUE.
        END IF
      END IF

      ! Skip the particle adding phase if this rank isn't adding the particle
      IF (skip_processor) THEN
        ! Calculate time for the next particle to be injected
        CALL read_injector_real(unit_t, next_time, injector)
        ! If there are no more particles to add, exit the loop, otherwise save
        ! the next injection time
        IF (injector%file_finished) EXIT
        injector%next_time = next_time
        CYCLE
      END IF

      ! Create the particle and assign properties
      CALL create_particle(new)
      new%part_pos(1) = x_in
      new%part_pos(2) = y_in
      new%part_p(1) = px_in
      new%part_p(2) = py_in
      new%part_p(3) = pz_in
#ifdef PER_PARTICLE_CHARGE_MASS
      new%charge = species_list(injector%species)%charge
      new%mass = mass
#endif
#ifndef PER_SPECIES_WEIGHT
      new%weight = w_in
#endif
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
      new%id = id_in
#endif

      ! Add the particle to our dummy particle list plist
      IF (no_particles_added) THEN
        ! Create an empty particle list the first time a particle is added
        CALL create_empty_partlist(plist)
        no_particles_added = .FALSE.
      END IF
      CALL add_particle_to_partlist(plist, new)

      ! Calculate the next time
      CALL read_injector_real(unit_t, next_time, injector)
      IF (injector%file_finished) THEN
        EXIT
      ELSE
        injector%next_time = next_time
      END IF
    END DO

    ! Append particles to the main particle list for this species on this
    ! processor if we have added particles
    IF (.NOT. no_particles_added) THEN
      CALL append_partlist(species_list(injector%species)%attached_list, plist)
    END IF

  END SUBROUTINE run_file_injection



  SUBROUTINE read_injector_real(unit_code, value, injector)

    ! Opens the file specified by the unit_code, reads a real variable and
    ! changes value to match. Use the unit codes defined in shared_data.F90 -
    ! e.g. unit_t for the t_data file. Also checks for the end of file condition

    TYPE(injector_block), POINTER :: injector
    REAL(num) :: value
    INTEGER :: unit_code, eof

    ! Read file, checking a value has been assigned to our variable, and that we
    ! have not reached the end of the file
    READ(custom_base_unit + (injector%custom_id-1)*custom_var_num + unit_code, &
        *, IOSTAT = eof) value

    ! File read successfully
    IF (eof == 0) RETURN

    ! Triggered if there are no more values to read
    IF (eof < 0) THEN
      injector%file_finished = .TRUE.

    ! Triggered if illegal value entered
    ELSE
      IF(rank == 0) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Illegal value found in read_injector_real'
        PRINT*,'Injected particle will not behave as expected'
        CALL abort_code(c_err_bad_value)
      END IF
    END IF

  END SUBROUTINE read_injector_real



#ifdef PARTICLE_ID
  SUBROUTINE read_injector_int(unit_code, value, injector)

    ! Opens the file specified by the unit_code, reads an integer variable and
    ! changes value to match. Use the unit codes defined in shared_data.F90 -
    ! e.g. unit_t for the t_data file. Also checks for the end of file condition

    TYPE(injector_block), POINTER :: injector
    INTEGER(i8) :: value
    INTEGER :: unit_code, eof

    ! Read file, checking a value has been assigned to our variable, and that we
    ! have not reached the end of the file
    READ(custom_base_unit + (injector%custom_id-1)*custom_var_num + unit_code, &
        *, IOSTAT = eof) value

    ! File read successfully
    IF (eof == 0) RETURN

    ! Triggered if there are no more values to read
    IF (eof < 0) THEN
      injector%file_finished = .TRUE.

    ! Triggered if illegal value entered
    ELSE
      IF(rank == 0) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Illegal value found in read_injector_int'
        PRINT*,'Injected particle will not behave as expected'
        CALL abort_code(c_err_bad_value)
      END IF
    END IF

  END SUBROUTINE read_injector_int
#endif



#ifdef PARTICLE_ID4
  SUBROUTINE read_injector_int4(unit_code, value, injector)

    ! Opens the file specified by the unit_code, reads an integer variable and
    ! changes value to match. Use the unit codes defined in shared_data.F90 -
    ! e.g. unit_t for the t_data file. Also checks for the end of file condition

    TYPE(injector_block), POINTER :: injector
    INTEGER :: value
    INTEGER :: unit_code, eof

    ! Read file, checking a value has been assigned to our variable, and that we
    ! have not reached the end of the file
    READ(custom_base_unit + (injector%custom_id-1)*custom_var_num + unit_code, &
        *, IOSTAT = eof) value

    ! File read successfully
    IF (eof == 0) RETURN

    ! Triggered if there are no more values to read
    IF (eof < 0) THEN
      injector%file_finished = .TRUE.

    ! Triggered if illegal value entered
    ELSE
      IF(rank == 0) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Illegal value found in read_injector_int4'
        PRINT*,'Injected particle will not behave as expected'
        CALL abort_code(c_err_bad_value)
      END IF
    END IF

  END SUBROUTINE read_injector_int4
#endif

END MODULE injectors
