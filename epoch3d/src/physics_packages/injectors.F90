! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2017 Chris Brady <C.S.Brady@warwick.ac.uk>
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

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_injector(boundary, injector)

    INTEGER, INTENT(IN) :: boundary
    TYPE(injector_block), INTENT(INOUT) :: injector

    injector%npart_per_cell = 0
    injector%species = -1
    injector%boundary = boundary
    injector%t_start = 0.0_num
    injector%t_end = t_end
    injector%has_t_end = .FALSE.
    injector%density_min = 0.0_num
    injector%density_max = HUGE(1.0_num)
    injector%use_flux_injector = .FALSE.

    IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
      ALLOCATE(injector%dt_inject(1-ng:ny+ng, 1-ng:nz+ng))
      ALLOCATE(injector%depth(1-ng:ny+ng, 1-ng:nz+ng))
    END IF

    IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
      ALLOCATE(injector%dt_inject(1-ng:nx+ng, 1-ng:nz+ng))
      ALLOCATE(injector%depth(1-ng:nx+ng, 1-ng:nz+ng))
    END IF

    IF (boundary == c_bd_z_min .OR. boundary == c_bd_z_max) THEN
      ALLOCATE(injector%dt_inject(1-ng:nx+ng, 1-ng:ny+ng))
      ALLOCATE(injector%depth(1-ng:nx+ng, 1-ng:ny+ng))
    END IF

    injector%depth = 1.0_num
    injector%dt_inject = -1.0_num
    NULLIFY(injector%next)

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
    ELSE IF (boundary == c_bd_z_min) THEN
      CALL attach_injector_to_list(injector_z_min, injector)
    ELSE IF (boundary == c_bd_z_max) THEN
      CALL attach_injector_to_list(injector_z_max, injector)
    END IF

  END SUBROUTINE attach_injector



  ! Actually does the attaching of the injector to the correct list
  SUBROUTINE attach_injector_to_list(list, injector)

    TYPE(injector_block), POINTER :: list
    TYPE(injector_block), POINTER :: injector
    TYPE(injector_block), POINTER :: current

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
    CALL deallocate_injector_list(injector_z_min)
    CALL deallocate_injector_list(injector_z_max)

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
      IF (ASSOCIATED(current%dt_inject)) DEALLOCATE(current%dt_inject)
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



  SUBROUTINE run_injectors

    TYPE(injector_block), POINTER :: current

    IF (x_min_boundary) THEN
      current => injector_x_min
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_x_min)
        current => current%next
      END DO
    END IF

    IF (x_max_boundary) THEN
      current => injector_x_max
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_x_max)
        current => current%next
      END DO
    END IF

    IF (y_min_boundary) THEN
      current => injector_y_min
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_y_min)
        current => current%next
      END DO
    END IF

    IF (y_max_boundary) THEN
      current => injector_y_max
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_y_max)
        current => current%next
      END DO
    END IF

    IF (z_min_boundary) THEN
      current => injector_z_min
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_z_min)
        current => current%next
      END DO
    END IF

    IF (z_max_boundary) THEN
      current => injector_z_max
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_z_max)
        current => current%next
      END DO
    END IF

  END SUBROUTINE run_injectors



  SUBROUTINE run_single_injector(injector, direction)

    TYPE(injector_block), POINTER :: injector
    INTEGER, INTENT(IN) :: direction
    REAL(num) :: bdy_pos, bdy_space
    TYPE(particle), POINTER :: new
    TYPE(particle_list) :: plist
    REAL(num) :: mass, typical_mc2, p_therm, p_inject_drift, density_grid
    REAL(num) :: gamma_mass, v_inject, density, vol
    REAL(num) :: npart_ideal, itemp, v_inject_s
    REAL(num), DIMENSION(3) :: temperature, drift
    INTEGER :: parts_this_time, ipart, idir, dir_index, ii, jj
    INTEGER, DIMENSION(c_ndims-1) :: perp_dir_index, nel
    REAL(num), DIMENSION(c_ndims-1) :: perp_cell_size, cur_cell
    TYPE(parameter_pack) :: parameters
    INTEGER, DIMENSION(2) :: i2d
    REAL(num), DIMENSION(3) :: dir_mult
    LOGICAL :: first_inject, flux_fn

    IF (time < injector%t_start .OR. time > injector%t_end) RETURN

    ! If you have a moving window that has started moving then unless you
    ! EXPLICITLY give a t_end value to the injector stop the injector
    IF (move_window .AND. window_started .AND. .NOT. injector%has_t_end) &
        RETURN

    flux_fn = .FALSE.
    dir_mult = 1.0_num

    IF (direction == c_bd_x_min) THEN
      parameters%pack_ix = 0
      nel = (/ny, nz/)
      perp_cell_size = (/dy, dz/)
      perp_dir_index = (/2, 3/)
      dir_index = 1
      bdy_pos = x_min
      bdy_space = -dx
      IF (injector%use_flux_injector) THEN
        flux_fn = .TRUE.
        dir_mult(dir_index) = 1.0_num
      END IF
    ELSE IF (direction == c_bd_x_max) THEN
      parameters%pack_ix = nx
      nel = (/ny, nz/)
      perp_cell_size = (/dy, dz/)
      perp_dir_index = (/2, 3/)
      dir_index = 1
      bdy_pos = x_max
      bdy_space = dx
      IF (injector%use_flux_injector) THEN
        flux_fn = .TRUE.
        dir_mult(dir_index) = -1.0_num
      END IF
    ELSE IF (direction == c_bd_y_min) THEN
      parameters%pack_iy = 0
      nel = (/nx, nz/)
      perp_cell_size = (/dx, dz/)
      perp_dir_index = (/1, 3/)
      dir_index = 2
      bdy_pos = y_min
      bdy_space = -dy
      IF (injector%use_flux_injector) THEN
        flux_fn = .TRUE.
        dir_mult(dir_index) = 1.0_num
      END IF
    ELSE IF (direction == c_bd_y_max) THEN
      parameters%pack_iy = ny
      nel = (/nx, nz/)
      perp_cell_size = (/dx, dz/)
      perp_dir_index = (/1, 3/)
      dir_index = 2
      bdy_pos = y_max
      bdy_space = dy
      IF (injector%use_flux_injector) THEN
        flux_fn = .TRUE.
        dir_mult(dir_index) = -1.0_num
      END IF
    ELSE IF (direction == c_bd_z_min) THEN
      parameters%pack_iz = 0
      nel = (/nx, ny/)
      perp_cell_size = (/dx, dy/)
      perp_dir_index = (/1, 2/)
      dir_index = 3
      bdy_pos = z_min
      bdy_space = -dz
      IF (injector%use_flux_injector) THEN
        flux_fn = .TRUE.
        dir_mult(dir_index) = 1.0_num
      END IF
    ELSE IF (direction == c_bd_z_max) THEN
      parameters%pack_iz = nz
      nel = (/nx, ny/)
      perp_cell_size = (/dx, dy/)
      perp_dir_index = (/1, 2/)
      dir_index = 3
      bdy_pos = z_max
      bdy_space = dz
      IF (injector%use_flux_injector) THEN
        flux_fn = .TRUE.
        dir_mult(dir_index) = -1.0_num
      END IF
    ELSE
      RETURN
    END IF

    vol = ABS(bdy_space)
    DO idir = 1, c_ndims-1
      vol = vol * perp_cell_size(idir)
    END DO

    mass = species_list(injector%species)%mass
    typical_mc2 = (mass * c)**2
    cur_cell = 0.0_num

    CALL create_empty_partlist(plist)
    DO ii = 1, nel(1)
      DO jj = 1, nel(2)
        i2d = (/ii, jj/)
        DO idir = 1, c_ndims-1
          IF (perp_dir_index(idir) == 1) cur_cell(idir) = x(i2d(idir))
          IF (perp_dir_index(idir) == 2) cur_cell(idir) = y(i2d(idir))
          IF (perp_dir_index(idir) == 3) cur_cell(idir) = z(i2d(idir))
        END DO

        parameters%use_grid_position = .TRUE.
        CALL assign_pack_value(parameters, perp_dir_index(1), ii)
        CALL assign_pack_value(parameters, perp_dir_index(2), jj)

        IF (injector%dt_inject(ii,jj) > 0.0_num) THEN
          npart_ideal = dt / injector%dt_inject(ii,jj)
          itemp = random_box_muller(0.5_num * SQRT(npart_ideal &
              * (1.0_num - npart_ideal / REAL(injector%npart_per_cell, num)))) &
              + npart_ideal
          injector%depth(ii,jj) = injector%depth(ii,jj) - itemp
          first_inject = .FALSE.

          IF (injector%depth(ii,jj) >= 0.0_num) CYCLE
        ELSE
          first_inject = .TRUE.
        END IF

        CALL populate_injector_properties(injector, parameters, density_grid, &
            temperature, drift)

        IF (density_grid < injector%density_min) CYCLE

        ! Assume agressive maximum thermal momentum, all components
        ! like hottest component
        p_therm = SQRT(mass * kb * MAXVAL(temperature))
        p_inject_drift = drift(dir_index)
        gamma_mass = SQRT((p_therm + p_inject_drift)**2 + typical_mc2) / c
        v_inject_s = p_inject_drift / gamma_mass
        v_inject = ABS(v_inject_s)

        injector%dt_inject(ii,jj) = ABS(bdy_space) &
            / MAX(injector%npart_per_cell * v_inject, c_tiny)
        IF (first_inject) THEN
          ! On the first run of the injectors it isn't possible to decrement
          ! the optical depth until this point
          npart_ideal = dt / injector%dt_inject(ii,jj)
          itemp = random_box_muller(0.5_num * SQRT(npart_ideal &
              * (1.0_num - npart_ideal / REAL(injector%npart_per_cell, num)))) &
              + npart_ideal
          injector%depth(ii,jj) = injector%depth(ii,jj) - itemp
        END IF

        parts_this_time = FLOOR(ABS(injector%depth(ii,jj) - 1.0_num))
        injector%depth(ii,jj) = injector%depth(ii,jj) &
            + REAL(parts_this_time, num)

        DO ipart = 1, parts_this_time
          CALL create_particle(new)

          new%part_pos = 0.0_num
          DO idir = 1, c_ndims-1
            new%part_pos(perp_dir_index(idir)) = &
                (random() - 0.5_num) * perp_cell_size(idir) + cur_cell(idir)
          END DO

          new%part_pos(dir_index) = bdy_pos + 0.5_num * bdy_space * png &
              - random() * v_inject_s * dt
          parameters%pack_pos = new%part_pos
          parameters%use_grid_position = .FALSE.

          CALL populate_injector_properties(injector, parameters, density, &
              temperature, drift)

          DO idir = 1, 3
            IF (flux_fn) THEN
              new%part_p(idir) = flux_momentum_from_temperature(mass, &
                  temperature(idir), drift(idir)) * dir_mult(idir)
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
          new%weight = vol * density / REAL(injector%npart_per_cell, num)
#endif
          CALL add_particle_to_partlist(plist, new)
        END DO
      END DO
    END DO

    CALL append_partlist(species_list(injector%species)%attached_list, plist)

  END SUBROUTINE run_single_injector



  SUBROUTINE populate_injector_properties(injector, parameters, density, &
      temperature, drift)

    TYPE(injector_block), POINTER :: injector
    TYPE(parameter_pack), INTENT(IN) :: parameters
    REAL(num), INTENT(OUT) :: density
    REAL(num), DIMENSION(3), INTENT(OUT) :: temperature, drift
    INTEGER :: errcode, i

    errcode = 0
    density = MAX(evaluate_with_parameters(injector%density_function, &
        parameters, errcode), 0.0_num)

    ! Stack can only be time varying if valid. Change if this isn't true
    DO i = 1, 3
      IF (injector%temperature_function(i)%init) THEN
        temperature(i) = &
            MAX(evaluate_with_parameters(injector%temperature_function(i), &
                parameters, errcode), 0.0_num)
      ELSE
        temperature(i) = 0.0_num
      END IF
      IF (injector%drift_function(i)%init) THEN
        drift(i) = &
            evaluate_with_parameters(injector%drift_function(i), &
                                     parameters, errcode)
      ELSE
        drift(i) = 0.0_num
      END IF
    END DO

    IF (errcode /= c_err_none) CALL abort_code(errcode)

  END SUBROUTINE populate_injector_properties



  SUBROUTINE assign_pack_value(parameters, dir_index, p_value)

    TYPE(parameter_pack), INTENT(INOUT) :: parameters
    INTEGER, INTENT(IN) :: dir_index
    INTEGER, INTENT(IN) :: p_value

    IF (dir_index == 1) THEN
      parameters%pack_ix = p_value
    ELSE IF (dir_index == 2) THEN
      parameters%pack_iy = p_value
    ELSE IF (dir_index == 3) THEN
      parameters%pack_iz = p_value
    END IF

  END SUBROUTINE assign_pack_value

END MODULE injectors
