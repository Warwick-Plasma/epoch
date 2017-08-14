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

  USE shared_data
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
    IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
      ALLOCATE(injector%dt_inject(1:ny))
      ALLOCATE(injector%depth(1:ny))
    ENDIF
    IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
      ALLOCATE(injector%dt_inject(1:nx))
      ALLOCATE(injector%depth(1:nx))
    ENDIF
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
    ENDIF

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
      ENDDO
      current%next => injector
    ELSE
      list => injector
    ENDIF

  END SUBROUTINE attach_injector_to_list



  SUBROUTINE run_injectors

    TYPE(injector_block), POINTER :: current

    IF (x_min_boundary) THEN
      current => injector_x_min
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_x_min, x_min, dx)
        current => current%next
      ENDDO
    ENDIF

    IF (x_max_boundary) THEN
      current => injector_x_max
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_x_max, x_max, dx)
        current => current%next
      ENDDO
    ENDIF

    IF (y_min_boundary) THEN
      current => injector_y_min
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_y_min, y_min, dy)
        current => current%next
      ENDDO
    ENDIF

    IF (y_max_boundary) THEN
      current => injector_y_max
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_bd_y_max, y_max, dy)
        current => current%next
      ENDDO
    ENDIF

  END SUBROUTINE run_injectors



  SUBROUTINE run_single_injector(injector, direction, bdy_pos, bdy_space)

    TYPE(injector_block), POINTER :: injector
    INTEGER, INTENT(IN) :: direction
    REAL(num), INTENT(IN) :: bdy_pos, bdy_space
    TYPE(particle), POINTER :: new
    TYPE(particle_list) :: plist
    REAL(num) :: mass, typical_mc2, p_therm, p_inject_drift
    REAL(num) :: gamma_mass, v_inject, n_eff, density, npart_par,vol
    REAL(num), DIMENSION(3) :: temperature, drift
    INTEGER :: parts_this_time, ipart, ct, idir, dir_index, ii, ipart_trans
    INTEGER, DIMENSION(c_ndims-1) :: perp_dir_index, nel, npart_perp
    REAL(num), DIMENSION(c_ndims-1) :: perp_cell_size, cur_cell
    TYPE(parameter_pack) :: parameters

    IF (time < injector%t_start .OR. time > injector%t_end) RETURN


    IF (direction == c_bd_x_min) THEN
      parameters%pack_ix = 0
      nel = (/ny/)
      perp_cell_size = (/dy/)
      perp_dir_index = (/2/)
      dir_index = 1
    ENDIF
    IF (direction == c_bd_x_max) THEN
      parameters%pack_ix = nx
      nel = (/ny/)
      perp_cell_size = (/dy/)
      perp_dir_index = (/2/)
      dir_index = 1
    ENDIF
    IF (direction == c_bd_y_min) THEN
      parameters%pack_iy = 0
      nel = (/nx/)
      perp_cell_size = (/dx/)
      perp_dir_index = (/1/)
      dir_index = 2
    ENDIF
    IF (direction == c_bd_y_max) THEN
      parameters%pack_ix = ny
      nel = (/nx/)
      perp_cell_size = (/dx/)
      perp_dir_index = (/1/)
      dir_index = 2
    ENDIF

    vol = bdy_space
    DO idir = 1, c_ndims-1
      vol = vol * perp_cell_size(idir)
    ENDDO

    mass = species_list(injector%species)%mass
    typical_mc2 = (mass * c)**2
    cur_cell = 0.0_num

    DO ii = 1, nel(1)
      CALL create_empty_partlist(plist)
      DO idir = 1, c_ndims-1
        IF (perp_dir_index(idir) == 1) cur_cell(idir) = x(ii)
        IF (perp_dir_index(idir) == 2) cur_cell(idir) = y(ii)
      ENDDO


      parameters%use_grid_position = .TRUE.
      CALL assign_pack_value(parameters, perp_dir_index(1), ii)

      IF (injector%dt_inject(ii) > 0.0_num) THEN
        injector%depth(ii) = injector%depth(ii) - random() * 2.0_num &
            * dt / injector%dt_inject(ii)
        IF (injector%depth(ii) >= 0.0_num) CYCLE
      ENDIF

      CALL populate_injector_properties(injector, parameters, density, &
          temperature, drift)

      !Fit particles
      !npart_par is REAL because you can meaningfully inject at non-integer
      !rate
      npart_par = MAX(SQRT(REAL(injector%npart_per_cell, num)),1.0_num)
      !npart_perp is integer because you are actually placing particles
      npart_perp = CEILING(npart_par)

      !Assume agressive maximum thermal momentum, all components
      !like hottest component
      p_therm = SQRT(mass * kb * MAXVAL(temperature))
      p_inject_drift = drift(dir_index)
      gamma_mass = SQRT((p_therm + p_inject_drift)**2 + typical_mc2) / c
      v_inject =  ABS(p_inject_drift / gamma_mass)

      injector%dt_inject(ii) = bdy_space/MAX(npart_par * v_inject,c_tiny)
      parts_this_time = FLOOR(ABS(injector%depth(ii) - 1.0_num))
      injector%depth(ii) = 1.0_num

      DO ipart = 1, parts_this_time
        DO ipart_trans = 1, npart_perp(1)
          CALL create_particle(new)
          new%part_pos = 0.0_num
          DO idir = 1, c_ndims-1
            new%part_pos(perp_dir_index(idir)) = (random() - 0.5_num)  &
                * perp_cell_size(idir) &
                + cur_cell(idir)
          ENDDO
          new%part_pos(dir_index) = bdy_pos - bdy_space*png/2.0_num
          parameters%pack_pos = new%part_pos
          parameters%use_grid_position = .FALSE.
          CALL populate_injector_properties(injector, parameters, density, &
            temperature, drift)
          DO idir = 1, 3
            new%part_p(idir) = momentum_from_temperature(mass, &
                temperature(idir), drift(idir))
          ENDDO
#ifdef PER_PARTICLE_CHARGE_MASS
          new%charge = species_list(injector%species)%charge
          new%mass = mass
#endif
          CALL add_particle_to_partlist(plist, new)
        ENDDO
      ENDDO
      new => plist%head
      v_inject = 0.0_num
      ct=0
      DO WHILE(ASSOCIATED(new))
        gamma_mass = SQRT(SUM(new%part_p**2) + typical_mc2) / c
        v_inject = v_inject + new%part_p(dir_index) / gamma_mass
        ct = ct + 1
        new => new%next
      ENDDO
      v_inject = ABS(v_inject) / MAX(REAL(ct, num), c_tiny)
      n_eff = bdy_space * parts_this_time /MAX(dt * v_inject, c_tiny)
      new => plist%head
      DO WHILE(ASSOCIATED(new))
        new%weight = vol * density/REAL(n_eff,num)
        new => new%next
      ENDDO
      CALL append_partlist(species_list(injector%species)%attached_list, plist)
    ENDDO

  END SUBROUTINE run_single_injector



  SUBROUTINE populate_injector_properties(injector, parameters, density, &
      temperature, drift)

    TYPE(injector_block), POINTER :: injector
    TYPE(parameter_pack), INTENT(IN) :: parameters
    REAL(num), INTENT(OUT) :: density
    REAL(num), DIMENSION(3), INTENT(OUT) :: temperature, drift
    INTEGER :: errcode, i

    errcode = 0
    density = evaluate_with_parameters(injector%density_function, &
        parameters, errcode)

    !Stack can only be time varying if valid. Change if this isn't true
    DO i = 1, 3
      IF (injector%temperature_function(i)%init) THEN
        temperature(i) = &
            evaluate_with_parameters(injector%temperature_function(i), &
            parameters, errcode)
      ELSE
        temperature(i) = 0.0_num
      ENDIF
      IF (injector%drift_function(i)%init) THEN
        drift(i) = &
            evaluate_with_parameters(injector%drift_function(i), &
            parameters, errcode)
      ELSE
        drift(i) = 0.0_num
      ENDIF
    ENDDO

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
    ENDIF

  END SUBROUTINE assign_pack_value

END MODULE injectors
