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

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_injector(boundary, injector)

    INTEGER, INTENT(IN) :: boundary
    TYPE(injector_block), INTENT(INOUT) :: injector

    injector%npart_per_cell = 0
    injector%species = -1
    injector%density = 0.0_num
    injector%drift = 0.0_num
    injector%temperature = 0.0_num
    injector%next_inject = 0.0_num
    injector%boundary = boundary
    injector%t_start = 0.0_num
    injector%t_end = t_end
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
        CALL run_single_injector(current, c_dir_x, x_min, dx)
        current => current%next
      ENDDO
    ENDIF

    IF (x_max_boundary) THEN
      current => injector_x_max
      DO WHILE(ASSOCIATED(current))
        CALL run_single_injector(current, c_dir_x, x_max, dx)
        current => current%next
      ENDDO
    ENDIF

  END SUBROUTINE run_injectors



  SUBROUTINE run_single_injector(injector, direction, bdy_pos, bdy_space)

    TYPE(injector_block), INTENT(INOUT) :: injector
    INTEGER, INTENT(IN) :: direction
    REAL(num), INTENT(IN) :: bdy_pos, bdy_space
    TYPE(particle), POINTER :: new
    TYPE(particle_list) :: plist
    REAL(num) :: mass, typical_mc2, p_therm, p_inject_drift
    REAL(num) :: gamma_mass, v_inject, dt_inject, n_eff
    INTEGER :: parts_this_time, ipart, ct, idir

    IF (time < injector%next_inject) RETURN
    IF (time < injector%t_start .OR. time > injector%t_end) RETURN
    mass = species_list(injector%species)%mass
    typical_mc2 = (mass * c)**2
    !Assume agressive maximum thermal momentum, all components
    !like hottest component
    p_therm = SQRT(mass * kb * MAXVAL(injector%temperature))
    p_inject_drift = injector%drift(direction)
    gamma_mass = SQRT((p_therm + p_inject_drift)**2 + typical_mc2) / c
    v_inject =  ABS(p_inject_drift / gamma_mass)

    dt_inject = bdy_space/(injector%npart_per_cell * v_inject)
    parts_this_time = MAX(FLOOR(dt/dt_inject),1)

    CALL create_empty_partlist(plist)
    DO ipart = 1, parts_this_time
      CALL create_particle(new)
      new%part_pos = bdy_pos - bdy_space*png/2.0_num
      DO idir = 1, 3
        new%part_p(idir) = momentum_from_temperature(mass, &
            injector%temperature(idir), injector%drift(idir))
#ifdef PER_PARTICLE_CHARGE_MASS
        new%charge = species_list(injector%species)%charge
        new%mass = mass
#endif
      ENDDO
      CALL add_particle_to_partlist(plist, new)
    ENDDO
    new => plist%head
    v_inject = 0.0_num
    ct=0
    DO WHILE(ASSOCIATED(new))
      gamma_mass = SQRT(SUM(new%part_p**2) + typical_mc2) / c
      v_inject = v_inject + new%part_p(direction) / gamma_mass
      ct = ct + 1
      new => new%next
    ENDDO
    v_inject = ABS(v_inject) / REAL(ct, num)
    n_eff = bdy_space/(dt / parts_this_time * v_inject)
    new => plist%head
    DO WHILE(ASSOCIATED(new))
      new%weight = bdy_space * injector%density/REAL(n_eff,num)
      new => new%next
    ENDDO
    CALL append_partlist(species_list(injector%species)%attached_list, plist)
    injector%next_inject = time + dt_inject

  END SUBROUTINE run_single_injector

END MODULE injectors
