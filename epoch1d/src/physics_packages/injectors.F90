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

SAVE

  TYPE(injector_block) :: block

CONTAINS

  SUBROUTINE initialise_injectors

    block%npart_per_cell = 128
    block%species = 1
    block%density = 1.0e23_num
    block%drift = 0.0_num
    block%drift(c_dir_x) = 1.0e-22_num
    block%temperature = 1.0e6_num
    block%next_inject = 0.0_num

  END SUBROUTINE initialise_injectors

  SUBROUTINE run_injector

    TYPE(particle), POINTER :: new
    TYPE(particle_list) :: plist
    REAL(num) :: mass, typical_mc2, p_therm, p_inject_drift
    REAL(num) :: gamma_mass, v_inject, dt_inject, n_eff
    INTEGER :: parts_this_time, ipart, ct, idir

    IF (x_min_boundary) THEN
        IF (time < block%next_inject) RETURN
        mass = species_list(block%species)%mass
        typical_mc2 = (mass * c)**2
        !Assume agressive maximum thermal momentum, all components
        !like hottest component
        p_therm = SQRT(mass * kb * MAXVAL(block%temperature))
        p_inject_drift = block%drift(c_dir_x)
        gamma_mass = SQRT((p_therm + p_inject_drift)**2 + typical_mc2) / c
        v_inject =  p_inject_drift / gamma_mass

        dt_inject = dx/(block%npart_per_cell * v_inject)
        parts_this_time = MAX(FLOOR(dt/dt_inject),1)

        CALL create_empty_partlist(plist)
        DO ipart = 1, parts_this_time
          CALL create_particle(new)
          new%part_pos = x_min - dx*png/2.0_num
          DO idir = 1, 3
            new%part_p(idir) = momentum_from_temperature(mass, &
                block%temperature(idir), block%drift(idir))
#ifdef PER_PARTICLE_CHARGE_MASS
            new%charge = species_list(block%species)%charge
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
          v_inject = v_inject + new%part_p(c_dir_x) / gamma_mass
          ct = ct + 1
          new => new%next
        ENDDO
        v_inject = v_inject / REAL(ct, num)
        n_eff = dx/(dt / parts_this_time * v_inject)
        new => plist%head
        DO WHILE(ASSOCIATED(new))
          new%weight = dx * block%density/REAL(n_eff,num)
          new => new%next
        ENDDO
        CALL append_partlist(species_list(block%species)%attached_list, plist)
        block%next_inject = block%next_inject + dt_inject
    ENDIF

  END SUBROUTINE run_injector

END MODULE injectors
