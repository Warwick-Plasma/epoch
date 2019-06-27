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

MODULE particle_pointer_advance

  USE shared_data

  IMPLICIT NONE

  SAVE
  TYPE(particle_species), POINTER :: current_head

CONTAINS

  SUBROUTINE start_particle_species(part_species, partlist, part)

    TYPE(particle_species), POINTER :: part_species
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: part

    part_species => species_list(1)

    CALL start_particle_list(part_species, partlist, part)

  END SUBROUTINE start_particle_species



  SUBROUTINE start_particle_species_only(part_species)

    TYPE(particle_species), POINTER :: part_species

    part_species => species_list(1)

  END SUBROUTINE start_particle_species_only



  SUBROUTINE start_particle_list(part_species, partlist, part)

    TYPE(particle_species), POINTER :: part_species
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: part

    IF (ASSOCIATED(part_species)) THEN
      partlist => part_species%attached_list
      IF (ASSOCIATED(partlist)) THEN
        part => partlist%head
      ELSE
        NULLIFY(part)
      END IF
    ELSE
      NULLIFY(partlist)
      NULLIFY(part)
    END IF

  END SUBROUTINE start_particle_list



  SUBROUTINE advance_particle_list(partlist, part)

    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: part

    partlist => partlist%next
    IF (ASSOCIATED(partlist)) THEN
      part => partlist%head
    ELSE
      NULLIFY(part)
    END IF

  END SUBROUTINE advance_particle_list



  SUBROUTINE advance_particle_species(part_species, partlist, part)

    TYPE(particle_species), POINTER :: part_species
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: part

    part_species => part_species%next
    CALL start_particle_list(part_species, partlist, part)

  END SUBROUTINE advance_particle_species



  SUBROUTINE advance_particle_species_only(part_species)

    TYPE(particle_species), POINTER :: part_species

    part_species => part_species%next

  END SUBROUTINE advance_particle_species_only

END MODULE particle_pointer_advance
