MODULE particle_pointer_advance

  USE shared_data

  IMPLICIT NONE

  SAVE
  TYPE(particle_family), POINTER :: current_head

CONTAINS

  SUBROUTINE start_particle_family(part_family, partlist, part)

    TYPE(particle_family), POINTER :: part_family
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: part

    part_family=>particle_species(1)

    CALL start_particle_list(part_family, partlist, part)

  END SUBROUTINE start_particle_family



  SUBROUTINE start_particle_family_only(part_family)

    TYPE(particle_family), POINTER :: part_family

    part_family=>particle_species(1)

  END SUBROUTINE start_particle_family_only



  SUBROUTINE start_particle_list(part_family, partlist, part)

    TYPE(particle_family), POINTER :: part_family
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: part

    IF (ASSOCIATED(part_family)) THEN
      partlist=>part_family%attached_list
      IF (ASSOCIATED(partlist)) THEN
        part=>partlist%head
      ELSE
        NULLIFY(part)
      ENDIF
    ELSE
      NULLIFY(partlist)
      NULLIFY(part)
    ENDIF

  END SUBROUTINE start_particle_list



  SUBROUTINE advance_particle_list(partlist, part)

    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: part

    partlist=>partlist%next
    IF (ASSOCIATED(partlist)) THEN
      part=>partlist%head
    ELSE
      NULLIFY(part)
    ENDIF

  END SUBROUTINE advance_particle_list



  SUBROUTINE advance_particle_family(part_family, partlist, part)

    TYPE(particle_family), POINTER :: part_family
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: part

    part_family=>part_family%next
    CALL start_particle_list(part_family, partlist, part)

  END SUBROUTINE advance_particle_family



  SUBROUTINE advance_particle_family_only(part_family)

    TYPE(particle_family), POINTER :: part_family

    part_family=>part_family%next

  END SUBROUTINE advance_particle_family_only

END MODULE particle_pointer_advance
