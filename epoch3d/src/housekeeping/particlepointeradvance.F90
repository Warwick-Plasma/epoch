MODULE particle_pointer_advance

  USE shared_data

  IMPLICIT NONE

  SAVE
  TYPE(particle_family), POINTER :: current_head

CONTAINS

  SUBROUTINE start_particle_family(part_family, partlist, part)

    TYPE(particle_family), POINTER, INTENT(INOUT) :: part_family
    TYPE(particle_list), POINTER, INTENT(INOUT) :: partlist
    TYPE(particle), POINTER, INTENT(INOUT) :: part

    part_family=>particle_species(1)

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

  END SUBROUTINE start_particle_family



  SUBROUTINE advance_particle_list(partlist, part)

    TYPE(particle_list), POINTER, INTENT(INOUT) :: partlist
    TYPE(particle), POINTER, INTENT(INOUT) :: part

    partlist=>partlist%next
    IF (ASSOCIATED(partlist)) THEN
      part=>partlist%head
    ELSE
      NULLIFY(part)
    ENDIF

  END SUBROUTINE advance_particle_list



  SUBROUTINE advance_particle_family(part_family, partlist, part)

    TYPE(particle_family), POINTER, INTENT(INOUT) :: part_family
    TYPE(particle_list), POINTER, INTENT(INOUT) :: partlist
    TYPE(particle), POINTER, INTENT(INOUT) :: part

    part_family=>part_family%next
    IF (ASSOCIATED(part_family)) THEN
      partlist=>part_family%attached_list
    ELSE
      NULLIFY(partlist)
    ENDIF

    IF (ASSOCIATED(partlist)) THEN
      part=>partlist%head
    ELSE
      NULLIFY(part)
    ENDIF

  END SUBROUTINE advance_particle_family

END MODULE particle_pointer_advance
