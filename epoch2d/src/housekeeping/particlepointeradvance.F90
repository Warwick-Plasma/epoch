MODULE ParticlePointerAdvance
  USE shared_data

  IMPLICIT NONE

SAVE
TYPE(ParticleFamily), POINTER, SAVE :: CurrentHead

CONTAINS

  SUBROUTINE Start_ParticleFamily(PartFamily,PartList,Part)

    TYPE(ParticleFamily), POINTER, INTENT(INOUT) :: PartFamily
    TYPE(ParticleList), POINTER, INTENT(INOUT) :: PartList
    TYPE(Particle), POINTER, INTENT(INOUT) :: Part

    PartFamily=>ParticleSpecies(1)

    IF (ASSOCIATED(PartFamily)) THEN
       PartList=>PartFamily%AttachedList
       IF (ASSOCIATED(PartList)) THEN
          Part=>PartList%Head
       ELSE
          NULLIFY(Part)
       ENDIF
    ELSE
       NULLIFY(PartList)
       NULLIFY(Part)
    ENDIF

  END SUBROUTINE Start_ParticleFamily

  SUBROUTINE Advance_ParticleList(PartList,Part)

    TYPE(ParticleList), POINTER, INTENT(INOUT) :: PartList
    TYPE(Particle), POINTER, INTENT(INOUT) :: Part

    PartList=>PartList%Next
    IF (ASSOCIATED(PartList)) THEN
       Part=>PartList%Head
    ELSE
       NULLIFY(Part)
    ENDIF


  END SUBROUTINE Advance_ParticleList

  SUBROUTINE Advance_ParticleFamily(PartFamily,PartList,Part)

    TYPE(ParticleFamily), POINTER, INTENT(INOUT) :: PartFamily
    TYPE(ParticleList), POINTER, INTENT(INOUT) :: PartList
    TYPE(Particle), POINTER, INTENT(INOUT) :: Part

    PartFamily=>PartFamily%Next
    IF (ASSOCIATED(PartFamily)) THEN
       PartList=>PartFamily%AttachedList
    ELSE
       NULLIFY(PartList)
    ENDIF

    IF (ASSOCIATED(PartList)) THEN
       Part=>PartList%Head
    ELSE
       NULLIFY(Part)
    ENDIF


  END SUBROUTINE Advance_ParticleFamily

END MODULE ParticlePointerAdvance
