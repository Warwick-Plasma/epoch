MODULE deck_eio_particle_probe_block

  USE probes
  USE shared_data
  USE strings_advanced
#ifndef PARTICLE_PROBES
CONTAINS

  SUBROUTINE Probe_Dummy
  END SUBROUTINE Probe_Dummy
#else

  SAVE
  TYPE(Particle_Probe),POINTER :: Working_Probe

CONTAINS

  SUBROUTINE Probe_Block_Start

    ALLOCATE(Working_Probe)
    CALL Init_Probe(Working_Probe)

  END SUBROUTINE Probe_Block_Start

  SUBROUTINE Probe_Block_End
    !Check whether or not the probe is valid
    REAL(num),DIMENSION(3) :: alpha, beta
    INTEGER :: iCorner, iDirection

    !The probe calculates the signed distance from a point to a plane using Hessian normal form
    alpha = working_probe%Corner(2,:)-working_probe%Corner(1,:)
    beta  = working_probe%Corner(3,:)-working_probe%Corner(1,:)
    !alpha (cross) beta
    working_probe%normal = (/alpha(2)*beta(3) - alpha(3)*beta(2), alpha(3)*beta(1) - alpha(1)*beta(3),&
         alpha(1)*beta(2) - alpha(2)*beta(1)/)
    IF (SUM(ABS(working_probe%normal)) .EQ. 0) THEN
       IF (rank .EQ. 0) PRINT*, "Points specified for the probe plane corners do not allow calculation of a normal to the plane. Probe ",TRIM(working_probe%name)," abandoned."
       DEALLOCATE(Working_probe)
       NULLIFY(working_probe)
       RETURN
    ENDIF
    !Normalise the normal (look, it could be worse OK)
    working_probe%normal=working_probe%normal/SQRT(SUM(working_probe%normal**2))

!!$    DO iCorner=1,4
!!$       DO iDirection=1,3
!!$          IF (working_Probe%Corner(iCorner,iDirection) .LT. working_Probe%Extents(iDirection,1))&
!!$               working_probe%Extents(iDirection,1)=Working_Probe%Corner(iCorner,iDirection)
!!$          IF (working_Probe%Corner(iCorner,iDirection) .GT. working_Probe%Extents(iDirection,2))&
!!$               working_probe%Extents(iDirection,2)=Working_Probe%Corner(iCorner,iDirection)
!!$       ENDDO
!!$    ENDDO

!!$    DO iDirection=1,3
!!$       IF (working_probe%Extents(iDirection,1) .GE. working_probe%Extents(iDirection,2)) THEN
!!$          IF (rank .EQ. 0) PRINT *,"Points specified for the probe plane corners collapse the probe plane to a lower dimension. This is not currently supported. Probe ",TRIM(working_probe%name)," abandoned."
!!$          DEALLOCATE(Working_probe)
!!$          NULLIFY(working_probe)
!!$          RETURN
!!$       ENDIF
!!$    ENDDO

    CALL Attach_Probe(Working_Probe)

  END SUBROUTINE Probe_Block_End

  FUNCTION HandleProbeDeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER :: HandleProbeDeck, iSpecies

    HandleProbeDeck=ERR_NONE

    IF (Element .EQ. blank .OR. Value .EQ. blank) RETURN

    ! get particle probe diagnostics (rolling total of all particles
    ! which pass through a given region of real space (defined by the line between two
    ! points in 2D).
    IF (StrCmp(Element,"dump")) THEN
       working_probe%dump=AsInteger(Value,HandleProbeDeck)
       RETURN
    ENDIF

    !Top left
    IF (StrCmp(Element,"x_tl")) THEN
       working_probe%Corner(1,1)=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"y_tl")) THEN
       working_probe%Corner(1,2)=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"z_tl")) THEN
       working_probe%Corner(1,3)=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF

    !Bottom right
    IF (StrCmp(Element,"x_br")) THEN
       working_probe%Corner(2,1)=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"y_br")) THEN
       working_probe%Corner(2,2)=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"z_br")) THEN
       working_probe%Corner(2,3)=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF

    !Top right
    IF (StrCmp(Element,"x_tr")) THEN
       working_probe%Corner(3,1)=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"y_tr")) THEN
       working_probe%Corner(3,2)=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"z_tr")) THEN
       working_probe%Corner(3,3)=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF

    !Bottom Left
    IF (StrCmp(Element,"x_bl")) THEN
       working_probe%Corner(4,1)=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"y_bl")) THEN
       working_probe%Corner(4,2)=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF
    IF (StrCmp(Element,"z_bl")) THEN
       working_probe%Corner(4,3)=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"probe_species")) THEN
       iSpecies=AsInteger(Value,HandleProbeDeck)
       IF (HandleProbeDeck .EQ. ERR_NONE) THEN
          IF (iSpecies .GT. 0 .AND. iSpecies .LE. nSpecies) THEN
             working_probe%probe_species=>ParticleSpecies(iSpecies)
          ELSE
             IF (rank .EQ. 0) PRINT *,"Unable to attach probe to non existant species ",iSpecies
             HandleProbeDeck=ERR_BAD_VALUE
          ENDIF
       ENDIF
       RETURN
    ENDIF

    IF (StrCmp(Element,"ek_min")) THEN
       working_probe%ek_min=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"ek_max")) THEN
       working_probe%ek_max=AsReal(Value,HandleProbeDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"name")) THEN
       working_probe%name=TRIM(Value)
       RETURN
    ENDIF

    HandleProbeDeck = ERR_UNKNOWN_ELEMENT

  END FUNCTION HandleProbeDeck

#endif

END MODULE deck_eio_particle_probe_block
