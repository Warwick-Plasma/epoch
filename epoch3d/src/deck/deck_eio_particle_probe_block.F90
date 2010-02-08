MODULE deck_eio_particle_probe_block

  USE probes
  USE shared_data
  USE strings_advanced
#ifndef PARTICLE_PROBES
CONTAINS

  SUBROUTINE probe_deck_dummy
  END SUBROUTINE probe_deck_dummy
#else

  SAVE
  TYPE(particle_probe),POINTER :: working_probe

CONTAINS

  SUBROUTINE probe_block_start

    ALLOCATE(working_probe)
    CALL init_probe(working_probe)

  END SUBROUTINE probe_block_start

  SUBROUTINE probe_block_end
    !Check whether or not the probe is valid
    REAL(num),DIMENSION(3) :: alpha, beta

    !The probe calculates the signed distance from a point to a plane using Hessian normal form
    alpha = working_probe%corner(2,:)-working_probe%corner(1,:)
    beta  = working_probe%corner(3,:)-working_probe%corner(1,:)
    !alpha (cross) beta
    working_probe%normal = (/alpha(2)*beta(3) - alpha(3)*beta(2), alpha(3)*beta(1) - alpha(1)*beta(3),&
         alpha(1)*beta(2) - alpha(2)*beta(1)/)
    IF (SUM(ABS(working_probe%normal)) .EQ. 0) THEN
       IF (rank .EQ. 0) PRINT*, "Points specified for the probe plane corners do not allow calculation of a normal to the plane. Probe ",TRIM(working_probe%name)," abandoned."
       DEALLOCATE(working_probe)
       NULLIFY(working_probe)
       RETURN
    ENDIF
    !Normalise the normal (look, it could be worse OK)
    working_probe%normal=working_probe%normal/SQRT(SUM(working_probe%normal**2))

!!$    DO iCorner=1,4
!!$       DO iDirection=1,3
!!$          IF (working_probe%corner(iCorner,iDirection) .LT. working_probe%Extents(iDirection,1))&
!!$               working_probe%Extents(iDirection,1)=working_probe%corner(iCorner,iDirection)
!!$          IF (working_probe%corner(iCorner,iDirection) .GT. working_probe%Extents(iDirection,2))&
!!$               working_probe%Extents(iDirection,2)=working_probe%corner(iCorner,iDirection)
!!$       ENDDO
!!$    ENDDO

!!$    DO iDirection=1,3
!!$       IF (working_probe%Extents(iDirection,1) .GE. working_probe%Extents(iDirection,2)) THEN
!!$          IF (rank .EQ. 0) PRINT *,"Points specified for the probe plane corners collapse the probe plane to a lower DIMENSION. This is not currently supported. probe ",TRIM(working_probe%name)," abandoned."
!!$          DEALLOCATE(working_probe)
!!$          NULLIFY(working_probe)
!!$          RETURN
!!$       ENDIF
!!$    ENDDO

    CALL attach_probe(working_probe)

  END SUBROUTINE probe_block_end

  FUNCTION handle_probe_deck(element,value)
    CHARACTER(*),INTENT(IN) :: element,value
    INTEGER :: handle_probe_deck, ispecies

    handle_probe_deck=c_err_none

    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    ! get particle probe diagnostics (rolling total of all particles
    ! which pass through a given region of real space (defined by the line between two
    ! points in 2D).
    IF (str_cmp(element,"dump")) THEN
       working_probe%dump=as_integer(value,handle_probe_deck)
       RETURN
    ENDIF

    !Top left
    IF (str_cmp(element,"x_tl")) THEN
       working_probe%corner(1,1)=as_real(value,handle_probe_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"y_tl")) THEN
       working_probe%corner(1,2)=as_real(value,handle_probe_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"z_tl")) THEN
       working_probe%corner(1,3)=as_real(value,handle_probe_deck)
       RETURN
    ENDIF

    !Bottom right
    IF (str_cmp(element,"x_br")) THEN
       working_probe%corner(2,1)=as_real(value,handle_probe_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"y_br")) THEN
       working_probe%corner(2,2)=as_real(value,handle_probe_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"z_br")) THEN
       working_probe%corner(2,3)=as_real(value,handle_probe_deck)
       RETURN
    ENDIF

    !Top right
    IF (str_cmp(element,"x_tr")) THEN
       working_probe%corner(3,1)=as_real(value,handle_probe_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"y_tr")) THEN
       working_probe%corner(3,2)=as_real(value,handle_probe_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"z_tr")) THEN
       working_probe%corner(3,3)=as_real(value,handle_probe_deck)
       RETURN
    ENDIF

    !Bottom Left
    IF (str_cmp(element,"x_bl")) THEN
       working_probe%corner(4,1)=as_real(value,handle_probe_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"y_bl")) THEN
       working_probe%corner(4,2)=as_real(value,handle_probe_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"z_bl")) THEN
       working_probe%corner(4,3)=as_real(value,handle_probe_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"probe_species")) THEN
       ispecies=as_integer(value,handle_probe_deck)
       IF (handle_probe_deck .EQ. c_err_none) THEN
          IF (ispecies .GT. 0 .AND. ispecies .LE. n_species) THEN
             working_probe%probe_species=>particle_species(ispecies)
          ELSE
             IF (rank .EQ. 0) PRINT *,"Unable to attach probe to non existant species ",ispecies
             handle_probe_deck=c_err_bad_value
          ENDIF
       ENDIF
       RETURN
    ENDIF

    IF (str_cmp(element,"ek_min")) THEN
       working_probe%ek_min=as_real(value,handle_probe_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"ek_max")) THEN
       working_probe%ek_max=as_real(value,handle_probe_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"name")) THEN
       working_probe%name=TRIM(value)
       RETURN
    ENDIF

    handle_probe_deck = c_err_unknown_element

  END FUNCTION handle_probe_deck

#endif

END MODULE deck_eio_particle_probe_block
