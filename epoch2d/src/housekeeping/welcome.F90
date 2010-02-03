MODULE welcome

  USE shared_data

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: welcome_message

CONTAINS

  SUBROUTINE welcome_message

    INTEGER,PARAMETER :: logo_x=39,logo_y=11
    INTEGER,DIMENSION(logo_x,logo_y) :: logo
    CHARACTER(logo_x*2+1) :: logo_string
    CHARACTER,DIMENSION(5) :: logo_els
    INTEGER :: ix,iy

    IF (rank .NE. 0) RETURN

    logo_els=(/' ','@'," "," "," "/)

    PRINT *,""
    PRINT *,""

    logo(:,1 )=(/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/)
    logo(:,2 )=(/3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4/)
    logo(:,3 )=(/3,0,1,1,1,1,0,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,1,0,0,1,0,0,1,1,1,1,0,1,1,1,0,0,4/)
    logo(:,4 )=(/3,0,1,0,0,0,0,1,0,0,1,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,1,0,1,0,0,4/)
    logo(:,5 )=(/3,0,1,0,0,0,0,1,0,0,1,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,4/)
    logo(:,6 )=(/3,0,1,1,1,1,0,1,1,1,0,0,1,0,0,1,0,1,0,0,0,0,1,1,1,1,0,0,1,1,1,1,0,1,0,0,1,0,4/)
    logo(:,7 )=(/3,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,1,0,0,1,0,4/)
    logo(:,8 )=(/3,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,1,0,1,0,0,4/)
    logo(:,9 )=(/3,0,1,1,1,1,0,1,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,0,0,1,0,0,1,1,1,1,0,1,1,1,0,0,4/)
    logo(:,10)=(/3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4/)
    logo(:,11)=(/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/)
    logo_string=" "
    DO iy=1,logo_y*2+1
       DO ix=1,logo_x
          logo_string(ix*2-1:ix*2-1)=logo_els(logo(ix,MAX(iy/2,1))+1)
          logo_string(ix*2:ix*2)=logo_els(logo(ix,MAX(iy/2,1))+1)
       ENDDO
       WRITE(*,*),logo_string
    ENDDO
    WRITE(*,*) ""
    WRITE(*,'("Welcome to EPOCH2D Version ",I1,".",I1)'),c_version,c_revision
    WRITE(*,*) ""

    CALL compiler_directives

  END SUBROUTINE welcome_message

  SUBROUTINE compiler_directives

    WRITE(*,*) "The code was compiled with the following compile time options"
    WRITE(*,*) "*************************************************************"
#ifdef PART_DEBUG
    WRITE(*,*) "Particle Debug information -DPART_DEBUG"
#endif
#ifdef FIELD_DEBUG
    WRITE(*,*) "Field Debug information -DFIELD_DEBUG"
#endif
#ifdef SPLINE_FOUR
   WRITE(*,*) "Fourth order spline interpolation -DSPLINE_FOUR"
#endif
#ifdef HIGH_ORDER_FIELDS
#ifdef ORDER_SIX
   WRITE(*,*) "6th order improved field solver -DHIGH_ORDER_FIELDS -DORDER_SIX"
#else
   WRITE(*,*) "4th order field solver -DHIGH_ORDER_FIELDS"
#endif
#endif
#ifdef HIGH_ORDER_SMOOTHING
  WRITE(*,*) "High order current smoothing (matches particle interpolation function) -DHIGH_ORDER_SMOOTHING"
#endif
#ifdef PARTICLE_CELL_DIVISION
    WRITE(*,*) "Particle/cell ordering -DPARTICLE_CELL_DIVISION"
#endif
#ifdef PER_PARTICLE_WEIGHT
    WRITE(*,*) "Per particle weighting -DPER_PARTICLE_WEIGHT"
#endif
#ifdef PARTICLE_COUNT_UPDATE
    WRITE(*,*) "Global particle counting -DPARTICLE_COUNT_UPDATE"
#endif
#ifdef TRACER_PARTICLES
    WRITE(*,*) "Tracer particle support -DTRACER_PARTICLES"
#endif
#ifdef PARTICLE_PROBES
    WRITE(*,*) "Particle probe support -DPARTICLE_PROBES"
#endif
#ifdef PER_PARTICLE_CHARGEMASS
    WRITE(*,*) "Per particle charge and mass -DPER_PARTICLE_CHARGEMASS"
#endif
#ifdef PART_IONISE
    WRITE(*,*) "Particle ionisation model -DPART_IONISE"
#endif
#ifdef NO_DECK
    WRITE(*,*) "Deactivated input deck support -DNO_DECK"
#endif
#ifdef NEWTONIAN
    WRITE(*,*) "Newtonian dynamics (no relativity) -DNEWTONIAN"
#endif
    WRITE(*,*) "*************************************************************"
    WRITE(*,*) ""

  END SUBROUTINE compiler_directives

END MODULE welcome
