MODULE welcome

  USE shared_data
  USE strings

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: welcome_message

CONTAINS

  SUBROUTINE welcome_message

    INTEGER,PARAMETER :: LOGOX=39,LOGOY=11
    INTEGER,DIMENSION(LOGOX,LOGOY) :: LOGO
    CHARACTER(LOGOX*2+1) :: LOGOSTRING
    CHARACTER,DIMENSION(5) :: LOGOELS
    INTEGER :: ix,iy
    CHARACTER(len=8) :: ver,rev

    IF (rank .NE. 0) RETURN

    LOGOELS=(/' ','@'," "," "," "/)

    PRINT *,""
    PRINT *,""

    LOGO(:,1 )=(/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/)
    LOGO(:,2 )=(/3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4/)
    LOGO(:,3 )=(/3,0,1,1,1,1,0,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,1,0,0,1,0,0,1,1,1,1,0,1,1,1,0,0,4/)
    LOGO(:,4 )=(/3,0,1,0,0,0,0,1,0,0,1,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,1,0,1,0,0,4/)
    LOGO(:,5 )=(/3,0,1,0,0,0,0,1,0,0,1,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,4/)
    LOGO(:,6 )=(/3,0,1,1,1,1,0,1,1,1,0,0,1,0,0,1,0,1,0,0,0,0,1,1,1,1,0,0,1,1,1,1,0,1,0,0,1,0,4/)
    LOGO(:,7 )=(/3,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,4/)
    LOGO(:,8 )=(/3,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,1,0,1,0,0,4/)
    LOGO(:,9 )=(/3,0,1,1,1,1,0,1,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,0,0,1,0,0,1,1,1,1,0,1,1,1,0,0,4/)
    LOGO(:,10)=(/3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4/)
    LOGO(:,11)=(/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/)
    LOGOSTRING=" "
    DO iy=1,LOGOY*2+1 
       DO ix=1,LOGOX
          LOGOSTRING(ix*2-1:ix*2-1)=LOGOELS(LOGO(ix,MAX(iy/2,1))+1)
          LOGOSTRING(ix*2:ix*2)=LOGOELS(LOGO(ix,MAX(iy/2,1))+1)
       ENDDO
       WRITE(*,*),LOGOSTRING
    ENDDO
    WRITE(*,*) ""
    CALL integer_as_string(version,ver)
    CALL integer_as_string(revision,rev)
    WRITE(*,*) "Welcome to EPOCH3D Version ",TRIM(ver),".",TRIM(ADJUSTL(rev)),"ALPHA"
    WRITE(*,*) ""

    CALL compiler_directives
    CALL mpi_status
    WRITE(*,*) ""

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
    WRITE(*,*) "*************************************************************"
    WRITE(*,*) ""

  END SUBROUTINE compiler_directives

  SUBROUTINE mpi_status
    CHARACTER(len=8) :: string

    CALL integer_as_string(nproc,string)

    WRITE(*,*) "Code is running on ",TRIM(string)," processing elements"
    WRITE(*,*) ""

  END SUBROUTINE mpi_status

END MODULE welcome
