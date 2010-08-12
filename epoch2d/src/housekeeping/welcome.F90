MODULE welcome

  USE strings
  USE version_data

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: welcome_message

CONTAINS

  SUBROUTINE welcome_message

    INTEGER, PARAMETER :: logo_x = 39, logo_y = 11
    INTEGER, DIMENSION(logo_x, logo_y) :: logo
    CHARACTER(logo_x*2+1) :: logo_string
    CHARACTER, DIMENSION(5) :: logo_els
    INTEGER :: ix, iy
    CHARACTER(LEN=11) :: ver, rev

    IF (rank .NE. 0) RETURN

    logo_els = (/' ', '@', " ", " ", " "/)

    PRINT *, ""
    PRINT *, ""

    logo(:,1 ) = (/2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2/)
    logo(:,2 ) = (/3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4/)
    logo(:,3 ) = (/3, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, &
        1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 4/)
    logo(:,4 ) = (/3, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, &
        0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 4/)
    logo(:,5 ) = (/3, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, &
        0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 4/)
    logo(:,6 ) = (/3, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, &
        0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 4/)
    logo(:,7 ) = (/3, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, &
        0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 4/)
    logo(:,8 ) = (/3, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, &
        0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 4/)
    logo(:,9 ) = (/3, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, &
        1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 4/)
    logo(:,10) = (/3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4/)
    logo(:,11) = (/2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2/)

    logo_string = " "
    DO iy = 1, logo_y*2+1
      DO ix = 1, logo_x
        logo_string(ix*2-1:ix*2-1) = logo_els(logo(ix, MAX(iy/2, 1))+1)
        logo_string(ix*2:ix*2) = logo_els(logo(ix, MAX(iy/2, 1))+1)
      ENDDO
      WRITE(*, *) logo_string
    ENDDO

    WRITE(*, *) ""
    CALL integer_as_string(c_version, ver)
    CALL integer_as_string(c_revision, rev)
    version_string = TRIM(ver) // "." // TRIM(ADJUSTL(rev))
    CALL integer_as_string(jobid%start_seconds, ver)
    CALL integer_as_string(jobid%start_milliseconds, rev)
    ascii_header = c_code_name // " v" // TRIM(version_string) // " " &
        // c_commit_id // " " // TRIM(ver) // "." // TRIM(ADJUSTL(rev))
    WRITE(*, *) "Welcome to EPOCH2D Version ", version_string
    WRITE(*, *) ""

    CALL compiler_directives
    CALL mpi_status
    WRITE(*, *) ""

  END SUBROUTINE welcome_message



  SUBROUTINE compiler_directives

    WRITE(*, *) "The code was compiled with the following compile time options"
    WRITE(*, *) "*************************************************************"
#ifdef PARTICLE_DEBUG
    defines = IOR(defines, c_def_particle_debug)
    WRITE(*, *) "Particle Debug information -DPARTICLE_DEBUG"
#endif
#ifdef FIELD_DEBUG
    defines = IOR(defines, c_def_field_debug)
    WRITE(*, *) "Field Debug information -DFIELD_DEBUG"
#endif
#ifdef PARSER_DEBUG
    defines = IOR(defines, c_def_parser_debug)
    WRITE(*, *) "Particle Debug information -DPARSER_DEBUG"
#endif
#ifdef PARTICLE_SHAPE_BSPLINE3
    defines = IOR(defines, c_def_particle_shape_bspline3)
    WRITE(*, *) "Third order B-spline particle shape -DPARTICLE_SHAPE_BSPLINE3"
#endif
#ifdef PARTICLE_SHAPE_TOPHAT
    defines = IOR(defines, c_def_particle_shape_tophat)
    WRITE(*, *) "Top-hat particle shape -DPARTICLE_SHAPE_TOPHAT"
#endif
#ifdef SPLIT_PARTICLES_AFTER_PUSH
    defines = IOR(defines, c_def_split_part_after_push)
    WRITE(*, *) "Particle/cell ordering -DSPLIT_PARTICLES_AFTER_PUSH"
#endif
#ifdef PER_PARTICLE_WEIGHT
    defines = IOR(defines, c_def_per_particle_weight)
    WRITE(*, *) "Per particle weighting -DPER_PARTICLE_WEIGHT"
#endif
#ifdef PARTICLE_COUNT_UPDATE
    defines = IOR(defines, c_def_particle_count_update)
    WRITE(*, *) "Global particle counting -DPARTICLE_COUNT_UPDATE"
#endif
#ifdef TRACER_PARTICLES
    defines = IOR(defines, c_def_tracer_particles)
    WRITE(*, *) "Tracer particle support -DTRACER_PARTICLES"
#endif
#ifdef PARTICLE_PROBES
    defines = IOR(defines, c_def_particle_probes)
    WRITE(*, *) "Particle probe support -DPARTICLE_PROBES"
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    defines = IOR(defines, c_def_per_particle_chargemass)
    WRITE(*, *) "Per particle charge and mass -DPER_PARTICLE_CHARGE_MASS"
#endif
#ifdef PARTICLE_IONISE
    defines = IOR(defines, c_def_particle_ionise)
    WRITE(*, *) "Particle ionisation model -DPARTICLE_IONISE"
#endif
#ifdef HIGH_ORDER_SMOOTHING
    defines = IOR(defines, c_def_high_order_smoothing)
    WRITE(*, *) "High order current smoothing (matches particle &
        &interpolation function) -DHIGH_ORDER_SMOOTHING"
#endif
    WRITE(*, *) "*************************************************************"
    WRITE(*, *) ""

  END SUBROUTINE compiler_directives



  SUBROUTINE mpi_status

    CHARACTER(LEN=8) :: string

    CALL integer_as_string(nproc, string)

    WRITE(*, *) "Code is running on ", TRIM(string), " processing elements"
    WRITE(*, *) ""

  END SUBROUTINE mpi_status

END MODULE welcome
