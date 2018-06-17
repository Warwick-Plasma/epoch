! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2015      Stephan Kuschel <stephan.kuschel@gmail.com>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
! Copyright (C) 2012      Martin Ramsay <M.G.Ramsay@warwick.ac.uk>
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

!******************************************************************************
! Welcome message routines
!******************************************************************************

MODULE welcome

  USE version_data
  USE shared_data
  USE terminal_controls

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: welcome_message, create_ascii_header

CONTAINS

  !****************************************************************************
  ! This routine prints the welcome message, MPI status
  !****************************************************************************

  SUBROUTINE welcome_message

    IF (rank /= 0) RETURN

    CALL set_term_attr(c_term_bold)
    CALL set_term_attr(c_term_blue)
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*)

    WRITE(*,'(A)') '        d########P  d########b        .######b          ' &
        // 'd#######  d##P      d##P'
    WRITE(*,'(A)') '       d########P  d###########    d###########     .###' &
        // '#######  d##P      d##P '
    CALL set_term_attr(c_term_cyan)
    WRITE(*,'(A)',ADVANCE='NO') '      ----        ----     ----  -----     ' &
        // '----   -----         ----      --'
    CALL set_term_attr(c_term_blue)
    WRITE(*,'(A)') ' P  '
    WRITE(*,'(A)') '     d########P  d####,,,####P ####.      .#### d###P   ' &
        // '       d############P   '
    WRITE(*,'(A)') '    d########P  d#########P   ####       .###P ####.    ' &
        // '      d############P    '
    WRITE(*,'(A)') '   d##P        d##P           ####     d####   ####.    ' &
        // '     d##P      d##P     '
    WRITE(*,'(A)') '  d########P  d##P            ###########P     #########' &
        // '#P  d##P      d##P      '
    WRITE(*,'(A)') ' d########P  d##P              d######P          #######' &
        // 'P  d##P      d##P       '
    WRITE(*,*)

    CALL create_ascii_header
    CALL set_term_attr(c_term_green)
    WRITE(*,*) 'Welcome to ', TRIM(c_code_name), ' version ', &
        TRIM(version_string) // '   (commit ' // TRIM(c_commit_id) // ')'
    WRITE(*,*)
    CALL set_term_attr(c_term_default_colour)
    CALL set_term_attr(c_term_reset_attributes)
    CALL compiler_directives
    CALL mpi_status_message

  END SUBROUTINE welcome_message



  SUBROUTINE compiler_directives

    LOGICAL :: found = .FALSE.

#ifdef PARTICLE_DEBUG
    found = .TRUE.
#endif
#ifdef PARSER_DEBUG
    found = .TRUE.
#endif
#ifdef PARSER_CHECKING
    found = .TRUE.
#endif
#ifdef PARTICLE_SHAPE_BSPLINE3
    found = .TRUE.
#endif
#ifdef PARTICLE_SHAPE_TOPHAT
    found = .TRUE.
#endif
#ifdef PER_SPECIES_WEIGHT
    found = .TRUE.
#endif
#ifdef NO_TRACER_PARTICLES
    found = .TRUE.
#endif
#ifdef NO_PARTICLE_PROBES
    found = .TRUE.
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    found = .TRUE.
#endif
#ifdef HIGH_ORDER_SMOOTHING
    found = .TRUE.
#endif
#ifdef PARTICLE_ID4
    found = .TRUE.
#endif
#ifdef PARTICLE_ID
    found = .TRUE.
#endif
#ifdef PHOTONS
    found = .TRUE.
#ifdef TRIDENT_PHOTONS
    found = .TRUE.
#endif
#endif
#ifdef PREFETCH
    found = .TRUE.
#endif
#ifdef MPI_DEBUG
    found = .TRUE.
#endif
#ifdef NO_IO
    found = .TRUE.
#endif

    IF (.NOT.found) THEN
      WRITE(*,*) '*************************************************************'
      WRITE(*,*) 'The code was compiled with no compile time options'
      WRITE(*,*) '*************************************************************'
      RETURN
    ENDIF

    WRITE(*,*) 'The code was compiled with the following compile time options'
    WRITE(*,*) '*************************************************************'
#ifdef PARTICLE_DEBUG
    defines = IOR(defines, c_def_particle_debug)
    WRITE(*,*) 'Particle Debug information -DPARTICLE_DEBUG'
#endif
#ifdef PARSER_DEBUG
    defines = IOR(defines, c_def_parser_debug)
    WRITE(*,*) 'Particle Debug information -DPARSER_DEBUG'
#endif
#ifdef PARSER_CHECKING
    defines = IOR(defines, c_def_parser_checking)
    WRITE(*,*) 'Parser FPE handling -DPARSER_CHECKING'
#endif
#ifdef PARTICLE_SHAPE_BSPLINE3
    defines = IOR(defines, c_def_particle_shape_bspline3)
    WRITE(*,*) 'Third order B-spline particle shape -DPARTICLE_SHAPE_BSPLINE3'
#endif
#ifdef PARTICLE_SHAPE_TOPHAT
    defines = IOR(defines, c_def_particle_shape_tophat)
    WRITE(*,*) 'Top-hat particle shape -DPARTICLE_SHAPE_TOPHAT'
#endif
#ifdef PER_SPECIES_WEIGHT
    WRITE(*,*) 'Per species weighting -DPER_SPECIES_WEIGHT'
#else
    defines = IOR(defines, c_def_per_particle_weight)
#endif
#ifdef NO_TRACER_PARTICLES
    WRITE(*,*) 'No tracer particle support -DNO_TRACER_PARTICLES'
#else
    defines = IOR(defines, c_def_tracer_particles)
#endif
#ifdef NO_PARTICLE_PROBES
    WRITE(*,*) 'No particle probe support -DNO_PARTICLE_PROBES'
#else
    defines = IOR(defines, c_def_particle_probes)
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    defines = IOR(defines, c_def_per_particle_chargemass)
    WRITE(*,*) 'Per particle charge and mass -DPER_PARTICLE_CHARGE_MASS'
#endif
#ifdef HIGH_ORDER_SMOOTHING
    defines = IOR(defines, c_def_high_order_smoothing)
    WRITE(*,*) 'High order current smoothing (matches particle ', &
        'interpolation function) -DHIGH_ORDER_SMOOTHING'
#endif
#ifdef PARTICLE_ID4
    defines = IOR(defines, c_def_particle_id4)
    WRITE(*,*) 'Particle ID tracking (4-bytes) -DPARTICLE_ID4'
#endif
#ifdef PARTICLE_ID
    defines = IOR(defines, c_def_particle_id)
    WRITE(*,*) 'Particle ID tracking (8-bytes) -DPARTICLE_ID'
#endif
#ifdef PHOTONS
    defines = IOR(defines, c_def_photons)
    WRITE(*,*) 'QED Effects -DPHOTONS'
#ifdef TRIDENT_PHOTONS
    defines = IOR(defines, c_def_trident_photons)
    WRITE(*,*) 'Pair production by Trident process -DTRIDENT_PHOTONS'
#endif
#endif
#ifdef PREFETCH
    defines = IOR(defines, c_def_prefetch)
    WRITE(*,*) 'Particle prefetching -DPREFETCH'
    WRITE(*,*) 'WARNING: sometimes causes errors'
#endif
#ifdef MPI_DEBUG
    defines = IOR(defines, c_def_mpi_debug)
    WRITE(*,*) 'MPI error handling -DMPI_DEBUG'
#endif
#ifdef NO_IO
    ! There is no need to add a c_def for this since no I/O occurs.
    WRITE(*,*) 'Perform no I/O -DNO_IO'
#endif
#ifdef DELTAF_METHOD
    defines = IOR(defines, c_def_deltaf_method)
    WRITE(*,*) 'Delta-f method -DDELTAF_METHOD'
#endif
#ifdef DELTAF_DEBUG
    defines = IOR(defines, c_def_deltaf_debug)
    WRITE(*,*) 'Delta-f debugging -DDELTAF_DEBUG'
#endif
#ifdef WORK_DONE_INTEGRATED
    defines = IOR(defines, c_def_work_done_integrated)
    WRITE(*,*) 'Work done on each particle -DWORK_DONE_INTEGRATED'
#endif
    WRITE(*,*) '*************************************************************'

  END SUBROUTINE compiler_directives


  !****************************************************************************
  ! This routine prints the mpi status information
  !****************************************************************************

  SUBROUTINE mpi_status_message

    CHARACTER(LEN=8) :: string

    CALL integer_as_string(nproc, string)

    WRITE(*,*) 'Code is running on ', TRIM(string), ' processing elements'
    WRITE(*,*)

  END SUBROUTINE mpi_status_message



  SUBROUTINE create_ascii_header

    CHARACTER(LEN=16) :: job1, job2
    CHARACTER(LEN=4) :: str
    INTEGER :: i, strmin, strmax, strlen

    ! Parse commit string to get version number
    ! Commit ID begins with the string v[0-9].[0-9].[0-9]-
    strlen = LEN_TRIM(c_commit_id)
    strmin = 2
    strmax = strmin + 4

    ! Version
    DO i = strmin, MIN(strmax,strlen)
      IF (c_commit_id(i:i) == '.') THEN
        str = c_commit_id(strmin:i-1)
        READ(str, '(i9)') c_version
        strmin = i + 1
        strmax = strmin + 4
        EXIT
      ENDIF
    ENDDO

    ! Revision
    DO i = strmin, MIN(strmax,strlen)
      IF (c_commit_id(i:i) == '.') THEN
        str = c_commit_id(strmin:i-1)
        READ(str, '(i9)') c_revision
        strmin = i + 1
        strmax = strmin + 4
        EXIT
      ENDIF
    ENDDO

    ! Minor revision
    DO i = strmin, MIN(strmax,strlen)
      IF (c_commit_id(i:i) == '-') THEN
        str = c_commit_id(strmin:i-1)
        READ(str, '(i9)') c_minor_rev
        strmax = i - 1
        EXIT
      ENDIF
    ENDDO

    version_string = c_commit_id(2:strmax)

    CALL integer_as_string(jobid%start_seconds, job1)
    CALL integer_as_string(jobid%start_milliseconds, job2)
    ascii_header = c_code_name // ' v' // TRIM(version_string) // '   ' &
        // c_commit_id // ' ' // TRIM(job1) // '.' // TRIM(ADJUSTL(job2))

  END SUBROUTINE create_ascii_header



  SUBROUTINE integer_as_string(int_in, string)

    INTEGER, INTENT(IN) :: int_in
    CHARACTER(LEN=*), INTENT(OUT) :: string

    INTEGER :: n_nums
    CHARACTER(LEN=9) :: numfmt

    IF (int_in == 0) THEN
      n_nums = 1
    ELSE
      n_nums = 1 + INT(LOG10(REAL(ABS(int_in), num)))
    ENDIF
    WRITE(numfmt, '(''(I'', I6.6, '')'')') n_nums
    WRITE(string, numfmt) int_in

  END SUBROUTINE integer_as_string

END MODULE welcome
