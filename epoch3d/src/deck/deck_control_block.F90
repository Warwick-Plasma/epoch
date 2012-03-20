MODULE deck_control_block

  USE mpi
  USE strings_advanced
  USE fields
  USE collisions

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: control_deck_initialise, control_deck_finalise
  PUBLIC :: control_block_start, control_block_end
  PUBLIC :: control_block_handle_element, control_block_check

  INTEGER, PARAMETER :: control_block_elements = 20 + 4 * c_ndims
  LOGICAL, DIMENSION(control_block_elements) :: control_block_done
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: coll_pairs_touched
  CHARACTER(LEN=string_length), DIMENSION(control_block_elements) :: &
      control_block_name = (/ &
          'nx                ', &
          'ny                ', &
          'nz                ', &
          'x_min             ', &
          'x_max             ', &
          'y_min             ', &
          'y_max             ', &
          'z_min             ', &
          'z_max             ', &
          'nprocx            ', &
          'nprocy            ', &
          'nprocz            ', &
          'npart             ', &
          'nsteps            ', &
          't_end             ', &
          'dt_multiplier     ', &
          'dlb_threshold     ', &
          'icfile            ', &
          'restart_snapshot  ', &
          'neutral_background', &
          'field_order       ', &
          'stdout_frequency  ', &
          'use_random_seed   ', &
          'smooth_currents   ', &
          'use_collisions    ', &
          'coulomb_log       ', &
          'collide           ', &
          'qed               ', &
          'qed_start_time    ', &
          'produce_photons   ', &
          'photon_energy_min ', &
          'produce_pairs     ' /)
  CHARACTER(LEN=string_length), DIMENSION(control_block_elements) :: &
      alternate_name = (/ &
          'nx                ', &
          'ny                ', &
          'nz                ', &
          'x_start           ', &
          'x_end             ', &
          'y_start           ', &
          'y_end             ', &
          'z_start           ', &
          'z_end             ', &
          'nprocx            ', &
          'nprocy            ', &
          'nprocz            ', &
          'npart             ', &
          'nsteps            ', &
          't_end             ', &
          'dt_multiplier     ', &
          'dlb_threshold     ', &
          'icfile            ', &
          'restart_snapshot  ', &
          'neutral_background', &
          'field_order       ', &
          'stdout_frequency  ', &
          'use_random_seed   ', &
          'smooth_currents   ', &
          'use_collisions    ', &
          'coulomb_log       ', &
          'collide           ', &
          'qed               ', &
          'qed_start_time    ', &
          'produce_photons   ', &
          'min_photon_energy ', &
          'produce_pairs     ' /)

CONTAINS

  SUBROUTINE control_deck_initialise

    IF (deck_state .EQ. c_ds_first) THEN
      control_block_done = .FALSE.
      use_collisions = .FALSE.
    ELSE
      ALLOCATE(coll_pairs_touched(1:n_species, 1:n_species))
      coll_pairs_touched = .FALSE.
      CALL setup_collisions
    ENDIF

  END SUBROUTINE control_deck_initialise



  SUBROUTINE control_deck_finalise

    INTEGER :: i, j

    IF (deck_state .EQ. c_ds_first) RETURN
    DEALLOCATE(coll_pairs_touched)

    IF (use_collisions) THEN
      use_collisions = .FALSE.
      DO j = 1, n_species
        DO i = 1, n_species
          IF (coll_pairs(i,j) .GT. 0) THEN
            use_collisions = .TRUE.
            EXIT
          ENDIF
        ENDDO
      ENDDO
      use_particle_lists = use_particle_lists .OR. use_collisions
    ENDIF

  END SUBROUTINE control_deck_finalise



  SUBROUTINE control_block_start

  END SUBROUTINE control_block_start



  SUBROUTINE control_block_end

  END SUBROUTINE control_block_end



  FUNCTION control_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    INTEGER :: loop, elementselected, field_order, ierr, io

    errcode = c_err_none

    IF (deck_state .NE. c_ds_first) THEN
      loop = 4*c_ndims + 15
      IF (str_cmp(element, TRIM(ADJUSTL(control_block_name(loop))))) THEN
        CALL set_collision_matrix(TRIM(ADJUSTL(value)), errcode)
      ENDIF
      RETURN
    ENDIF

    errcode = c_err_unknown_element
    elementselected = 0

    DO loop = 1, control_block_elements
      IF (str_cmp(element, TRIM(ADJUSTL(control_block_name(loop)))) &
          .OR. str_cmp(element, TRIM(ADJUSTL(alternate_name(loop))))) THEN
        elementselected = loop
        EXIT
      ENDIF
    ENDDO

    IF (elementselected .EQ. 0) RETURN

    IF (control_block_done(elementselected)) THEN
      errcode = c_err_preset_element
      RETURN
    ENDIF

    control_block_done(elementselected) = .TRUE.
    control_block_done(4*c_ndims+15) = .FALSE.
    errcode = c_err_none

    SELECT CASE (elementselected)
    CASE(1)
      nx_global = as_integer(value, errcode)
    CASE(2)
      ny_global = as_integer(value, errcode)
    CASE(3)
      nz_global = as_integer(value, errcode)
    CASE(c_ndims+1)
      x_min = as_real(value, errcode)
    CASE(c_ndims+2)
      x_max = as_real(value, errcode)
    CASE(c_ndims+3)
      y_min = as_real(value, errcode)
    CASE(c_ndims+4)
      y_max = as_real(value, errcode)
    CASE(c_ndims+5)
      z_min = as_real(value, errcode)
    CASE(c_ndims+6)
      z_max = as_real(value, errcode)
    CASE(3*c_ndims+1)
      nprocx = as_integer(value, errcode)
    CASE(3*c_ndims+2)
      nprocy = as_integer(value, errcode)
    CASE(3*c_ndims+3)
      nprocz = as_integer(value, errcode)
    CASE(4*c_ndims+1)
      npart_global = as_long_integer(value, errcode)
    CASE(4*c_ndims+2)
      nsteps = as_integer(value, errcode)
    CASE(4*c_ndims+3)
      t_end = as_real(value, errcode)
    CASE(4*c_ndims+4)
      dt_multiplier = as_real(value, errcode)
    CASE(4*c_ndims+5)
      dlb_threshold = as_real(value, errcode)
      dlb = .TRUE.
    CASE(4*c_ndims+6)
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The "icfile" option is no longer supported.'
          WRITE(io,*) 'Please use the "import" directive instead'
        ENDDO
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    CASE(4*c_ndims+7)
      restart_snapshot = as_integer(value, errcode)
      ic_from_restart = .TRUE.
    CASE(4*c_ndims+8)
      neutral_background = as_logical(value, errcode)
    CASE(4*c_ndims+9)
      field_order = as_integer(value, errcode)
      IF (field_order .NE. 2 .AND. field_order .NE. 4 &
          .AND. field_order .NE. 6) THEN
        errcode = c_err_bad_value
      ELSE
        CALL set_field_order(field_order)
      ENDIF
    CASE(4*c_ndims+10)
      stdout_frequency = as_integer(value, errcode)
    CASE(4*c_ndims+11)
      use_random_seed = as_logical(value, errcode)
    CASE(4*c_ndims+12)
      smooth_currents = as_logical(value, errcode)
    CASE(4*c_ndims+13)
      use_collisions = as_logical(value, errcode)
    CASE(4*c_ndims+14)
      IF (str_cmp(value, 'auto')) THEN
        coulomb_log_auto = .TRUE.
      ELSE
        coulomb_log_auto = .FALSE.
        coulomb_log = as_real(value, errcode)
      ENDIF
    ! 4*c_ndims+15 is the special case of the collision matrix
    CASE(4*c_ndims+16)
#ifdef PHOTONS
      use_qed = as_logical(value, errcode)
#else
      extended_error_string = '-DPHOTONS'
      errcode = c_err_pp_options_wrong
#endif
    CASE(4*c_ndims+17)
#ifdef PHOTONS
      qed_start_time = as_real(value, errcode)
#else
      extended_error_string = '-DPHOTONS'
      errcode = c_err_pp_options_wrong
#endif
    CASE(4*c_ndims+18)
#ifdef PHOTONS
      produce_photons = as_logical(value, errcode)
#else
      extended_error_string = '-DPHOTONS'
      errcode = c_err_pp_options_wrong
#endif
    CASE(4*c_ndims+19)
#ifdef PHOTONS
      photon_energy_min = as_real(value, errcode)
#else
      extended_error_string = '-DPHOTONS'
      errcode = c_err_pp_options_wrong
#endif
    CASE(4*c_ndims+20)
#ifdef PHOTONS
      produce_pairs = as_logical(value, errcode)
#else
      extended_error_string = '-DPHOTONS'
      errcode = c_err_pp_options_wrong
#endif
    END SELECT

  END FUNCTION control_block_handle_element



  FUNCTION control_block_check() RESULT(errcode)

    INTEGER :: errcode, index, io

    errcode = c_err_none

    ! nprocx/y/z and npart are optional
    control_block_done(3*c_ndims+1:4*c_ndims+1) = .TRUE.

    ! Only one of nsteps or t_end need be specified
    IF (control_block_done(4*c_ndims+2)) &
        control_block_done(4*c_ndims+3) = .TRUE.
    IF (control_block_done(4*c_ndims+3)) &
        control_block_done(4*c_ndims+2) = .TRUE.

    ! All entries after t_end are optional
    control_block_done(4*c_ndims+4:) = .TRUE.

    ! QED stuff is incorrect if not compiled with the correct options
#ifndef PHOTONS
    DO index = 4*c_ndims+16, 4*c_ndims+20
      IF (control_block_done(index)) THEN
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Cannot turn on QED option ', &
                TRIM(ADJUSTL(control_block_name(index))), &
                ' unless code is compiled with the correct ', &
                'preprocessor option (-DPHOTONS)'
          ENDDO
        ENDIF
        extended_error_string = '-DPHOTONS'
        errcode = c_err_pp_options_wrong
      ENDIF
    ENDDO
#endif

    DO index = 1, control_block_elements
      IF (.NOT. control_block_done(index)) THEN
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Required control block element ', &
                TRIM(ADJUSTL(control_block_name(index))), &
                ' absent. Please create this entry in the input deck'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF
    ENDDO

    IF (.NOT. neutral_background) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The option "neutral_background=F" is not supported', &
              ' in this version of EPOCH.'
        ENDDO
      ENDIF
      errcode = c_err_terminate
    ENDIF

  END FUNCTION control_block_check



! The following code is all about reading the coll_pairs from the input deck

  SUBROUTINE get_token(str_in, str_out, token_out, err)

    CHARACTER(*), INTENT(IN) :: str_in
    CHARACTER(*), INTENT(OUT) :: str_out
    CHARACTER(*), INTENT(OUT) :: token_out
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: str_len, char, pos
    CHARACTER(1) :: c

    str_len = LEN(str_in)
    pos = str_len

    DO char = 1, str_len
      c = str_in(char:char)
      IF (c .EQ. ' ')  THEN
        pos = char
        EXIT
      ENDIF
    ENDDO

    IF (pos .LT. str_len) THEN
      str_out = TRIM(ADJUSTL(str_in(pos+1:str_len)))
    ELSE
      str_out = ''
    ENDIF

    token_out = TRIM(str_in(1:pos))

  END SUBROUTINE get_token



  SUBROUTINE set_collision_matrix(str_in, errcode)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: errcode
    CHARACTER(LEN=string_length) :: tstr1, tstr2
    CHARACTER(LEN=string_length) :: species1, species2
    REAL(num) :: collstate
    INTEGER :: io, sp1, sp2

    IF (deck_state .NE. c_ds_last) RETURN

    IF (str_cmp(TRIM(str_in), 'all')) THEN
      coll_pairs = 1.0_num
      RETURN
    ENDIF

    IF (str_cmp(TRIM(str_in), 'none')) THEN
      coll_pairs = -1.0_num
      RETURN
    ENDIF

    CALL get_token(str_in, tstr1, species1, errcode)
    IF (errcode .NE. 0) RETURN

    sp1 = as_integer(species1, errcode)
    IF (errcode .NE. 0) RETURN

    CALL get_token(tstr1, tstr2, species2, errcode)
    IF (errcode .NE. 0) RETURN

    sp2 = as_integer(species2, errcode)
    IF (errcode .NE. 0) RETURN

    collstate = 1.0_num
    IF (str_cmp(TRIM(tstr2), 'on') .OR. str_cmp(TRIM(tstr2), '')) THEN
      collstate = 1.0_num
    ELSEIF (str_cmp(TRIM(tstr2), 'off')) THEN
      collstate = -1.0_num
    ELSE
      collstate = as_real(tstr2, errcode)
      IF (errcode .NE. 0) RETURN
    ENDIF

    IF (coll_pairs_touched(sp1, sp2) .AND. rank .EQ. 0) THEN
      DO io = stdout, du, du - stdout ! Print to stdout and to file
        WRITE(io,*)
        WRITE(io,*) '*** WARNING ***'
        WRITE(io,*) 'The collide parameter for ' // TRIM(species1) // ' <-> ' &
            // TRIM(species2)
        WRITE(io,*) 'has been set multiple times!'
        WRITE(io,*) 'Collisions will only be carried out once per species pair.'
        WRITE(io,*) 'Later specifications will always override earlier ones.'
        WRITE(io,*)
      ENDDO
    ENDIF

    coll_pairs(sp1, sp2) = collstate
    coll_pairs(sp2, sp1) = collstate
    coll_pairs_touched(sp1, sp2) = .TRUE.
    coll_pairs_touched(sp2, sp1) = .TRUE.

  END SUBROUTINE set_collision_matrix

END MODULE deck_control_block
