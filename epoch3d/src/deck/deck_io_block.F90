MODULE deck_io_block

  USE mpi
  USE strings_advanced

  IMPLICIT NONE

  SAVE

  INTEGER, PARAMETER :: n_var_special = 8
  INTEGER, PARAMETER :: io_block_elements = n_var_special + num_vars_to_dump
  LOGICAL, DIMENSION(io_block_elements) :: io_block_done = .FALSE.
  LOGICAL :: need_dt = .FALSE.
  INTEGER, PARAMETER :: c_dump_part_grid         = 1
  INTEGER, PARAMETER :: c_dump_grid              = 2
  INTEGER, PARAMETER :: c_dump_part_species      = 3
  INTEGER, PARAMETER :: c_dump_part_weight       = 4
  INTEGER, PARAMETER :: c_dump_part_px           = 5
  INTEGER, PARAMETER :: c_dump_part_py           = 6
  INTEGER, PARAMETER :: c_dump_part_pz           = 7
  INTEGER, PARAMETER :: c_dump_part_vx           = 8
  INTEGER, PARAMETER :: c_dump_part_vy           = 9
  INTEGER, PARAMETER :: c_dump_part_vz           = 10
  INTEGER, PARAMETER :: c_dump_part_charge       = 11
  INTEGER, PARAMETER :: c_dump_part_mass         = 12
  INTEGER, PARAMETER :: c_dump_ex                = 13
  INTEGER, PARAMETER :: c_dump_ey                = 14
  INTEGER, PARAMETER :: c_dump_ez                = 15
  INTEGER, PARAMETER :: c_dump_bx                = 16
  INTEGER, PARAMETER :: c_dump_by                = 17
  INTEGER, PARAMETER :: c_dump_bz                = 18
  INTEGER, PARAMETER :: c_dump_jx                = 19
  INTEGER, PARAMETER :: c_dump_jy                = 20
  INTEGER, PARAMETER :: c_dump_jz                = 21
  INTEGER, PARAMETER :: c_dump_ekbar             = 22
  INTEGER, PARAMETER :: c_dump_mass_density      = 23
  INTEGER, PARAMETER :: c_dump_charge_density    = 24
  INTEGER, PARAMETER :: c_dump_number_density    = 25
  INTEGER, PARAMETER :: c_dump_temperature       = 26
  INTEGER, PARAMETER :: c_dump_dist_fns          = 27
  INTEGER, PARAMETER :: c_dump_probes            = 28
  INTEGER, PARAMETER :: c_dump_ejected_particles = 29
  CHARACTER(LEN=string_length), DIMENSION(io_block_elements) :: &
      io_block_name = (/ &
          "dt_snapshot                  ", & ! s1
          "full_dump_every              ", & ! s2
          "restart_dump_every           ", & ! s3
          "force_final_to_be_restartable", & ! s4
          "use_offset_grid              ", & ! s5
          "extended_io_file             ", & ! s6
          "averaging_period             ", & ! s7
          "min_cycles_per_average       ", & ! s8
          "particles                    ", & ! 1
          "grid                         ", & ! 2
          "species_id                   ", & ! 3
          "particle_weight              ", & ! 4
          "px                           ", & ! 5
          "py                           ", & ! 6
          "pz                           ", & ! 7
          "vx                           ", & ! 8
          "vy                           ", & ! 9
          "vz                           ", & ! 10
          "charge                       ", & ! 11
          "mass                         ", & ! 12
          "ex                           ", & ! 13
          "ey                           ", & ! 14
          "ez                           ", & ! 15
          "bx                           ", & ! 16
          "by                           ", & ! 17
          "bz                           ", & ! 18
          "jx                           ", & ! 19
          "jy                           ", & ! 20
          "jz                           ", & ! 21
          "ekbar                        ", & ! 22
          "mass_density                 ", & ! 23
          "charge_density               ", & ! 24
          "number_density               ", & ! 25
          "temperature                  ", & ! 26
          "distribution_functions       ", & ! 27
          "particle_probes              ", & ! 28
          "ejected_particles            " /) ! 29

CONTAINS

  FUNCTION handle_io_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_io_deck
    INTEGER :: loop, elementselected, mask, mask_element, ierr
    LOGICAL :: bad

    handle_io_deck = c_err_unknown_element

    elementselected = 0

    DO loop = 1, io_block_elements
      IF (str_cmp(element, TRIM(ADJUSTL(io_block_name(loop))))) THEN
        elementselected = loop
        EXIT
      ENDIF
    ENDDO

    IF (elementselected .EQ. 0) RETURN
    IF (io_block_done(elementselected)) THEN
      handle_io_deck = c_err_preset_element
      RETURN
    ENDIF
    io_block_done(elementselected) = .TRUE.
    handle_io_deck = c_err_none

    SELECT CASE (elementselected)
    CASE(1)
      dt_snapshots = as_real(value, handle_io_deck)
    CASE(2)
      full_dump_every = as_integer(value, handle_io_deck)
    CASE(3)
      restart_dump_every = as_integer(value, handle_io_deck)
    CASE(4)
      force_final_to_be_restartable = as_logical(value, handle_io_deck)
    CASE(5)
      use_offset_grid = as_logical(value, handle_io_deck)
    CASE(6)
      IF (rank .EQ. 0) THEN
        WRITE(*, *) '***ERROR***'
        WRITE(*, *) 'The "extended_io_file" option is no longer supported.'
        WRITE(*, *) 'Please use the "import" directive instead'
        WRITE(40,*) '***ERROR***'
        WRITE(40,*) 'The "extended_io_file" option is no longer supported.'
        WRITE(40,*) 'Please use the "import" directive instead'
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    CASE(7)
      average_time = as_real(value, handle_io_deck)
    CASE(8)
      min_cycles_per_average = as_integer(value, handle_io_deck)
    END SELECT

    IF (elementselected .LE. n_var_special) RETURN

    need_dt = .TRUE.

    mask = as_integer(value, handle_io_deck)
    mask_element = elementselected - n_var_special

    ! If setting dumpmask for particle probes then report if the code wasn't
    ! compiled for particle probes
#ifndef PARTICLE_PROBES
    IF (mask_element .EQ. c_dump_probes .AND. mask .NE. c_io_never) THEN
      handle_io_deck = c_err_pp_options_wrong
      extended_error_string = "-DPARTICLE_PROBES"
      mask = c_io_never
    ENDIF
#endif

    ! Setting some flags like species and average 
    ! wastes memory if the parameters make no sense. Do sanity checking.

    IF (IAND(mask, c_io_species) .NE. 0) THEN
      bad = .TRUE.
      ! Check for sensible per species variables
      IF (mask_element .EQ. c_dump_ekbar) bad = .FALSE.
      IF (mask_element .EQ. c_dump_mass_density) bad = .FALSE.
      IF (mask_element .EQ. c_dump_charge_density) bad = .FALSE.
      IF (mask_element .EQ. c_dump_number_density) bad = .FALSE.
      IF (mask_element .EQ. c_dump_temperature) bad = .FALSE.
      IF (bad) THEN
        IF (rank .EQ. 0) THEN
          PRINT*, '*** WARNING ***'
          PRINT*, 'Attempting to set per species property for "' &
              // TRIM(element) // '" which'
          PRINT*, 'does not support this property. Ignoring.'
        ENDIF
        mask = IAND(mask, NOT(c_io_species))
      ENDIF
    ENDIF

    IF (IAND(mask, c_io_averaged) .NE. 0) THEN
      bad = .TRUE.
      ! Check for sensible averaged variables
      IF (mask_element .EQ. c_dump_ex) bad = .FALSE.
      IF (mask_element .EQ. c_dump_ey) bad = .FALSE.
      IF (mask_element .EQ. c_dump_ez) bad = .FALSE.
      IF (mask_element .EQ. c_dump_bx) bad = .FALSE.
      IF (mask_element .EQ. c_dump_by) bad = .FALSE.
      IF (mask_element .EQ. c_dump_bz) bad = .FALSE.
      IF (mask_element .EQ. c_dump_jx) bad = .FALSE.
      IF (mask_element .EQ. c_dump_jy) bad = .FALSE.
      IF (mask_element .EQ. c_dump_jz) bad = .FALSE.
      IF (mask_element .EQ. c_dump_ekbar) bad = .FALSE.
      IF (mask_element .EQ. c_dump_mass_density) bad = .FALSE.
      IF (mask_element .EQ. c_dump_charge_density) bad = .FALSE.
      IF (mask_element .EQ. c_dump_number_density) bad = .FALSE.
      IF (mask_element .EQ. c_dump_temperature) bad = .FALSE.
      IF (bad) THEN
        IF (rank .EQ. 0) THEN
          WRITE(*,*) '*** WARNING ***'
          WRITE(*,*) 'Attempting to set average property for "' &
              // TRIM(element) // '" which'
          WRITE(*,*) 'does not support this property. Ignoring.'
        ENDIF
        mask = IAND(mask, NOT(c_io_averaged))
      ELSE
        any_average = .TRUE.
      ENDIF
    ENDIF

    dumpmask(mask_element) = mask

  END FUNCTION handle_io_deck



  FUNCTION check_io_block()

    INTEGER :: check_io_block, index

    ! Just assume that anything not included except for the compulsory
    ! elements is not wanted
    check_io_block = c_err_none

    IF (.NOT. need_dt) io_block_done(1) = .TRUE.
    ! Other control parameters are optional
    io_block_done(2:6) = .TRUE.
    ! Averaging info not compulsory unless averaged variable selected
    IF (.NOT. any_average) io_block_done(7:8) = .TRUE.

    IF (dt_snapshots .LT. average_time) THEN
      IF (rank .EQ. 0) THEN
        WRITE(*,*) '*** WARNING ***'
        WRITE(*,*) 'Averaging time is longer than dt_snapshot, will set', &
            ' dt_snapshot equal'
        WRITE(*,*) 'to averaging time.'
      ENDIF
      dt_snapshots = average_time
    ENDIF

    ! Particles
    dumpmask(c_dump_part_grid) = &
        IOR(dumpmask(c_dump_part_grid), c_io_restartable)
    dumpmask(c_dump_part_species) = &
        IOR(dumpmask(c_dump_part_species), c_io_restartable)
    dumpmask(c_dump_part_weight) = &
        IOR(dumpmask(c_dump_part_weight), c_io_restartable)
    dumpmask(c_dump_part_px) = IOR(dumpmask(c_dump_part_px), c_io_restartable)
    dumpmask(c_dump_part_py) = IOR(dumpmask(c_dump_part_py), c_io_restartable)
    dumpmask(c_dump_part_pz) = IOR(dumpmask(c_dump_part_pz), c_io_restartable)
    ! Fields
    dumpmask(c_dump_grid) = IOR(dumpmask(c_dump_grid), c_io_restartable)
    dumpmask(c_dump_ex) = IOR(dumpmask(c_dump_ex), c_io_restartable)
    dumpmask(c_dump_ey) = IOR(dumpmask(c_dump_ey), c_io_restartable)
    dumpmask(c_dump_ez) = IOR(dumpmask(c_dump_ez), c_io_restartable)
    dumpmask(c_dump_bx) = IOR(dumpmask(c_dump_bx), c_io_restartable)
    dumpmask(c_dump_by) = IOR(dumpmask(c_dump_by), c_io_restartable)
    dumpmask(c_dump_bz) = IOR(dumpmask(c_dump_bz), c_io_restartable)
    dumpmask(c_dump_jx) = IOR(dumpmask(c_dump_jx), c_io_restartable)
    dumpmask(c_dump_jy) = IOR(dumpmask(c_dump_jy), c_io_restartable)
    dumpmask(c_dump_jz) = IOR(dumpmask(c_dump_jz), c_io_restartable)

    DO index = 1, n_var_special
      IF (.NOT. io_block_done(index)) THEN
        IF (rank .EQ. 0) THEN
          PRINT *, "***ERROR***"
          PRINT *, "Required output block element ", &
              TRIM(ADJUSTL(io_block_name(index))), &
              " absent. Please create this entry in the input deck"
          WRITE(40, *) ""
          WRITE(40, *) "***ERROR***"
          WRITE(40, *) "Required output block element ", &
              TRIM(ADJUSTL(io_block_name(index))), &
              " absent. Please create this entry in the input deck"
        ENDIF
        check_io_block = c_err_missing_elements
      ENDIF
    ENDDO

  END FUNCTION check_io_block

END MODULE deck_io_block
