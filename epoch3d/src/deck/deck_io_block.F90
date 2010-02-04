MODULE deck_io_block

  USE shared_data
  USE strings
  USE strings_advanced

  IMPLICIT NONE

  SAVE

  INTEGER, PARAMETER :: n_var_special = 7
  INTEGER, PARAMETER :: io_block_elements = n_var_special+num_vars_to_dump
  LOGICAL, DIMENSION(io_block_elements)  :: io_block_done
  CHARACTER(LEN=string_length), DIMENSION(io_block_elements) :: io_block_name = (/"dt_snapshot", "full_dump_every", "restart_dump_every", "force_final_to_be_restartable", "use_offset_grid", "use_extended_io", "extended_io_file", "particles", "grid", "px", "py", "pz", "vx", "vy", "vz", "ex", "ey", "ez", "bx", "by", "bz", "jx", "jy", "jz", "charge", "mass", "ekbar", "mass_density", "charge_density", "number_density", "particle_weight", "species_id", "distribution_functions", "particle_probes", "temperature"/)

CONTAINS

  FUNCTION handle_io_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_io_deck
    INTEGER :: loop, elementselected

    handle_io_deck = c_err_unknown_element

    elementselected = 0

    DO loop = 1, io_block_elements
      IF(str_cmp(element, TRIM(ADJUSTL(io_block_name(loop))))) THEN
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
      use_extended_io = as_logical(value, handle_io_deck)
    CASE(7)
      extended_io_file = TRIM(value)
    END SELECT

    IF (elementselected .LE. n_var_special) RETURN
    IF (elementselected .GT. n_var_special) dumpmask(elementselected-n_var_special) = as_real(value, handle_io_deck)

    ! If setting dumpmask for particle probes then report if the code wasn't compiled for particle probes
#ifndef PARTICLE_PROBES
    IF (elementselected-n_var_special .EQ. 27) THEN
      handle_io_deck = c_err_pp_options_wrong
      extended_error_string = "-DPARTICLE_PROBES"
    ENDIF
#endif

  END FUNCTION handle_io_deck



  FUNCTION check_io_block()

    INTEGER :: check_io_block, index

    ! Just assume that anything not included except for the compulsory elements is not wanted
    check_io_block = c_err_none

    ! If not requesting extended io then don't check for extended_io_file
    IF (.NOT. io_block_done(6) .OR. .NOT. use_extended_io) io_block_done(6:7) = .TRUE.

    ! Particle Positions
    dumpmask(1:5) = IOR(dumpmask(1:5), c_io_restartable)
    ! Fields
    dumpmask(9:14) = IOR(dumpmask(9:14), c_io_restartable)
    ! Weight and species info
    dumpmask(24:25) = IOR(dumpmask(24:25), c_io_restartable)

    DO index = 1, n_var_special
      IF (.NOT. io_block_done(index)) THEN
        IF (rank .EQ. 0) THEN
          PRINT *, "***ERROR***"
          PRINT *, "Required output block element ", TRIM(ADJUSTL(io_block_name(index))), " absent. Please create this entry in the input deck"
          WRITE(40, *) ""
          WRITE(40, *) "***ERROR***"
          WRITE(40, *) "Required output block element ", TRIM(ADJUSTL(io_block_name(index))), " absent. Please create this entry in the input deck"
        ENDIF
        check_io_block = c_err_missing_elements
      ENDIF
    ENDDO

  END FUNCTION check_io_block

END MODULE deck_io_block
