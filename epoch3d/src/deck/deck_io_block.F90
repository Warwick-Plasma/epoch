MODULE deck_io_block

  USE mpi
  USE strings_advanced

  IMPLICIT NONE

  SAVE

  INTEGER, PARAMETER :: io_block_elements = num_vars_to_dump + 11
  LOGICAL, DIMENSION(io_block_elements) :: io_block_done = .FALSE.
  CHARACTER(LEN=string_length), DIMENSION(io_block_elements) :: io_block_name

CONTAINS

  SUBROUTINE io_deck_initialise

    INTEGER :: i

    io_block_name(c_dump_part_grid        ) = 'particles'
    io_block_name(c_dump_grid             ) = 'grid'
    io_block_name(c_dump_part_species     ) = 'species_id'
    io_block_name(c_dump_part_weight      ) = 'particle_weight'
    io_block_name(c_dump_part_px          ) = 'px'
    io_block_name(c_dump_part_py          ) = 'py'
    io_block_name(c_dump_part_pz          ) = 'pz'
    io_block_name(c_dump_part_vx          ) = 'vx'
    io_block_name(c_dump_part_vy          ) = 'vy'
    io_block_name(c_dump_part_vz          ) = 'vz'
    io_block_name(c_dump_part_charge      ) = 'charge'
    io_block_name(c_dump_part_mass        ) = 'mass'
    io_block_name(c_dump_ex               ) = 'ex'
    io_block_name(c_dump_ey               ) = 'ey'
    io_block_name(c_dump_ez               ) = 'ez'
    io_block_name(c_dump_bx               ) = 'bx'
    io_block_name(c_dump_by               ) = 'by'
    io_block_name(c_dump_bz               ) = 'bz'
    io_block_name(c_dump_jx               ) = 'jx'
    io_block_name(c_dump_jy               ) = 'jy'
    io_block_name(c_dump_jz               ) = 'jz'
    io_block_name(c_dump_ekbar            ) = 'ekbar'
    io_block_name(c_dump_mass_density     ) = 'mass_density'
    io_block_name(c_dump_charge_density   ) = 'charge_density'
    io_block_name(c_dump_number_density   ) = 'number_density'
    io_block_name(c_dump_temperature      ) = 'temperature'
    io_block_name(c_dump_dist_fns         ) = 'distribution_functions'
    io_block_name(c_dump_probes           ) = 'particle_probes'
    io_block_name(c_dump_ejected_particles) = 'ejected_particles'

    i = num_vars_to_dump
    io_block_name(i+1 ) = 'dt_snapshot'
    io_block_name(i+2 ) = 'full_dump_every'
    io_block_name(i+3 ) = 'restart_dump_every'
    io_block_name(i+4 ) = 'force_final_to_be_restartable'
    io_block_name(i+5 ) = 'use_offset_grid'
    io_block_name(i+6 ) = 'extended_io_file'
    io_block_name(i+7 ) = 'averaging_period'
    io_block_name(i+8 ) = 'min_cycles_per_average'
    io_block_name(i+9 ) = 'nstep_snapshot'
    io_block_name(i+10) = 'dump_source_code'
    io_block_name(i+11) = 'dump_input_decks'

  END SUBROUTINE io_deck_initialise



  SUBROUTINE io_deck_finalise

  END SUBROUTINE io_deck_finalise



  SUBROUTINE io_block_start

  END SUBROUTINE io_block_start



  SUBROUTINE io_block_end

  END SUBROUTINE io_block_end



  FUNCTION io_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    INTEGER :: loop, elementselected, mask, mask_element, ierr, io
    LOGICAL :: bad

    errcode = c_err_unknown_element

    elementselected = 0

    DO loop = 1, io_block_elements
      IF (str_cmp(element, TRIM(ADJUSTL(io_block_name(loop))))) THEN
        elementselected = loop
        EXIT
      ENDIF
    ENDDO

    IF (elementselected .EQ. 0) RETURN
    IF (io_block_done(elementselected)) THEN
      errcode = c_err_preset_element
      RETURN
    ENDIF
    io_block_done(elementselected) = .TRUE.
    errcode = c_err_none

    SELECT CASE (elementselected-num_vars_to_dump)
    CASE(1)
      dt_snapshot = as_real(value, errcode)
      IF (dt_snapshot .LT. 0.0_num) dt_snapshot = 0.0_num
    CASE(2)
      full_dump_every = as_integer(value, errcode)
    CASE(3)
      restart_dump_every = as_integer(value, errcode)
    CASE(4)
      force_final_to_be_restartable = as_logical(value, errcode)
    CASE(5)
      use_offset_grid = as_logical(value, errcode)
    CASE(6)
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The "extended_io_file" option is no longer supported.'
          WRITE(io,*) 'Please use the "import" directive instead'
        ENDDO
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    CASE(7)
      average_time = as_real(value, errcode)
    CASE(8)
      min_cycles_per_average = as_integer(value, errcode)
    CASE(9)
      nstep_snapshot = as_integer(value, errcode)
      IF (nstep_snapshot .LT. 0) nstep_snapshot = 0
    CASE(10)
      dump_source_code = as_logical(value, errcode)
    CASE(11)
      dump_input_decks = as_logical(value, errcode)
    END SELECT

    IF (elementselected .GT. num_vars_to_dump) RETURN

    mask = as_integer(value, errcode)
    mask_element = elementselected

    ! If setting dumpmask for particle probes then report if the code wasn't
    ! compiled for particle probes
#ifndef PARTICLE_PROBES
    IF (mask_element .EQ. c_dump_probes .AND. mask .NE. c_io_never) THEN
      errcode = c_err_pp_options_wrong
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
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Attempting to set per species property for "' &
                // TRIM(element) // '" which'
            WRITE(io,*) 'does not support this property. Ignoring.'
          ENDDO
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
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Attempting to set average property for "' &
                // TRIM(element) // '" which'
            WRITE(io,*) 'does not support this property. Ignoring.'
          ENDDO
        ENDIF
        mask = IAND(mask, NOT(c_io_averaged))
      ELSE
        any_average = .TRUE.
      ENDIF
    ENDIF

    dumpmask(mask_element) = mask

  END FUNCTION io_block_handle_element



  FUNCTION io_block_check() RESULT(errcode)

    INTEGER :: errcode, index, io, i

    ! Just assume that anything not included except for the compulsory
    ! elements is not wanted
    errcode = c_err_none

    ! Other control parameters are optional
    i = num_vars_to_dump
    io_block_done(i+1:i+6) = .TRUE.
    io_block_done(i+9:io_block_elements) = .TRUE.
    ! Averaging info not compulsory unless averaged variable selected
    IF (.NOT. any_average) io_block_done(i+7:i+8) = .TRUE.

    DO index = i+7, i+8
      IF (.NOT. io_block_done(index)) THEN
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Required output block element ', &
                TRIM(ADJUSTL(io_block_name(index))), &
                ' absent. Please create this entry in the input deck'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF
    ENDDO

    IF (dt_snapshot .LT. average_time) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Averaging time is longer than dt_snapshot, will set', &
              ' dt_snapshot equal'
          WRITE(io,*) 'to averaging time.'
        ENDDO
      ENDIF
      dt_snapshot = average_time
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

  END FUNCTION io_block_check

END MODULE deck_io_block
