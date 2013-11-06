MODULE deck_qed_block

  USE strings_advanced

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: qed_deck_initialise, qed_deck_finalise
  PUBLIC :: qed_block_start, qed_block_end
  PUBLIC :: qed_block_handle_element, qed_block_check

CONTAINS

  SUBROUTINE qed_deck_initialise

#ifdef PHOTONS
    IF (deck_state .EQ. c_ds_first) THEN
      qed_table_location = 'src/physics_packages/TABLES'
      use_radiation_reaction = .TRUE.
    ENDIF
#endif

  END SUBROUTINE qed_deck_initialise



  SUBROUTINE qed_deck_finalise

    INTEGER :: io, iu, ierr
#ifdef PHOTONS
    LOGICAL :: exists

    IF (deck_state .EQ. c_ds_first) RETURN

    IF (rank .EQ. 0 .AND. use_qed) THEN
      INQUIRE(file=TRIM(qed_table_location)//'/hsokolov.table', exist=exists)
      IF (.NOT.exists) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to find QED tables in the ', &
              'directory "' // TRIM(qed_table_location) // '"'
        ENDDO
        CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      ENDIF
    ENDIF

    IF (use_qed) need_random_state = .TRUE.
#else
    IF (use_qed) THEN
      IF (rank .EQ. 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to set "use_qed=T" in the "qed" block.'
          WRITE(io,*) 'Please recompile with the -DPHOTONS preprocessor flag.'
        ENDDO
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    ENDIF
#endif

  END SUBROUTINE qed_deck_finalise



  SUBROUTINE qed_block_start

  END SUBROUTINE qed_block_start



  SUBROUTINE qed_block_end

  END SUBROUTINE qed_block_end



  FUNCTION qed_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    errcode = c_err_none
    IF (deck_state .EQ. c_ds_first) RETURN
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element, 'use_qed') .OR. str_cmp(element, 'qed')) THEN
      use_qed = as_logical_print(value, element, errcode)
      RETURN
    ENDIF

#ifdef PHOTONS
    IF (str_cmp(element, 'qed_start_time')) THEN
      qed_start_time = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'produce_photons')) THEN
      produce_photons = as_logical_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'use_radiation_reaction')) THEN
      use_radiation_reaction = as_logical_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'photon_energy_min') &
        .OR. str_cmp(element, 'min_photon_energy')) THEN
      photon_energy_min = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'produce_pairs')) THEN
      produce_pairs = as_logical_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'qed_table_location')) THEN
      qed_table_location = TRIM(ADJUSTL(value))
      RETURN
    ENDIF

    IF (str_cmp(element, 'photon_dynamics')) THEN
      photon_dynamics = as_logical_print(value, element, errcode)
      RETURN
    ENDIF

    errcode = c_err_unknown_element
#endif

  END FUNCTION qed_block_handle_element



  FUNCTION qed_block_check() RESULT(errcode)

    INTEGER :: errcode
#ifdef PHOTONS
    INTEGER :: io, iu
#endif

    errcode = c_err_none

#ifdef PHOTONS
    IF (produce_pairs .AND. .NOT. photon_dynamics) THEN
      IF (rank .EQ. 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'You cannot set photon_dynamics=F when ', &
            'produce_pairs=T. Without ', 'photon motion, pair ', &
            'creation will be incorrect.'
          WRITE(io,*) 'Code will terminate.'
        ENDDO
      ENDIF
      errcode = c_err_bad_value + c_err_terminate
    ENDIF
#endif

  END FUNCTION qed_block_check

END MODULE deck_qed_block
