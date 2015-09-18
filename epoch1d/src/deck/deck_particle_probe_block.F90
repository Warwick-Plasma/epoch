MODULE deck_particle_probe_block

#ifdef NO_PARTICLE_PROBES
  IMPLICIT NONE

CONTAINS

  SUBROUTINE probe_deck_dummy

  END SUBROUTINE probe_deck_dummy

#else
  USE strings_advanced
  USE probes

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: probe_deck_initialise, probe_deck_finalise
  PUBLIC :: probe_block_start, probe_block_end
  PUBLIC :: probe_block_handle_element, probe_block_check

  TYPE(particle_probe), POINTER :: working_probe
  LOGICAL :: got_name, got_point, got_normal

CONTAINS

  SUBROUTINE probe_deck_initialise

  END SUBROUTINE probe_deck_initialise



  SUBROUTINE probe_deck_finalise

  END SUBROUTINE probe_deck_finalise



  SUBROUTINE probe_block_start

    IF (deck_state == c_ds_first) RETURN

    ALLOCATE(working_probe)
    CALL init_probe(working_probe)
    got_name = .FALSE.
    got_point = .FALSE.
    got_normal = .FALSE.

  END SUBROUTINE probe_block_start



  SUBROUTINE probe_block_end

    LOGICAL :: discard
    INTEGER :: io, iu, ierr

    IF (deck_state == c_ds_first) RETURN

    IF (.NOT.got_name) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) '"probe" block does not have a "name" entry.'
        ENDDO
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, c_err_required_element_not_set, ierr)
    ENDIF

    discard = .NOT.(got_point .AND. got_normal)

    IF (discard) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Position of probe "' // TRIM(working_probe%name) &
              // '" ', 'not fully specified. ', 'It will be discarded.'
          WRITE(io,*) 'Both "point" and "normal" are required.'
        ENDDO
      ENDIF

      DEALLOCATE(working_probe)
      NULLIFY(working_probe)
    ELSE
      ! Normalise the normal. Not really necessary but doesn't hurt.
      working_probe%normal = SIGN(1.0_num, working_probe%normal)

      CALL attach_probe(working_probe)
    ENDIF

  END SUBROUTINE probe_block_end



  FUNCTION probe_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode, ispecies, io, iu

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    ! get particle probe diagnostics (rolling total of all particles which
    ! pass through a given region of real space (defined by a point on a plane
    ! and the normal to that plane.
    IF (str_cmp(element, 'dumpmask') .OR. str_cmp(element, 'dump')) THEN
      working_probe%dumpmask = as_integer_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'point') .OR. str_cmp(element, 'probe_point')) THEN
      got_point = .TRUE.
      working_probe%point = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'normal')) THEN
      got_normal = .TRUE.
      working_probe%normal = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'left_to_right')) THEN
      got_normal = as_logical_print(value, element, errcode)
      IF (got_normal) THEN
        working_probe%normal = 1.0_num
      ELSE
        working_probe%normal = -1.0_num
      ENDIF
      got_normal = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'include_species') &
        .OR. str_cmp(element, 'probe_species')) THEN
      ispecies = as_integer_print(value, element, errcode)
      IF (errcode == c_err_none) THEN
        IF (ispecies > 0 .AND. ispecies <= n_species) THEN
          working_probe%use_species(ispecies) = .TRUE.
        ELSE
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Unable to attach probe to non existant species ', &
                  ispecies
            ENDDO
          ENDIF
          errcode = c_err_bad_value
        ENDIF
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'ek_min')) THEN
      working_probe%ek_min = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'ek_max')) THEN
      working_probe%ek_max = as_real_print(value, element, errcode)
      IF (working_probe%ek_max < 0) working_probe%ek_max = HUGE(1.0_num)
      RETURN
    ENDIF

    IF (str_cmp(element, 'name')) THEN
      got_name = .TRUE.
      working_probe%name = TRIM(value)
      RETURN
    ENDIF

    errcode = c_err_unknown_element

  END FUNCTION probe_block_handle_element



  FUNCTION probe_block_check() RESULT(errcode)

    INTEGER :: errcode
    errcode = c_err_none

  END FUNCTION probe_block_check
#endif

END MODULE deck_particle_probe_block
