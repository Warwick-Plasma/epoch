MODULE deck_subset_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: subset_deck_initialise, subset_deck_finalise
  PUBLIC :: subset_block_start, subset_block_end
  PUBLIC :: subset_block_handle_element, subset_block_check

  INTEGER :: subset_id, current_block
  INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0
  CHARACTER(LEN=string_length), DIMENSION(:), POINTER :: subset_names
  INTEGER, DIMENSION(:), POINTER :: subset_blocks
  LOGICAL :: got_name
  INTEGER :: check_block = c_err_none

CONTAINS

  SUBROUTINE subset_deck_initialise

    INTEGER :: i

    current_block = 0
    IF (deck_state .EQ. c_ds_first) THEN
      n_subsets = 0
      ALLOCATE(subset_names(4))
      ALLOCATE(subset_blocks(4))
    ELSE
      DO i = 1, n_subsets
        ALLOCATE(subset_list(i)%use_species(n_species))
        subset_list(i)%use_species = .FALSE.
      ENDDO
    ENDIF

  END SUBROUTINE subset_deck_initialise



  SUBROUTINE subset_deck_finalise

    INTEGER :: i

    IF (deck_state .EQ. c_ds_first) THEN
      CALL setup_subsets

      DO i = 1, n_subsets
        subset_list(i)%name = subset_names(i)
      ENDDO
      DEALLOCATE(subset_names)
    ELSE
      DEALLOCATE(subset_blocks)
    ENDIF

  END SUBROUTINE subset_deck_finalise



  SUBROUTINE subset_block_start

    current_block = current_block + 1
    got_name = .FALSE.
    IF (deck_state .EQ. c_ds_first) RETURN
    subset_id = subset_blocks(current_block)
    offset = 0

  END SUBROUTINE subset_block_start



  SUBROUTINE subset_block_end

    CHARACTER(LEN=8) :: id_string
    INTEGER :: io

    IF (.NOT.got_name) THEN
      IF (rank .EQ. 0) THEN
        CALL integer_as_string(current_block, id_string)
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Restriction block number ', TRIM(id_string), &
              ' has no "name" element.'
        ENDDO
      ENDIF

      check_block = c_err_missing_elements
    ENDIF

  END SUBROUTINE subset_block_end



  FUNCTION subset_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    INTEGER :: io, ispecies

    errcode = c_err_none
    IF (value .EQ. blank .OR. element .EQ. blank) RETURN

    IF (str_cmp(element, 'name')) THEN
      IF (got_name) THEN
        errcode = c_err_preset_element
        RETURN
      ENDIF
      got_name = .TRUE.
      IF (deck_state .NE. c_ds_first) RETURN
      CALL grow_array(subset_blocks, current_block)
      subset_blocks(current_block) = subset_number_from_name(value)
      RETURN
    ENDIF

    IF (deck_state .EQ. c_ds_first) RETURN

    IF (str_cmp(element, 'random_fraction')) THEN
      subset_list(subset_id)%random_fraction = as_real(value, errcode)
      subset_list(subset_id)%use_random = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'gamma_min')) THEN
      subset_list(subset_id)%gamma_min = as_real(value, errcode)
      subset_list(subset_id)%use_gamma_min = .TRUE.
      subset_list(subset_id)%use_gamma = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'gamma_max')) THEN
      subset_list(subset_id)%gamma_max = as_real(value, errcode)
      subset_list(subset_id)%use_gamma_max = .TRUE.
      subset_list(subset_id)%use_gamma = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'x_min')) THEN
      subset_list(subset_id)%x_min = as_real(value, errcode)
      subset_list(subset_id)%use_x_min = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'x_max')) THEN
      subset_list(subset_id)%x_max = as_real(value, errcode)
      subset_list(subset_id)%use_x_max = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'y_min')) THEN
      subset_list(subset_id)%y_min = as_real(value, errcode)
      subset_list(subset_id)%use_y_min = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'y_max')) THEN
      subset_list(subset_id)%y_max = as_real(value, errcode)
      subset_list(subset_id)%use_y_max = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'px_min')) THEN
      subset_list(subset_id)%px_min = as_real(value, errcode)
      subset_list(subset_id)%use_px_min = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'px_max')) THEN
      subset_list(subset_id)%px_max = as_real(value, errcode)
      subset_list(subset_id)%use_px_max = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'py_min')) THEN
      subset_list(subset_id)%py_min = as_real(value, errcode)
      subset_list(subset_id)%use_py_min = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'py_max')) THEN
      subset_list(subset_id)%py_max = as_real(value, errcode)
      subset_list(subset_id)%use_py_max = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'pz_min')) THEN
      subset_list(subset_id)%pz_min = as_real(value, errcode)
      subset_list(subset_id)%use_pz_min = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'pz_max')) THEN
      subset_list(subset_id)%pz_max = as_real(value, errcode)
      subset_list(subset_id)%use_pz_max = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'weight_min')) THEN
      subset_list(subset_id)%weight_min = as_real(value, errcode)
      subset_list(subset_id)%use_weight_min = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'weight_max')) THEN
      subset_list(subset_id)%weight_max = as_real(value, errcode)
      subset_list(subset_id)%use_weight_max = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'charge_min')) THEN
      subset_list(subset_id)%charge_min = as_real(value, errcode)
      subset_list(subset_id)%use_charge_min = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'charge_max')) THEN
      subset_list(subset_id)%charge_max = as_real(value, errcode)
      subset_list(subset_id)%use_charge_max = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'mass_min')) THEN
      subset_list(subset_id)%mass_min = as_real(value, errcode)
      subset_list(subset_id)%use_mass_min = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'mass_max')) THEN
      subset_list(subset_id)%mass_max = as_real(value, errcode)
      subset_list(subset_id)%use_mass_max = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'id_min')) THEN
      subset_list(subset_id)%id_min = as_integer(value, errcode)
      subset_list(subset_id)%use_id_min = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'id_max')) THEN
      subset_list(subset_id)%id_max = as_integer(value, errcode)
      subset_list(subset_id)%use_id_max = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'dumpmask')) THEN
      subset_list(subset_id)%mask = as_integer(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'include_species')) THEN
      ispecies = as_integer(value, errcode)
      IF (errcode .EQ. c_err_none) THEN
        IF (ispecies .GT. 0 .AND. ispecies .LE. n_species) THEN
          subset_list(subset_id)%use_species(ispecies) = .TRUE.
        ELSE
          IF (rank .EQ. 0) THEN
            DO io = stdout, du, du - stdout ! Print to stdout and to file
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Unable to apply subset to non existant ', &
                  'species ', ispecies
            ENDDO
          ENDIF
          errcode = c_err_bad_value
        ENDIF
      ENDIF
      RETURN
    ENDIF

    errcode = c_err_unknown_element

  END FUNCTION subset_block_handle_element



  FUNCTION subset_block_check() RESULT(errcode)

    INTEGER :: errcode

    errcode = check_block

  END FUNCTION subset_block_check



  FUNCTION subset_number_from_name(name)

    CHARACTER(*), INTENT(IN) :: name
    INTEGER :: subset_number_from_name
    INTEGER :: i

    DO i = 1, n_subsets
      IF (str_cmp(name, subset_names(i))) THEN
        subset_number_from_name = i
        RETURN
      ENDIF
    ENDDO
    n_subsets = n_subsets + 1
    CALL grow_array(subset_names, n_subsets)
    subset_names(n_subsets) = TRIM(name)
    subset_number_from_name = n_subsets

  END FUNCTION subset_number_from_name



  SUBROUTINE setup_subsets

    INTEGER :: i

    ALLOCATE(subset_list(n_subsets))

    DO i = 1, n_subsets
      subset_list(i)%name = blank
      subset_list(i)%use_random     = .FALSE.
      subset_list(i)%use_gamma      = .FALSE.
      subset_list(i)%use_gamma_min  = .FALSE.
      subset_list(i)%use_gamma_max  = .FALSE.
      subset_list(i)%use_x_min      = .FALSE.
      subset_list(i)%use_x_max      = .FALSE.
      subset_list(i)%use_y_min      = .FALSE.
      subset_list(i)%use_y_max      = .FALSE.
      subset_list(i)%use_px_min     = .FALSE.
      subset_list(i)%use_px_max     = .FALSE.
      subset_list(i)%use_py_min     = .FALSE.
      subset_list(i)%use_py_max     = .FALSE.
      subset_list(i)%use_pz_min     = .FALSE.
      subset_list(i)%use_pz_max     = .FALSE.
      subset_list(i)%use_weight_min = .FALSE.
      subset_list(i)%use_weight_max = .FALSE.
      subset_list(i)%use_charge_min = .FALSE.
      subset_list(i)%use_charge_max = .FALSE.
      subset_list(i)%use_mass_min   = .FALSE.
      subset_list(i)%use_mass_max   = .FALSE.
      subset_list(i)%random_fraction = 0.0_num
      subset_list(i)%gamma_min  = -HUGE(1.0_num)
      subset_list(i)%gamma_max  =  HUGE(1.0_num)
      subset_list(i)%x_min      = -HUGE(1.0_num)
      subset_list(i)%x_max      =  HUGE(1.0_num)
      subset_list(i)%y_min      = -HUGE(1.0_num)
      subset_list(i)%y_max      =  HUGE(1.0_num)
      subset_list(i)%px_min     = -HUGE(1.0_num)
      subset_list(i)%px_max     =  HUGE(1.0_num)
      subset_list(i)%py_min     = -HUGE(1.0_num)
      subset_list(i)%py_max     =  HUGE(1.0_num)
      subset_list(i)%pz_min     = -HUGE(1.0_num)
      subset_list(i)%pz_max     =  HUGE(1.0_num)
      subset_list(i)%weight_min = -HUGE(1.0_num)
      subset_list(i)%weight_max =  HUGE(1.0_num)
      subset_list(i)%charge_min = -HUGE(1.0_num)
      subset_list(i)%charge_max =  HUGE(1.0_num)
      subset_list(i)%mass_min   = -HUGE(1.0_num)
      subset_list(i)%mass_max   =  HUGE(1.0_num)
      subset_list(i)%mask = c_io_always
      ALLOCATE(subset_list(i)%dumpmask(n_io_blocks,num_vars_to_dump))
      subset_list(i)%dumpmask = c_io_never
    ENDDO

  END SUBROUTINE setup_subsets

END MODULE deck_subset_block
