MODULE antennae

  USE shared_data
  USE evaluator
  USE shunt
  IMPLICIT NONE

  LOGICAL :: any_antennae = .FALSE.

  TYPE :: antenna
    TYPE(primitive_stack) :: jx_expression, jy_expression, jz_expression
    TYPE(primitive_stack) :: ranges
    TYPE(primitive_stack) :: omega
    REAL(num) :: start_time, stop_time
    REAL(num) :: omega_value
    REAL(num) :: phase_history
    LOGICAL :: active = .FALSE.
  END TYPE antenna

  TYPE(antenna), DIMENSION(:), ALLOCATABLE :: antenna_list

  CONTAINS

  !> Initialise an antenna
  SUBROUTINE initialise_antenna(antenna_in)
    TYPE(antenna), INTENT(INOUT) :: antenna_in
    antenna_in%active = .TRUE.

    IF (antenna_in%jx_expression%init) &
        CALL deallocate_stack(antenna_in%jx_expression)
    IF (antenna_in%jy_expression%init) &
        CALL deallocate_stack(antenna_in%jy_expression)
    IF (antenna_in%jz_expression%init) &
        CALL deallocate_stack(antenna_in%jz_expression)

    IF (antenna_in%ranges%init) &
        CALL deallocate_stack(antenna_in%ranges)

    antenna_in%phase_history = 0.0_num
    antenna_in%start_time = -1.0_num
    antenna_in%stop_time = c_largest_number

  END SUBROUTINE initialise_antenna



  !> Add an antenna to the list of antennae
  SUBROUTINE add_antenna(antenna_in)
    TYPE(antenna), INTENT(IN) :: antenna_in
    TYPE(antenna), DIMENSION(:), ALLOCATABLE :: temp
    LOGICAL :: copyback
    INTEGER :: sz

    any_antennae = .TRUE.

    copyback = ALLOCATED(antenna_list)
    sz = 0
    IF (copyback) THEN
      sz = SIZE(antenna_list)
      ALLOCATE(temp(sz))
      temp = antenna_list
      DEALLOCATE(antenna_list)
    END IF

    sz = sz + 1
    ALLOCATE(antenna_list(sz))
    antenna_list(sz) = antenna_in

    IF (copyback) THEN
      antenna_list(1:sz-1) = temp
      DEALLOCATE(temp)
    END IF

  END SUBROUTINE add_antenna



  !> Get the currents from the antennae
  SUBROUTINE generate_antennae_currents

    TYPE(parameter_pack) :: parameters
    INTEGER :: iant, ix, iy, sz, err, nels
    REAL(num), DIMENSION(:), POINTER :: ranges
    REAL(num) :: oscil_dat
    LOGICAL :: use_ranges

    IF (.NOT. ALLOCATED(antenna_list)) RETURN
    sz = SIZE(antenna_list)

    err = c_err_none

    ranges => NULL()
    DO iant = 1, sz
      IF (.NOT. antenna_list(iant)%active) CYCLE
      IF (time < antenna_list(iant)%start_time &
          .OR. time > antenna_list(iant)%stop_time) CYCLE
      use_ranges = antenna_list(iant)%ranges%init
      IF (antenna_list(iant)%ranges%init) THEN
        CALL evaluate_and_return_all(antenna_list(iant)%ranges, nels, ranges, &
            err)
        IF (err /= c_err_none .OR. nels /= c_ndims * 2) CYCLE
      ELSE
        ALLOCATE(ranges(c_ndims*2))
      END IF
      IF (antenna_list(iant)%omega%init) THEN
        IF (antenna_list(iant)%omega%is_time_varying) THEN
          antenna_list(iant)%phase_history = &
              antenna_list(iant)%phase_history + &
              evaluate(antenna_list(iant)%omega, err) * dt
        ELSE
          antenna_list(iant)%phase_history = antenna_list(iant)%omega_value &
              * time
        END IF
        oscil_dat = SIN(antenna_list(iant)%phase_history)
      ELSE
        oscil_dat = 1.0_num
      END IF
      DO iy = 1-ng, ny+ng
        IF ((y(iy) < ranges(c_dir_y * 2 - 1) .OR. &
            y(iy) > ranges(c_dir_y * 2)) .AND. use_ranges) CYCLE
        parameters%pack_iy = iy
        DO ix = 1-ng, nx+ng
          IF ((x(ix) < ranges(c_dir_x * 2 - 1) .OR. &
              x(ix) > ranges(c_dir_x * 2)) .AND. use_ranges) CYCLE
          parameters%pack_ix = ix
          parameters%pack_iy = iy
          IF(antenna_list(iant)%jx_expression%init) THEN
            jx(ix,iy) = jx(ix,iy) + evaluate_with_parameters(&
                antenna_list(iant)%jx_expression, &
                parameters, err) * oscil_dat
          END IF
          IF(antenna_list(iant)%jy_expression%init) THEN
            jy(ix,iy) = jy(ix,iy) + evaluate_with_parameters(&
                antenna_list(iant)%jy_expression, &
                parameters, err) * oscil_dat
          END IF
          IF(antenna_list(iant)%jz_expression%init) THEN
            jz(ix,iy) = jz(ix,iy) + evaluate_with_parameters(&
                antenna_list(iant)%jz_expression, &
                parameters, err) * oscil_dat
          END IF
        END DO
      END DO
    END DO

    DEALLOCATE(ranges)

  END SUBROUTINE generate_antennae_currents



  !> Deallocate antennae
  SUBROUTINE deallocate_antennae
    IF (ALLOCATED(antenna_list)) DEALLOCATE(antenna_list)
  END SUBROUTINE deallocate_antennae

END MODULE antennae

