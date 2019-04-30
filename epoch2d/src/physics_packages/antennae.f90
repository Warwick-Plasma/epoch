MODULE antennae

  USE shared_data
  USE evaluator
  USE shunt

  IMPLICIT NONE

  LOGICAL :: any_antennae = .FALSE.

  TYPE antenna
    TYPE(primitive_stack) :: jx_expression, jy_expression, jz_expression
    TYPE(primitive_stack) :: ranges
    TYPE(primitive_stack) :: omega
    REAL(num) :: start_time, stop_time
    REAL(num) :: omega_value
    REAL(num) :: phase_history
    LOGICAL :: active = .FALSE.
  END TYPE antenna

  TYPE antenna_holder
    TYPE(antenna), POINTER :: contents
  END TYPE antenna_holder

  TYPE(antenna_holder), DIMENSION(:), ALLOCATABLE :: antenna_list

CONTAINS

  !> Initialise an antenna
  SUBROUTINE initialise_antenna(antenna_in)

    TYPE(antenna), INTENT(INOUT) :: antenna_in

    antenna_in%active = .TRUE.
    antenna_in%phase_history = 0.0_num
    antenna_in%start_time = -1.0_num
    antenna_in%stop_time = c_largest_number

  END SUBROUTINE initialise_antenna



  !> Routine to clean up an antenna after use
  SUBROUTINE finalise_antenna(antenna_in)

    TYPE(antenna), INTENT(INOUT) :: antenna_in

    antenna_in%active = .FALSE.

    IF (antenna_in%jx_expression%init) &
        CALL deallocate_stack(antenna_in%jx_expression)
    IF (antenna_in%jy_expression%init) &
        CALL deallocate_stack(antenna_in%jy_expression)
    IF (antenna_in%jz_expression%init) &
        CALL deallocate_stack(antenna_in%jz_expression)

    IF (antenna_in%ranges%init) &
        CALL deallocate_stack(antenna_in%ranges)
    IF (antenna_in%omega%init) &
        CALL deallocate_stack(antenna_in%omega)

  END SUBROUTINE finalise_antenna



  !> Add an antenna to the list of antennae
  SUBROUTINE add_antenna(antenna_in)

    TYPE(antenna), INTENT(IN), POINTER :: antenna_in
    TYPE(antenna_holder), DIMENSION(:), ALLOCATABLE :: temp
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
    antenna_list(sz)%contents => antenna_in

    IF (copyback) THEN
      antenna_list(1:sz-1) = temp
      DEALLOCATE(temp)
    END IF

  END SUBROUTINE add_antenna



  !> Get the currents from the antennae
  SUBROUTINE generate_antennae_currents

    TYPE(parameter_pack) :: parameters
    INTEGER :: iant, ix, iy, sz, err, nels
    REAL(num), DIMENSION(:), ALLOCATABLE :: ranges
    REAL(num) :: oscil_dat
    LOGICAL :: use_ranges
    TYPE(antenna), POINTER :: current_antenna

    IF (.NOT. ALLOCATED(antenna_list)) RETURN
    sz = SIZE(antenna_list)

    err = c_err_none

    DO iant = 1, sz
      current_antenna => antenna_list(iant)%contents
      IF (.NOT. current_antenna%active) CYCLE
      IF (time < current_antenna%start_time &
          .OR. time > current_antenna%stop_time) CYCLE
      use_ranges = current_antenna%ranges%init
      IF (current_antenna%ranges%init) THEN
        CALL evaluate_and_return_all(current_antenna%ranges, nels, ranges, err)
        IF (err /= c_err_none .OR. nels /= c_ndims * 2) CYCLE
      ELSE
        ALLOCATE(ranges(c_ndims*2))
      END IF
      IF (current_antenna%omega%init) THEN
        IF (current_antenna%omega%is_time_varying) THEN
          current_antenna%phase_history = current_antenna%phase_history &
              + evaluate(current_antenna%omega, err) * dt
        ELSE
          current_antenna%phase_history = current_antenna%omega_value * time
        END IF
        oscil_dat = SIN(current_antenna%phase_history)
      ELSE
        oscil_dat = 1.0_num
      END IF
      DO iy = 1-ng, ny+ng
        IF ((y(iy) < ranges(c_dir_y * 2 - 1) &
            .OR. y(iy) > ranges(c_dir_y * 2)) .AND. use_ranges) CYCLE
        parameters%pack_iy = iy
        DO ix = 1-ng, nx+ng
          IF ((x(ix) < ranges(c_dir_x * 2 - 1) &
              .OR. x(ix) > ranges(c_dir_x * 2)) .AND. use_ranges) CYCLE
          parameters%pack_ix = ix
          IF (current_antenna%jx_expression%init) THEN
            jx(ix,iy) = jx(ix,iy) + evaluate_with_parameters(&
                current_antenna%jx_expression, &
                parameters, err) * oscil_dat
          END IF
          IF (current_antenna%jy_expression%init) THEN
            jy(ix,iy) = jy(ix,iy) + evaluate_with_parameters(&
                current_antenna%jy_expression, &
                parameters, err) * oscil_dat
          END IF
          IF (current_antenna%jz_expression%init) THEN
            jz(ix,iy) = jz(ix,iy) + evaluate_with_parameters(&
                current_antenna%jz_expression, &
                parameters, err) * oscil_dat
          END IF
        END DO
      END DO
    END DO

    DEALLOCATE(ranges)

  END SUBROUTINE generate_antennae_currents



  !> Deallocate antennae
  SUBROUTINE deallocate_antennae

    INTEGER :: iant

    IF (.NOT. ALLOCATED(antenna_list)) RETURN

    DO iant = 1, SIZE(antenna_list)
      CALL finalise_antenna(antenna_list(iant)%contents)
      DEALLOCATE(antenna_list(iant)%contents)
    END DO

    DEALLOCATE(antenna_list)

  END SUBROUTINE deallocate_antennae

END MODULE antennae
