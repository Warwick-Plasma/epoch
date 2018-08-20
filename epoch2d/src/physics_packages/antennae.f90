MODULE antennae

  USE shared_data
  USE evaluator
  IMPLICIT NONE

  LOGICAL :: any_antennae = .FALSE.

  TYPE :: antenna
    TYPE(primitive_stack) :: jx_expression, jy_expression, jz_expression
    TYPE(primitive_stack) :: ranges
    LOGICAL :: active = .FALSE.
  END TYPE antenna

  TYPE(antenna), DIMENSION(:), ALLOCATABLE :: antenna_list

  CONTAINS

  !> Initialise an antenna
  SUBROUTINE initialise_antenna(antenna_in)
    TYPE(antenna), INTENT(INOUT) :: antenna_in
    antenna_in%active = .TRUE.
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
  SUBROUTINE generate_antennae_currents(jx, jy, jz)

    REAL(num), DIMENSION(1-ng:, 1-ng:), INTENT(INOUT) :: jx, jy, jz
    TYPE(parameter_pack) :: parameters
    INTEGER :: iant, ix, iy, sz, err, nels
    REAL(num), DIMENSION(:), POINTER :: ranges

    IF (.NOT. ALLOCATED(antenna_list)) RETURN
    sz = SIZE(antenna_list)

    err = 0

    ranges => NULL()
    DO iant = 1, sz
      IF (.NOT. antenna_list(iant)%active) CYCLE
      CALL evaluate_and_return_all(antenna_list(iant)%ranges, nels, ranges, err)
      IF (err /= c_err_none .OR. nels /= c_ndims * 2) CYCLE
      DO iy = 1-ng, ny+ng
        IF (y(iy) < ranges(c_dir_y * 2 - 1) .OR. &
            y(iy) > ranges(c_dir_y * 2)) CYCLE
        DO ix = 1-ng, nx+ng
          IF (x(ix) < ranges(c_dir_x * 2 - 1) .OR. &
              x(ix) > ranges(c_dir_x * 2)) CYCLE
          parameters%pack_ix = ix
          parameters%pack_iy = iy
          IF(antenna_list(iant)%jx_expression%init) THEN
            jx(ix,iy) = jx(ix,iy) + evaluate_with_parameters(&
                antenna_list(iant)%jx_expression, &
                parameters, err)
          END IF
          IF(antenna_list(iant)%jy_expression%init) THEN
            jy(ix,iy) = jy(ix,iy) + evaluate_with_parameters(&
                antenna_list(iant)%jy_expression, &
                parameters, err)
          END IF
          IF(antenna_list(iant)%jz_expression%init) THEN
            jz(ix,iy) = jz(ix,iy) + evaluate_with_parameters(&
                antenna_list(iant)%jz_expression, &
                parameters, err)
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

