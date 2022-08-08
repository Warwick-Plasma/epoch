! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE deck_k_alpha_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: k_alpha_deck_initialise, k_alpha_deck_finalise
  PUBLIC :: k_alpha_block_start, k_alpha_block_end
  PUBLIC :: k_alpha_block_handle_element, k_alpha_block_check

CONTAINS

  SUBROUTINE k_alpha_deck_initialise

#ifdef K_ALPHA
    IF (deck_state /= c_ds_first) RETURN
    use_k_alpha_recoil = .TRUE.
    use_k_alpha = .FALSE.
    k_alpha_photon_species = -1
#ifndef PHOTONS
#ifndef BREMSSTRAHLUNG
    photon_species = -1
#endif
#endif
    photon_energy_min_k_alpha = EPSILON(1.0_num)
    k_alpha_start_time = 0.0_num
    k_alpha_weight = 1.0_num
    produce_k_alpha_photons = .FALSE.
    k_alpha_photon_dynamics = .FALSE.
#endif

  END SUBROUTINE k_alpha_deck_initialise



  SUBROUTINE k_alpha_deck_finalise

    INTEGER :: io, iu
#ifdef K_ALPHA
    LOGICAL :: exists

    IF (deck_state == c_ds_first) RETURN

    IF (use_k_alpha) need_random_state = .TRUE.
#else
    IF (use_k_alpha) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to set "use_k_alpha=T" in the ', &
              '"k_alpha" block.'
          WRITE(io,*) 'Please recompile with the -DK_ALPHA ', &
              'preprocessor flag.'
        END DO
      END IF
      CALL abort_code(c_err_pp_options_missing)
    END IF
#endif

  END SUBROUTINE k_alpha_deck_finalise



  SUBROUTINE k_alpha_block_start

#ifdef K_ALPHA
    ! Enable by default if any k_alpha block is defined
    use_k_alpha = .TRUE.
#endif

  END SUBROUTINE k_alpha_block_start



  SUBROUTINE k_alpha_block_end

  END SUBROUTINE k_alpha_block_end



  FUNCTION k_alpha_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'use_k_alpha') &
        .OR. str_cmp(element, 'k_alpha') &
        .OR. str_cmp(element, 'enable')) THEN
      use_k_alpha = as_logical_print(value, element, errcode)
      RETURN
    END IF

#ifdef K_ALPHA
    IF (str_cmp(element, 'k_alpha_start_time') &
        .OR. str_cmp(element, 'start_time')) THEN
      k_alpha_start_time = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'produce_k_alpha_photons') &
        .OR. str_cmp(element, 'produce_photons')) THEN
      produce_k_alpha_photons = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_k_alpha_recoil') &
        .OR. str_cmp(element, 'use_radiation_reaction')) THEN
      use_k_alpha_recoil = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'photon_energy_min') &
        .OR. str_cmp(element, 'min_photon_energy') &
        .OR. str_cmp(element, 'photon_energy_min_k_alpha')) THEN
      photon_energy_min_k_alpha = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'auger_frac')) THEN
      auger_frac = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'photon_weight') &
        .OR. str_cmp(element, 'photon_weight_multiplier')) THEN
      k_alpha_weight = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'photon_dynamics') &
        .OR. str_cmp(element, 'k_alpha_photon_dynamics')) THEN
      k_alpha_photon_dynamics = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'sig_ka1_0')) THEN
      sig_ka1_0 = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'sig_ka1_1')) THEN
      sig_ka1_1 = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'sig_ka1_2')) THEN
      sig_ka1_2 = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'sig_ka1_3')) THEN
      sig_ka1_3 = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'sig_ka1_4')) THEN
      sig_ka1_4 = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'sig_ka2_0')) THEN
      sig_ka2_0 = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'sig_ka2_1')) THEN
      sig_ka2_1 = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'sig_ka2_2')) THEN
      sig_ka2_2 = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'sig_ka2_3')) THEN
      sig_ka2_3 = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'sig_ka2_4')) THEN
      sig_ka2_4 = as_real_print(value, element, errcode)
      RETURN
    END IF
    
    errcode = c_err_unknown_element
#endif

  END FUNCTION k_alpha_block_handle_element



  FUNCTION k_alpha_block_check() RESULT(errcode)

    INTEGER :: errcode
#ifdef K_ALPHA
    INTEGER :: io, iu
#endif

    errcode = c_err_none

#ifdef K_ALPHA
    IF (k_alpha_weight <= 0.0_num) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'You cannot set the photon_weight to less than or ', &
              'equal to zero. To prevent k_alpha photons from being ', &
              'emitted, set use_k_alpha = F.'
          WRITE(io,*) 'Code will terminate.'
        END DO
      END IF
      errcode = c_err_bad_value + c_err_terminate
    END IF

    IF (photon_weight > 1.0_num) THEN
      photon_weight = 1.0_num
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'You cannot set photon_weight > 1.0. This variable ', &
              'has been truncated to 1.0.'
        END DO
      END IF
    END IF
#endif

  END FUNCTION k_alpha_block_check

END MODULE deck_k_alpha_block
