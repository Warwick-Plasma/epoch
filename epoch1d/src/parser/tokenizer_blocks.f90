! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2014      Stephan Kuschel <stephan.kuschel@gmail.com>
! Copyright (C) 2009      Chris Brady <C.S.Brady@warwick.ac.uk>
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

MODULE tokenizer_blocks

  USE strings

  IMPLICIT NONE

  SAVE

  INTEGER, PARAMETER :: c_max_new_func = 256
  INTEGER :: n_new_func = 0
  TYPE(string_type), DIMENSION(c_max_new_func) :: new_func_name
  INTEGER, DIMENSION(c_max_new_func) :: new_func_code

  INTEGER, PARAMETER :: c_max_new_const = 256
  INTEGER :: n_new_constant = 0
  TYPE(string_type), DIMENSION(c_max_new_const) :: new_constant_name
  INTEGER, DIMENSION(c_max_new_const) :: new_constant_code
  INTEGER :: last_block_type, tokenize_stagger

CONTAINS

  ! Functions to register new functions and constants
  FUNCTION register_function(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: register_function

    IF (n_new_func == c_max_new_func) THEN
      register_function = -1
      RETURN
    ENDIF
    n_new_func = n_new_func+1
    new_func_name(n_new_func)%value = name
    new_func_code(n_new_func) = c_func_custom_lowbound+n_new_func
    register_function = c_func_custom_lowbound+n_new_func

  END FUNCTION register_function



  FUNCTION register_constant(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: register_constant

    IF (n_new_constant == c_max_new_const) THEN
      register_constant = -1
      RETURN
    ENDIF
    n_new_constant = n_new_constant+1
    new_constant_name(n_new_constant)%value = name
    new_constant_code(n_new_constant) = c_const_custom_lowbound+n_new_constant
    register_constant = c_const_custom_lowbound+n_new_constant

  END FUNCTION register_constant



  FUNCTION as_species(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_species
    INTEGER :: i

    as_species = c_prc_not_this_type
    IF (.NOT. ASSOCIATED(species_list)) RETURN

    DO i = 1, n_species
      IF (str_cmp(name, species_list(i)%name)) THEN
        as_species = i
        RETURN
      ENDIF
    ENDDO

  END FUNCTION as_species



  FUNCTION as_subset(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_subset
    INTEGER :: i

    as_subset = c_prc_not_this_type
    IF (.NOT. ASSOCIATED(subset_list)) RETURN

    DO i = 1, n_subsets
      IF (str_cmp(name, subset_list(i)%name)) THEN
        as_subset = i
        RETURN
      ENDIF
    ENDDO

  END FUNCTION as_subset



  FUNCTION as_constant(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_constant
    INTEGER :: i

    as_constant = c_prc_not_this_type

    ! Numeric constants
    IF (str_cmp(name, 'pi')) as_constant = c_const_pi
    IF (str_cmp(name, 'kb')) as_constant = c_const_kb
    IF (str_cmp(name, 'me')) as_constant = c_const_me
    IF (str_cmp(name, 'qe')) as_constant = c_const_qe
    IF (str_cmp(name, 'c')) as_constant = c_const_c
    IF (str_cmp(name, 'eps0') .OR. str_cmp(name, 'epsilon0') &
        .OR. str_cmp(name, 'epsilonnought')) as_constant = c_const_eps0
    IF (str_cmp(name, 'mu0') &
        .OR. str_cmp(name, 'munought')) as_constant = c_const_mu0
    IF (str_cmp(name, 'ev')) as_constant = c_const_ev
    IF (str_cmp(name, 'kev')) as_constant = c_const_kev
    IF (str_cmp(name, 'mev')) as_constant = c_const_mev
    IF (str_cmp(name, 'milli')) as_constant = c_const_milli
    IF (str_cmp(name, 'micro') &
        .OR. str_cmp(name, 'micron') &
        .OR. str_cmp(name, 'cm3') &
        .OR. str_cmp(name, 'cc')) as_constant = c_const_micro
    IF (str_cmp(name, 'nano')) as_constant = c_const_nano
    IF (str_cmp(name, 'pico')) as_constant = c_const_pico
    IF (str_cmp(name, 'femto')) as_constant = c_const_femto
    IF (str_cmp(name, 'atto')) as_constant = c_const_atto
    IF (str_cmp(name, 'time')) as_constant = c_const_time
    IF (str_cmp(name, 'x'))  as_constant = c_const_x
    IF (str_cmp(name, 'xb')) as_constant = c_const_xb
    IF (str_cmp(name, 'ix')) as_constant = c_const_ix
    IF (str_cmp(name, 'nx')) as_constant = c_const_nx
    IF (str_cmp(name, 'dx')) as_constant = c_const_dx
    IF (str_cmp(name, 'x_min') .OR. str_cmp(name, 'x_start')) &
        as_constant = c_const_x_min
    IF (str_cmp(name, 'x_max') .OR. str_cmp(name, 'x_end')) &
        as_constant = c_const_x_max
    IF (str_cmp(name, 'lengthx') .OR. str_cmp(name, 'length_x')) &
        as_constant = c_const_lx
    IF (str_cmp(name, 'nprocx')) as_constant = c_const_nprocx
    IF (str_cmp(name, 'nsteps')) as_constant = c_const_nsteps
    IF (str_cmp(name, 't_end')) as_constant = c_const_t_end
    IF (str_cmp(name, 'ndims')) as_constant = c_const_ndims

    ! Dumpmask constants
    IF (str_cmp(name, 'never')) as_constant = c_const_io_never
    IF (str_cmp(name, 'always')) as_constant = c_const_io_always
    IF (str_cmp(name, 'full')) as_constant = c_const_io_full
    IF (str_cmp(name, 'restartable')) as_constant = c_const_io_restartable
    IF (str_cmp(name, 'restart')) as_constant = c_const_io_restartable
    IF (str_cmp(name, 'average')) as_constant = c_const_io_average
    IF (str_cmp(name, 'snapshot')) as_constant = c_const_io_snapshot
    IF (str_cmp(name, 'species')) as_constant = c_const_io_species
    IF (str_cmp(name, 'no_sum')) as_constant = c_const_io_no_sum
    IF (str_cmp(name, 'single')) as_constant = c_const_io_dump_single
    IF (str_cmp(name, 'average_single')) as_constant = c_const_io_average_single

    ! Distribution function constants
    IF (str_cmp(name, 'dir_x')) as_constant = c_const_dir_x
    IF (str_cmp(name, 'dir_px')) as_constant = c_const_dir_px
    IF (str_cmp(name, 'dir_py')) as_constant = c_const_dir_py
    IF (str_cmp(name, 'dir_pz')) as_constant = c_const_dir_pz
    IF (str_cmp(name, 'dir_en')) as_constant = c_const_dir_en
    IF (str_cmp(name, 'dir_gamma_m1')) as_constant = c_const_dir_gamma_m1
    IF (str_cmp(name, 'dir_xy_angle')) as_constant = c_const_dir_xy_angle
    IF (str_cmp(name, 'dir_yz_angle')) as_constant = c_const_dir_yz_angle
    IF (str_cmp(name, 'dir_zx_angle')) as_constant = c_const_dir_zx_angle

    IF (as_constant == c_const_x .AND. stagger(c_dir_x,tokenize_stagger)) &
        as_constant = c_const_xb

    IF (str_cmp(name, 'yee')) as_constant = c_const_maxwell_solver_yee
    IF (str_cmp(name, 'cowan')) as_constant = c_const_maxwell_solver_cowan
    IF (str_cmp(name, 'pukhov')) as_constant = c_const_maxwell_solver_pukhov
    IF (str_cmp(name, 'lehe_x')) as_constant = c_const_maxwell_solver_lehe
    IF (str_cmp(name, 'custom')) as_constant = c_const_maxwell_solver_custom

    ! User submitted constant using 'Register'
    DO i = 1, n_new_constant
      IF (str_cmp(TRIM(name), TRIM(new_constant_name(i)%value))) &
          as_constant = new_constant_code(i)
    ENDDO

  END FUNCTION as_constant



  FUNCTION as_default_constant(name) RESULT(as_constant)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_constant
    INTEGER :: io, iu
    LOGICAL, SAVE :: warn = .TRUE.

    as_constant = c_prc_not_this_type

    IF (str_cmp(name, 'y'))  as_constant = c_const_y
    IF (str_cmp(name, 'z'))  as_constant = c_const_z
    IF (str_cmp(name, 'yb')) as_constant = c_const_yb
    IF (str_cmp(name, 'zb')) as_constant = c_const_zb
    IF (str_cmp(name, 'iy')) as_constant = c_const_iy
    IF (str_cmp(name, 'iz')) as_constant = c_const_iz
    IF (str_cmp(name, 'ny')) as_constant = c_const_ny
    IF (str_cmp(name, 'nz')) as_constant = c_const_nz
    IF (str_cmp(name, 'dy')) as_constant = c_const_dy
    IF (str_cmp(name, 'dz')) as_constant = c_const_dz
    IF (str_cmp(name, 'y_min') .OR. str_cmp(name, 'y_start')) &
        as_constant = c_const_y_min
    IF (str_cmp(name, 'z_min') .OR. str_cmp(name, 'z_start')) &
        as_constant = c_const_z_min
    IF (str_cmp(name, 'y_max') .OR. str_cmp(name, 'y_end')) &
        as_constant = c_const_y_max
    IF (str_cmp(name, 'z_max') .OR. str_cmp(name, 'z_end')) &
        as_constant = c_const_z_max
    IF (str_cmp(name, 'lengthy') .OR. str_cmp(name, 'length_y')) &
        as_constant = c_const_ly
    IF (str_cmp(name, 'lengthz') .OR. str_cmp(name, 'length_z')) &
        as_constant = c_const_lz
    IF (str_cmp(name, 'r_xy')) as_constant = c_const_r_xy
    IF (str_cmp(name, 'r_yz')) as_constant = c_const_r_yz
    IF (str_cmp(name, 'r_xz')) as_constant = c_const_r_xz
    IF (str_cmp(name, 'r_xyz')) as_constant = c_const_r_xyz
    IF (str_cmp(name, 'nprocy')) as_constant = c_const_nprocy
    IF (str_cmp(name, 'nprocz')) as_constant = c_const_nprocz
    IF (str_cmp(name, 'dir_y')) as_constant = c_const_dir_y
    IF (str_cmp(name, 'dir_z')) as_constant = c_const_dir_z

    IF (as_constant == c_const_y .AND. stagger(c_dir_y,tokenize_stagger)) &
        as_constant = c_const_yb
    IF (as_constant == c_const_z .AND. stagger(c_dir_z,tokenize_stagger)) &
        as_constant = c_const_zb

    IF (warn .AND. as_constant /= c_prc_not_this_type) THEN
      warn = .FALSE.
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'A default value (set y and z to 0) was used for the ', &
              'constant ', '"' // TRIM(name) // '"'
          WRITE(io,*)
        ENDDO
      ENDIF
    ENDIF

  END FUNCTION as_default_constant



  FUNCTION as_deck_constant(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_deck_constant
    INTEGER :: i

    as_deck_constant = 0

    DO i = 1, n_deck_constants
      IF (str_cmp(TRIM(name), TRIM(deck_constant_list(i)%name))) THEN
        as_deck_constant = i
        RETURN
      ENDIF
    ENDDO

  END FUNCTION as_deck_constant



  FUNCTION as_function(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_function
    INTEGER :: i

    as_function = c_prc_not_this_type

    IF (str_cmp(name, '-'))     as_function = c_func_neg
    IF (str_cmp(name, 'abs'))   as_function = c_func_abs
    IF (str_cmp(name, 'floor')) as_function = c_func_floor
    IF (str_cmp(name, 'ceil'))  as_function = c_func_ceil
    IF (str_cmp(name, 'nint'))  as_function = c_func_nint
    IF (str_cmp(name, 'sqrt'))  as_function = c_func_sqrt
    IF (str_cmp(name, 'sin'))   as_function = c_func_sine
    IF (str_cmp(name, 'cos'))   as_function = c_func_cosine
    IF (str_cmp(name, 'tan'))   as_function = c_func_tan
    IF (str_cmp(name, 'asin'))  as_function = c_func_arcsine
    IF (str_cmp(name, 'acos'))  as_function = c_func_arccosine
    IF (str_cmp(name, 'atan'))  as_function = c_func_arctan
    IF (str_cmp(name, 'sinh'))  as_function = c_func_sinh
    IF (str_cmp(name, 'cosh'))  as_function = c_func_cosh
    IF (str_cmp(name, 'tanh'))  as_function = c_func_tanh
    IF (str_cmp(name, 'exp'))   as_function = c_func_exp
    IF (str_cmp(name, 'loge'))  as_function = c_func_loge
    IF (str_cmp(name, 'log10')) as_function = c_func_log10
    IF (str_cmp(name, 'log_base')) as_function = c_func_log_base
    IF (str_cmp(name, 'gauss')) as_function = c_func_gauss
    IF (str_cmp(name, 'semigauss')) as_function = c_func_semigauss
    IF (str_cmp(name, 'supergauss')) as_function = c_func_supergauss
    IF (str_cmp(name, 'interpolate')) as_function = c_func_interpolate
    IF (str_cmp(name, 'if')) as_function = c_func_if

    IF (str_cmp(name, 'density') .OR. str_cmp(name, 'rho')) &
        as_function = c_func_rho
    IF (str_cmp(name, 'temp_x') .OR. str_cmp(name, 'temp_x_k')) &
        as_function = c_func_tempx
    IF (str_cmp(name, 'temp_y') .OR. str_cmp(name, 'temp_y_k')) &
        as_function = c_func_tempy
    IF (str_cmp(name, 'temp_z') .OR. str_cmp(name, 'temp_z_k')) &
        as_function = c_func_tempz
    IF (str_cmp(name, 'temp_x_ev')) as_function = c_func_tempx_ev
    IF (str_cmp(name, 'temp_y_ev')) as_function = c_func_tempy_ev
    IF (str_cmp(name, 'temp_z_ev')) as_function = c_func_tempz_ev
    IF (str_cmp(name, 'ex')) as_function = c_func_ex
    IF (str_cmp(name, 'ey')) as_function = c_func_ey
    IF (str_cmp(name, 'ez')) as_function = c_func_ez
    IF (str_cmp(name, 'bx')) as_function = c_func_bx
    IF (str_cmp(name, 'by')) as_function = c_func_by
    IF (str_cmp(name, 'bz')) as_function = c_func_bz
    IF (str_cmp(name, 'critical')) as_function = c_func_crit

    DO i = 1, n_new_func
      IF (str_cmp(TRIM(name), TRIM(new_func_name(i)%value))) THEN
        as_function = new_func_code(i)
      ENDIF
    ENDDO

  END FUNCTION as_function



  FUNCTION as_operator(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_operator

    as_operator = c_prc_not_this_type

    IF (str_cmp(name, '+')) THEN
      IF (last_block_type == c_pt_variable &
          .OR. last_block_type == c_pt_constant &
          .OR. last_block_type == c_pt_default_constant &
          .OR. last_block_type == c_pt_deck_constant &
          .OR. last_block_type == c_pt_species &
          .OR. last_block_type == c_pt_subset) THEN
        as_operator = c_opcode_plus
      ELSE
        as_operator = c_opcode_unary_plus
      ENDIF
    ENDIF
    IF (str_cmp(name, '-'))  THEN
      IF (last_block_type == c_pt_variable &
          .OR. last_block_type == c_pt_constant &
          .OR. last_block_type == c_pt_default_constant &
          .OR. last_block_type == c_pt_deck_constant &
          .OR. last_block_type == c_pt_species &
          .OR. last_block_type == c_pt_subset) THEN
        as_operator = c_opcode_minus
      ELSE
        as_operator = c_opcode_unary_minus
      ENDIF
    ENDIF
    IF (str_cmp(name, '*')) THEN
      as_operator = c_opcode_times
    ENDIF
    IF (str_cmp(name, '/')) THEN
      as_operator = c_opcode_divide
    ENDIF
    IF (str_cmp(name, '^')) THEN
      as_operator = c_opcode_power
    ENDIF
    IF (str_cmp(name, 'e')) THEN
      as_operator = c_opcode_expo
    ENDIF
    IF (str_cmp(name, 'lt')) as_operator = c_opcode_lt
    IF (str_cmp(name, 'gt')) as_operator = c_opcode_gt
    IF (str_cmp(name, 'eq')) as_operator = c_opcode_eq
    IF (str_cmp(name, 'and')) as_operator = c_opcode_and
    IF (str_cmp(name, 'or'))  as_operator = c_opcode_or

  END FUNCTION as_operator



  SUBROUTINE check_deprecated(name)

    CHARACTER(LEN=*), INTENT(IN) :: name

    IF (str_cmp(name, 'ln')) THEN
      WRITE(*,*) '"' // TRIM(name) // '" is deprecated.'
      WRITE(*,*) 'Use "loge" instead.'
      RETURN
    ENDIF

  END SUBROUTINE check_deprecated

END MODULE tokenizer_blocks
