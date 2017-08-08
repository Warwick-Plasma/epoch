! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2014-2015 Stephan Kuschel <stephan.kuschel@gmail.com>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
! Copyright (C) 2012      Martin Ramsay <M.G.Ramsay@warwick.ac.uk>
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

! ****************************************************************
! All global variables defined here (cf F77 COMMON block).
! ****************************************************************

MODULE constants

  USE sdf, ONLY : c_stagger_edge_x, c_stagger_edge_y, c_stagger_edge_z, &
      c_stagger_face_x, c_stagger_face_y, c_stagger_face_z, &
      c_stagger_cell_centre, c_stagger_vertex, c_id_length

  IMPLICIT NONE

  INTEGER, PARAMETER :: num = KIND(1.d0)
  INTEGER, PARAMETER :: dbl = KIND(1.d0)
  INTEGER, PARAMETER :: r4  = SELECTED_REAL_KIND(r=30)
  INTEGER, PARAMETER :: i4  = SELECTED_INT_KIND(9)  ! 4-byte 2^31 ~ 10^9
  INTEGER, PARAMETER :: i8  = SELECTED_INT_KIND(18) ! 8-byte 2^63 ~ 10^18
  REAL(num), PARAMETER :: c_tiny = TINY(1.0_num)
  REAL(num), PARAMETER :: c_largest_number = HUGE(1.0_num)
  REAL(num), PARAMETER :: c_maxexponent = MAXEXPONENT(1.0_num)
  REAL(num), PARAMETER :: c_log2 = 0.69314718055994530941723212145817657_num
  REAL(num), PARAMETER :: c_largest_exp = c_maxexponent * c_log2
  REAL(num), PARAMETER :: &
      c_smallest_exp = (MINEXPONENT(1.0_num) - 1.0_num) * c_log2

  INTEGER, PARAMETER :: c_ndims = 1

  ! File unit for deck parser diagnostic output
  INTEGER, PARAMETER :: stdout = 6
  INTEGER, PARAMETER :: du = 40
  INTEGER, PARAMETER :: lu = 41
#ifdef NO_IO
  INTEGER, PARAMETER :: io_units(1) = (/ stdout /)
#else
  INTEGER, PARAMETER :: io_units(2) = (/ stdout,du /)
#endif
  INTEGER, PARAMETER :: nio_units = SIZE(io_units)

  ! Boundary type codes
  INTEGER, PARAMETER :: c_bc_null = -1
  INTEGER, PARAMETER :: c_bc_periodic = 1
  INTEGER, PARAMETER :: c_bc_other = 2
  INTEGER, PARAMETER :: c_bc_simple_laser = 3
  INTEGER, PARAMETER :: c_bc_simple_outflow = 4
  INTEGER, PARAMETER :: c_bc_open = 5
  INTEGER, PARAMETER :: c_bc_dump = 6
  INTEGER, PARAMETER :: c_bc_zero_gradient = 7
  INTEGER, PARAMETER :: c_bc_clamp = 8
  INTEGER, PARAMETER :: c_bc_reflect = 9
  INTEGER, PARAMETER :: c_bc_conduct = 10
  INTEGER, PARAMETER :: c_bc_thermal = 11
  INTEGER, PARAMETER :: c_bc_cpml_laser = 12
  INTEGER, PARAMETER :: c_bc_cpml_outflow = 13

  ! Boundary location codes
  INTEGER, PARAMETER :: c_bd_x_min = 1
  INTEGER, PARAMETER :: c_bd_x_max = 2
  INTEGER, PARAMETER :: c_bd_y_min = 3
  INTEGER, PARAMETER :: c_bd_y_max = 4
  INTEGER, PARAMETER :: c_bd_z_min = 5
  INTEGER, PARAMETER :: c_bd_z_max = 6

  ! Frequency function type codes
  INTEGER, PARAMETER :: c_of_omega = 1
  INTEGER, PARAMETER :: c_of_freq = 2
  INTEGER, PARAMETER :: c_of_lambda = 3

  ! Error codes
  INTEGER, PARAMETER :: c_err_none = 0
  INTEGER, PARAMETER :: c_err_unknown_block = 2**0
  INTEGER, PARAMETER :: c_err_unknown_element = 2**1
  INTEGER, PARAMETER :: c_err_preset_element = 2**2
  INTEGER, PARAMETER :: c_err_preset_element_use_later = 2**3
  INTEGER, PARAMETER :: c_err_bad_value = 2**4
  INTEGER, PARAMETER :: c_err_missing_elements = 2**5
  INTEGER, PARAMETER :: c_err_terminate = 2**6
  INTEGER, PARAMETER :: c_err_required_element_not_set = 2**7
  INTEGER, PARAMETER :: c_err_pp_options_missing = 2**8
  INTEGER, PARAMETER :: c_err_bad_array_length = 2**9
  INTEGER, PARAMETER :: c_err_other = 2**10
  INTEGER, PARAMETER :: c_err_warn_bad_value = 2**11
  INTEGER, PARAMETER :: c_err_generic_warning = 2**12
  INTEGER, PARAMETER :: c_err_generic_error = 2**13
  INTEGER, PARAMETER :: c_err_pp_options_wrong = 2**14
  INTEGER, PARAMETER :: c_err_io_error = 2**15
  INTEGER, PARAMETER :: c_err_bad_setup = 2**16

  INTEGER, PARAMETER :: c_ds_first = 1
  INTEGER, PARAMETER :: c_ds_last = 2

  ! IO codes
  INTEGER, PARAMETER :: c_io_none = 0
  INTEGER, PARAMETER :: c_io_always = 2**0
  INTEGER, PARAMETER :: c_io_full = 2**1
  INTEGER, PARAMETER :: c_io_restartable = 2**2
  INTEGER, PARAMETER :: c_io_species = 2**3
  INTEGER, PARAMETER :: c_io_no_sum = 2**4
  INTEGER, PARAMETER :: c_io_averaged = 2**5
  INTEGER, PARAMETER :: c_io_snapshot = 2**6
  INTEGER, PARAMETER :: c_io_field = 2**7
  INTEGER, PARAMETER :: c_io_dump_single = 2**8
  INTEGER, PARAMETER :: c_io_average_single = 2**9
  INTEGER, PARAMETER :: c_io_never = 2**10

  ! Maxwell Solvers
  INTEGER, PARAMETER :: c_maxwell_solver_yee = 0
  INTEGER, PARAMETER :: c_maxwell_solver_lehe = 1
  INTEGER, PARAMETER :: c_maxwell_solver_cowan = 2
  INTEGER, PARAMETER :: c_maxwell_solver_pukhov = 3

  ! domain codes
  INTEGER, PARAMETER :: c_do_full = 0
  INTEGER, PARAMETER :: c_do_decomposed = 1

  ! Load balance codes
  INTEGER, PARAMETER :: c_lb_x = 1
  INTEGER, PARAMETER :: c_lb_y = 2
  INTEGER, PARAMETER :: c_lb_z = 4
  INTEGER, PARAMETER :: c_lb_all = c_lb_x + c_lb_y + c_lb_z
  INTEGER, PARAMETER :: c_lb_auto = c_lb_all + 1

  ! Taken from http://physics.nist.gov/cuu/Constants (05/07/2012)
  REAL(num), PARAMETER :: pi = 3.141592653589793238462643383279503_num
  REAL(num), PARAMETER :: q0 = 1.602176565e-19_num ! C (+/- 3.5e-27)
  REAL(num), PARAMETER :: m0 = 9.10938291e-31_num ! kg (+/- 4e-38)
  REAL(num), PARAMETER :: c  = 2.99792458e8_num   ! m/s^2 (exact)
  REAL(num), PARAMETER :: kb = 1.3806488e-23_num  ! J/K (+/- 1.3e-29)
  REAL(num), PARAMETER :: mu0 = 4.e-7_num * pi ! N/A^2 (exact)
  ! epsilon0 = 1.0_num / mu0 / c**2 ! F/m (exact)
  REAL(num), PARAMETER :: epsilon0 = 8.854187817620389850536563031710750e-12_num
  REAL(num), PARAMETER :: h_planck = 6.62606957e-34_num ! J s (+/- 2.9e-41)
  REAL(num), PARAMETER :: ev = q0 ! J
  ! Derived physical parameters used in ionisation
  ! h_bar = h_planck / 2.0_num / pi
  REAL(num), PARAMETER :: h_bar = 1.054571725336289397963133257349698e-34_num
  ! Bohr radius
  ! a0 = 4.0_num * pi * epsilon0 * (h_bar / q0)**2.0_num / m0
  REAL(num), PARAMETER :: a0 = 5.291772101121111395947216558438464e-11_num
  ! Hartree energy = 2 * Rydberg energy
  ! hartree = (h_bar / a0)**2.0_num / m0
  REAL(num), PARAMETER :: hartree = 4.359744350823120007758594450644308e-18_num
  ! Fine structure constant, alpha = h_bar / a0 / m0 / c
  REAL(num), PARAMETER :: alpha = 7.2973525755230202568508027295271584628e-3_num
  ! atomic_time = h_bar / hartree
  REAL(num), PARAMETER :: &
      atomic_time = 2.418884320905619591809404261549867e-17_num
  ! atomic_electric_field = hartree / q0 / a0
  REAL(num), PARAMETER :: &
      atomic_electric_field = 5.142206538736485312185213306837419e11_num
  ! m0 * c
  REAL(num), PARAMETER :: mc0 = 2.73092429345209278e-22_num

  ! Constants used in pair production
#ifdef PHOTONS
  ! b_s = mc0**2 / (h_bar * q0)
  REAL(num), PARAMETER :: b_s = 4.414005028109566589829741352306303e9_num
  ! e_s = b_s * c
  REAL(num), PARAMETER :: e_s = 1.323285417001326061279735961512150e18_num
  ! alpha_f = q0**2 / (2.0_num * epsilon0 * h_planck * c)
  REAL(num), PARAMETER :: alpha_f = 7.297352575523020256850802729527158e-3_num
  ! tau_c = h_bar / (m0 * c**2)
  REAL(num), PARAMETER :: tau_c = 1.288088667367242662108649212042082e-21_num
#endif

  ! define special particle IDs
  INTEGER, PARAMETER :: c_species_id_generic = 0
  INTEGER, PARAMETER :: c_species_id_photon = 1
  INTEGER, PARAMETER :: c_species_id_electron = 2
  INTEGER, PARAMETER :: c_species_id_positron = 3
  INTEGER, PARAMETER :: c_species_id_proton = 4

  ! direction parameters
  INTEGER, PARAMETER :: c_dir_x = 1
  INTEGER, PARAMETER :: c_dir_y = 2
  INTEGER, PARAMETER :: c_dir_z = 3
  INTEGER, PARAMETER :: c_dir_px = c_ndims + 1
  INTEGER, PARAMETER :: c_dir_py = c_ndims + 2
  INTEGER, PARAMETER :: c_dir_pz = c_ndims + 3
  INTEGER, PARAMETER :: c_dir_en = c_ndims + 4
  INTEGER, PARAMETER :: c_dir_gamma_m1 = c_ndims + 5
  INTEGER, PARAMETER :: c_dir_xy_angle = c_ndims + 6
  INTEGER, PARAMETER :: c_dir_yz_angle = c_ndims + 7
  INTEGER, PARAMETER :: c_dir_zx_angle = c_ndims + 8
  INTEGER, PARAMETER :: c_dir_mod_p = c_ndims + 9

  ! constants defining the maximum number of dimensions and directions
  ! in a distribution function
  INTEGER, PARAMETER :: c_df_maxdirs = c_dir_mod_p
  INTEGER, PARAMETER :: c_df_maxdims = 3

  ! define flags
  INTEGER(i8), PARAMETER :: c_def_particle_debug = 2**0
  INTEGER(i8), PARAMETER :: c_def_field_debug = 2**1
  INTEGER(i8), PARAMETER :: c_def_particle_shape_bspline3 = 2**2
  INTEGER(i8), PARAMETER :: c_def_split_part_after_push = 2**3
  INTEGER(i8), PARAMETER :: c_def_per_particle_weight = 2**4
  INTEGER(i8), PARAMETER :: c_def_particle_count_update = 2**5
  INTEGER(i8), PARAMETER :: c_def_tracer_particles = 2**6
  INTEGER(i8), PARAMETER :: c_def_particle_probes = 2**7
  INTEGER(i8), PARAMETER :: c_def_per_particle_chargemass = 2**8
  INTEGER(i8), PARAMETER :: c_def_particle_ionise = 2**9
  INTEGER(i8), PARAMETER :: c_def_high_order_smoothing = 2**10
  INTEGER(i8), PARAMETER :: c_def_particle_shape_tophat = 2**11
  INTEGER(i8), PARAMETER :: c_def_parser_debug = 2**12
  INTEGER(i8), PARAMETER :: c_def_particle_id4 = 2**13
  INTEGER(i8), PARAMETER :: c_def_particle_id  = 2**14
  INTEGER(i8), PARAMETER :: c_def_photons = 2**15
  INTEGER(i8), PARAMETER :: c_def_trident_photons = 2**16
  INTEGER(i8), PARAMETER :: c_def_prefetch = 2**17
  INTEGER(i8), PARAMETER :: c_def_mpi_debug = 2**18
  INTEGER(i8), PARAMETER :: c_def_parser_checking = 2**19

  ! Stagger types
  INTEGER, PARAMETER :: c_stagger_ex = c_stagger_face_x
  INTEGER, PARAMETER :: c_stagger_ey = c_stagger_face_y
  INTEGER, PARAMETER :: c_stagger_ez = c_stagger_face_z
  INTEGER, PARAMETER :: c_stagger_bx = c_stagger_edge_x
  INTEGER, PARAMETER :: c_stagger_by = c_stagger_edge_y
  INTEGER, PARAMETER :: c_stagger_bz = c_stagger_edge_z
  INTEGER, PARAMETER :: c_stagger_centre = c_stagger_cell_centre
  INTEGER, PARAMETER :: c_stagger_max = c_stagger_vertex
  INTEGER, PARAMETER :: c_stagger_jx = c_stagger_ex
  INTEGER, PARAMETER :: c_stagger_jy = c_stagger_ey
  INTEGER, PARAMETER :: c_stagger_jz = c_stagger_ez

  ! Length of a standard string
  INTEGER, PARAMETER :: string_length = 256

  INTEGER, PARAMETER :: stat_unit = 20

END MODULE constants



MODULE shared_parser_data

  USE constants

  IMPLICIT NONE

  INTEGER, PARAMETER :: c_char_numeric = 1
  INTEGER, PARAMETER :: c_char_alpha = 2
  INTEGER, PARAMETER :: c_char_delimiter = 3
  INTEGER, PARAMETER :: c_char_space = 4
  INTEGER, PARAMETER :: c_char_opcode = 5
  INTEGER, PARAMETER :: c_char_unknown = 1024

  INTEGER, PARAMETER :: c_prc_not_this_type = 0

  ! block type constants
  INTEGER, PARAMETER :: c_pt_variable = 1
  INTEGER, PARAMETER :: c_pt_constant = 2
  INTEGER, PARAMETER :: c_pt_operator = 3
  INTEGER, PARAMETER :: c_pt_function = 4
  INTEGER, PARAMETER :: c_pt_parenthesis = 5
  INTEGER, PARAMETER :: c_pt_separator = 6
  INTEGER, PARAMETER :: c_pt_character = 7
  INTEGER, PARAMETER :: c_pt_deck_constant = 8
  INTEGER, PARAMETER :: c_pt_species = 9
  INTEGER, PARAMETER :: c_pt_subset = 10
  INTEGER, PARAMETER :: c_pt_default_constant = 11
  INTEGER, PARAMETER :: c_pt_bad = 1024
  INTEGER, PARAMETER :: c_pt_null = 1025

  ! Opcode constants
  INTEGER, PARAMETER :: c_opcode_plus = 1
  INTEGER, PARAMETER :: c_opcode_minus = 2
  INTEGER, PARAMETER :: c_opcode_times = 3
  INTEGER, PARAMETER :: c_opcode_divide = 4
  INTEGER, PARAMETER :: c_opcode_power = 5
  INTEGER, PARAMETER :: c_opcode_expo = 6
  INTEGER, PARAMETER :: c_opcode_lt = 7
  INTEGER, PARAMETER :: c_opcode_gt = 8
  INTEGER, PARAMETER :: c_opcode_eq = 9
  INTEGER, PARAMETER :: c_opcode_and = 10
  INTEGER, PARAMETER :: c_opcode_or = 11
  INTEGER, PARAMETER :: c_opcode_unary_minus = 12
  INTEGER, PARAMETER :: c_opcode_unary_plus = 13
  INTEGER, PARAMETER :: c_num_ops = 13

  ! Associativity constants
  INTEGER, PARAMETER :: c_assoc_a = 1
  INTEGER, PARAMETER :: c_assoc_la = 2
  INTEGER, PARAMETER :: c_assoc_ra = 3

  INTEGER, DIMENSION(c_num_ops), PARAMETER :: &
      opcode_precedence = (/2, 2, 3, 3, 4, 4, 1, 1, 1, 0, 0, 4, 4/)
  INTEGER, DIMENSION(c_num_ops), PARAMETER :: &
      opcode_assoc = (/c_assoc_a, c_assoc_la, c_assoc_a, c_assoc_la, &
          c_assoc_ra, c_assoc_ra, c_assoc_la, c_assoc_la, c_assoc_la, &
          c_assoc_la, c_assoc_la, c_assoc_ra, c_assoc_ra/)

  INTEGER, PARAMETER :: c_paren_left_bracket = 1
  INTEGER, PARAMETER :: c_paren_right_bracket = 2

  ! Actual constants
  INTEGER, PARAMETER :: c_const_pi = 1
  INTEGER, PARAMETER :: c_const_kb = 2
  INTEGER, PARAMETER :: c_const_me = 3
  INTEGER, PARAMETER :: c_const_qe = 4
  INTEGER, PARAMETER :: c_const_eps0 = 5
  INTEGER, PARAMETER :: c_const_mu0 = 6
  INTEGER, PARAMETER :: c_const_c = 7
  INTEGER, PARAMETER :: c_const_ev = 8
  INTEGER, PARAMETER :: c_const_kev = 9
  INTEGER, PARAMETER :: c_const_mev = 10
  INTEGER, PARAMETER :: c_const_milli = 11
  INTEGER, PARAMETER :: c_const_micro = 12
  INTEGER, PARAMETER :: c_const_nano = 13
  INTEGER, PARAMETER :: c_const_pico = 14
  INTEGER, PARAMETER :: c_const_femto = 15
  INTEGER, PARAMETER :: c_const_atto = 16

  ! Constants refering to grid properties
  INTEGER, PARAMETER :: c_const_xb = 22
  INTEGER, PARAMETER :: c_const_yb = 23
  INTEGER, PARAMETER :: c_const_zb = 24
  INTEGER, PARAMETER :: c_const_x = 25
  INTEGER, PARAMETER :: c_const_y = 26
  INTEGER, PARAMETER :: c_const_z = 27
  INTEGER, PARAMETER :: c_const_lx = 28
  INTEGER, PARAMETER :: c_const_ly = 29
  INTEGER, PARAMETER :: c_const_lz = 30
  INTEGER, PARAMETER :: c_const_dx = 31
  INTEGER, PARAMETER :: c_const_dy = 32
  INTEGER, PARAMETER :: c_const_dz = 33
  INTEGER, PARAMETER :: c_const_x_min = 34
  INTEGER, PARAMETER :: c_const_y_min = 35
  INTEGER, PARAMETER :: c_const_z_min = 36
  INTEGER, PARAMETER :: c_const_x_max = 37
  INTEGER, PARAMETER :: c_const_y_max = 38
  INTEGER, PARAMETER :: c_const_z_max = 39
  INTEGER, PARAMETER :: c_const_ix = 40
  INTEGER, PARAMETER :: c_const_iy = 41
  INTEGER, PARAMETER :: c_const_iz = 42
  INTEGER, PARAMETER :: c_const_nx = 43
  INTEGER, PARAMETER :: c_const_ny = 44
  INTEGER, PARAMETER :: c_const_nz = 45
  INTEGER, PARAMETER :: c_const_time = 46
  INTEGER, PARAMETER :: c_const_r_xy = 47
  INTEGER, PARAMETER :: c_const_r_yz = 48
  INTEGER, PARAMETER :: c_const_r_xz = 49
  INTEGER, PARAMETER :: c_const_nprocx = 50
  INTEGER, PARAMETER :: c_const_nprocy = 51
  INTEGER, PARAMETER :: c_const_nprocz = 52
  INTEGER, PARAMETER :: c_const_nsteps = 53
  INTEGER, PARAMETER :: c_const_t_end = 54
  INTEGER, PARAMETER :: c_const_ndims = 55
  INTEGER, PARAMETER :: c_const_r_xyz = 56

  INTEGER, PARAMETER :: c_const_io_never = 60
  INTEGER, PARAMETER :: c_const_io_always = 61
  INTEGER, PARAMETER :: c_const_io_full = 62
  INTEGER, PARAMETER :: c_const_io_restartable = 63
  INTEGER, PARAMETER :: c_const_io_species = 64
  INTEGER, PARAMETER :: c_const_io_no_sum = 65
  INTEGER, PARAMETER :: c_const_io_average = 66
  INTEGER, PARAMETER :: c_const_io_snapshot = 67
  INTEGER, PARAMETER :: c_const_io_dump_single = 68
  INTEGER, PARAMETER :: c_const_io_average_single = 69

  INTEGER, PARAMETER :: c_const_dir_x = 80
  INTEGER, PARAMETER :: c_const_dir_y = 81
  INTEGER, PARAMETER :: c_const_dir_z = 82
  INTEGER, PARAMETER :: c_const_dir_px = 83
  INTEGER, PARAMETER :: c_const_dir_py = 84
  INTEGER, PARAMETER :: c_const_dir_pz = 85
  INTEGER, PARAMETER :: c_const_dir_en = 86
  INTEGER, PARAMETER :: c_const_dir_gamma_m1 = 87
  INTEGER, PARAMETER :: c_const_dir_xy_angle = 88
  INTEGER, PARAMETER :: c_const_dir_yz_angle = 89
  INTEGER, PARAMETER :: c_const_dir_zx_angle = 90
  INTEGER, PARAMETER :: c_const_dir_mod_p = 91

  INTEGER, PARAMETER :: c_const_maxwell_solver_yee = 100
  INTEGER, PARAMETER :: c_const_maxwell_solver_lehe = 101
  INTEGER, PARAMETER :: c_const_maxwell_solver_cowan = 102
  INTEGER, PARAMETER :: c_const_maxwell_solver_pukhov = 103

  ! Custom constants
  INTEGER, PARAMETER :: c_const_deck_lowbound = 4096
  INTEGER, PARAMETER :: c_const_custom_lowbound = 8192

  INTEGER, PARAMETER :: c_func_sine = 1
  INTEGER, PARAMETER :: c_func_cosine = 2
  INTEGER, PARAMETER :: c_func_tan = 3
  INTEGER, PARAMETER :: c_func_exp = 4
  INTEGER, PARAMETER :: c_func_arcsine = 10
  INTEGER, PARAMETER :: c_func_arccosine = 11
  INTEGER, PARAMETER :: c_func_arctan = 12
  INTEGER, PARAMETER :: c_func_neg = 13
  INTEGER, PARAMETER :: c_func_if = 14
  INTEGER, PARAMETER :: c_func_floor = 15
  INTEGER, PARAMETER :: c_func_ceil = 16
  INTEGER, PARAMETER :: c_func_nint = 17
  INTEGER, PARAMETER :: c_func_rho = 18
  INTEGER, PARAMETER :: c_func_tempx = 19
  INTEGER, PARAMETER :: c_func_tempy = 20
  INTEGER, PARAMETER :: c_func_tempz = 21
  INTEGER, PARAMETER :: c_func_interpolate = 22
  INTEGER, PARAMETER :: c_func_tanh = 23
  INTEGER, PARAMETER :: c_func_sinh = 24
  INTEGER, PARAMETER :: c_func_cosh = 25
  INTEGER, PARAMETER :: c_func_ex = 26
  INTEGER, PARAMETER :: c_func_ey = 27
  INTEGER, PARAMETER :: c_func_ez = 28
  INTEGER, PARAMETER :: c_func_bx = 29
  INTEGER, PARAMETER :: c_func_by = 30
  INTEGER, PARAMETER :: c_func_bz = 31
  INTEGER, PARAMETER :: c_func_sqrt = 32
  INTEGER, PARAMETER :: c_func_gauss = 33
  INTEGER, PARAMETER :: c_func_semigauss = 34
  INTEGER, PARAMETER :: c_func_crit = 35
  INTEGER, PARAMETER :: c_func_abs = 36
  INTEGER, PARAMETER :: c_func_loge = 37
  INTEGER, PARAMETER :: c_func_log10 = 38
  INTEGER, PARAMETER :: c_func_log_base = 39
  INTEGER, PARAMETER :: c_func_supergauss = 40
  INTEGER, PARAMETER :: c_func_tempx_ev = 41
  INTEGER, PARAMETER :: c_func_tempy_ev = 42
  INTEGER, PARAMETER :: c_func_tempz_ev = 43
  INTEGER, PARAMETER :: c_func_driftx = 44
  INTEGER, PARAMETER :: c_func_drifty = 45
  INTEGER, PARAMETER :: c_func_driftz = 46

  INTEGER, PARAMETER :: c_func_custom_lowbound = 4096

  ! This type represents parameters given to the parser.
  ! It can be extended by a developer freely
  ! It is the responsibility of the developer to ensure that a parameter is
  ! specified when needed

  ! If you set the use_grid_position parameter to .FALSE. then the deck parser
  ! will evaluate position x, y, z as being at the location pack_pos(1,2,3)
  ! rather than x(pack%ix), y(pack%iy), z(pack%iz). It is essential that the
  ! ix, parameters are still set to match, because other functions
  ! will still use them
  TYPE parameter_pack
    LOGICAL :: use_grid_position = .TRUE.
    INTEGER :: pack_ix = 1
    REAL(num) :: pack_pos = 0.0_num
  END TYPE parameter_pack

  TYPE stack_element
    INTEGER :: ptype
    INTEGER :: value
    REAL(num) :: numerical_data
#ifdef PARSER_DEBUG
    CHARACTER(LEN=string_length) :: text
#endif
  END TYPE stack_element

  TYPE primitive_stack
    TYPE(stack_element), POINTER :: entries(:)
    INTEGER :: stack_point, stack_size
    LOGICAL :: init = .FALSE.
    LOGICAL :: is_time_varying = .FALSE.
    LOGICAL :: should_simplify
  END TYPE primitive_stack

  TYPE deck_constant
    CHARACTER(LEN=string_length) :: name
    TYPE(primitive_stack) :: execution_stream
  END TYPE deck_constant

  INTEGER :: n_deck_constants = 0
  TYPE(deck_constant), DIMENSION(:), ALLOCATABLE :: deck_constant_list

END MODULE shared_parser_data



MODULE shared_data

  USE mpi
  USE shared_parser_data
  USE sdf_job_info

  IMPLICIT NONE

  !----------------------------------------------------------------------------
  ! string handling
  !----------------------------------------------------------------------------
  CHARACTER(LEN=string_length) :: blank
  TYPE string_type
    CHARACTER(string_length) :: value
  END TYPE string_type
  CHARACTER(LEN=string_length) :: extended_error_string

  !----------------------------------------------------------------------------
  ! Particles
  !----------------------------------------------------------------------------

  ! Time to start the particle push - 0 by default, can be set in the control
  ! block of the deck using 'particle_tstart'.
  REAL(num) :: particle_push_start_time = 0.0_num

  ! The order for the spline interpolation used as a particle representation.
  ! png is the number of ghost cells needed by the particles
#ifdef PARTICLE_SHAPE_BSPLINE3
  INTEGER, PARAMETER :: sf_min = -2
  INTEGER, PARAMETER :: sf_max =  2
  INTEGER, PARAMETER :: png =  4
#elif  PARTICLE_SHAPE_TOPHAT
  INTEGER, PARAMETER :: sf_min =  0
  INTEGER, PARAMETER :: sf_max =  1
  INTEGER, PARAMETER :: png =  2
#else
  INTEGER, PARAMETER :: sf_min = -1
  INTEGER, PARAMETER :: sf_max =  1
  INTEGER, PARAMETER :: png =  3
#endif

  ! Object representing a particle
  ! If you add or remove from this section then you *must* update the
  ! particle pack and unpack routines
  TYPE particle
    REAL(num), DIMENSION(3) :: part_p
    REAL(num) :: part_pos
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS)
    REAL(num) :: weight
#endif
#ifdef DELTAF_METHOD
    REAL(num) :: pvol
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    REAL(num) :: charge
    REAL(num) :: mass
#endif
    TYPE(particle), POINTER :: next, prev
#ifdef PARTICLE_DEBUG
    INTEGER :: processor
    INTEGER :: processor_at_t0
#endif
#ifdef PARTICLE_ID4
    INTEGER :: id
#elif PARTICLE_ID
    INTEGER(i8) :: id
#endif
#ifdef COLLISIONS_TEST
    INTEGER :: coll_count
#endif
#ifdef PHOTONS
    REAL(num) :: optical_depth
    REAL(num) :: particle_energy
#ifdef TRIDENT_PHOTONS
    REAL(num) :: optical_depth_tri
#endif
#endif
    REAL(num) :: force_multiplier = 1.0_num
  END TYPE particle

  ! Data for migration between species
  TYPE particle_species_migration
    LOGICAL :: this_species, fluid, done
    LOGICAL :: promoteable, demoteable
    INTEGER :: promote_to_species, demote_to_species
    REAL(num) :: promotion_energy_factor, demotion_energy_factor
    REAL(num) :: promotion_density, demotion_density
    REAL(num), DIMENSION(:), POINTER :: fluid_energy, fluid_density
  END TYPE particle_species_migration

  LOGICAL :: use_particle_migration = .FALSE.
  INTEGER :: particle_migration_interval = 1

  ! Object representing a collection of particles
  ! Used internally by the MPI particle transfer code
  TYPE particle_list
    TYPE(particle), POINTER :: head
    TYPE(particle), POINTER :: tail
    INTEGER(i8) :: count
    INTEGER :: id_update
    ! Pointer is safe if the particles in it are all unambiguously linked
    LOGICAL :: safe

    TYPE(particle_list), POINTER :: next, prev
  END TYPE particle_list

  ! Represents the initial conditions of a species
  TYPE initial_condition_block
    REAL(num), DIMENSION(:), POINTER :: density
    REAL(num), DIMENSION(:,:), POINTER :: temp
    REAL(num), DIMENSION(:,:), POINTER :: drift

    REAL(num) :: density_min
    REAL(num) :: density_max

    ! Delta-f settings
    REAL(num) :: temp_back(3)
    REAL(num) :: drift_back(3)
    REAL(num) :: density_back
  END TYPE initial_condition_block

  ! Object representing a particle species
  TYPE particle_species
    ! Core properties
    CHARACTER(string_length) :: name
    TYPE(particle_species), POINTER :: next, prev
    INTEGER :: id
    INTEGER :: dumpmask
    INTEGER :: count_update_step

    REAL(num) :: charge
    REAL(num) :: mass
    REAL(num) :: weight
    INTEGER(i8) :: count
    TYPE(particle_list) :: attached_list
    LOGICAL :: immobile

#ifndef NO_TRACER_PARTICLES
    LOGICAL :: tracer
#endif

    ! ID code which identifies if a species is of a special type
    INTEGER :: species_type

    ! particle cell division
    INTEGER(i8) :: global_count
    LOGICAL :: split
    INTEGER(i8) :: npart_max
    ! Secondary list
    TYPE(particle_list), DIMENSION(:), POINTER :: secondary_list

    ! Injection of particles
    REAL(num) :: npart_per_cell
    TYPE(primitive_stack) :: density_function, temperature_function(3)
    TYPE(primitive_stack) :: drift_function(3)

    ! Thermal boundaries
    REAL(num), DIMENSION(:), POINTER :: ext_temp_x_min, ext_temp_x_max

    ! Species_ionisation
    LOGICAL :: electron
    LOGICAL :: ionise
    INTEGER :: ionise_to_species
    INTEGER :: release_species
    INTEGER :: n
    INTEGER :: l
    REAL(num) :: ionisation_energy

    ! Attached probes for this species
#ifndef NO_PARTICLE_PROBES
    TYPE(particle_probe), POINTER :: attached_probes
#endif

    ! Particle migration
    TYPE(particle_species_migration) :: migrate

    ! Initial conditions
    TYPE(initial_condition_block) :: initial_conditions

    !Per species boundary conditions
    INTEGER, DIMENSION(2 * c_ndims) :: bc_particle
  END TYPE particle_species

  !----------------------------------------------------------------------------
  ! file handling
  !----------------------------------------------------------------------------
  INTEGER :: deck_state
  INTEGER :: subtype_field, subtype_field_r4
  INTEGER :: subarray_field, subarray_field_r4
  INTEGER :: subarray_field_big, subarray_field_big_r4
  INTEGER(KIND=MPI_OFFSET_KIND) :: initialdisp
  INTEGER :: full_dump_every, restart_dump_every
  LOGICAL :: force_first_to_be_restartable
  LOGICAL :: force_final_to_be_restartable
  LOGICAL :: use_offset_grid
  INTEGER :: n_zeros_control, n_zeros = 4
  INTEGER, PARAMETER :: c_max_zeros = 9
  INTEGER, PARAMETER :: c_dump_part_grid         = 1
  INTEGER, PARAMETER :: c_dump_grid              = 2
  INTEGER, PARAMETER :: c_dump_part_species      = 3
  INTEGER, PARAMETER :: c_dump_part_weight       = 4
  INTEGER, PARAMETER :: c_dump_part_px           = 5
  INTEGER, PARAMETER :: c_dump_part_py           = 6
  INTEGER, PARAMETER :: c_dump_part_pz           = 7
  INTEGER, PARAMETER :: c_dump_part_vx           = 8
  INTEGER, PARAMETER :: c_dump_part_vy           = 9
  INTEGER, PARAMETER :: c_dump_part_vz           = 10
  INTEGER, PARAMETER :: c_dump_part_charge       = 11
  INTEGER, PARAMETER :: c_dump_part_mass         = 12
  INTEGER, PARAMETER :: c_dump_part_id           = 13
  INTEGER, PARAMETER :: c_dump_ex                = 14
  INTEGER, PARAMETER :: c_dump_ey                = 15
  INTEGER, PARAMETER :: c_dump_ez                = 16
  INTEGER, PARAMETER :: c_dump_bx                = 17
  INTEGER, PARAMETER :: c_dump_by                = 18
  INTEGER, PARAMETER :: c_dump_bz                = 19
  INTEGER, PARAMETER :: c_dump_jx                = 20
  INTEGER, PARAMETER :: c_dump_jy                = 21
  INTEGER, PARAMETER :: c_dump_jz                = 22
  INTEGER, PARAMETER :: c_dump_ekbar             = 23
  INTEGER, PARAMETER :: c_dump_mass_density      = 24
  INTEGER, PARAMETER :: c_dump_charge_density    = 25
  INTEGER, PARAMETER :: c_dump_number_density    = 26
  INTEGER, PARAMETER :: c_dump_temperature       = 27
  INTEGER, PARAMETER :: c_dump_dist_fns          = 28
  INTEGER, PARAMETER :: c_dump_probes            = 29
  INTEGER, PARAMETER :: c_dump_ejected_particles = 30
  INTEGER, PARAMETER :: c_dump_ekflux            = 31
  INTEGER, PARAMETER :: c_dump_poynt_flux        = 32
  INTEGER, PARAMETER :: c_dump_cpml_psi_eyx      = 33
  INTEGER, PARAMETER :: c_dump_cpml_psi_ezx      = 34
  INTEGER, PARAMETER :: c_dump_cpml_psi_byx      = 35
  INTEGER, PARAMETER :: c_dump_cpml_psi_bzx      = 36
  INTEGER, PARAMETER :: c_dump_absorption        = 37
  INTEGER, PARAMETER :: c_dump_part_ek           = 38
  INTEGER, PARAMETER :: c_dump_part_opdepth      = 39
  INTEGER, PARAMETER :: c_dump_part_qed_energy   = 40
  INTEGER, PARAMETER :: c_dump_part_opdepth_tri  = 41
  INTEGER, PARAMETER :: c_dump_total_energy      = 42
  INTEGER, PARAMETER :: c_dump_total_energy_sum  = 43
  INTEGER, PARAMETER :: c_dump_part_rel_mass     = 44
  INTEGER, PARAMETER :: c_dump_part_gamma        = 45
  INTEGER, PARAMETER :: c_dump_part_proc         = 46
  INTEGER, PARAMETER :: c_dump_part_proc0        = 47
  INTEGER, PARAMETER :: num_vars_to_dump         = 47
  INTEGER, DIMENSION(num_vars_to_dump) :: dumpmask

  !----------------------------------------------------------------------------
  ! Time averaged IO
  !----------------------------------------------------------------------------
  TYPE averaged_data_block
    REAL(num), DIMENSION(:,:), POINTER :: array
    REAL(r4), DIMENSION(:,:), POINTER :: r4array
    REAL(num) :: real_time
    INTEGER :: species_sum, n_species
    LOGICAL :: started, dump_single
  END TYPE averaged_data_block
  LOGICAL :: any_average = .FALSE.

  TYPE io_block_type
    CHARACTER(LEN=string_length) :: name
    REAL(num) :: dt_snapshot, time_prev, time_first
    REAL(num) :: dt_average, dt_min_average, average_time, average_time_start
    REAL(num) :: time_start, time_stop
    REAL(num), POINTER :: dump_at_times(:)
    INTEGER, POINTER :: dump_at_nsteps(:)
    INTEGER :: nstep_snapshot, nstep_prev, nstep_first, nstep_average
    INTEGER :: nstep_start, nstep_stop, dump_cycle, prefix_index
    INTEGER :: dump_cycle_first_index
    LOGICAL :: restart, dump, any_average, dump_first, dump_last
    LOGICAL :: dump_source_code, dump_input_decks, rolling_restart
    LOGICAL :: dump_first_after_restart
    LOGICAL :: disabled
    INTEGER, DIMENSION(num_vars_to_dump) :: dumpmask
    TYPE(averaged_data_block), DIMENSION(num_vars_to_dump) :: averaged_data
  END TYPE io_block_type

  TYPE(io_block_type), POINTER :: io_block_list(:)
  INTEGER :: n_io_blocks
  LOGICAL :: track_ejected_particles, new_style_io_block
  INTEGER, DIMENSION(num_vars_to_dump) :: averaged_var_block
  REAL(num) :: time_start, time_stop
  INTEGER :: nstep_start, nstep_stop
  CHARACTER(LEN=c_id_length), ALLOCATABLE :: file_prefixes(:)
  INTEGER, ALLOCATABLE :: file_numbers(:)
  INTEGER(i8) :: sdf_buffer_size

  !----------------------------------------------------------------------------
  ! Extended IO information
  !----------------------------------------------------------------------------

  ! Represents a 2 or 3D distribution
  TYPE distribution_function_block
    CHARACTER(LEN=string_length) :: name

    ! The number of dimensions left
    INTEGER :: ndims

    ! The dumpmask for the distribution function
    INTEGER :: dumpmask

    ! The variables which define the ranges and resolutions of the
    ! distribution function
    INTEGER, DIMENSION(c_df_maxdims) :: directions
    REAL(num), DIMENSION(2,c_df_maxdims) :: ranges
    INTEGER, DIMENSION(c_df_maxdims) :: resolution
    LOGICAL, DIMENSION(:), POINTER :: use_species
    REAL(num), DIMENSION(2,c_df_maxdirs) :: restrictions
    LOGICAL, DIMENSION(c_df_maxdirs) :: use_restrictions

    ! Whether to output deltaf or totalf
    LOGICAL :: output_deltaf

    ! Pointer to next distribution function
    TYPE(distribution_function_block), POINTER :: next
  END TYPE distribution_function_block
  TYPE(distribution_function_block), POINTER :: dist_fns

  ! Represents a data subset
  TYPE subset
    CHARACTER(LEN=string_length) :: name

    ! The dumpmask for the subset
    INTEGER :: mask
    INTEGER, DIMENSION(:,:), POINTER :: dumpmask
    LOGICAL, DIMENSION(:), POINTER :: use_species
    LOGICAL :: use_gamma, use_gamma_min, use_gamma_max, use_random
    LOGICAL :: use_x_min, use_x_max
    LOGICAL :: use_px_min, use_px_max
    LOGICAL :: use_py_min, use_py_max
    LOGICAL :: use_pz_min, use_pz_max
    LOGICAL :: use_weight_min, use_weight_max
    LOGICAL :: use_charge_min, use_charge_max
    LOGICAL :: use_mass_min, use_mass_max
    LOGICAL :: use_id_min, use_id_max
    LOGICAL :: space_restrictions
    LOGICAL :: skip, dump_field_grid
    REAL(num) :: gamma_min, gamma_max, random_fraction
    REAL(num) :: x_min, x_max
    REAL(num) :: px_min, px_max, py_min, py_max, pz_min, pz_max
    REAL(num) :: weight_min, weight_max
    REAL(num) :: charge_min, charge_max
    REAL(num) :: mass_min, mass_max
    INTEGER(i8) :: id_min, id_max
    INTEGER :: subtype, subarray, subtype_r4, subarray_r4
    INTEGER, DIMENSION(c_ndims) :: skip_dir, n_local, n_global, n_start

    ! Pointer to next subset
    TYPE(subset), POINTER :: next
  END TYPE subset
  TYPE(subset), DIMENSION(:), POINTER :: subset_list
  INTEGER :: n_subsets

#ifndef NO_PARTICLE_PROBES
  TYPE particle_probe
    ! Arbitrary point on the plane
    REAL(num) :: point
    ! The normal to the plane
    REAL(num) :: normal
    REAL(num) :: ek_min, ek_max
    CHARACTER(LEN=string_length) :: name

    LOGICAL, DIMENSION(:), POINTER :: use_species
    TYPE(particle_list) :: sampled_particles
    TYPE(particle_probe), POINTER :: next
    INTEGER :: dumpmask
  END TYPE particle_probe
#endif

  INTEGER :: cpml_thickness
  INTEGER :: cpml_x_min_start, cpml_x_min_end, cpml_x_min_offset
  INTEGER :: cpml_x_max_start, cpml_x_max_end, cpml_x_max_offset
  ! Indicate that we have a boundary on the current processor
  LOGICAL :: cpml_x_min = .FALSE., cpml_x_max = .FALSE.
  LOGICAL :: cpml_boundaries
  ! Indicate that the laser injection is located on the current processor
  INTEGER :: cpml_x_min_laser_idx, cpml_x_max_laser_idx
  REAL(num) :: cpml_kappa_max, cpml_a_max, cpml_sigma_max
  REAL(num), ALLOCATABLE, DIMENSION(:) :: cpml_kappa_ex, cpml_kappa_bx
  REAL(num), ALLOCATABLE, DIMENSION(:) :: cpml_a_ex, cpml_sigma_ex
  REAL(num), ALLOCATABLE, DIMENSION(:) :: cpml_a_bx, cpml_sigma_bx
  REAL(num), ALLOCATABLE, DIMENSION(:) :: cpml_psi_eyx, cpml_psi_ezx
  REAL(num), ALLOCATABLE, DIMENSION(:) :: cpml_psi_byx, cpml_psi_bzx

  !----------------------------------------------------------------------------
  ! Core code
  !----------------------------------------------------------------------------
  INTEGER :: mpireal = MPI_DOUBLE_PRECISION
  INTEGER :: realsize

  ! ng is the number of ghost cells allocated in the arrays
  ! fng is the number of ghost cells needed by the field solver
  ! jng is the number of ghost cells needed by the current arrays
  INTEGER, PARAMETER :: ng = png * 2
  INTEGER, PARAMETER :: jng =  MAX(ng,png)
  INTEGER :: fng, nx
  INTEGER :: nx_global
  INTEGER(i8) :: npart_global, particles_max_id
  INTEGER :: nsteps, n_species = -1
  LOGICAL :: smooth_currents
  REAL(num), ALLOCATABLE, DIMENSION(:) :: ex, ey, ez, bx, by, bz, jx, jy, jz
  REAL(num), ALLOCATABLE, DIMENSION(:) :: wk_array

  REAL(num) :: ex_x_min, ex_x_max
  REAL(num) :: ey_x_min, ey_x_max
  REAL(num) :: ez_x_min, ez_x_max
  REAL(num) :: bx_x_min, bx_x_max
  REAL(num) :: by_x_min, by_x_max
  REAL(num) :: bz_x_min, bz_x_max

  REAL(num) :: initial_jx, initial_jy, initial_jz

  TYPE(particle_species), DIMENSION(:), POINTER :: species_list
  TYPE(particle_species), DIMENSION(:), POINTER :: ejected_list
  TYPE(particle_species), DIMENSION(:), POINTER :: io_list, io_list_data

  REAL(num), ALLOCATABLE, DIMENSION(:) :: x, xb

  INTEGER, PARAMETER :: c_max_string_length = 64
  INTEGER, PARAMETER :: c_max_prefix = 16
  ! Maximum path length on Linux machines
  INTEGER, PARAMETER :: c_max_path_length = 4096 + c_max_prefix
  CHARACTER(LEN=c_max_path_length) :: data_dir
  CHARACTER(LEN=c_max_prefix) :: filesystem

  LOGICAL :: neutral_background = .TRUE.
  LOGICAL :: use_random_seed = .FALSE.
  LOGICAL :: use_particle_lists = .FALSE.

  REAL(num) :: dt, t_end, time, dt_multiplier, dt_laser, dt_plasma_frequency
  REAL(num) :: dt_from_restart
  REAL(num) :: dt_min_average, cfl
  ! x_min is the left-hand edge of the simulation domain as specified in
  ! the input deck.
  ! x_grid_min is the location of x(1). Since the grid is cell-centred,
  ! this is usually at x_min + dx/2.
  ! If CPML boundaries are used then the whole grid is shifted along by
  ! cpml_thicknes cells and then x(1) (and also x_grid_min) is at
  ! the location x_min + dx*(1/2-cpml_thickness)
  REAL(num) :: length_x, dx, x_grid_min, x_grid_max, x_min, x_max
  REAL(num) :: x_grid_min_local, x_grid_max_local, x_min_local, x_max_local
  REAL(num), DIMENSION(:), ALLOCATABLE :: x_grid_mins, x_grid_maxs

  REAL(num) :: total_ohmic_heating = 0.0_num

  LOGICAL :: ic_from_restart = .FALSE.
  LOGICAL :: need_random_state
  LOGICAL :: use_exact_restart
  LOGICAL :: allow_cpu_reduce
  LOGICAL :: simplify_deck
  LOGICAL :: print_deck_constants
  LOGICAL :: allow_missing_restart
  LOGICAL :: done_mpi_initialise = .FALSE.
  LOGICAL :: use_current_correction
  INTEGER, DIMENSION(2*c_ndims) :: bc_field, bc_particle
  INTEGER :: restart_number, step
  CHARACTER(LEN=c_max_path_length) :: full_restart_filename, restart_filename

  TYPE particle_sort_element
    TYPE(particle), POINTER :: particle
  END TYPE particle_sort_element

  TYPE(particle_sort_element), POINTER, DIMENSION(:) :: coll_sort_array
  INTEGER :: coll_sort_array_size = 0

  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: coll_pairs
  REAL(num) :: coulomb_log
  LOGICAL :: coulomb_log_auto, use_collisions

  LOGICAL :: use_field_ionisation, use_collisional_ionisation
  LOGICAL :: use_multiphoton, use_bsi

  INTEGER :: maxwell_solver = c_maxwell_solver_yee

  !----------------------------------------------------------------------------
  ! Moving window
  !----------------------------------------------------------------------------
  LOGICAL :: move_window, inject_particles
  REAL(num) :: window_v_x
  REAL(num) :: window_start_time
  INTEGER :: bc_x_min_after_move
  INTEGER :: bc_x_max_after_move
  REAL(num) :: window_shift

#ifdef PHOTONS
  !----------------------------------------------------------------------------
  ! QED - Written by C. P. Ridgers
  !----------------------------------------------------------------------------
  REAL(num), ALLOCATABLE :: log_chi2(:), epsilon_split(:), p_energy(:,:)
  REAL(num), ALLOCATABLE :: log_hsokolov(:,:), p_photon_energy(:,:)
  REAL(num), ALLOCATABLE :: log_eta(:), log_chi(:,:)
  REAL(num), ALLOCATABLE :: log_tpair(:,:), chimin_table(:), log_omegahat(:,:)
  INTEGER :: n_photon, n_pair
  INTEGER :: n_sample_epsilon, n_sample_chi2, n_sample_h
  INTEGER :: n_sample_eta, n_sample_chi, n_sample_t

  ! These track which species should be the species used by the QED routines
  INTEGER :: photon_species = -1, trident_electron_species = -1
  INTEGER :: breit_wheeler_electron_species = -1
  INTEGER :: trident_positron_species = -1, breit_wheeler_positron_species = -1

  REAL(num) :: photon_energy_min = EPSILON(1.0_num)
  REAL(num) :: qed_start_time = 0.0_num
  LOGICAL :: produce_pairs = .FALSE., use_radiation_reaction = .TRUE.
  LOGICAL :: produce_photons = .FALSE., photon_dynamics = .FALSE.
  CHARACTER(LEN=string_length) :: qed_table_location
#endif
  LOGICAL :: use_qed = .FALSE.

  !----------------------------------------------------------------------------
  ! MPI data
  !----------------------------------------------------------------------------
  INTEGER :: coordinates(c_ndims), neighbour(-1:1)
  INTEGER :: x_coords, proc_x_min, proc_x_max
  INTEGER :: errcode, comm, tag, rank
  INTEGER :: nproc, nprocx
  INTEGER :: nprocdir(c_ndims)
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nx_each_rank
  INTEGER(i8), ALLOCATABLE, DIMENSION(:) :: npart_each_rank
  LOGICAL :: x_min_boundary, x_max_boundary

  !----------------------------------------------------------------------------
  ! domain and loadbalancing
  !----------------------------------------------------------------------------
  LOGICAL :: use_balance
  REAL(num) :: dlb_threshold
  INTEGER(i8), PARAMETER :: npart_per_it = 1000000
  REAL(num), DIMENSION(:), ALLOCATABLE :: x_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: xb_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: xb_offset_global
  ! The location of the processors
  INTEGER, DIMENSION(:), ALLOCATABLE :: cell_x_min, cell_x_max
  INTEGER, DIMENSION(:), ALLOCATABLE :: old_x_max
  INTEGER :: nx_global_min, nx_global_max
  INTEGER :: n_global_min(c_ndims), n_global_max(c_ndims)
  INTEGER :: balance_mode
  LOGICAL :: debug_mode

  !----------------------------------------------------------------------------
  ! Particle injectors
  !----------------------------------------------------------------------------
  TYPE injector_block

    INTEGER :: boundary
    INTEGER :: id
    INTEGER :: species
    INTEGER(i8) :: npart_per_cell

    TYPE(primitive_stack) :: density_function
    TYPE(primitive_stack) :: temperature_function(3)
    TYPE(primitive_stack) :: drift_function(3)

    REAL(num) :: t_start, t_end
    REAL(num) :: depth, dt_inject

    TYPE(injector_block), POINTER :: next
  END TYPE injector_block

  TYPE(injector_block), POINTER :: injector_x_min, injector_x_max

  !----------------------------------------------------------------------------
  ! laser boundaries
  !----------------------------------------------------------------------------
  TYPE laser_block
    ! Boundary to which laser is attached
    INTEGER :: boundary
    ! A unique id number for the laser (not used directly by EPOCH)
    ! Only used if hard coding time profiles
    INTEGER :: id
    REAL(num) :: profile
    REAL(num) :: phase
    REAL(num) :: current_integral_phase

    LOGICAL :: use_time_function, use_phase_function, use_profile_function
    LOGICAL :: use_omega_function
    TYPE(primitive_stack) :: time_function, phase_function, profile_function
    TYPE(primitive_stack) :: omega_function

    REAL(num) :: amp, omega, pol_angle, t_start, t_end
    INTEGER :: omega_func_type

    TYPE(laser_block), POINTER :: next
  END TYPE laser_block

  TYPE(laser_block), POINTER :: laser_x_min, laser_x_max
  INTEGER :: n_laser_x_min = 0, n_laser_x_max = 0
  LOGICAL, DIMENSION(2*c_ndims) :: add_laser = .FALSE.

  TYPE(jobid_type) :: jobid

  INTEGER(i4) :: run_date
  INTEGER(i8) :: defines

  REAL(num) :: walltime_start, real_walltime_start
  REAL(num) :: stop_at_walltime
  INTEGER :: stdout_frequency, check_stop_frequency
  LOGICAL :: check_walltime, print_eta_string

  LOGICAL, DIMENSION(c_dir_x:c_dir_z,0:c_stagger_max) :: stagger
  INTEGER(i8) :: push_per_field = 5

  ! Absorption diagnostic
  REAL(num) :: laser_inject_local = 0.0_num
  REAL(num) :: laser_absorb_local = 0.0_num
  REAL(num) :: laser_injected = 0.0_num
  REAL(num) :: laser_absorbed = 0.0_num
  LOGICAL :: dump_absorption = .FALSE.

  REAL(num) :: total_particle_energy = 0.0_num
  REAL(num) :: total_field_energy = 0.0_num

  !----------------------------------------------------------------------------
  ! custom particle loading - written by MP Tooley
  !----------------------------------------------------------------------------

  INTEGER :: n_custom_loaders = 0
  INTEGER, PARAMETER :: c_loader_chunk_size = 131072 ! 1MB/variable data chunks

  TYPE custom_particle_loader
    INTEGER :: species_id

    ! Position Data
    CHARACTER(LEN=string_length) :: x_data
    INTEGER(KIND=MPI_OFFSET_KIND) :: x_data_offset

    ! Momentum Data
    CHARACTER(LEN=string_length) :: px_data
    CHARACTER(LEN=string_length) :: py_data
    CHARACTER(LEN=string_length) :: pz_data
    INTEGER(KIND=MPI_OFFSET_KIND) :: px_data_offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: py_data_offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: pz_data_offset
    LOGICAL :: px_data_given
    LOGICAL :: py_data_given
    LOGICAL :: pz_data_given
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS)
    ! Weight data
    CHARACTER(LEN=string_length) :: w_data
    INTEGER(KIND=MPI_OFFSET_KIND) :: w_data_offset
#endif
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
    ! ID data
    CHARACTER(LEN=string_length) :: id_data
    INTEGER(KIND=MPI_OFFSET_KIND) :: id_data_offset
    LOGICAL :: id_data_given
    LOGICAL :: id_data_4byte
#endif
  END TYPE custom_particle_loader

  TYPE(custom_particle_loader), DIMENSION(:), POINTER :: custom_loaders_list

END MODULE shared_data
