! Copyright (C) 2009-2018 University of Warwick
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

  INTEGER, PARAMETER :: c_ndims = 3
  INTEGER, PARAMETER :: c_ndirs = 3 ! Number of directions for velocity

  ! File unit for deck parser diagnostic output
  INTEGER, PARAMETER :: stdout = 6
  INTEGER, PARAMETER :: du = 40
  INTEGER, PARAMETER :: lu = 41
  INTEGER, PARAMETER :: duc = 42
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
  INTEGER, PARAMETER :: c_bc_mixed = 14

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
  INTEGER, PARAMETER :: c_err_window = 2**17

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
  INTEGER, PARAMETER :: c_maxwell_solver_custom = -1
  INTEGER, PARAMETER :: c_maxwell_solver_yee = 0
  INTEGER, PARAMETER :: c_maxwell_solver_lehe = 1
  INTEGER, PARAMETER :: c_maxwell_solver_lehe_x = 2
  INTEGER, PARAMETER :: c_maxwell_solver_lehe_y = 3
  INTEGER, PARAMETER :: c_maxwell_solver_lehe_z = 4
  INTEGER, PARAMETER :: c_maxwell_solver_cowan = 5
  INTEGER, PARAMETER :: c_maxwell_solver_pukhov = 6

  ! Particle distribution type codes
  INTEGER, PARAMETER :: c_ic_df_thermal = 1
  INTEGER, PARAMETER :: c_ic_df_relativistic_thermal = 2

  ! domain codes
  INTEGER, PARAMETER :: c_do_full = 0
  INTEGER, PARAMETER :: c_do_decomposed = 1

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
  INTEGER(i8), PARAMETER :: c_def_deltaf_method = 2**20
  INTEGER(i8), PARAMETER :: c_def_deltaf_debug = 2**21
  INTEGER(i8), PARAMETER :: c_def_work_done_integrated = 2**22
  INTEGER(i8), PARAMETER :: c_def_hc_push = 2**23

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

#ifdef PARTICLE_ID4
  INTEGER, PARAMETER :: idkind = i4
#else
  INTEGER, PARAMETER :: idkind = i8
#endif

END MODULE constants
