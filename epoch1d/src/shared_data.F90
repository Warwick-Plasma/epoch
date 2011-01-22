! ****************************************************************
! All global variables defined here (cf F77 COMMON block).
! All the names in here are public provided the MODULE is USE'd
! ****************************************************************

MODULE constants

  USE sdf

  IMPLICIT NONE

  INTEGER, PARAMETER :: num = KIND(1.d0)
  INTEGER, PARAMETER :: dbl = KIND(1.d0)
  REAL(num), PARAMETER :: c_non_zero = TINY(1.0_num)
  REAL(num), PARAMETER :: largest_number = HUGE(1.0_num)

  INTEGER, PARAMETER :: c_ndims = 1

  ! Boundary type codes
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

  ! Boundary location codes
  INTEGER, PARAMETER :: c_bd_x_min = 1
  INTEGER, PARAMETER :: c_bd_x_max = 2
  INTEGER, PARAMETER :: c_bd_y_min = 3
  INTEGER, PARAMETER :: c_bd_y_max = 4
  INTEGER, PARAMETER :: c_bd_z_min = 5
  INTEGER, PARAMETER :: c_bd_z_max = 6

  ! Error codes
  INTEGER, PARAMETER :: c_err_none = 0
  INTEGER, PARAMETER :: c_err_unknown_block = 1
  INTEGER, PARAMETER :: c_err_unknown_element = 2
  INTEGER, PARAMETER :: c_err_preset_element = 4
  INTEGER, PARAMETER :: c_err_preset_element_use_later = 8
  INTEGER, PARAMETER :: c_err_bad_value = 16
  INTEGER, PARAMETER :: c_err_missing_elements = 32
  INTEGER, PARAMETER :: c_err_terminate = 64
  INTEGER, PARAMETER :: c_err_required_element_not_set = 128
  INTEGER, PARAMETER :: c_err_pp_options_wrong = 256
  INTEGER, PARAMETER :: c_err_bad_array_length = 512
  INTEGER, PARAMETER :: c_err_other = 1024

  INTEGER, PARAMETER :: c_ds_deck = 1
  INTEGER, PARAMETER :: c_ds_ic = 2
  INTEGER, PARAMETER :: c_ds_eio = 3

  ! IO codes
  INTEGER, PARAMETER :: c_io_never = 0
  INTEGER, PARAMETER :: c_io_always = 1
  INTEGER, PARAMETER :: c_io_full = 2
  INTEGER, PARAMETER :: c_io_restartable = 4
  INTEGER, PARAMETER :: c_io_species = 8
  INTEGER, PARAMETER :: c_io_no_intrinsic = 16
  INTEGER, PARAMETER :: c_io_averaged = 32
  INTEGER, PARAMETER :: c_io_snapshot = 64
  ! domain codes
  INTEGER, PARAMETER :: c_do_full = 0
  INTEGER, PARAMETER :: c_do_decomposed = 1

  ! Load balance codes
  INTEGER, PARAMETER :: c_lb_x = 1
  INTEGER, PARAMETER :: c_lb_all = c_lb_x
  INTEGER, PARAMETER :: c_lb_auto = c_lb_all + 1

  ! Taken from http://physics.nist.gov/cuu/Constants
  REAL(num), PARAMETER :: pi = 3.141592653589793238462643383279502884_num
  REAL(num), PARAMETER :: q0 = 1.602176487e-19_num ! C (+/- 4e-27)
  REAL(num), PARAMETER :: m0 = 9.10938215e-31_num ! kg (+/- 4.5e-38)
  REAL(num), PARAMETER :: c  = 2.99792458e8_num   ! m/s^2 (exact)
  REAL(num), PARAMETER :: kb = 1.3806504e-23_num  ! J/K (+/- 2.4e-29)
  REAL(num), PARAMETER :: mu0 = 4.e-7_num * pi ! N/A^2 (exact)
  REAL(num), PARAMETER :: epsilon0 = 1.0_num / mu0 / c**2 ! F/m (exact)
  REAL(num), PARAMETER :: h_planck = 6.62606896e-34_num ! J s (+/- 3.3e-41)
  REAL(num), PARAMETER :: ev = q0 ! J

  ! direction parameters
  INTEGER, PARAMETER :: c_dir_x = 1
  INTEGER, PARAMETER :: c_dir_y = 2
  INTEGER, PARAMETER :: c_dir_z = 4
  INTEGER, PARAMETER :: c_dir_px = 8
  INTEGER, PARAMETER :: c_dir_py = 16
  INTEGER, PARAMETER :: c_dir_pz = 32
  INTEGER, PARAMETER :: c_dir_en = 64
  INTEGER, PARAMETER :: c_dir_gamma_m1 = 128

  ! define flags
  INTEGER(4), PARAMETER :: c_def_particle_debug = 1
  INTEGER(4), PARAMETER :: c_def_field_debug = 2
  INTEGER(4), PARAMETER :: c_def_particle_shape_bspline3 = 4
  INTEGER(4), PARAMETER :: c_def_split_part_after_push = 8
  INTEGER(4), PARAMETER :: c_def_per_particle_weight = 16
  INTEGER(4), PARAMETER :: c_def_particle_count_update = 32
  INTEGER(4), PARAMETER :: c_def_tracer_particles = 64
  INTEGER(4), PARAMETER :: c_def_particle_probes = 128
  INTEGER(4), PARAMETER :: c_def_per_particle_chargemass = 256
  INTEGER(4), PARAMETER :: c_def_particle_ionise = 512
  INTEGER(4), PARAMETER :: c_def_high_order_smoothing = 1024
  INTEGER(4), PARAMETER :: c_def_particle_shape_tophat = 2048
  INTEGER(4), PARAMETER :: c_def_parser_debug = 4096

  ! constants defining the maximum number of dimensions and directions
  ! in a distribution function
  INTEGER, PARAMETER :: c_df_maxdirs = 5 + c_ndims
  INTEGER, PARAMETER :: c_df_maxdims = 3

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
  INTEGER, PARAMETER :: string_length = 128

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

  ! Constants refering to grid properties
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
  INTEGER, PARAMETER :: c_const_time = 43

  INTEGER, PARAMETER :: c_const_io_never = 44
  INTEGER, PARAMETER :: c_const_io_always = 45
  INTEGER, PARAMETER :: c_const_io_full = 46
  INTEGER, PARAMETER :: c_const_io_restartable = 47
  INTEGER, PARAMETER :: c_const_io_species = 48
  INTEGER, PARAMETER :: c_const_io_average = 49
  INTEGER, PARAMETER :: c_const_io_snapshot = 50
  INTEGER, PARAMETER :: c_const_dir_x = 51
  INTEGER, PARAMETER :: c_const_dir_y = 52
  INTEGER, PARAMETER :: c_const_dir_z = 53
  INTEGER, PARAMETER :: c_const_dir_px = 54
  INTEGER, PARAMETER :: c_const_dir_py = 55
  INTEGER, PARAMETER :: c_const_dir_pz = 56
  INTEGER, PARAMETER :: c_const_dir_en = 57
  INTEGER, PARAMETER :: c_const_dir_gamma_m1 = 58
  INTEGER, PARAMETER :: c_const_nx = 59

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

  INTEGER, PARAMETER :: c_func_custom_lowbound = 4096

  ! Associativity constants
  INTEGER, PARAMETER :: c_assoc_a = 1
  INTEGER, PARAMETER :: c_assoc_la = 2
  INTEGER, PARAMETER :: c_assoc_ra = 3

  INTEGER, PARAMETER :: num_ops = 12
  INTEGER, DIMENSION(num_ops), PARAMETER :: &
      opcode_precedence = (/1, 1, 2, 2, 3, 4, 1, 1, 1, 2, 2, 5/)
  INTEGER, DIMENSION(num_ops), PARAMETER :: &
      opcode_assoc = (/c_assoc_a, c_assoc_la, c_assoc_a, c_assoc_la, &
          c_assoc_la, c_assoc_a, c_assoc_a, c_assoc_a, c_assoc_a, c_assoc_a, &
          c_assoc_a, c_assoc_ra/)

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
  END TYPE primitive_stack

  TYPE deck_constant
    CHARACTER(LEN=string_length) :: name
    TYPE(primitive_stack) :: execution_stream
  END TYPE deck_constant

  INTEGER :: n_deck_constants = 0
  TYPE(deck_constant), DIMENSION(:), ALLOCATABLE :: deck_constant_list

END MODULE shared_parser_data



MODULE shared_data

  USE sdf_job_info
  USE constants
  USE shared_parser_data
  USE mpi

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

  ! The order for the spline interpolation used as a particle representation.
#ifdef PARTICLE_SHAPE_BSPLINE3
  INTEGER, PARAMETER :: sf_min = -2
  INTEGER, PARAMETER :: sf_max =  2
#elif  PARTICLE_SHAPE_TOPHAT
  INTEGER, PARAMETER :: sf_min =  0
  INTEGER, PARAMETER :: sf_max =  1
#else
  INTEGER, PARAMETER :: sf_min = -1
  INTEGER, PARAMETER :: sf_max =  1
#endif

  ! Object representing a particle
  ! If you add or remove from this section then you *must* update the
  ! particle pack and unpack routines
  TYPE particle
    REAL(num), DIMENSION(3) :: part_p
    REAL(num) :: part_pos
#ifdef PER_PARTICLE_WEIGHT
    REAL(num) :: weight
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
  END TYPE particle

  ! Object representing a collection of particles
  ! Used internally by the MPI particle transfer code
  TYPE particle_list
    TYPE(particle), POINTER :: head
    TYPE(particle), POINTER :: tail
    INTEGER(KIND=8) :: count
    ! Pointer is safe if the particles in it are all unambiguously linked
    LOGICAL :: safe

    TYPE(particle_list), POINTER :: next, prev
  END TYPE particle_list

  ! Object representing a particle species
  TYPE particle_family
    ! Core properties
    CHARACTER(string_length) :: name
    TYPE(particle_family), POINTER :: next, prev
    INTEGER :: id
    LOGICAL :: dump

    REAL(num) :: charge
    REAL(num) :: mass
    INTEGER(KIND=8) :: count
    TYPE(particle_list) :: attached_list

#ifdef TRACER_PARTICLES
    LOGICAL :: tracer
#endif

    ! particle cell division
#ifdef SPLIT_PARTICLES_AFTER_PUSH
    INTEGER(KIND=8) :: global_count
    LOGICAL :: split
    INTEGER(KIND=8) :: npart_max
    ! Secondary list
    TYPE(particle_list), DIMENSION(:), POINTER :: secondary_list
#endif

    ! Injection of particles
    INTEGER(KIND=8) :: npart_per_cell
    REAL(num) :: density
    REAL(num), DIMENSION(3) :: temperature

    ! Species_ionisation
#ifdef PARTICLE_IONISE
    LOGICAL :: ionise
    INTEGER :: ionise_to_species
    INTEGER :: release_species
    REAL(num) :: critical_field
    REAL(num) :: ionisation_energy
#endif
    ! Attached probes for this species
#ifdef PARTICLE_PROBES
    TYPE(particle_probe), POINTER :: attached_probes
#endif
  END TYPE particle_family

  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  ! Represents the initial conditions of a species
  TYPE initial_condition_block
    REAL(num), DIMENSION(:), POINTER :: rho
    REAL(num), DIMENSION(:,:), POINTER :: temp
    REAL(num), DIMENSION(:,:), POINTER :: drift

    REAL(num) :: minrho
    REAL(num) :: maxrho
  END TYPE initial_condition_block

  INTEGER :: deck_state
  TYPE(initial_condition_block), DIMENSION(:), ALLOCATABLE :: initial_conditions

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

    ! Pointer to next distribution function
    TYPE(distribution_function_block), POINTER :: next
  END TYPE distribution_function_block
  TYPE(distribution_function_block), POINTER :: dist_fns

#ifdef PARTICLE_PROBES
  TYPE particle_probe
    REAL(num) :: probe_point
    LOGICAL :: left_to_right
    REAL(num) :: ek_min, ek_max
    CHARACTER(LEN=string_length) :: name

    TYPE(particle_family), POINTER :: probe_species
    TYPE(particle_list) :: sampled_particles
    TYPE(particle_probe), POINTER :: next
    INTEGER :: dump
  END TYPE particle_probe
#endif

  !----------------------------------------------------------------------------
  ! Core code
  !----------------------------------------------------------------------------
  INTEGER :: mpireal = MPI_DOUBLE_PRECISION
  INTEGER :: realsize

  INTEGER :: nx
  INTEGER :: nx_global
  INTEGER(KIND=8) :: npart_global
  INTEGER :: nprocx
  INTEGER :: nsteps, n_species = -1
  REAL(num), ALLOCATABLE, DIMENSION(:) :: ex, ey, ez, bx, by, bz, jx, jy, jz
  REAL(num), ALLOCATABLE, DIMENSION(:) :: wk_array

  TYPE(particle_family), DIMENSION(:), POINTER :: particle_species

  REAL(num), ALLOCATABLE, DIMENSION(:) :: x

  INTEGER, PARAMETER :: data_dir_max_length = 64
  CHARACTER(LEN=data_dir_max_length) :: data_dir
  INTEGER, PARAMETER :: c_max_string_length = 64

  LOGICAL :: neutral_background = .TRUE.
  LOGICAL :: use_random_seed = .FALSE.

  REAL(num) :: dt, t_end, time, dt_multiplier, dt_laser, dt_plasma_frequency
  REAL(num) :: dt_snapshots, dt_min_average
  REAL(num) :: length_x, dx, x_min, x_max
  REAL(num) :: x_min_local, x_max_local, length_x_local
  REAL(num), DIMENSION(:), ALLOCATABLE :: x_mins, x_maxs

  REAL(num) :: total_ohmic_heating = 0.0_num
  REAL(num) :: weight

  LOGICAL :: ic_from_restart = .FALSE.
  INTEGER, DIMENSION(2*c_ndims) :: bc_field, bc_particle
  INTEGER :: restart_snapshot

  !----------------------------------------------------------------------------
  ! Moving window
  !----------------------------------------------------------------------------
  LOGICAL :: move_window, inject_particles
  LOGICAL :: window_started
  REAL(num) :: window_shift_fraction
  REAL(num) :: window_v_x
  REAL(num) :: window_start_time
  INTEGER :: bc_x_min_after_move
  INTEGER :: bc_x_max_after_move
  REAL(num) :: window_shift
  TYPE(particle_list) :: ejected_particles

  !----------------------------------------------------------------------------
  ! MPI data
  !----------------------------------------------------------------------------
  INTEGER :: rank, proc_x_min, proc_x_max
  INTEGER :: coordinates(1), neighbour(-1:1)
  INTEGER :: errcode, comm, tag, nproc, icycle_max = 1000000
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nx_each_rank
  INTEGER(KIND=8), ALLOCATABLE, DIMENSION(:) :: npart_each_rank

  !----------------------------------------------------------------------------
  ! domain and loadbalancing
  !----------------------------------------------------------------------------
  LOGICAL :: dlb
  REAL(num) :: dlb_threshold
  INTEGER(KIND=8), PARAMETER :: npart_per_it = 1000000
  REAL(num), DIMENSION(:), ALLOCATABLE :: x_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: xb_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: xb_offset_global
  ! The location of the processors
  INTEGER, DIMENSION(:), ALLOCATABLE :: cell_x_min, cell_x_max
  INTEGER :: nx_global_min, nx_global_max
  INTEGER :: balance_mode
  LOGICAL :: debug_mode

  !----------------------------------------------------------------------------
  ! file handling
  !----------------------------------------------------------------------------
  INTEGER :: subtype_field, subarray_field
  INTEGER(KIND=MPI_OFFSET_KIND) :: initialdisp
  INTEGER :: full_dump_every, restart_dump_every
  INTEGER, PARAMETER :: num_vars_to_dump = 29
  INTEGER, DIMENSION(num_vars_to_dump) :: dumpmask
  INTEGER :: output_file
  LOGICAL :: force_final_to_be_restartable
  LOGICAL :: use_offset_grid
  INTEGER :: n_zeros = 4

  !----------------------------------------------------------------------------
  ! Time averaged IO
  !----------------------------------------------------------------------------
  TYPE averaged_data_block
    REAL(num), DIMENSION(:,:), POINTER :: array
    REAL(num) :: average_over_real_time
    REAL(num) :: real_time_after_average
  END TYPE averaged_data_block
  TYPE(averaged_data_block), DIMENSION(num_vars_to_dump), SAVE :: averaged_data
  INTEGER :: min_cycles_per_average = -1
  LOGICAL :: any_average = .FALSE.
  REAL(num) :: average_time = -1.0_num

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

    LOGICAL :: use_time_function
    TYPE(primitive_stack) :: time_function

    REAL(num) :: amp, freq, k, pol_angle, angle, t_start, t_end

    TYPE(laser_block), POINTER :: next
  END TYPE laser_block

  TYPE(laser_block), POINTER :: laser_x_min, laser_x_max
  INTEGER :: n_laser_x_min, n_laser_x_max

  TYPE(jobid_type) :: jobid

  INTEGER(4) :: run_date
  INTEGER(8) :: defines

  REAL(num) :: walltime_start
  INTEGER :: stdout_frequency
  INTEGER(KIND=MPI_OFFSET_KIND), DIMENSION(:), ALLOCATABLE :: &
      particle_file_lengths, particle_file_offsets

  INTEGER, DIMENSION(3,0:c_stagger_max) :: stagger

END MODULE shared_data
