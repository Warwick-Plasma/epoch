!****************************************************************
! All global variables defined here (cf F77 COMMON block).
! All the names in here are public provided the MODULE is USE'd
!**************************************************************** 

MODULE constants
  IMPLICIT NONE
  INTEGER, PARAMETER :: num = KIND(1.d0)
  INTEGER, PARAMETER :: dbl = KIND(1.d0)
  REAL(num), PARAMETER :: pi = 3.14159265358979323_num
  REAL(num), PARAMETER :: c_non_zero = TINY(1.0_num)
  REAL(num), PARAMETER :: largest_number = HUGE(1.0_num)

  !Boundary type codes
  INTEGER, PARAMETER :: c_bc_periodic = 1, c_bc_other = 2, c_bc_simple_laser = 3, c_bc_simple_outflow = 4
  INTEGER, PARAMETER :: c_bc_open = 5, c_bc_dump=6, c_bc_zero_gradient=7, c_bc_clamp=8, c_bc_reflect=9

  !Boundary location codes
  INTEGER, PARAMETER :: c_bd_left=1, c_bd_right=2, c_bd_up=3, c_bd_down=4, c_bd_front=5, c_bd_back=6

  INTEGER, PARAMETER :: c_version = 1, c_revision = 3

  INTEGER,PARAMETER ::c_err_none=0,c_err_unknown_block=1,c_err_unknown_element=2
  INTEGER,PARAMETER ::c_err_preset_element=4,c_err_preset_element_use_later=8
  INTEGER,PARAMETER ::c_err_bad_value=16,c_err_missing_elements=32,c_err_terminate=64
  INTEGER,PARAMETER ::c_err_required_element_not_set=128, c_err_pp_options_wrong=256
  INTEGER,PARAMETER ::c_err_bad_array_length=512
  INTEGER,PARAMETER ::c_err_other=1024

  !IC codes
  !This is a bitmask, remember that
  INTEGER, PARAMETER :: c_ic_early_internal=1, c_ic_late_internal=2,  c_ic_external=4, c_ic_manual=8, c_ic_restart=16
  INTEGER, PARAMETER :: c_ic_autoload = c_ic_early_internal + c_ic_late_internal + c_ic_external

  INTEGER, PARAMETER :: c_ds_deck=1, c_ds_ic=2, c_ds_eio=3

  !IO codes
  INTEGER, PARAMETER :: c_io_never=0, c_io_always=1, c_io_full=2, c_io_restartable=4, c_io_species=8, c_io_no_intrinsic=16
  !domain codes
  INTEGER, PARAMETER :: c_do_full=0,c_do_decomposed=1

  !Load balance codes
  INTEGER,PARAMETER :: c_lb_x=1, c_lb_y=2, c_lb_z=4, c_lb_auto=8, c_lb_all=c_lb_x+c_lb_y+c_lb_z

  REAL(num),PARAMETER :: q0 = 1.60217646e-19_num !c
  REAL(num),PARAMETER :: m0 = 9.10938188e-31_num !kg
  REAL(num),PARAMETER :: c  = 2.99792458e8_num   !ms^(-2)
  REAL(num),PARAMETER :: kb = 1.3806503e-23_num  !m^2kgs(-2)K^(-1)
  REAL(num),PARAMETER :: epsilon0 = 8.85418782e-12_num
  REAL(num),PARAMETER :: mu0 = 1.0_num/(c**2*epsilon0)

  !direction parameters
  INTEGER, PARAMETER :: c_dir_x=1, c_dir_y=2, c_dir_z=4, c_dir_px=8, c_dir_py=16, c_dir_pz=32

  !Length of a standard string
  INTEGER,PARAMETER :: string_length=400

END MODULE constants

MODULE shared_parser_data

  USE constants

  IMPLICIT NONE

  INTEGER,PARAMETER :: c_char_numeric=1, c_char_alpha=2, c_char_delimiter=3, c_char_space = 4
  INTEGER,PARAMETER :: c_char_opcode=5,c_char_unknown=1024

  INTEGER,PARAMETER :: c_prc_not_this_type=0

  !block type constants
  INTEGER,PARAMETER :: c_pt_variable=1, c_pt_constant=2, c_pt_operator=3, c_pt_function=4, c_pt_parenthesis=5, c_pt_separator=6, c_pt_character=7
  INTEGER,PARAMETER :: c_pt_deferred_execution_object=8
  INTEGER,PARAMETER :: c_pt_bad=1024, c_pt_null=1025

  !Opcode constants
  INTEGER,PARAMETER :: c_opcode_plus=1, c_opcode_minus=2, c_opcode_times=3
  INTEGER,PARAMETER :: c_opcode_divide=4, c_opcode_power=5, c_opcode_expo=6
  INTEGER,PARAMETER :: c_opcode_lt=7, c_opcode_gt=8, c_opcode_eq=9
  INTEGER,PARAMETER :: c_opcode_and=10, c_opcode_or=11, c_opcode_unary_minus=12

  INTEGER,PARAMETER :: c_paren_left_bracket=1, c_paren_right_bracket=2

  !Actual constants
  INTEGER,PARAMETER :: c_const_pi=1, c_const_kb=2, c_const_me=3, c_const_qe=4
  INTEGER,PARAMETER :: c_const_eps0=5, c_const_mu0=6, c_const_c=7

  !Constants refering to grid properties
  INTEGER,PARAMETER :: c_const_x=25, c_const_y=26, c_const_z=27
  INTEGER,PARAMETER :: c_const_lx=28, c_const_ly=29, c_const_lz=30
  INTEGER,PARAMETER :: c_const_dx=31, c_const_dy=32, c_const_dz=33
  INTEGER,PARAMETER :: c_const_start_x=34, c_const_start_y=35, c_const_start_z=36
  INTEGER,PARAMETER :: c_const_end_x=37, c_const_end_y=38, c_const_end_z=39
  INTEGER,PARAMETER :: c_const_ix=40, c_const_iy=41, c_const_iz=42
  INTEGER,PARAMETER :: c_const_time=43

  INTEGER,PARAMETER :: c_const_io_never=44, c_const_io_always=45, c_const_io_full=46, c_const_io_restartable=47, c_const_io_species=48
  INTEGER,PARAMETER :: c_const_dir_x=49, c_const_dir_y=50, c_const_dir_z=51, c_const_dir_px=52
  INTEGER,PARAMETER :: c_const_dir_py=53, c_const_dir_pz=54, c_const_r_xy=55, c_const_r_yz=56, c_const_r_xz=57

  !Constants for initial conditions
  INTEGER,PARAMETER :: c_const_autoearly=20, c_const_autolate=21, c_const_external=22
  INTEGER,PARAMETER :: c_const_manual=23, c_const_restart=24

  !Custom constants
  INTEGER,PARAMETER :: c_const_deck_lowbound=4096
  INTEGER,PARAMETER :: c_const_custom_lowbound=8192

  INTEGER,PARAMETER :: c_func_sine=1, c_func_cosine=2, c_func_tan=3, c_func_exp=4
  INTEGER,PARAMETER :: c_func_arcsine=10, c_func_arccosine=11, c_func_arctan=12
  INTEGER,PARAMETER :: c_func_neg=13, c_func_if=14, c_func_floor=15, c_func_ceil=16
  INTEGER,PARAMETER :: c_func_nint=17, c_func_rho=18, c_func_tempx=19, c_func_tempy=20
  INTEGER,PARAMETER :: c_func_tempz=21, c_func_interpolate=22, c_func_tanh=23
  INTEGER,PARAMETER :: c_func_sinh=24, c_func_cosh=25, c_func_ex=26
  INTEGER,PARAMETER :: c_func_ey=27, c_func_ez=28, c_func_bx=29
  INTEGER,PARAMETER :: c_func_by=30, c_func_bz=31, c_func_sqrt=32
  INTEGER,PARAMETER :: c_func_gauss=33, c_func_abs=34

  INTEGER,PARAMETER :: c_func_custom_lowbound=4096

  !Associativity constants
  INTEGER,PARAMETER :: c_assoc_a=1, c_assoc_la=2, c_assoc_ra=3

  INTEGER,PARAMETER :: num_ops=12
  INTEGER, DIMENSION(num_ops), PARAMETER :: opcode_precedence = (/1,1,2,2,3,4,1,1,1,2,2,5/)
  INTEGER, DIMENSION(num_ops), PARAMETER :: opcode_assoc = (/c_assoc_a,c_assoc_la,c_assoc_a,c_assoc_la,c_assoc_la,&
       c_assoc_a,c_assoc_a,c_assoc_a,c_assoc_a,c_assoc_a,c_assoc_a,c_assoc_ra/)

  INTEGER,PARAMETER :: stack_size=10000

  TYPE :: stack_element
     INTEGER :: ptype
     INTEGER :: data
     REAL(num) :: numerical_data
#ifdef PARSER_DEBUG
	  CHARACTER(len=string_length) :: text
#endif
  END TYPE stack_element

  TYPE :: primitive_stack
     TYPE(stack_element),DIMENSION(stack_size) :: data
     INTEGER :: stack_point
  END TYPE primitive_stack

  TYPE :: deck_constant
     REAL(num) :: value
     CHARACTER(len=string_length) :: name
  END TYPE deck_constant

  TYPE :: deferred_execution_object
     CHARACTER(len=string_length) :: name
     TYPE(primitive_stack) :: execution_stream
  END TYPE deferred_execution_object

  INTEGER :: n_deck_constants=0
  INTEGER :: n_deferred_execution_objects=0
  TYPE(deck_constant),DIMENSION(:),ALLOCATABLE :: deck_constant_list
  TYPE(deferred_execution_object),DIMENSION(:),ALLOCATABLE :: deferred_objects

END MODULE shared_parser_data


MODULE shared_data
  USE constants
  USE shared_parser_data
  IMPLICIT NONE

  !---------------------------------------------------------------------------------------
  !string handling
  !---------------------------------------------------------------------------------------
  INCLUDE 'mpif.h'
  CHARACTER(len=string_length) :: blank
  TYPE :: string_type
     CHARACTER(string_length) :: value
  END TYPE string_type
  CHARACTER(len=string_length) :: extended_error_string

  !---------------------------------------------------------------------------------------
  !Particles
  !---------------------------------------------------------------------------------------

  !The order for the spline interpolation used as a particle representation.
#ifdef SPLINE_FOUR
  INTEGER,PARAMETER :: sf_order=2
#else
  INTEGER,PARAMETER :: sf_order=1
#endif

  !Object representing a particle
  TYPE :: particle
     REAL(num), DIMENSION(3) :: part_p
     REAL(num),DIMENSION(3) :: part_pos
#ifdef PER_PARTICLE_WEIGHT
     REAL(num) :: weight
#endif
#ifdef PER_PARTICLE_CHARGEMASS
     REAL(num) :: charge
     REAL(num) :: mass
#endif
     TYPE(particle),POINTER :: next, prev
#ifdef PART_DEBUG
     INTEGER :: processor
     INTEGER :: processor_at_t0
#endif
  END TYPE particle

  !Object representing a collection of particles
  !Used internally by the MPI particle transfer code
  TYPE :: particle_list
     TYPE(particle),POINTER :: head
     TYPE(particle),POINTER :: tail
     INTEGER(KIND=8) :: count
     !Pointer is safe if the particles in it are all unambiguously linked
     LOGICAL :: safe

     TYPE(particle_list), POINTER :: next, prev
  END TYPE particle_list

  !Object representing a particle species
  TYPE :: particle_family
     !Core properties
     CHARACTER(string_length) :: name
     TYPE(particle_family),POINTER :: next,prev
     INTEGER :: id
     LOGICAL :: dump

     REAL(num) :: charge
     REAL(num) :: mass
     INTEGER(KIND=8) :: count
     TYPE(particle_list) :: attached_list

#ifdef TRACER_PARTICLES
     LOGICAL :: tracer
#endif

     !particle cell division
#ifdef PARTICLE_CELL_DIVISION
     INTEGER(KIND=8) :: global_count
     LOGICAL :: split
     INTEGER(KIND=8) :: npart_max
#endif
     !Secondary list
#ifdef USE_SECONDARY_LIST
     TYPE(particle_list),DIMENSION(:,:,:),POINTER :: secondary_list
#endif

     !Injection of particles
     INTEGER(KIND=8) :: window_npart_per_cell
     REAL(num),DIMENSION(:,:),ALLOCATABLE :: window_density
     REAL(num),DIMENSION(:,:,:),ALLOCATABLE :: window_temperature

     !Species_ionisation
#ifdef PART_IONISE
     LOGICAL :: ionise
     INTEGER :: ionise_to_species
     INTEGER :: release_species
     REAL(num) :: critical_field
     REAL(num) :: ionisation_energy
#endif
     !Attached probes for this species
#ifdef PARTICLE_PROBES
     TYPE(particle_probe),POINTER :: attached_probes
#endif
  END TYPE particle_family

  !---------------------------------------------------------------------------------------
  !Initial conditions
  !---------------------------------------------------------------------------------------
  !Represents the initial conditions of a species
  TYPE :: initial_condition_block

     REAL(num),DIMENSION(:,:,:),POINTER :: rho
     REAL(num),DIMENSION(:,:,:,:),POINTER :: temp
     REAL(num),DIMENSION(3) :: drift

     REAL(num) :: minrho
     REAL(num) :: maxrho

  END TYPE initial_condition_block

  INTEGER :: deck_state
  TYPE(initial_condition_block),DIMENSION(:),ALLOCATABLE :: initial_conditions
  INTEGER :: ictype
  TYPE(string_type) :: icfile

  !---------------------------------------------------------------------------------------
  !Extended IO information
  !---------------------------------------------------------------------------------------
  LOGICAL :: use_extended_io
  CHARACTER(len=string_length) :: extended_io_file
  !Represents a 2 or 3D distribution
  TYPE :: distribution_function_block

     CHARACTER(len=string_length) :: name

     !The number of dimensions left
     INTEGER :: ndims

     !The dumpmask for the distribution function
     INTEGER :: dumpmask

     !Whether or not to store the range returned from the distribtution function code
     !This allows the code to store auto determined ranges(experimental)
     LOGICAL :: store_ranges

     INTEGER,DIMENSION(3) :: directions
     REAL(num),DIMENSION(3,2) :: ranges
     INTEGER,DIMENSION(3) :: resolution
     LOGICAL,DIMENSION(:),ALLOCATABLE :: use_species
     REAL(num),DIMENSION(6,2) :: restrictions
     LOGICAL,DIMENSION(6) :: use_restrictions

     TYPE(distribution_function_block),POINTER :: next

  END TYPE distribution_function_block
  TYPE(distribution_function_block),POINTER :: dist_fns

#ifdef PARTICLE_PROBES
  TYPE :: particle_probe
     !The corner points
     REAL(num), DIMENSION(4,3) :: corner
     !The normal to the plane
     REAL(num), DIMENSION(3) :: normal
     REAL(num) :: ek_min, ek_max
     CHARACTER(len=string_length) :: name

     TYPE(particle_family), POINTER :: probe_species
     TYPE(particle_list) :: sampled_particles
     TYPE(particle_probe),POINTER :: next
     INTEGER :: dump
  END TYPE particle_probe
#endif

  !---------------------------------------------------------------------------------------
  !Core code
  !---------------------------------------------------------------------------------------
  INTEGER :: mpireal = MPI_DOUBLE_PRECISION

  INTEGER :: nx,ny,nz
  INTEGER :: nx_global, ny_global, nz_global
  INTEGER(8) :: npart_global
  INTEGER :: nprocx,nprocy,nprocz
  INTEGER :: nsteps,n_species=-1
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:)     :: ex,ey,ez,bx,by,bz,jx,jy,jz

  TYPE(particle_family),DIMENSION(:),POINTER :: particle_species
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:,:)   :: ekbar !temperature per species
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:,:)   :: ekbar_sum,sum_sq,ct !Used to calculate temperature

  REAL(num), ALLOCATABLE, DIMENSION(:)   ::  x,y,z

  INTEGER,PARAMETER :: data_dir_max_length=64
  CHARACTER(len=data_dir_max_length) :: data_dir

  LOGICAL :: neutral_background= .TRUE.

  REAL(num) :: dt, dt_plasma, t_end, time,dt_multiplier, dt_laser
  REAL(num) :: dt_snapshots
  REAL(num) :: length_x, dx, x_start, x_end, x_start_local, x_end_local, length_x_local
  REAL(num) :: length_y, dy, y_start, y_end, y_start_local, y_end_local, length_y_local
  REAL(num) :: length_z, dz, z_start, z_end, z_start_local, z_end_local, length_z_local
  REAL(num), DIMENSION(:), ALLOCATABLE :: x_starts, x_ends, y_starts, y_ends, z_starts, z_ends

  REAL(num) :: total_ohmic_heating = 0.0_num
  REAL(num) :: weight

  LOGICAL :: SAVE, restart,deckfile
  INTEGER :: xbc_right, xbc_left,ybc_up,ybc_down,zbc_front,zbc_back
  INTEGER :: xbc_right_particle, xbc_left_particle&
       ,ybc_up_particle,ybc_down_particle&
       ,zbc_front_particle,zbc_back_particle
  INTEGER :: xbc_right_field, xbc_left_field&
       ,ybc_up_field,ybc_down_field&
       ,zbc_front_field,zbc_back_field
  INTEGER :: restart_snapshot
  INTEGER(8) :: ix,iy,iz,ipart

  !---------------------------------------------------------------------------------------
  !Moving window
  !---------------------------------------------------------------------------------------
  LOGICAL :: move_window, inject_particles
  LOGICAL :: window_started
  REAL(num) :: window_shift_fraction
  REAL(num) :: window_v_x
  REAL(num) :: window_start_time
  INTEGER :: xbc_left_after_move
  INTEGER :: xbc_right_after_move
  REAL(num),DIMENSION(3) :: window_shift
  LOGICAL :: any_open
  TYPE(particle_list) :: ejected_particles

  !---------------------------------------------------------------------------------------
  ! MPI data
  !---------------------------------------------------------------------------------------
  INTEGER :: rank, left, right, up, down, front, back
  INTEGER :: coordinates(3),neighbour(-1:1,-1:1,-1:1)
  INTEGER :: errcode, comm, tag, nproc,icycle_max=1000000
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nx_each_rank, ny_each_rank, nz_each_rank
  INTEGER(8), ALLOCATABLE, DIMENSION(:) :: npart_each_rank
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: start_each_rank,end_each_rank

  !---------------------------------------------------------------------------------------
  !domain and loadbalancing
  !---------------------------------------------------------------------------------------
  LOGICAL :: dlb 
  REAL(num) :: dlb_threshold
  INTEGER(KIND=8), PARAMETER :: npart_per_it = 1000000
  REAL(num),DIMENSION(:),ALLOCATABLE :: x_global,y_global,z_global
  REAL(num),DIMENSION(:),ALLOCATABLE :: x_offset_global, y_offset_global,z_offset_global
  !The location of the processors
  INTEGER,DIMENSION(:),ALLOCATABLE :: cell_x_start, cell_x_end, cell_y_start, cell_y_end&
       ,cell_z_start,cell_z_end
  INTEGER :: balance_mode
  LOGICAL :: debug_mode

  !---------------------------------------------------------------------------------------
  ! file handling
  !---------------------------------------------------------------------------------------
  INTEGER :: subtype_field,subtype_particle_var,subtype_particle,subtype_particle_int
  INTEGER(KIND=MPI_OFFSET_KIND) :: initialdisp 
  INTEGER :: full_dump_every,restart_dump_every
  INTEGER, PARAMETER :: num_vars_to_dump = 28
  INTEGER, DIMENSION(num_vars_to_dump) :: dumpmask 
  INTEGER :: output_file
  LOGICAL :: force_final_to_be_restartable
  LOGICAL :: use_offset_grid
  INTEGER :: n_zeros=4

  !---------------------------------------------------------------------------------------
  !laser boundaries
  !---------------------------------------------------------------------------------------
  TYPE :: laser_block
     !Boundary to which laser is attached
     INTEGER :: direction
     !A unique id number for the laser (not used directly by EPOCH)
     !Only used if hard coding time profiles
     INTEGER :: id
     REAL(num), DIMENSION(:,:), POINTER :: profile
     REAL(num), DIMENSION(:,:), POINTER :: phase

     LOGICAL :: use_time_function
     TYPE(primitive_stack) :: time_function

     REAL(num) :: amp=0.0_num,freq=1.0_num, k=1.0_num
     REAL(num) :: pol=0.0_num, angle=0.0_num
     REAL(num) :: t_start=0.0_num, t_end=0.0_num

     TYPE(laser_block),POINTER :: next
  END TYPE laser_block
  TYPE(laser_block),POINTER :: laser_left, laser_right
  TYPE(laser_block),POINTER :: laser_up, laser_down
  TYPE(laser_block),POINTER :: laser_front, laser_back
  INTEGER :: n_laser_left, n_laser_right, n_laser_up, n_laser_down

END MODULE shared_data






