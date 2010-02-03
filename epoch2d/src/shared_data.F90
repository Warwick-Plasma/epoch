!****************************************************************
! All global variables defined here (cf F77 COMMON block).
! All the names in here are public provided the MODULE is USE'd
!**************************************************************** 


MODULE constants
  !These tell the code to setup secondary particles lists for each particle species
#ifdef PART_IONISE
#define USE_SECONDARY_LIST
#endif

#ifdef PARTICLE_CELL_DIVISION
#define USE_SECONDARY_LIST
#endif

#ifdef PART_IONISE_FULL
#define PART_IONISE
#endif
  IMPLICIT NONE
  INTEGER, PARAMETER :: num = KIND(1.D0)
  INTEGER, PARAMETER :: dbl = KIND(1.D0)
  REAL(num), PARAMETER :: pi = 3.14159265358979323_num
  REAL(num), PARAMETER :: none_zero = TINY(1.0_num)
  REAL(num), PARAMETER :: largest_number = HUGE(1.0_num)

  !Boundary type codes
  INTEGER, PARAMETER :: BC_PERIODIC = 1, BC_OTHER = 2, BC_SIMPLE_LASER = 3, BC_SIMPLE_OUTFLOW = 4
  INTEGER, PARAMETER :: BC_OPEN = 5, BC_DUMP=6, BC_ZERO_GRADIENT=7, BC_CLAMP=8, BC_REFLECT=9

  !Boundary location codes
  INTEGER, PARAMETER :: BD_LEFT=1, BD_RIGHT=2, BD_UP=3, BD_DOWN=4, BD_FRONT=5, BC_BACK=6

  INTEGER, PARAMETER :: Version = 1, Revision = 3

  !Error codes
  INTEGER,PARAMETER ::ERR_NONE=0,ERR_UNKNOWN_BLOCK=1,ERR_UNKNOWN_ELEMENT=2
  INTEGER,PARAMETER ::ERR_PRESET_ELEMENT=4,ERR_PRESET_ELEMENT_USE_LATER=8
  INTEGER,PARAMETER ::ERR_BAD_VALUE=16,ERR_MISSING_ELEMENTS=32,ERR_TERMINATE=64
  INTEGER,PARAMETER ::ERR_REQUIRED_ELEMENT_NOT_SET=128, ERR_PP_OPTIONS_WRONG=256
  INTEGER,PARAMETER ::ERR_BAD_ARRAY_LENGTH=512
  INTEGER,PARAMETER ::ERR_OTHER=2048

  !IC codes
  !This is a bitmask, remember that
  INTEGER, PARAMETER :: IC_EARLY_INTERNAL=1, IC_LATE_INTERNAL=2,  IC_EXTERNAL=4, IC_MANUAL=8, IC_RESTART=16
  INTEGER, PARAMETER :: IC_AUTOLOAD = IC_EARLY_INTERNAL + IC_LATE_INTERNAL + IC_EXTERNAL

  INTEGER, PARAMETER :: DS_DECK=1, DS_IC=2, DS_EIO=3

  !IO codes
  INTEGER, PARAMETER :: IO_NEVER=0, IO_ALWAYS=1, IO_FULL=2, IO_RESTARTABLE=4, IO_SPECIES=8, IO_AVERAGED=16
  INTEGER, PARAMETER :: IO_NO_INTRINSIC=32
  !Domain codes
  INTEGER, PARAMETER :: DO_FULL=0,DO_DECOMPOSED=1

  !Load balance codes
  INTEGER,PARAMETER :: LB_X=1, LB_Y=2, LB_AUTO=4, LB_BOTH=LB_X+LB_Y

  REAL(num),PARAMETER :: Q0 = 1.60217646e-19_num !C
  REAL(num),PARAMETER :: M0 = 9.10938188e-31_num !kg
  REAL(num),PARAMETER :: C  = 2.99792458e8_num   !ms^(-2)
  REAL(num),PARAMETER :: kb = 1.3806503e-23_num  !m^2kgs(-2)K^(-1)
  REAL(num),PARAMETER :: epsilon0 = 8.85418782e-12_num
  REAL(num),PARAMETER :: mu0 = 1.0_num/(C**2*epsilon0)
  REAL(num),PARAMETER :: h_planck=6.626068e-34_num
  REAL(num),PARAMETER :: ev = Q0 !J

  !Direction parameters
  INTEGER, PARAMETER :: DIR_X=1, DIR_Y=2, DIR_Z=4, DIR_PX=8, DIR_PY=16, DIR_PZ=32

  !Length of a standard string
  INTEGER,PARAMETER :: EntryLength=128

END MODULE constants

MODULE shared_parser_data

  USE constants

  IMPLICIT NONE

  INTEGER,PARAMETER :: CHAR_NUMERIC=1, CHAR_ALPHA=2, CHAR_DELIMITER=3, CHAR_SPACE = 4
  INTEGER,PARAMETER :: CHAR_OPCODE=5,CHAR_UNKNOWN=1024

  INTEGER,PARAMETER :: PRC_NOT_THIS_TYPE=0

  !Block type constants
  INTEGER,PARAMETER :: PT_VARIABLE=1, PT_CONSTANT=2, PT_OPERATOR=3, PT_FUNCTION=4, PT_PARENTHESIS=5, PT_SEPARATOR=6, PT_CHARACTER=7
  INTEGER,PARAMETER :: PT_DEFERRED_EXECUTION_OBJECT=8
  INTEGER,PARAMETER :: PT_BAD=1024, PT_NULL=1025

  !Opcode constants
  INTEGER,PARAMETER :: OPCODE_PLUS=1, OPCODE_MINUS=2, OPCODE_TIMES=3
  INTEGER,PARAMETER :: OPCODE_DIVIDE=4, OPCODE_POWER=5, OPCODE_EXPO=6
  INTEGER,PARAMETER :: OPCODE_LT=7, OPCODE_GT=8, OPCODE_EQ=9
  INTEGER,PARAMETER :: OPCODE_AND=10, OPCODE_OR=11, OPCODE_UNARY_MINUS=12

  INTEGER,PARAMETER :: PAREN_LEFT_BRACKET=1, PAREN_RIGHT_BRACKET=2

  !Actual constants
  INTEGER,PARAMETER :: CONST_PI=1, CONST_KB=2, CONST_ME=3, CONST_QE=4
  INTEGER,PARAMETER :: CONST_EPS0=5, CONST_MU0=6, CONST_C=7, CONST_EV=8
  INTEGER,PARAMETER :: CONST_KEV=9, CONST_MEV=10

  !Constants refering to grid properties
  INTEGER,PARAMETER :: CONST_X=25, CONST_Y=26, CONST_Z=27
  INTEGER,PARAMETER :: CONST_LX=28, CONST_LY=29, CONST_LZ=30
  INTEGER,PARAMETER :: CONST_DX=31, CONST_DY=32, CONST_DZ=33
  INTEGER,PARAMETER :: CONST_START_X=34, CONST_START_Y=35, CONST_START_Z=36
  INTEGER,PARAMETER :: CONST_END_X=37, CONST_END_Y=38, CONST_END_Z=39
  INTEGER,PARAMETER :: CONST_IX=40, CONST_IY=41, CONST_IZ=42
  INTEGER,PARAMETER :: CONST_TIME=43

  INTEGER,PARAMETER :: CONST_IO_NEVER=44, CONST_IO_ALWAYS=45, CONST_IO_FULL=46, CONST_IO_RESTARTABLE=47, CONST_IO_SPECIES=48
  INTEGER,PARAMETER :: CONST_DIR_X=49, CONST_DIR_Y=50, CONST_DIR_Z=51, CONST_DIR_PX=52
  INTEGER,PARAMETER :: CONST_DIR_PY=53, CONST_DIR_PZ=54, CONST_VAR=55

  !Constants for initial conditions
  INTEGER,PARAMETER :: CONST_AUTOEARLY=20, CONST_AUTOLATE=21, CONST_EXTERNAL=22
  INTEGER,PARAMETER :: CONST_MANUAL=23, CONST_RESTART=24

  !Custom constants
  INTEGER,PARAMETER :: CONSTANT_DECK_LOWBOUND=4096
  INTEGER,PARAMETER :: CONSTANT_CUSTOM_LOWBOUND=8192

  INTEGER,PARAMETER :: FUNC_SINE=1, FUNC_COSINE=2, FUNC_TAN=3, FUNC_EXP=4
  INTEGER,PARAMETER :: FUNC_ARCSINE=10, FUNC_ARCCOSINE=11, FUNC_ARCTAN=12
  INTEGER,PARAMETER :: FUNC_NEG=13, FUNC_IF=14, FUNC_FLOOR=15, FUNC_CEIL=16
  INTEGER,PARAMETER :: FUNC_NINT=17, FUNC_RHO=18, FUNC_TEMPX=19, FUNC_TEMPY=20
  INTEGER,PARAMETER :: FUNC_TEMPZ=21, FUNC_INTERPOLATE=22, FUNC_TANH=23
  INTEGER,PARAMETER :: FUNC_SINH=24, FUNC_COSH=25, FUNC_EX=26
  INTEGER,PARAMETER :: FUNC_EY=27, FUNC_EZ=28, FUNC_BX=29
  INTEGER,PARAMETER :: FUNC_BY=30, FUNC_BZ=31, FUNC_SQRT=32
  INTEGER,PARAMETER :: FUNC_GAUSS=33, FUNC_SEMIGAUSS=34, FUNC_CRIT=35
  INTEGER,PARAMETER :: FUNC_ABS=36

  INTEGER,PARAMETER :: FUNC_CUSTOM_LOWBOUND=4096

  !Associativity constants
  INTEGER,PARAMETER :: ASSOC_A=1, ASSOC_LA=2, ASSOC_RA=3

  INTEGER,PARAMETER :: Num_Ops=12
  INTEGER, DIMENSION(Num_Ops), PARAMETER :: OPCODE_PRECEDENCE = (/1,1,2,2,3,4,1,1,1,2,2,5/)
  INTEGER, DIMENSION(Num_Ops), PARAMETER :: OPCODE_ASSOC = (/ASSOC_A,ASSOC_LA,ASSOC_A,ASSOC_LA,ASSOC_LA,&
       ASSOC_A,ASSOC_A,ASSOC_A,ASSOC_A,ASSOC_A,ASSOC_A,ASSOC_RA/)

  INTEGER,PARAMETER :: StackSize=10000

  TYPE :: StackElement
     INTEGER :: type
     INTEGER :: Data
     REAL(num) :: NumericalData
#ifdef PARSER_DEBUG
	  CHARACTER(LEN=ENTRYLENGTH) :: Text
#endif
  END TYPE StackElement

  TYPE :: PrimitiveStack
     TYPE(StackElement),DIMENSION(StackSize) :: Data
     INTEGER :: StackPoint
  END TYPE PrimitiveStack

  TYPE :: Deck_Constant
     REAL(num) :: Value
     CHARACTER(LEN=ENTRYLENGTH) :: Name
  END TYPE Deck_Constant

  TYPE :: Deferred_Execution_Object
     CHARACTER(LEN=ENTRYLENGTH) :: Name
     TYPE(PrimitiveStack) :: Execution_Stream
  END TYPE Deferred_Execution_Object

  INTEGER :: n_Deck_Constants=0
  INTEGER :: n_Deferred_Execution_Objects=0
  TYPE(Deck_Constant),DIMENSION(:),ALLOCATABLE :: Deck_Constant_List
  TYPE(Deferred_Execution_Object),DIMENSION(:),ALLOCATABLE :: Deferred_Objects

  REAL(num) :: context_variable !A context dependant variable which is used in various places in the code

END MODULE shared_parser_data


MODULE shared_data
  USE constants
  USE shared_parser_data
  IMPLICIT NONE

  !---------------------------------------------------------------------------------------
  !String handling
  !---------------------------------------------------------------------------------------
  INCLUDE 'mpif.h'
  CHARACTER(LEN=EntryLength) :: Blank
  TYPE :: Entry
     CHARACTER(EntryLength) :: Value
  END TYPE Entry
  CHARACTER(LEN=EntryLength) :: Extended_Error_String

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
  TYPE :: Particle
     REAL(num), DIMENSION(3) :: Part_P
     REAL(num),DIMENSION(2) :: Part_pos
#ifdef PER_PARTICLE_WEIGHT
     REAL(num) :: weight
#endif
#ifdef PER_PARTICLE_CHARGEMASS
     REAL(num) :: charge
     REAL(num) :: mass
#endif
     TYPE(Particle),POINTER :: Next, Prev
#ifdef PART_DEBUG
     INTEGER :: Processor
     INTEGER :: Processor_at_t0
#endif
  END TYPE Particle

  !Object representing a collection of particles
  !Used internally by the MPI particle transfer code
  TYPE :: ParticleList
     TYPE(Particle),POINTER :: Head
     TYPE(Particle),POINTER :: Tail
     INTEGER(KIND=8) :: Count
     !Pointer is safe if the particles in it are all unambiguously linked
     LOGICAL :: Safe

     TYPE(ParticleList), POINTER :: Next, Prev
  END TYPE ParticleList

  !Object representing a particle species
  TYPE :: ParticleFamily
     !Core properties
     CHARACTER(EntryLength) :: Name
     TYPE(ParticleFamily),POINTER :: Next,Prev
     INTEGER :: ID
     LOGICAL :: Dump

     REAL(num) :: Charge
     REAL(num) :: Mass
     INTEGER(KIND=8) :: Count
     TYPE(ParticleList) :: AttachedList

#ifdef TRACER_PARTICLES
     LOGICAL :: Tracer
#endif

     !Particle cell division
#ifdef PARTICLE_CELL_DIVISION
     INTEGER(KIND=8) :: GlobalCount
     LOGICAL :: Split
     INTEGER(KIND=8) :: npart_max
#endif
     !Secondary list
#ifdef USE_SECONDARY_LIST
     TYPE(ParticleList),DIMENSION(:,:),POINTER :: SecondaryList
#endif

     !Injection of particles
     INTEGER(KIND=8) :: npart_per_cell
     REAL(num),DIMENSION(:),ALLOCATABLE :: Density
     REAL(num),DIMENSION(:,:),ALLOCATABLE :: Temperature

     !Species_ionisation
#ifdef PART_IONISE
     LOGICAL :: ionise
     INTEGER :: ionise_to_species
     INTEGER :: release_species
     REAL(num) :: critical_field
     REAL(num) :: ionisation_energy
#endif
     !Attached Probes for this species
#ifdef PARTICLE_PROBES
     TYPE(Particle_Probe),POINTER :: AttachedProbes
#endif
  END TYPE ParticleFamily

  !---------------------------------------------------------------------------------------
  !Initial conditions
  !---------------------------------------------------------------------------------------
  !Represents the initial conditions of a species
  TYPE :: Initial_Condition_Block

     REAL(num),DIMENSION(:,:),POINTER :: Rho
     REAL(num),DIMENSION(:,:,:),POINTER :: Temp
	  REAL(num),DIMENSION(3) :: Drift

     REAL(num) :: minrho
     REAL(num) :: maxrho

  END TYPE Initial_Condition_Block

  INTEGER :: Deck_State
  TYPE(Initial_Condition_Block),DIMENSION(:),ALLOCATABLE :: InitialConditions
  INTEGER :: ictype
  TYPE(Entry) :: icfile

  !---------------------------------------------------------------------------------------
  !Extended IO information
  !---------------------------------------------------------------------------------------
  LOGICAL :: use_extended_io
  CHARACTER(LEN=EntryLength) :: extended_io_file
  !Represents a 2 or 3D distribution
  TYPE :: Distribution_Function_Block

     CHARACTER(LEN=ENTRYLENGTH) :: Name

     !The number of dimensions left
     INTEGER :: nDims

     !The dumpmask for the distribution function
     INTEGER :: DumpMask

     !Whether or not to store the range returned from the distribtution function code
     !This allows the code to store auto determined ranges(experimental)
     LOGICAL :: Store_Ranges

     !The variables which define the ranges and resolutions of the distribution function
     INTEGER,DIMENSION(3) :: Directions
     REAL(num),DIMENSION(3,2) :: Ranges
     INTEGER,DIMENSION(3) :: Resolution
     LOGICAL,DIMENSION(:),ALLOCATABLE :: Use_Species
     REAL(num),DIMENSION(5,2) :: Restrictions
     LOGICAL,DIMENSION(5) :: Use_Restrictions

     !Pointer to next distribution function
     TYPE(Distribution_Function_Block),POINTER :: Next

  END TYPE Distribution_Function_Block
  TYPE(Distribution_Function_Block),POINTER :: Dist_Fns

#ifdef PARTICLE_PROBES
  TYPE :: Particle_Probe
     REAL(num), DIMENSION(2) :: vertex_top, vertex_bottom
     REAL(num) :: ek_min, ek_max
     CHARACTER(LEN=ENTRYLENGTH) :: Name

     TYPE(ParticleFamily), POINTER :: probe_species
     TYPE(ParticleList) :: sampled_particles
     TYPE(particle_probe),POINTER :: Next
     INTEGER :: Dump
  END TYPE Particle_Probe
#endif

  !---------------------------------------------------------------------------------------
  !Core code
  !---------------------------------------------------------------------------------------
  INTEGER :: mpireal = MPI_DOUBLE_PRECISION

#ifdef NEWTONIAN
  !In a Newtonian code, can't guarantee that particles won't exceed lightspeed
  !Therefore record the fastest particle speed here
  REAL(num) :: Max_Part_V
#endif
  INTEGER :: nx,ny
  INTEGER :: nx_global, ny_global
  INTEGER(8) :: npart_global
  INTEGER :: nprocx,nprocy
  INTEGER :: nsteps,nspecies=-1
  LOGICAL :: Smooth_Currents
  REAL(num), ALLOCATABLE, DIMENSION(:,:)     :: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz
  REAL(num), ALLOCATABLE, DIMENSION(:,:)     :: wk_array

  TYPE(ParticleFamily),DIMENSION(:),POINTER :: ParticleSpecies
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:)   :: ekbar !Temperature per species
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:)   :: ekbar_sum,sum_sq,ct !Used to calculate temperature

  REAL(num), ALLOCATABLE, DIMENSION(:)   ::  x,y

  INTEGER,PARAMETER :: data_dir_max_length=64
  CHARACTER(LEN=data_dir_max_length) :: data_dir

  LOGICAL :: Neutral_Background= .TRUE.

  REAL(num) :: dt, t_end, time,dt_multiplier, dt_laser, dt_plasma_frequency
  REAL(num) :: dt_snapshots
  REAL(num) :: length_x, dx, x_start, x_end, x_start_local, x_end_local, length_x_local
  REAL(num) :: length_y, dy, y_start, y_end, y_start_local, y_end_local, length_y_local
  REAL(num), DIMENSION(:), ALLOCATABLE :: x_starts, x_ends, y_starts, y_ends

  REAL(num) :: total_ohmic_heating = 0.0_num
  REAL(num) :: weight

  LOGICAL :: save, restart,deckfile
  INTEGER :: xbc_right, xbc_left,ybc_up,ybc_down
  INTEGER :: xbc_right_particle, xbc_left_particle&
       ,ybc_up_particle,ybc_down_particle
  INTEGER :: xbc_right_field, xbc_left_field&
       ,ybc_up_field,ybc_down_field
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
  REAL(num),DIMENSION(3) :: Window_Shift
  LOGICAL :: AnyOpen
  TYPE(ParticleList) :: Ejected_Particles

  !---------------------------------------------------------------------------------------
  ! MPI data
  !---------------------------------------------------------------------------------------
  INTEGER :: rank, left, right, up, down, coordinates(2),neighbour(-1:1,-1:1)
  INTEGER :: errcode, comm, tag, nproc,icycle_max=1000000
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nx_each_rank, ny_each_rank
  INTEGER(8), ALLOCATABLE, DIMENSION(:) :: npart_each_rank
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: start_each_rank,end_each_rank

  !---------------------------------------------------------------------------------------
  !Domain and loadbalancing
  !---------------------------------------------------------------------------------------
  LOGICAL :: DLB 
  REAL(num) :: DLB_Threshold
  INTEGER(KIND=8), PARAMETER :: npart_per_it = 1000000
  REAL(num),DIMENSION(:),ALLOCATABLE :: x_global,y_global
  REAL(num),DIMENSION(:),ALLOCATABLE :: x_offset_global, y_offset_global
  !The location of the processors
  INTEGER,DIMENSION(:),ALLOCATABLE :: cell_x_start, cell_x_end, cell_y_start, cell_y_end
  INTEGER :: Balance_Mode
  LOGICAL :: DEBUG_MODE

  !---------------------------------------------------------------------------------------
  ! file handling
  !---------------------------------------------------------------------------------------
  INTEGER :: subtype_field,subtype_particle_var,subtype_particle,subtype_particle_int
  INTEGER(KIND=MPI_OFFSET_KIND) :: initialdisp 
  INTEGER :: Full_Dump_Every,Restart_Dump_Every
  INTEGER, PARAMETER :: num_vars_to_dump = 29
  INTEGER, DIMENSION(num_vars_to_dump) :: DumpMask 
  INTEGER :: output_file
  LOGICAL :: force_final_to_be_restartable
  LOGICAL :: Use_Offset_Grid
  INTEGER :: n_zeros=4

!---------------------------------------------------------------------------------------
!Time averaged IO
!---------------------------------------------------------------------------------------

TYPE :: averaged_data_block
	INTEGER :: AverageType=0
	REAL(num),DIMENSION(:,:),POINTER :: data
	INTEGER :: average_over_iterations=-1
	REAL(num) :: average_over_real_time
	INTEGER :: number_of_iterations
END TYPE averaged_data_block
TYPE(averaged_data_block), DIMENSION(num_vars_to_dump), SAVE :: averaged_data

  !---------------------------------------------------------------------------------------
  !Laser boundaries
  !---------------------------------------------------------------------------------------
  TYPE :: Laser_Block
     !Boundary to which laser is attached
     INTEGER :: Direction
     !A unique ID number for the laser (not used directly by EPOCH)
     !Only used if hard coding time profiles
     INTEGER :: ID
     REAL(num), DIMENSION(:), POINTER :: Profile
     REAL(num), DIMENSION(:), POINTER :: Phase

     LOGICAL :: UseTimeFunction
     TYPE(PrimitiveStack) :: TimeFunction

     REAL(num) :: amp=0.0_num,freq=1.0_num, k=1.0_num
     REAL(num) :: pol=0.0_num, angle=0.0_num
     REAL(num) :: t_start=0.0_num, t_end=0.0_num

     TYPE(Laser_Block),POINTER :: Next
  END TYPE Laser_Block
  TYPE(Laser_Block),POINTER :: Laser_Left, Laser_Right
  TYPE(Laser_Block),POINTER :: Laser_Up, Laser_Down
  INTEGER :: n_Laser_Left, n_Laser_Right, n_Laser_up, n_Laser_Down

END MODULE shared_data






