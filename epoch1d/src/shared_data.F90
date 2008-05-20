!****************************************************************
! All global variables defined here (cf F77 COMMON block).
! All the names in here are public provided the MODULE is USE'd
!**************************************************************** 

MODULE constants
  IMPLICIT NONE
  INTEGER, PARAMETER :: num = KIND(1.D0)
  INTEGER, PARAMETER :: dbl = KIND(1.D0)
  REAL(num), PARAMETER :: pi = 3.14159265358979323_num
  REAL(num), PARAMETER :: none_zero = TINY(1.0_num)
  REAL(num), PARAMETER :: largest_number = HUGE(1.0_num)

  !Boundary type codes
  INTEGER, PARAMETER :: BC_PERIODIC = 1, BC_OTHER = 2, BC_SIMPLE_LASER = 3

  INTEGER, PARAMETER :: Version = 1, Revision = 1

  !Input deck error codes
  INTEGER,PARAMETER ::ERR_NONE=0,ERR_UNKNOWN_BLOCK=1,ERR_UNKNOWN_ELEMENT=2
  INTEGER,PARAMETER ::ERR_PRESET_ELEMENT=4,ERR_PRESET_ELEMENT_USE_LATER=8
  INTEGER,PARAMETER ::ERR_BAD_VALUE=16,ERR_MISSING_ELEMENTS=32,ERR_TERMINATE=64
  INTEGER,PARAMETER ::ERR_OTHER=1024

  !IO codes
  INTEGER, PARAMETER :: IO_NEVER=0, IO_ALWAYS=1, IO_FULL=2, IO_RESTARTABLE=4
  !Domain codes
  INTEGER, PARAMETER :: DO_FULL=0,DO_DECOMPOSED=1

  !Physical constants
  REAL(num),PARAMETER :: Q0 = 1.60217646e-19_num !C Charge of electron
  REAL(num),PARAMETER :: M0 = 9.10938188e-31_num !kg Mass of electron
  REAL(num),PARAMETER :: kb = 1.3806503e-23_num  !m^2kgs(-2)K^(-1)
  REAL(num),PARAMETER :: epsilon0 = 8.85418782e-12_num
  REAL(num),PARAMETER :: c = 2.99792458e8_num

END MODULE constants


MODULE shared_data
  USE constants
  IMPLICIT NONE

  INCLUDE 'mpif.h'
  TYPE :: Entry
     character(30) :: Value
  END TYPE Entry


  TYPE :: ParticlePointer
     TYPE(Particle),POINTER :: Head
     TYPE(Particle),POINTER :: Tail
     INTEGER(KIND=8) :: Count
     !Pointer is safe if the particles in it are all unambiguously linked
     LOGICAL :: Safe
  END TYPE ParticlePointer

  !Object representing a particle
  TYPE :: Particle
     REAL(num), DIMENSION(3) :: Part_P
     REAL(num) :: Part_pos
     INTEGER :: part_Species
#ifdef PER_PARTICLE_WEIGHT
     REAL(num) :: weight
#endif
#ifdef PER_PARTICLE_CHARGEMASS
     REAL(num) :: charge
     REAL(num) :: mass
#endif
     TYPE(Particle),POINTER :: Next, Prev
  END TYPE Particle

  INTEGER :: mpireal = MPI_DOUBLE_PRECISION

  INTEGER :: nx,nx_global
  INTEGER(8) :: npart_global,npart

  INTEGER :: nprocx
  INTEGER :: nsteps,nspecies=-1
  REAL(num), ALLOCATABLE, DIMENSION(:)     :: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,jx_sum,jx_sum2
  TYPE(ParticlePointer),ALLOCATABLE, DIMENSION(:) :: ParticlesByCell
  TYPE(ParticlePointer) :: MainRoot
  REAL(num), ALLOCATABLE, DIMENSION(:,:)   :: Ek_Bar !Temperature per species
  REAL(num), ALLOCATABLE, DIMENSION(:,:)   :: ekbar_sum,ct !Used to calculate temperature
  REAL(num), ALLOCATABLE, DIMENSION(:,:)   :: Species !(Species Number,Unique ID type) (1=charge, 2=mass, 3=number fraction)
  TYPE(Entry), ALLOCATABLE, DIMENSION(:)   :: Species_Name !Name of species n

  REAL(num), ALLOCATABLE, DIMENSION(:)   ::  x,y

  INTEGER,PARAMETER :: data_dir_max_length=64
  INTEGER,PARAMETER :: n_zeros = 4
  CHARACTER(LEN=data_dir_max_length) :: data_dir

  REAL(num) :: dt, t_end, time,dt_multiplier
  REAL(num) :: dt_snapshots
  REAL(num) :: length_x, dx, x_start, x_end, x_start_local, x_end_local
  REAL(num), DIMENSION(:), ALLOCATABLE :: x_starts, x_ends

  REAL(num) :: total_ohmic_heating = 0.0_num
  REAL(num) :: weight

  LOGICAL :: save, restart,deckfile
  INTEGER :: xbc_right, xbc_left,xbc_right_particle, xbc_left_particle
  INTEGER :: restart_snapshot
  INTEGER(8) :: ix,iy,iz,ipart

  !Input deck stuff
  LOGICAL :: code_setup

  ! MPI data

  INTEGER :: rank, left, right, up, down, coordinates(1)
  INTEGER :: errcode, comm, tag, nproc,icycle_max=100000
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nx_each_rank
  INTEGER(8), ALLOCATABLE, DIMENSION(:) :: npart_each_rank
  REAL(num), ALLOCATABLE, DIMENSION(:) :: x_start_each_rank,x_end_each_rank

  !Domain and loadbalancing
  INTEGER :: domain
  LOGICAL :: DLB 
  REAL(num) :: DLB_Threshold
  INTEGER, PARAMETER :: npart_per_it = 10000
  REAL(num),DIMENSION(:),ALLOCATABLE :: x_global

  ! file handling
  INTEGER :: subtype_field,subtype_particle,subtype_particle_var,subtype_particle_int
  INTEGER(KIND=MPI_OFFSET_KIND) :: initialdisp 
  INTEGER :: Full_Dump_Every
  INTEGER, PARAMETER :: num_vars_to_dump = 25
  INTEGER, DIMENSION(num_vars_to_dump) :: DumpMask

  !Termination control
  LOGICAL :: halt = .FALSE.

  !Laser parameters
  REAL(num) :: laser_amp=1.0e8_num,laser_omega=2.0_num * pi /1.0e-8_num

END MODULE shared_data



