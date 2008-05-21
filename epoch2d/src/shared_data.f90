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
  INTEGER, PARAMETER :: periodic = 1, other = 2
  INTEGER, PARAMETER :: open = 3

  INTEGER, PARAMETER :: Version = 1, Revision = 0

  INTEGER,PARAMETER ::ERR_NONE=0,ERR_UNKNOWN_BLOCK=1,ERR_UNKNOWN_ELEMENT=2
  INTEGER,PARAMETER ::ERR_PRESET_ELEMENT=4,ERR_PRESET_ELEMENT_USE_LATER=8
  INTEGER,PARAMETER ::ERR_BAD_VALUE=16,ERR_MISSING_ELEMENTS=32,ERR_TERMINATE=64
  INTEGER,PARAMETER ::ERR_OTHER=1024

  INTEGER, PARAMETER :: IO_NEVER=0, IO_ALWAYS=1, IO_FULL=2

  REAL(num),PARAMETER :: Q0 = 1.60217646e-19_num !C
  REAL(num),PARAMETER :: M0 = 9.10938188e-31_num !kg
  REAL(num),PARAMETER :: C = 2.98e8_num         !ms^(-2)
  REAL(num),PARAMETER :: kb = 1.3806503e-23_num  !m^2kgs(-2)K^(-1)
  REAL(num),PARAMETER :: epsilon0 = 8.85418782e-12_num


END MODULE constants


MODULE shared_data
  USE constants
  IMPLICIT NONE

  INCLUDE 'mpif.h'
  TYPE :: Entry
     character(30) :: Value
  END TYPE Entry

  INTEGER :: mpireal = MPI_DOUBLE_PRECISION

  INTEGER :: nx,ny
  INTEGER(8) :: npart_global,npart
  INTEGER :: nprocx,nprocy
  INTEGER :: nsteps,nspecies=-1
  REAL(num), ALLOCATABLE, DIMENSION(:,:)     :: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz
  REAL(num), ALLOCATABLE, DIMENSION(:,:)   :: Part_pos, Part_P !(Particle_number,Direction)
  INTEGER,   ALLOCATABLE, DIMENSION(:)     :: Part_species 
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:)   :: Temperature !Temperature per species
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:)   :: sum,sum_sq,ct !Used to calculate temperature
  REAL(num), ALLOCATABLE, DIMENSION(:,:)   :: Species !(Species Number,Unique ID type)
  TYPE(Entry), ALLOCATABLE, DIMENSION(:)   :: Species_Name

  REAL(num), ALLOCATABLE, DIMENSION(:)   ::  x,y

  INTEGER,PARAMETER :: data_dir_max_length=64
  CHARACTER(LEN=data_dir_max_length) :: data_dir

  REAL(num) :: dt, t_end, time,dt_multiplier
  REAL(num) :: dt_snapshots
  REAL(num) :: length_x, dx, x_start, x_end
  REAL(num) :: length_y, dy, y_start, y_end

  REAL(num) :: total_ohmic_heating = 0.0_num
  REAL(num) :: weight

  LOGICAL :: save, restart,deckfile
  INTEGER :: xbc_right, xbc_left,ybc_up,ybc_down
  INTEGER :: restart_snapshot
  INTEGER(8) :: ix,iy,iz,ipart

  ! MPI data

  INTEGER :: rank, left, right, up, down, coordinates(2)
  INTEGER :: errcode, comm, tag, nproc,icycle_max=100000
  INTEGER :: status(MPI_STATUS_SIZE)

  ! file handling
  INTEGER :: subtype,subtype_particle_mesh,subtype_int
  INTEGER(KIND=MPI_OFFSET_KIND) :: initialdisp 
  INTEGER :: Full_Dump_Every
  INTEGER, PARAMETER :: num_vars_to_dump = 24
  INTEGER, DIMENSION(num_vars_to_dump) :: DumpMask 
  INTEGER :: output_file

  REAL(num) :: las_amp,las_omega
  REAL(num) :: t_drop

END MODULE shared_data



