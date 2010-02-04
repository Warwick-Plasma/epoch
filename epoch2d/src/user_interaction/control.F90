MODULE control

  USE shared_data
  USE setup

  IMPLICIT NONE

CONTAINS
  !****************************************************
  !Equivalent to "control" block in input deck
  !****************************************************
  SUBROUTINE setup_control_block

    !Number of gridpoints in x direction
    nx_global = 256

    !Number of gridpoints in y direction
    ny_global = 50

    !Number of particles of all species(at t=0)
    npart_global=8e5

    !maximum number of iterations
    !set to -1 to run until finished
    nsteps = -1
    !Final runtime of simulation in seconds (NOT walltime)
    t_end = 3.0e-12

    !size of domain
    x_start=-350.0e-6_num
    x_end=150.0e-6_num
    y_start=-50.0e-6_num
    y_end=50.0e-6_num

    !CFL multiplier
    dt_multiplier=0.8_num

    !Dynamic load balancing
    dlb=.FALSE.
    dlb_threshold=1.0_num

    !Initial conditions
    !Bitmask, combine with IOR (or addition)
    !IC_EARLY_INTERNAL
    !IC_EXTERNAL
    !IC_LATE_INTERNAL
    !IC_MANUAL
    !IC_RESTART
    ictype=IC_EARLY_INTERNAL
    icfile%value="ic.deck"
    restart_snapshot=0
    
    !Neutralising background
    neutral_background=.TRUE.

    !Setup the moving window
    move_window=.TRUE.
    window_v_x=2.8e8_num
    window_start_time=7.0e-13_num
    xbc_left_after_move=BC_SIMPLE_OUTFLOW
    xbc_right_after_move=BC_SIMPLE_OUTFLOW


  END SUBROUTINE setup_control_block
  !****************************************************
  !Equivalent to "boundaries" block in input deck
  !****************************************************
  SUBROUTINE setup_boundaries_block

    !Boundary values, can be
    !BC_OTHER - Reflecting boundary, easily changed
    !BC_PERIODIC - Periodic boundary
    !BC_SIMPLE_OUTFLOW - Outflow particle and wave boundary
    !BC_SIMPLE_LASER - Outflow particle and wave boundary with superimposed
    ! inward propagating laser
    xbc_left=BC_SIMPLE_LASER
    xbc_right=BC_SIMPLE_OUTFLOW
    ybc_down=BC_SIMPLE_OUTFLOW
    ybc_up=BC_SIMPLE_OUTFLOW

  END SUBROUTINE setup_boundaries_block
  !****************************************************
  !Equivalent to "species" block in input deck
  !****************************************************
  SUBROUTINE setup_species_block
    !Must set n_species before the call to setup_species
    n_species=2
    CALL setup_species

    !m0 is the mass of an electron
    !q0 is the charge of an electron
    particle_species(1)%name="Electron"
    particle_species(1)%mass=1.0_num*m0
    particle_species(1)%charge=-1.0_num*q0
    particle_species(1)%count=0.5_num * npart_global
    particle_species(1)%dump=.TRUE.

    particle_species(2)%name="Ions"
    particle_species(2)%mass=1836.0_num*m0
    particle_species(2)%charge=1.0_num*q0
    particle_species(2)%count=0.5_num * npart_global
    particle_species(2)%dump=.TRUE.


  END SUBROUTINE setup_species_block
  !****************************************************
  !Equivalent to "output" block in input deck
  !****************************************************
  SUBROUTINE setup_output_block

    dt_snapshots=3.0e-14_num
    full_dump_every=1
    restart_dump_every=-1
    force_final_to_be_restartable=.TRUE.
    use_offset_grid=.TRUE.
    use_extended_io=.FALSE.

    dumpmask=IO_NEVER

    !dumpmask is a bitmasked variable which determines whether or not to
    !dump a given variable at a given type of output dump. The possible values are
    !IO_NEVER - Default value, never dump
    !IO_ALWAYS - dump at every output
    !IO_FULL - dump at only full outputs
    !IO_SPECIES - When ORed with another value, tells the code to output 
    !             some additional per species information

    !The numbers for each variable are
    !1   - The particle positions
    !2   - The cartesian grid
    !3   - particle px
    !4   - particle py
    !5   - particle Pz
    !6   - particle vx
    !7   - particle vy
    !8   - particle vz
    !9   - ex
    !10  - ey
    !11  - ez
    !12  - bx
    !13  - by
    !14  - bz
    !15  - jx
    !16  - jy
    !17  - jz
    !18  - particle charge
    !19  - particle mass
    !20  - Kinetic energy (on grid)
    !21  - mass density (Can have IO_SPECIES)
    !22  - charge density (Can have IO_SPECIES)
    !23  - number density (Can have IO_SPECIES)
    !24  - particle weighting value
    !25  - particle species information
    !26  - Distribution functions
    !27  - particle probes
    !28  - temperature
    !29  - Ejected particles

    dumpmask(2)=IO_ALWAYS
    dumpmask(9:10)=IO_ALWAYS
    dumpmask(14)=IO_ALWAYS
    dumpmask(15:16)=IO_ALWAYS
    dumpmask(20)=IO_ALWAYS
    dumpmask(21)=IOR(IO_ALWAYS,IO_SPECIES)
    dumpmask(22)=IO_ALWAYS
    dumpmask(23)=IOR(IO_ALWAYS,IO_SPECIES)

  END SUBROUTINE  setup_output_block

END MODULE control
