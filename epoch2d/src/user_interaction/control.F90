MODULE Control

  USE shared_data
  USE setup

  IMPLICIT NONE

CONTAINS
  !****************************************************
  !Equivalent to "control" block in input deck
  !****************************************************
  SUBROUTINE Setup_Control_Block

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
    icfile%Value="ic.deck"
    restart_snapshot=0
    
    !Neutralising background
    Neutral_Background=.TRUE.

    !Setup the moving window
    move_window=.TRUE.
    window_v_x=2.8e8_num
    window_start_time=7.0e-13_num
    xbc_left_after_move=BC_SIMPLE_OUTFLOW
    xbc_right_after_move=BC_SIMPLE_OUTFLOW


  END SUBROUTINE Setup_Control_Block
  !****************************************************
  !Equivalent to "boundaries" block in input deck
  !****************************************************
  SUBROUTINE Setup_Boundaries_Block

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

  END SUBROUTINE Setup_Boundaries_Block
  !****************************************************
  !Equivalent to "species" block in input deck
  !****************************************************
  SUBROUTINE Setup_Species_Block
    !Must set nSpecies before the call to Setup_Species
    nSpecies=2
    CALL Setup_Species

    !M0 is the mass of an electron
    !Q0 is the charge of an electron
    ParticleSpecies(1)%Name="Electron"
    ParticleSpecies(1)%Mass=1.0_num*M0
    ParticleSpecies(1)%Charge=-1.0_num*Q0
    ParticleSpecies(1)%Count=0.5_num * npart_global
    ParticleSpecies(1)%dump=.TRUE.

    ParticleSpecies(2)%Name="Ions"
    ParticleSpecies(2)%Mass=1836.0_num*M0
    ParticleSpecies(2)%Charge=1.0_num*Q0
    ParticleSpecies(2)%Count=0.5_num * npart_global
    ParticleSpecies(2)%dump=.TRUE.


  END SUBROUTINE Setup_Species_Block
  !****************************************************
  !Equivalent to "output" block in input deck
  !****************************************************
  SUBROUTINE Setup_Output_Block

    dt_snapshots=3.0e-14_num
    full_dump_every=1
    restart_dump_every=-1
    force_final_to_be_restartable=.TRUE.
    use_offset_grid=.TRUE.
    use_extended_io=.FALSE.

    DumpMask=IO_NEVER

    !Dumpmask is a bitmasked variable which determines whether or not to
    !dump a given variable at a given type of output dump. The possible values are
    !IO_NEVER - Default value, never dump
    !IO_ALWAYS - Dump at every output
    !IO_FULL - Dump at only full outputs
    !IO_SPECIES - When ORed with another value, tells the code to output 
    !             some additional per species information

    !The numbers for each variable are
    !1   - The particle positions
    !2   - The cartesian grid
    !3   - Particle px
    !4   - Particle py
    !5   - Particle Pz
    !6   - Particle vx
    !7   - Particle vy
    !8   - Particle vz
    !9   - Ex
    !10  - Ey
    !11  - Ez
    !12  - Bx
    !13  - By
    !14  - Bz
    !15  - Jx
    !16  - Jy
    !17  - Jz
    !18  - Particle Charge
    !19  - Particle Mass
    !20  - Kinetic energy (on grid)
    !21  - Mass density (Can have IO_SPECIES)
    !22  - charge density (Can have IO_SPECIES)
    !23  - number density (Can have IO_SPECIES)
    !24  - Particle weighting value
    !25  - Particle species information
    !26  - Distribution functions
    !27  - Particle Probes
    !28  - Temperature
    !29  - Ejected particles

    DumpMask(2)=IO_ALWAYS
    DumpMask(9:10)=IO_ALWAYS
    DumpMask(14)=IO_ALWAYS
    DumpMask(15:16)=IO_ALWAYS
    DumpMask(20)=IO_ALWAYS
    DumpMask(21)=IOR(IO_ALWAYS,IO_SPECIES)
    DumpMask(22)=IO_ALWAYS
    DumpMask(23)=IOR(IO_ALWAYS,IO_SPECIES)

  END SUBROUTINE  Setup_Output_Block

END MODULE Control
