MODULE particles
  USE shared_data
  USE boundary

  IMPLICIT NONE

CONTAINS

  SUBROUTINE push_particles

    !2nd order accurate particle pusher using parabolic weighting
    !on and off the grid. The calculation of J looks rather odd
    !Since it works by solving d(rho)/dt=div(J) and doing a 1st order
    !Estimate of rho(t+1.5*dt) rather than calculating J directly
    !This gives exact charge conservation on the grid

    !Contains the integer cell position of the particle in x,y,z
    INTEGER :: cell_x1,cell_x2,cell_x3,cell_y1,cell_y2,cell_y3
    INTEGER :: cell_z1,cell_z2,cell_z3


    !Xi (space factor see page 38 in manual)
    REAL(num),ALLOCATABLE,DIMENSION(:) :: Xi0x, Xi0y,Xi0z
    REAL(num),ALLOCATABLE,DIMENSION(:) :: Xi1x, Xi1y,Xi1z
    !J from a given particle, can be spread over up to 3 cells in 
    !Each direction due to parabolic weighting. We allocate 4 or 5
    !Cells because the position of the particle at t=t+1.5dt is not
    !known until later. This part of the algorithm could probably be
    !Improved, but at the moment, this is just a straight copy of
    !The core of the PSC algorithm
    REAL(num),ALLOCATABLE,DIMENSION(:,:,:) :: jxh,jyh,jzh

    !temp is used in the MPI
    REAL(num),ALLOCATABLE,DIMENSION(:,:,:) :: temp
    REAL(num),ALLOCATABLE,DIMENSION(:,:,:,:) :: temp2

    !Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x,part_y,part_z,part_px,part_py,part_pz,part_q,part_m
    REAL(num) :: root,part_vx,part_vy,part_vz
    INTEGER :: part_species

    !Used for particle probes (to see of probe conditions are satisfied)
#ifdef PARTICLE_PROBES
    REAL(num) :: init_part_x, init_part_y, init_part_z
    REAL(num) :: final_part_x, final_part_y, final_part_z
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: d_init, d_final
    REAL(num) :: probe_temp,probe_energy
#endif

    !Contains the floating point version of the cell number (never actually used)
    REAL(num) :: cell_x_r,cell_y_r,cell_z_r

    !The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x,cell_frac_y,cell_frac_z

    !Weighting factors as Eqn 4.77 page 25 of manual
    !Eqn 4.77 would be written as
    !F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    !Defined at the particle position
    REAL(num) :: gmx,gmy,gmz,g0x,g0y,g0z,gpx,gpy,gpz


    !Weighting factors as Eqn 4.77 page 25 of manual
    !Eqn 4.77 would be written as
    !F(j-1) * hmx + F(j) * h0x + F(j+1) * hpx
    !Defined at the particle position - 0.5 grid cell in each direction
    !This is to deal with the grid stagger
    REAL(num) :: hmx,hmy,hmz,h0x,h0y,h0z,hpx,hpy,hpz

    !Fields at particle location
    REAL(num) :: Ex_part,Ey_part,Ez_part,Bx_part,By_part,Bz_part

    !P+ and P- from Page27 of manual
    REAL(num) :: pxp,pxm,pyp,pym,pzp,pzm

    !Charge to mass ratio modified by normalisation
    REAL(num) :: cmratio

    !Tau variables from Page27 of manual
    REAL(num) :: tau,taux,tauy,tauz 

    !Used by J update
    INTEGER :: xmin,xmax,ymin,ymax,zmin,zmax
    REAL(num) :: wx,wy,wz,cell_y_r0,third,part_weight

    !Temporary variables
    REAL(num) :: sum_local,sum_local_sqr,mean,jmx
    INTEGER :: ispecies, fail
    INTEGER(KIND=8) :: ipart

    TYPE(particle),POINTER :: Current

    ALLOCATE(Xi0x(-2:2), Xi0y(-2:2), Xi0z(-2:2))
    ALLOCATE(Xi1x(-2:2), Xi1y(-2:2), Xi1z(-2:2))

    ALLOCATE(jxh(-3:2,-2:2,-2:2))
    ALLOCATE(jyh(-2:2,-3:2,-2:2))
    ALLOCATE(jzh(-2:2,-2:2,-3:2))

    Jx=0.0_num
    Jy=0.0_num
    Jz=0.0_num

    jxh=0.0_num
    jyh=0.0_num
    jzh=0.0_num

    ekbar_sum=0.0_num
    ct=0.0_num
    third=1.0_num/3.0_num

        !RETURN

    DO iSpecies=1,nSpecies
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO ipart=1,ParticleSpecies(iSpecies)%AttachedList%Count
          !Set the weighting functions to zero for each new particle
          Xi0x=0.0_num
          Xi1x=0.0_num
          Xi0y=0.0_num
          Xi1y=0.0_num
          Xi0z=0.0_num
          Xi1z=0.0_num

          !Copy the particle properties out for speed
          part_x  = Current%Part_pos(1) - x_start_local
          part_y  = Current%Part_pos(2) - y_start_local
          part_z  = Current%Part_pos(3) - z_start_local
          part_px = Current%Part_P(1)
          part_py = Current%Part_P(2)
          part_pz = Current%Part_P(3)
          part_species = iSpecies
          !Use a lookup table for charge and mass to save memory
          !No reason not to do this (I think), check properly later
#ifdef PER_PARTICLE_CHARGEMASS
          part_q=cur%Charge
          part_m=cur%Mass
#else
          part_q  = ParticleSpecies(iSpecies)%Charge
          part_m  = ParticleSpecies(iSpecies)%Mass
#endif

#ifdef PER_PARTICLE_WEIGHT
          part_weight=Current%Weight
#else
          part_weight=weight
#endif


#ifdef PARTICLE_PROBES
          init_part_x = Current%Part_pos(1)
          init_part_y = Current%Part_pos(2)
          init_part_z = Current%Part_pos(3)
#endif
          !Calculate v(t+0.5dt) from p(t)
          !See PSC manual page (25-27)
          root=1.0_num/SQRT(part_m**2 + (part_px**2 + part_py**2 + part_pz**2)/c**2)
          part_vx = part_px * root
          part_vy = part_py * root
          part_vz = part_pz * root

          !Move particles to half timestep position to first order
          part_x=part_x + part_vx*dt/2.0_num
          part_y=part_y + part_vy*dt/2.0_num
          part_z=part_z + part_vz*dt/2.0_num


          !Work out number of grid cells in the particle is
          !Not in general an integer
          cell_x_r = part_x/dx
          !Round cell position to nearest cell
          cell_x1=NINT(cell_x_r)
          !Calculate fraction of cell between nearest cell boundary and particle
          cell_frac_x = REAL(cell_x1,num) - cell_x_r
          cell_x1=cell_x1+1

          !Work out number of grid cells in the particle is
          !Not in general an integer
          cell_y_r = part_y/dy
          !Round cell position to nearest cell
          cell_y1=NINT(cell_y_r)
          !Calculate fraction of cell between nearest cell boundary and particle
          cell_frac_y = REAL(cell_y1,num) - cell_y_r
          cell_y1=cell_y1+1

          !Work out number of grid cells in the particle is
          !Not in general an integer
          cell_z_r = part_z/dz
          !Round cell position to nearest cell
          cell_z1=NINT(cell_z_r)
          !Calculate fraction of cell between nearest cell boundary and particle
          cell_frac_z = REAL(cell_z1,num) - cell_z_r
          cell_z1=cell_z1+1

          fail=0
          IF (cell_x1 .LT. 0) fail=fail+1
          IF (cell_x1 .GT. nx+1) fail=fail+2
          IF (cell_y1 .LT. 0) fail=fail+4
          IF (cell_y1 .GT. ny+1) fail=fail+8
          IF (cell_z1 .LT. 0) fail=fail+16
          IF (cell_z1 .GT. nz+1) fail=fail+32

          IF (fail .NE. 0) THEN
             PRINT *,part_px,part_py,part_pz,part_m
             !          CALL MPI_ABORT(comm,errcode)
          ENDIF

          !Grid weighting factors in 2D (2D analogue of equation 4.77 page 25 of manual)
          !These weight grid properties onto particles
          gmx=0.5_num * (0.5_num + cell_frac_x)**2
          g0x=0.75_num - cell_frac_x**2
          gpx=0.5_num * (0.5_num - cell_frac_x)**2

          gmy=0.5_num * (0.5_num + cell_frac_y)**2
          g0y=0.75_num - cell_frac_y**2
          gpy=0.5_num * (0.5_num - cell_frac_y)**2

          gmz=0.5_num * (0.5_num + cell_frac_z)**2
          g0z=0.75_num - cell_frac_z**2
          gpz=0.5_num * (0.5_num - cell_frac_z)**2

          sum_local=SQRT(((part_px*part_weight)**2+(part_py*part_weight)**2+(part_pz*part_weight)**2)*c**2 + (part_m*part_weight)**2*c**4) - (part_m*part_weight)*c**2
          !Calculate the sum of the particle velocities
          ekbar_sum(cell_x1-1,cell_y1-1,cell_z1-1,part_species)=ekbar_sum(cell_x1-1,cell_y1-1,cell_z1-1,part_species)+sum_local*gmx*gmy*gmz
          ekbar_sum(cell_x1,cell_y1-1,cell_z1-1,part_species)=ekbar_sum(cell_x1,cell_y1-1,cell_z1-1,part_species)+sum_local*g0x*gmy*gmz
          ekbar_sum(cell_x1+1,cell_y1-1,cell_z1-1,part_species)=ekbar_sum(cell_x1+1,cell_y1-1,cell_z1-1,part_species)+sum_local*gpx*gmy*gmz

          ekbar_sum(cell_x1-1,cell_y1,cell_z1-1,part_species)=ekbar_sum(cell_x1-1,cell_y1-1,cell_z1-1,part_species)+sum_local*gmx*g0y*gmz
          ekbar_sum(cell_x1,cell_y1,cell_z1-1,part_species)=ekbar_sum(cell_x1,cell_y1-1,cell_z1-1,part_species)+sum_local*g0x*g0y*gmz
          ekbar_sum(cell_x1+1,cell_y1,cell_z1-1,part_species)=ekbar_sum(cell_x1+1,cell_y1-1,cell_z1-1,part_species)+sum_local*gpx*g0y*gmz

          ekbar_sum(cell_x1-1,cell_y1+1,cell_z1-1,part_species)=ekbar_sum(cell_x1-1,cell_y1-1,cell_z1-1,part_species)+sum_local*gmx*gpy*gmz
          ekbar_sum(cell_x1,cell_y1+1,cell_z1-1,part_species)=ekbar_sum(cell_x1,cell_y1-1,cell_z1-1,part_species)+sum_local*g0x*gpy*gmz
          ekbar_sum(cell_x1+1,cell_y1+1,cell_z1-1,part_species)=ekbar_sum(cell_x1+1,cell_y1-1,cell_z1-1,part_species)+sum_local*gpx*gpy*gmz

          ekbar_sum(cell_x1-1,cell_y1-1,cell_z1,part_species)=ekbar_sum(cell_x1-1,cell_y1-1,cell_z1,part_species)+sum_local*gmx*gmy*g0z
          ekbar_sum(cell_x1,cell_y1-1,cell_z1,part_species)=ekbar_sum(cell_x1,cell_y1-1,cell_z1,part_species)+sum_local*g0x*gmy*g0z
          ekbar_sum(cell_x1+1,cell_y1-1,cell_z1,part_species)=ekbar_sum(cell_x1+1,cell_y1-1,cell_z1,part_species)+sum_local*gpx*gmy*g0z

          ekbar_sum(cell_x1-1,cell_y1,cell_z1,part_species)=ekbar_sum(cell_x1-1,cell_y1-1,cell_z1,part_species)+sum_local*gmx*g0y*g0z
          ekbar_sum(cell_x1,cell_y1,cell_z1,part_species)=ekbar_sum(cell_x1,cell_y1-1,cell_z1,part_species)+sum_local*g0x*g0y*g0z
          ekbar_sum(cell_x1+1,cell_y1,cell_z1,part_species)=ekbar_sum(cell_x1+1,cell_y1-1,cell_z1,part_species)+sum_local*gpx*g0y*g0z

          ekbar_sum(cell_x1-1,cell_y1+1,cell_z1,part_species)=ekbar_sum(cell_x1-1,cell_y1-1,cell_z1,part_species)+sum_local*gmx*gpy*g0z
          ekbar_sum(cell_x1,cell_y1+1,cell_z1,part_species)=ekbar_sum(cell_x1,cell_y1-1,cell_z1,part_species)+sum_local*g0x*gpy*g0z
          ekbar_sum(cell_x1+1,cell_y1+1,cell_z1,part_species)=ekbar_sum(cell_x1+1,cell_y1-1,cell_z1,part_species)+sum_local*gpx*gpy*g0z

          ekbar_sum(cell_x1-1,cell_y1-1,cell_z1+1,part_species)=ekbar_sum(cell_x1-1,cell_y1-1,cell_z1+1,part_species)+sum_local*gmx*gmy*gpz
          ekbar_sum(cell_x1,cell_y1-1,cell_z1+1,part_species)=ekbar_sum(cell_x1,cell_y1-1,cell_z1+1,part_species)+sum_local*g0x*gmy*gpz
          ekbar_sum(cell_x1+1,cell_y1-1,cell_z1+1,part_species)=ekbar_sum(cell_x1+1,cell_y1-1,cell_z1+1,part_species)+sum_local*gpx*gmy*gpz

          ekbar_sum(cell_x1-1,cell_y1,cell_z1+1,part_species)=ekbar_sum(cell_x1-1,cell_y1-1,cell_z1+1,part_species)+sum_local*gmx*g0y*gpz
          ekbar_sum(cell_x1,cell_y1,cell_z1+1,part_species)=ekbar_sum(cell_x1,cell_y1-1,cell_z1+1,part_species)+sum_local*g0x*g0y*gpz
          ekbar_sum(cell_x1+1,cell_y1,cell_z1+1,part_species)=ekbar_sum(cell_x1+1,cell_y1-1,cell_z1+1,part_species)+sum_local*gpx*g0y*gpz

          ekbar_sum(cell_x1-1,cell_y1+1,cell_z1+1,part_species)=ekbar_sum(cell_x1-1,cell_y1-1,cell_z1+1,part_species)+sum_local*gmx*gpy*gpz
          ekbar_sum(cell_x1,cell_y1+1,cell_z1+1,part_species)=ekbar_sum(cell_x1,cell_y1-1,cell_z1+1,part_species)+sum_local*g0x*gpy*gpz
          ekbar_sum(cell_x1+1,cell_y1+1,cell_z1+1,part_species)=ekbar_sum(cell_x1+1,cell_y1-1,cell_z1+1,part_species)+sum_local*gpx*gpy*gpz

          !Calculate the particle weights on the grid
          ct(cell_x1-1,cell_y1-1,cell_z1-1,part_species)=ct(cell_x1-1,cell_y1-1,cell_z1-1,part_species)+ gmx*gmy*gmz
          ct(cell_x1,cell_y1-1,cell_z1-1,part_species)=ct(cell_x1,cell_y1-1,cell_z1-1,part_species)+ g0x*gmy*gmz
          ct(cell_x1+1,cell_y1-1,cell_z1-1,part_species)=ct(cell_x1+1,cell_y1-1,cell_z1-1,part_species)+gpx*gmy*gmz

          ct(cell_x1-1,cell_y1,cell_z1-1,part_species)=ct(cell_x1-1,cell_y1-1,cell_z1-1,part_species)+gmx*g0y*gmz
          ct(cell_x1,cell_y1,cell_z1-1,part_species)=ct(cell_x1,cell_y1-1,cell_z1-1,part_species)+g0x*g0y*gmz
          ct(cell_x1+1,cell_y1,cell_z1-1,part_species)=ct(cell_x1+1,cell_y1-1,cell_z1-1,part_species)+gpx*g0y*gmz

          ct(cell_x1-1,cell_y1+1,cell_z1-1,part_species)=ct(cell_x1-1,cell_y1-1,cell_z1-1,part_species)+gmx*gpy*gmz
          ct(cell_x1,cell_y1+1,cell_z1-1,part_species)=ct(cell_x1,cell_y1-1,cell_z1-1,part_species)+g0x*gpy*gmz
          ct(cell_x1+1,cell_y1+1,cell_z1-1,part_species)=ct(cell_x1+1,cell_y1-1,cell_z1-1,part_species)+gpx*gpy*gmz

          ct(cell_x1-1,cell_y1-1,cell_z1,part_species)=ct(cell_x1-1,cell_y1-1,cell_z1,part_species)+gmx*gmy*g0z
          ct(cell_x1,cell_y1-1,cell_z1,part_species)=ct(cell_x1,cell_y1-1,cell_z1,part_species)+g0x*gmy*g0z
          ct(cell_x1+1,cell_y1-1,cell_z1,part_species)=ct(cell_x1+1,cell_y1-1,cell_z1,part_species)+gpx*gmy*g0z

          ct(cell_x1-1,cell_y1,cell_z1,part_species)=ct(cell_x1-1,cell_y1-1,cell_z1,part_species)+gmx*g0y*g0z
          ct(cell_x1,cell_y1,cell_z1,part_species)=ct(cell_x1,cell_y1-1,cell_z1,part_species)+g0x*g0y*g0z
          ct(cell_x1+1,cell_y1,cell_z1,part_species)=ct(cell_x1+1,cell_y1-1,cell_z1,part_species)+gpx*g0y*g0z

          ct(cell_x1-1,cell_y1+1,cell_z1,part_species)=ct(cell_x1-1,cell_y1-1,cell_z1,part_species)+gmx*gpy*g0z
          ct(cell_x1,cell_y1+1,cell_z1,part_species)=ct(cell_x1,cell_y1-1,cell_z1,part_species)+g0x*gpy*g0z
          ct(cell_x1+1,cell_y1+1,cell_z1,part_species)=ct(cell_x1+1,cell_y1-1,cell_z1,part_species)+gpx*gpy*g0z

          ct(cell_x1-1,cell_y1-1,cell_z1+1,part_species)=ct(cell_x1-1,cell_y1-1,cell_z1+1,part_species)+gmx*gmy*gpz
          ct(cell_x1,cell_y1-1,cell_z1+1,part_species)=ct(cell_x1,cell_y1-1,cell_z1+1,part_species)+g0x*gmy*gpz
          ct(cell_x1+1,cell_y1-1,cell_z1+1,part_species)=ct(cell_x1+1,cell_y1-1,cell_z1+1,part_species)+gpx*gmy*gpz

          ct(cell_x1-1,cell_y1,cell_z1+1,part_species)=ct(cell_x1-1,cell_y1-1,cell_z1+1,part_species)+gmx*g0y*gpz
          ct(cell_x1,cell_y1,cell_z1+1,part_species)=ct(cell_x1,cell_y1-1,cell_z1+1,part_species)+g0x*g0y*gpz
          ct(cell_x1+1,cell_y1,cell_z1+1,part_species)=ct(cell_x1+1,cell_y1-1,cell_z1+1,part_species)+gpx*g0y*gpz

          ct(cell_x1-1,cell_y1+1,cell_z1+1,part_species)=ct(cell_x1-1,cell_y1-1,cell_z1+1,part_species)+gmx*gpy*gpz
          ct(cell_x1,cell_y1+1,cell_z1+1,part_species)=ct(cell_x1,cell_y1-1,cell_z1+1,part_species)+g0x*gpy*gpz
          ct(cell_x1+1,cell_y1+1,cell_z1+1,part_species)=ct(cell_x1+1,cell_y1-1,cell_z1+1,part_species)+gpx*gpy*gpz

          !Particle weighting factors in 2D (2D analogue of 4.140 page 38 of manual)
          !These wieght particle properties onto grid
          !This is used later to calculate J
          Xi0x(-1)=0.5_num * (1.5_num - ABS(cell_frac_x-1.0_num))**2
          Xi0x(+0)=0.75_num - ABS(cell_frac_x)**2
          Xi0x(+1)=0.5_num * (1.5_num - ABS(cell_frac_x + 1.0_num))**2
          Xi0y(-1)=0.5_num * (1.5_num - ABS(cell_frac_y-1.0_num))**2
          Xi0y(+0)=0.75_num - ABS(cell_frac_y)**2
          Xi0y(+1)=0.5_num * (1.5_num - ABS(cell_frac_y + 1.0_num))**2
          Xi0z(-1)=0.5_num * (1.5_num - ABS(cell_frac_z-1.0_num))**2
          Xi0z(+0)=0.75_num - ABS(cell_frac_z)**2
          Xi0z(+1)=0.5_num * (1.5_num - ABS(cell_frac_z + 1.0_num))**2

          !Now redo shifted by half a cell due to grid stagger.
          !Use shifted version for Ex in X, Ey in Y, Ez in Z
          !And in Y&Z for Bx, X&Z for By, X&Y for Bz
          cell_x_r = part_x/dx - 0.5_num
          cell_x2  = NINT(cell_x_r)
          cell_frac_x = REAL(cell_x2,num) - cell_x_r
          cell_x2=cell_x2+1

          cell_y_r = part_y/dy - 0.5_num
          cell_y2  = NINT(cell_y_r)
          cell_frac_y = REAL(cell_y2,num) - cell_y_r
          cell_y2=cell_y2+1

          cell_z_r = part_z/dz - 0.5_num
          cell_z2  = NINT(cell_z_r)
          cell_frac_z = REAL(cell_z2,num) - cell_z_r
          cell_z2 = cell_z2 + 1

!!$

          !Grid weighting factors in 3D (3D analogue of equation 4.77 page 25 of manual)
          !These weight grid properties onto particles
          hmx=0.5_num * (0.5_num + cell_frac_x)**2
          h0x=0.75_num - cell_frac_x**2
          hpx=0.5_num * (0.5_num - cell_frac_x)**2

          hmy=0.5_num * (0.5_num + cell_frac_y)**2
          h0y=0.75_num - cell_frac_y**2
          hpy=0.5_num * (0.5_num - cell_frac_y)**2

          hmz=0.5_num * (0.5_num + cell_frac_z)**2
          h0z=0.75_num - cell_frac_z**2
          hpz=0.5_num * (0.5_num - cell_frac_z)**2

          ex_part=0.0_num
          ey_part=0.0_num
          ez_part=0.0_num
          bx_part=0.0_num
          by_part=0.0_num
          bz_part=0.0_num

          !These are the electric an magnetic fields interpolated to the
          !Particle position. They have been checked and are correct.
          !Actually checking this is messy.
          ex_part=gmz * (gmy * (hmx*ex(cell_x2-1,cell_y1-1,cell_z1-1) + h0x*ex(cell_x2,cell_y1-1,cell_z1-1) + hpx*ex(cell_x2+1,cell_y1-1,cell_z1-1))&
               +g0y * (hmx*ex(cell_x2-1,cell_y1,cell_z1-1) + h0x*ex(cell_x2,cell_y1,cell_z1-1) + hpx*ex(cell_x2+1,cell_y1,cell_z1-1))&
               +gpy * (hmx*ex(cell_x2-1,cell_y1+1,cell_z1-1) + h0x*ex(cell_x2,cell_y1+1,cell_z1-1) + hpx*ex(cell_x2+1,cell_y1+1,cell_z1-1)))&
               +g0z   * (gmy * (hmx*ex(cell_x2-1,cell_y1-1,cell_z1) + h0x*ex(cell_x2,cell_y1-1,cell_z1) + hpx*ex(cell_x2+1,cell_y1-1,cell_z1))&
               +g0y * (hmx*ex(cell_x2-1,cell_y1,cell_z1) + h0x*ex(cell_x2,cell_y1,cell_z1) + hpx*ex(cell_x2+1,cell_y1,cell_z1))&
               +gpy * (hmx*ex(cell_x2-1,cell_y1+1,cell_z1) + h0x*ex(cell_x2,cell_y1+1,cell_z1) + hpx*ex(cell_x2+1,cell_y1+1,cell_z1)))+&
               gpz    * (gmy * (hmx*ex(cell_x2-1,cell_y1-1,cell_z1+1) + h0x*ex(cell_x2,cell_y1-1,cell_z1+1) + hpx*ex(cell_x2+1,cell_y1-1,cell_z1+1))&
               +g0y * (hmx*ex(cell_x2-1,cell_y1,cell_z1+1) + h0x*ex(cell_x2,cell_y1,cell_z1+1) + hpx*ex(cell_x2+1,cell_y1,cell_z1+1))&
               +gpy * (hmx*ex(cell_x2-1,cell_y1+1,cell_z1+1) + h0x*ex(cell_x2,cell_y1+1,cell_z1+1) + hpx*ex(cell_x2+1,cell_y1+1,cell_z1+1)))

          ey_part=gmz * (hmy * (gmx*ey(cell_x2-1,cell_y1-1,cell_z1-1) + g0x*ey(cell_x2,cell_y1-1,cell_z1-1) + gpx*ey(cell_x2+1,cell_y1-1,cell_z1-1))&
               +h0y * (gmx*ey(cell_x2-1,cell_y1,cell_z1-1) + g0x*ey(cell_x2,cell_y1,cell_z1-1) + gpx*ey(cell_x2+1,cell_y1,cell_z1-1))&
               +hpy * (gmx*ey(cell_x2-1,cell_y1+1,cell_z1-1) + g0x*ey(cell_x2,cell_y1+1,cell_z1-1) + gpx*ey(cell_x2+1,cell_y1+1,cell_z1-1)))&
               +g0z   * (hmy * (gmx*ey(cell_x2-1,cell_y1-1,cell_z1) + g0x*ey(cell_x2,cell_y1-1,cell_z1) + gpx*ey(cell_x2+1,cell_y1-1,cell_z1))&
               +h0y * (gmx*ey(cell_x2-1,cell_y1,cell_z1) + g0x*ey(cell_x2,cell_y1,cell_z1) + gpx*ey(cell_x2+1,cell_y1,cell_z1))&
               +hpy * (gmx*ey(cell_x2-1,cell_y1+1,cell_z1) + g0x*ey(cell_x2,cell_y1+1,cell_z1) + gpx*ey(cell_x2+1,cell_y1+1,cell_z1)))+&
               gpz    * (hmy * (gmx*ey(cell_x2-1,cell_y1-1,cell_z1+1) + g0x*ey(cell_x2,cell_y1-1,cell_z1+1) + gpx*ey(cell_x2+1,cell_y1-1,cell_z1+1))&
               +h0y * (gmx*ey(cell_x2-1,cell_y1,cell_z1+1) + g0x*ey(cell_x2,cell_y1,cell_z1+1) + gpx*ey(cell_x2+1,cell_y1,cell_z1+1))&
               +hpy * (gmx*ey(cell_x2-1,cell_y1+1,cell_z1+1) + g0x*ey(cell_x2,cell_y1+1,cell_z1+1) + gpx*ey(cell_x2+1,cell_y1+1,cell_z1+1)))

          ez_part=hmz * (gmy * (gmx*ez(cell_x2-1,cell_y1-1,cell_z1-1) + g0x*ez(cell_x2,cell_y1-1,cell_z1-1) + gpx*ez(cell_x2+1,cell_y1-1,cell_z1-1))&
               +g0y * (gmx*ez(cell_x2-1,cell_y1,cell_z1-1) + g0x*ez(cell_x2,cell_y1,cell_z1-1) + gpx*ez(cell_x2+1,cell_y1,cell_z1-1))&
               +gpy * (gmx*ez(cell_x2-1,cell_y1+1,cell_z1-1) + g0x*ez(cell_x2,cell_y1+1,cell_z1-1) + gpx*ez(cell_x2+1,cell_y1+1,cell_z1-1)))&
               +h0z   * (gmy * (gmx*ez(cell_x2-1,cell_y1-1,cell_z1) + g0x*ez(cell_x2,cell_y1-1,cell_z1) + gpx*ez(cell_x2+1,cell_y1-1,cell_z1))&
               +g0y * (gmx*ez(cell_x2-1,cell_y1,cell_z1) + g0x*ez(cell_x2,cell_y1,cell_z1) + gpx*ez(cell_x2+1,cell_y1,cell_z1))&
               +gpy * (gmx*ez(cell_x2-1,cell_y1+1,cell_z1) + g0x*ez(cell_x2,cell_y1+1,cell_z1) + gpx*ez(cell_x2+1,cell_y1+1,cell_z1)))+&
               hpz    * (gmy * (gmx*ez(cell_x2-1,cell_y1-1,cell_z1+1) + g0x*ez(cell_x2,cell_y1-1,cell_z1+1) + gpx*ez(cell_x2+1,cell_y1-1,cell_z1+1))&
               +g0y * (gmx*ez(cell_x2-1,cell_y1,cell_z1+1) + g0x*ez(cell_x2,cell_y1,cell_z1+1) + gpx*ez(cell_x2+1,cell_y1,cell_z1+1))&
               +gpy * (gmx*ez(cell_x2-1,cell_y1+1,cell_z1+1) + g0x*ez(cell_x2,cell_y1+1,cell_z1+1) + gpx*ez(cell_x2+1,cell_y1+1,cell_z1+1)))

          bx_part=hmz * (hmy * (gmx*bx(cell_x2-1,cell_y1-1,cell_z1-1) + g0x*bx(cell_x2,cell_y1-1,cell_z1-1) + gpx*bx(cell_x2+1,cell_y1-1,cell_z1-1))&
               +h0y * (gmx*bx(cell_x2-1,cell_y1,cell_z1-1) + g0x*bx(cell_x2,cell_y1,cell_z1-1) + gpx*bx(cell_x2+1,cell_y1,cell_z1-1))&
               +hpy * (gmx*bx(cell_x2-1,cell_y1+1,cell_z1-1) + g0x*bx(cell_x2,cell_y1+1,cell_z1-1) + gpx*bx(cell_x2+1,cell_y1+1,cell_z1-1)))&
               +h0z   * (hmy * (gmx*bx(cell_x2-1,cell_y1-1,cell_z1) + g0x*bx(cell_x2,cell_y1-1,cell_z1) + gpx*bx(cell_x2+1,cell_y1-1,cell_z1))&
               +h0y * (gmx*bx(cell_x2-1,cell_y1,cell_z1) + g0x*bx(cell_x2,cell_y1,cell_z1) + gpx*bx(cell_x2+1,cell_y1,cell_z1))&
               +hpy * (gmx*bx(cell_x2-1,cell_y1+1,cell_z1) + g0x*bx(cell_x2,cell_y1+1,cell_z1) + gpx*bx(cell_x2+1,cell_y1+1,cell_z1)))+&
               hpz    * (hmy * (gmx*bx(cell_x2-1,cell_y1-1,cell_z1+1) + g0x*bx(cell_x2,cell_y1-1,cell_z1+1) + gpx*bx(cell_x2+1,cell_y1-1,cell_z1+1))&
               +h0y * (gmx*bx(cell_x2-1,cell_y1,cell_z1+1) + g0x*bx(cell_x2,cell_y1,cell_z1+1) + gpx*bx(cell_x2+1,cell_y1,cell_z1+1))&
               +hpy * (gmx*bx(cell_x2-1,cell_y1+1,cell_z1+1) + g0x*bx(cell_x2,cell_y1+1,cell_z1+1) + gpx*bx(cell_x2+1,cell_y1+1,cell_z1+1)))

          by_part=hmz * (gmy * (hmx*by(cell_x2-1,cell_y1-1,cell_z1-1) + h0x*by(cell_x2,cell_y1-1,cell_z1-1) + hpx*by(cell_x2+1,cell_y1-1,cell_z1-1))&
               +g0y * (hmx*by(cell_x2-1,cell_y1,cell_z1-1) + h0x*by(cell_x2,cell_y1,cell_z1-1) + hpx*by(cell_x2+1,cell_y1,cell_z1-1))&
               +gpy * (hmx*by(cell_x2-1,cell_y1+1,cell_z1-1) + h0x*by(cell_x2,cell_y1+1,cell_z1-1) + hpx*by(cell_x2+1,cell_y1+1,cell_z1-1)))&
               +h0z   * (gmy * (gmx*by(cell_x2-1,cell_y1-1,cell_z1) + g0x*by(cell_x2,cell_y1-1,cell_z1) + gpx*by(cell_x2+1,cell_y1-1,cell_z1))&
               +g0y * (hmx*by(cell_x2-1,cell_y1,cell_z1) + h0x*by(cell_x2,cell_y1,cell_z1) + hpx*by(cell_x2+1,cell_y1,cell_z1))&
               +gpy * (gmx*by(cell_x2-1,cell_y1+1,cell_z1) + g0x*by(cell_x2,cell_y1+1,cell_z1) + gpx*by(cell_x2+1,cell_y1+1,cell_z1)))+&
               hpz    * (gmy * (hmx*by(cell_x2-1,cell_y1-1,cell_z1+1) + h0x*by(cell_x2,cell_y1-1,cell_z1+1) + hpx*by(cell_x2+1,cell_y1-1,cell_z1+1))&
               +g0y * (hmx*by(cell_x2-1,cell_y1,cell_z1+1) + h0x*by(cell_x2,cell_y1,cell_z1+1) + hpx*by(cell_x2+1,cell_y1,cell_z1+1))&
               +gpy * (hmx*by(cell_x2-1,cell_y1+1,cell_z1+1) + h0x*by(cell_x2,cell_y1+1,cell_z1+1) + hpx*by(cell_x2+1,cell_y1+1,cell_z1+1)))

          bz_part=gmz * (hmy * (hmx*bz(cell_x2-1,cell_y1-1,cell_z1-1) + h0x*bz(cell_x2,cell_y1-1,cell_z1-1) + hpx*bz(cell_x2+1,cell_y1-1,cell_z1-1))&
               +h0y * (hmx*bz(cell_x2-1,cell_y1,cell_z1-1) + h0x*bz(cell_x2,cell_y1,cell_z1-1) + hpx*bz(cell_x2+1,cell_y1,cell_z1-1))&
               +hpy * (hmx*bz(cell_x2-1,cell_y1+1,cell_z1-1) + h0x*bz(cell_x2,cell_y1+1,cell_z1-1) + hpx*bz(cell_x2+1,cell_y1+1,cell_z1-1)))&
               +g0z   * (hmy * (hmx*bz(cell_x2-1,cell_y1-1,cell_z1) + h0x*bz(cell_x2,cell_y1-1,cell_z1) + hpx*bz(cell_x2+1,cell_y1-1,cell_z1))&
               +h0y * (hmx*bz(cell_x2-1,cell_y1,cell_z1) + h0x*bz(cell_x2,cell_y1,cell_z1) + hpx*bz(cell_x2+1,cell_y1,cell_z1))&
               +hpy * (hmx*bz(cell_x2-1,cell_y1+1,cell_z1) + h0x*bz(cell_x2,cell_y1+1,cell_z1) + hpx*bz(cell_x2+1,cell_y1+1,cell_z1)))+&
               gpz    * (hmy * (gmx*bz(cell_x2-1,cell_y1-1,cell_z1+1) + h0x*bz(cell_x2,cell_y1-1,cell_z1+1) + hpx*bz(cell_x2+1,cell_y1-1,cell_z1+1))&
               +h0y * (hmx*bz(cell_x2-1,cell_y1,cell_z1+1) + h0x*bz(cell_x2,cell_y1,cell_z1+1) + hpx*bz(cell_x2+1,cell_y1,cell_z1+1))&
               +hpy * (hmx*bz(cell_x2-1,cell_y1+1,cell_z1+1) + h0x*bz(cell_x2,cell_y1+1,cell_z1+1) + hpx*bz(cell_x2+1,cell_y1+1,cell_z1+1)))

          !update particle momenta using weighted fields
          cmratio = part_q * 0.5_num * dt
          pxm = part_px + cmratio * ex_part
          pym = part_py + cmratio * ey_part
          pzm = part_pz + cmratio * ez_part

          !Half timestep,then use Boris1970 rotation, see Birdsall and Langdon

          root = cmratio / SQRT(part_m**2 + (pxm**2 + pym**2 + pzm**2)/c**2)
          taux = bx_part * root
          tauy = by_part * root
          tauz = bz_part * root

          tau=1.0_num / (1.0_num + taux**2 + tauy**2 + tauz**2)
          pxp=((1.0_num+taux*taux-tauy*tauy-tauz*tauz)*pxm&
               +(2.0_num*taux*tauy+2.0_num*tauz)*pym&
               +(2.0_num*taux*tauz-2.0_num*tauy)*pzm)*tau
          pyp=((2.0_num*taux*tauy-2.0_num*tauz)*pxm&
               +(1.0_num-taux*taux+tauy*tauy-tauz*tauz)*pym&
               +(2.0_num*tauy*tauz+2.0_num*taux)*pzm)*tau
          pzp=((2.0_num*taux*tauz+2.0_num*tauy)*pxm&
               +(2.0_num*tauy*tauz-2.0_num*taux)*pym&
               +(1.0_num-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau

          !Rotation over, go to full timestep
          part_px = pxp + cmratio * ex_part
          part_py = pyp + cmratio * ey_part
          part_pz = pzp + cmratio * ez_part

          !Calculate particle velocity from particle momentum
          root = 1.0_num/SQRT(part_m**2 + (part_px**2 + part_py**2 + part_pz**2)/c**2)
          part_vx=part_px * root
          part_vy=part_py * root
          part_vz=part_pz * root

          !Move particles to end of time step at 2nd order accuracy
          part_x = part_x + part_vx * dt/2.0_num
          part_y = part_y + part_vy * dt/2.0_num
          part_z = part_z + part_vz * dt/2.0_num

          !Particle has now finished move to end of timestep, so copy back into particle array
          Current%Part_pos(1) = part_x + x_start_local
          Current%Part_pos(2) = part_y + y_start_local
          Current%Part_pos(3) = part_z + z_start_local
          Current%Part_p  (1) = part_px
          Current%Part_p  (2) = part_py
          Current%Part_p  (3) = part_pz


#ifdef PARTICLE_PROBES
          final_part_x = Current%Part_pos(1)
          final_part_y = Current%Part_pos(2)
          final_part_z = Current%Part_pos(3)
#endif

          !Original code calculates densities of electrons, ions and neutrals here
          !This has been removed to reduce memory footprint


#ifdef TRACER_PARTICLES
          IF (.NOT. ParticleSpecies(iSpecies)%Tracer) THEN
#endif
             !Now advance to t+1.5dt to calculate current. This is detailed in the manual
             !Between pages 37 and 41. The version coded up looks completely different to that
             !In the manual, but is equivalent
             !Use t+1.5 dt so that can update J to t+dt at 2nd order
             part_x = part_x + part_vx * dt/2.0_num
             part_y = part_y + part_vy * dt/2.0_num
             part_z = part_z + part_vz * dt/2.0_num

             cell_x_r = part_x / dx
             cell_x3  = NINT(cell_x_r)
             cell_frac_x = REAL(cell_x3,num) - cell_x_r
             cell_x3=cell_x3+1

             cell_y_r = part_y / dy
             cell_y3  = NINT(cell_y_r)
             cell_frac_y = REAL(cell_y3,num) - cell_y_r
             cell_y3=cell_y3+1

             cell_z_r = part_z / dz
             cell_z3  = NINT(cell_z_r)
             cell_frac_z = REAL(cell_z3,num) - cell_z_r
             cell_z3 = cell_z3 + 1


!!$

             Xi1x(cell_x3 - cell_x1 - 1) = 0.5_num * (1.5_num - ABS(cell_frac_x - 1.0_num))**2
             Xi1x(cell_x3 - cell_x1 + 0) = 0.75_num - ABS(cell_frac_x)**2
             Xi1x(cell_x3 - cell_x1 + 1) = 0.5_num * (1.5_num - ABS(cell_frac_x + 1.0_num))**2

             Xi1y(cell_y3 - cell_y1 - 1) = 0.5_num * (1.5_num - ABS(cell_frac_y - 1.0_num))**2
             Xi1y(cell_y3 - cell_y1 + 0) = 0.75_num - ABS(cell_frac_y)**2
             Xi1y(cell_y3 - cell_y1 + 1) = 0.5_num * (1.5_num - ABS(cell_frac_y + 1.0_num))**2

             Xi1z(cell_z3 - cell_z1 - 1) = 0.5_num * (1.5_num - ABS(cell_frac_z - 1.0_num))**2
             Xi1z(cell_z3 - cell_z1 + 0) = 0.75_num - ABS(cell_frac_z)**2
             Xi1z(cell_z3 - cell_z1 + 1) = 0.5_num * (1.5_num - ABS(cell_frac_z + 1.0_num))**2

             !Now change Xi1* to be Xi1*-Xi0*. This makes the representation of the current update much simpler
             Xi1x = Xi1x - Xi0x
             Xi1y = Xi1y - Xi0y
             Xi1z = Xi1z - Xi0z


             !Remember that due to CFL condition particle can never cross more than one gridcell
             !In one timestep
             IF (cell_x3 == cell_x1) THEN !Particle is still in same cell at t+1.5dt as at t+0.5dt
                xmin = -1
                xmax = +1
             ELSE IF (cell_x3 == cell_x1 - 1) THEN !Particle has moved one cell to left
                xmin = -2
                xmax = +1
             ELSE IF (cell_x3 == cell_x1 + 1) THEN !Particle has moved one cell to right
                xmin=-1
                xmax=+2
             ENDIF

             IF (cell_y3 == cell_y1) THEN !Particle is still in same cell at t+1.5dt as at t+0.5dt
                ymin = -1
                ymax = +1
             ELSE IF (cell_y3 == cell_y1 - 1) THEN !Particle has moved one cell to left
                ymin = -2
                ymax = +1
             ELSE IF (cell_y3 == cell_y1 + 1) THEN !Particle has moved one cell to right
                ymin=-1
                ymax=+2
             ENDIF

             IF (cell_z3 == cell_z1) THEN !Particle is still in same cell at t+1.5dt as at t+0.5dt
                zmin = -1
                zmax = +1
             ELSE IF (cell_z3 == cell_z1 - 1) THEN !Particle has moved one cell to left
                zmin = -2
                zmax = +1
             ELSE IF (cell_z3 == cell_z1 + 1) THEN !Particle has moved one cell to right
                zmin=-1
                zmax=+2
             ENDIF

             !Set these to zero due to diffential inside loop
             jxh=0.0_num
             jyh=0.0_num
             jzh=0.0_num

             !IF (part_vz .NE. 0) PRINT *,"Moving PARTICLE ERROR"

             !      jmx=MAXVAL(ABS(jz))

             DO iz=zmin,zmax
                DO iy=ymin,ymax
                   DO ix=xmin,xmax
                      wx = Xi1x(ix) * (Xi0y(iy) * Xi0z(iz) +&
                           0.5_num * Xi1y(iy) * Xi0z(iz)+&
                           0.5_num * Xi0y(iy) * Xi1z(iz)+&
                           third * Xi1y(iy) * Xi1z(iz))
                      wy = Xi1y(iy) * (Xi0x(ix) * Xi0z(iz) +&
                           0.5_num * Xi1x(ix) * Xi0z(iz)+&
                           0.5_num * Xi0x(ix) * Xi1z(iz)+&
                           third * Xi1y(iy) * Xi1z(iz))
                      wz = Xi1z(iz) * (Xi0y(iy) * Xi0x(ix) +&
                           0.5_num * Xi1y(iy) * Xi0x(ix)+&
                           0.5_num * Xi0y(iy) * Xi1x(ix)+&
                           third * Xi1y(iy) * Xi1x(ix))

                      !This is the bit that actually solves d(rho)/dt=-div(J)
                      jxh(ix,iy,iz)=jxh(ix-1,iy,iz) - Part_q * wx * 1.0_num/dt * part_weight/(dy*dz)
                      jyh(ix,iy,iz)=jyh(ix,iy-1,iz) - Part_q * wy * 1.0_num/dt * part_weight/(dx*dz)
                      jzh(ix,iy,iz)=jzh(ix,iy,iz-1) - Part_q * wz * 1.0_num/dt * part_weight/(dx*dy)

                      Jx(cell_x1+ix,cell_y1+iy,cell_z1+iz)=Jx(cell_x1+ix,cell_y1+iy,cell_z1+iz)&
                           +jxh(ix,iy,iz)
                      Jy(cell_x1+ix,cell_y1+iy,cell_z1+iz)=Jy(cell_x1+ix,cell_y1+iy,cell_z1+iz)&
                           +jyh(ix,iy,iz)
                      Jz(cell_x1+ix,cell_y1+iy,cell_z1+iz)=Jz(cell_x1+ix,cell_y1+iy,cell_z1+iz)&
                           +jzh(ix,iy,iz)

                   ENDDO
                ENDDO
             ENDDO
#ifdef TRACER_PARTICLES
          ENDIF
#endif
          !       IF (MAXVAL(ABS(jz)) .NE. jmx) PRINT *,"J Change",jmx, wz,Part_vz
#ifdef PARTICLE_PROBES
          ! Compare the current particle with the parameters of any probes in the system. 
          ! These particles are copied into a separate part of the output file.

          Current_probe=>ParticleSpecies(iSpecies)%AttachedProbes

          ! Cycle through probes
          DO WHILE(ASSOCIATED(Current_probe))
             !Note that this is the energy of a single REAL particle in the pseudoparticle, NOT the energy of the pseudoparticle
             probe_energy=(SQRT(1.0_num + (part_px**2 + part_py**2 + part_pz**2)/(part_m * c)**2) - 1.0_num)&
                  * (part_m * c**2)

             ! right energy? (in J)
             IF(probe_energy.GT.current_probe%ek_min)THEN
                IF((probe_energy.LT.current_probe%ek_max).OR.(current_probe%ek_max.LT.0.0_num)) THEN

                   d_init=SUM(Current_Probe%Normal * (Current_Probe%Corner(1,:)-(/init_part_x,init_part_y,init_part_z/)))
                   d_final=SUM(Current_Probe%Normal * (Current_Probe%Corner(1,:)- (/final_part_x,final_part_y,final_part_z/)))
                   IF (SIGN(1.0_num,d_init)*SIGN(1.0_num,d_final) .LE. 0.0_num) THEN
                      ! this particle is wanted so copy it to the list associated with this probe
                      ALLOCATE(particle_copy)
                      particle_copy = current
                      CALL add_Particle_To_PartList(current_probe%sampled_particles,particle_copy)
                      NULLIFY(particle_copy)
                   ENDIF

                ENDIF
             ENDIF
             current_probe => current_probe%next
          ENDDO

#endif
          Current=>Current%Next
!!$
       ENDDO
    ENDDO

    CALL Processor_Summation_BCS(Jx)
    CALL Processor_Summation_BCS(Jy)
    CALL Processor_Summation_BCS(Jz)

    DO iSpecies=1,nSpecies
       CALL Processor_Summation_BCS(ekbar_sum(:,:,:,iSpecies))
       CALL Processor_Summation_BCS(ct(:,:,:,iSpecies))
    ENDDO

    DO ipart=1,nspecies
       DO iz=1,nz
          DO iy=1,ny
             DO ix=1,nx
                IF (ct(ix,iy,iz,ipart) .GT. 0) THEN
                   mean=ekbar_sum(ix,iy,iz,ipart)/ct(ix,iy,iz,ipart)
                ELSE
                   mean=0.0_num
                ENDIF
                ekbar(ix,iy,iz,ipart)=mean
             ENDDO
          ENDDO
       ENDDO
    ENDDO


    DEALLOCATE(Xi0x)
    DEALLOCATE(Xi1x)
    DEALLOCATE(Xi0y)
    DEALLOCATE(Xi1y)
    DEALLOCATE(Xi0z)
    DEALLOCATE(Xi1z)
    DEALLOCATE(jxh)
    DEALLOCATE(jyh)
    DEALLOCATE(jzh)

    !    jx=0.0_num
    !    jy=0.0_num
    !    jz=0.0_num

    CALL Particle_bcs

  END SUBROUTINE push_particles

END MODULE particles














