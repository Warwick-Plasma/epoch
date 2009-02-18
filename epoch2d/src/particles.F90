MODULE particles
  USE shared_data
  USE boundary
  USE partlist

  IMPLICIT NONE

CONTAINS

  SUBROUTINE push_particles

    IMPLICIT NONE

    !2nd order accurate particle pusher using parabolic weighting
    !on and off the grid. The calculation of J looks rather odd
    !Since it works by solving d(rho)/dt=div(J) and doing a 1st order
    !Estimate of rho(t+1.5*dt) rather than calculating J directly
    !This gives exact charge conservation on the grid

    !Contains the integer cell position of the particle in x,y,z
    INTEGER :: cell_x1,cell_x2,cell_x3,cell_y1,cell_y2,cell_y3


    !Xi (space factor see page 38 in manual)
    REAL(num),ALLOCATABLE,DIMENSION(:) :: Xi0x, Xi0y
    REAL(num),ALLOCATABLE,DIMENSION(:) :: Xi1x, Xi1y
    !J from a given particle, can be spread over up to 3 cells in 
    !Each direction due to parabolic weighting. We allocate 4 or 5
    !Cells because the position of the particle at t=t+1.5dt is not
    !known until later. This part of the algorithm could probably be
    !Improved, but at the moment, this is just a straight copy of
    !The core of the PSC algorithm
    REAL(num),ALLOCATABLE,DIMENSION(:,:) :: jxh,jyh,jzh

    !J_temp is used in the MPI
    REAL(num),ALLOCATABLE,DIMENSION(:,:) :: J_temp

    !Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x,part_y,part_px,part_py,part_pz,part_q,part_m
    REAL(num) :: root,part_vx,part_vy,part_vz,part_weight
    INTEGER :: part_species

    !Used for particle probes (to see of probe conditions are satisfied)
#ifdef PARTICLE_PROBES
    REAL(num) :: init_part_x, init_part_y, final_part_x, final_part_y
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: probe_x1, probe_y1, probe_x2, probe_y2, probe_a, probe_b
    REAL(num) :: probe_temp,probe_energy
#endif

    !Contains the floating point version of the cell number (never actually used)
    REAL(num) :: cell_x_r,cell_y_r

    !The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x,cell_frac_y

    !Weighting factors as Eqn 4.77 page 25 of manual
    !Eqn 4.77 would be written as
    !F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    !Defined at the particle position
    REAL(num) :: gmx,gmy,g0x,g0y,gpx,gpy


    !Weighting factors as Eqn 4.77 page 25 of manual
    !Eqn 4.77 would be written as
    !F(j-1) * hmx + F(j) * h0x + F(j+1) * hpx
    !Defined at the particle position - 0.5 grid cell in each direction
    !This is to deal with the grid stagger
    REAL(num) :: hmx,hmy,h0x,h0y,hpx,hpy

    !Fields at particle location
    REAL(num) :: Ex_part,Ey_part,Ez_part,Bx_part,By_part,Bz_part,e_part

    !P+ and P- from Page27 of manual
    REAL(num) :: pxp,pxm,pyp,pym,pzp,pzm

    !Charge to mass ratio modified by normalisation
    REAL(num) :: cmratio

    !Tau variables from Page27 of manual
    REAL(num) :: tau,taux,tauy,tauz 

    !Used by J update
    INTEGER :: xmin,xmax,ymin,ymax
    REAL(num) :: wx,wy,wz,cell_y_r0

    !Temporary variables
    REAL(num) :: sum_local,sum_local_sqr,mean,jmx
    INTEGER :: iSpecies,ixt,iyt,next_species

    TYPE(Particle),POINTER :: Current,Next,New_Part

    ALLOCATE(Xi0x(-2:2), Xi0y(-2:2))
    ALLOCATE(Xi1x(-2:2), Xi1y(-2:2))

    ALLOCATE(jxh(-3:2,-2:2))
    ALLOCATE(jyh(-2:2,-3:2))
    ALLOCATE(jzh(-2:2,-2:2))

    Jx=0.0_num
    Jy=0.0_num
    Jz=0.0_num

    jxh=0.0_num
    jyh=0.0_num
    jzh=0.0_num

    EKBAR_SUM=0.0_num
    ct=0.0_num

    part_weight=weight

    DO iSpecies=1,nspecies
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO ipart=1,ParticleSpecies(iSpecies)%AttachedList%Count
          Next=>Current%Next
#ifdef PER_PARTICLE_WEIGHT
          part_weight=Current%Weight
#endif
          !Set the weighting functions to zero for each new particle
          Xi0x=0.0_num
          Xi1x=0.0_num
          Xi0y=0.0_num
          Xi1y=0.0_num

          !Copy the particle properties out for speed
          part_x  = Current%Part_pos(1) - x_start_local
          part_y  = Current%Part_pos(2) - y_start_local
          part_px = Current%Part_P(1)
          part_py = Current%Part_P(2)
          part_pz = Current%Part_P(3)
          part_species=ParticleSpecies(iSpecies)%ID
          !Use a lookup table for charge and mass to save memory
          !No reason not to do this (I think), check properly later
#ifdef PER_PARTICLE_CHARGEMASS
          part_q=current%Charge
          part_m=current%Mass
#else
          part_q  = ParticleSpecies(iSpecies)%Charge
          part_m  = ParticleSpecies(iSpecies)%Mass
#endif

#ifdef PARTICLE_PROBES
          init_part_x = Current%Part_pos(1)
          init_part_y = Current%Part_pos(2)
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

          !Grid weighting factors in 2D (2D analogue of equation 4.77 page 25 of manual)
          !These weight grid properties onto particles


          gmx = 0.5_num * (1.5_num - ABS(cell_frac_x - 1.0_num))**2
          g0x = 0.75_num - ABS(cell_frac_x)**2
          gpx = 0.5_num * (1.5_num - ABS(cell_frac_x + 1.0_num))**2

          gmy = 0.5_num * (1.5_num - ABS(cell_frac_y - 1.0_num))**2
          g0y = 0.75_num - ABS(cell_frac_y)**2
          gpy = 0.5_num * (1.5_num - ABS(cell_frac_y + 1.0_num))**2

          sum_local=SQRT(((part_px*part_weight)**2+(part_py*part_weight)**2+(part_pz*part_weight)**2)*c**2 + (part_m*part_weight)**2*c**4) - (part_m*part_weight)*c**2
          !Calculate the sum of the particle velocities
          ekbar_sum(cell_x1-1,cell_y1-1,part_species)=ekbar_sum(cell_x1-1,cell_y1-1,part_species)+sum_local*gmx*gmy
          ekbar_sum(cell_x1,cell_y1-1,part_species)=ekbar_sum(cell_x1,cell_y1-1,part_species)+sum_local*g0x*gmy
          ekbar_sum(cell_x1+1,cell_y1-1,part_species)=ekbar_sum(cell_x1+1,cell_y1-1,part_species)+sum_local*gpx*gmy

          ekbar_sum(cell_x1-1,cell_y1,part_species)=ekbar_sum(cell_x1-1,cell_y1,part_species)+sum_local*gmx*g0y
          ekbar_sum(cell_x1,cell_y1,part_species)=ekbar_sum(cell_x1,cell_y1,part_species)+sum_local*g0x*g0y
          ekbar_sum(cell_x1+1,cell_y1,part_species)=ekbar_sum(cell_x1+1,cell_y1,part_species)+sum_local*gpx*g0y

          ekbar_sum(cell_x1-1,cell_y1+1,part_species)=ekbar_sum(cell_x1-1,cell_y1+1,part_species)+sum_local*gmx*gpy
          ekbar_sum(cell_x1,cell_y1+1,part_species)=ekbar_sum(cell_x1,cell_y1+1,part_species)+sum_local*g0x*gpy
          ekbar_sum(cell_x1+1,cell_y1+1,part_species)=ekbar_sum(cell_x1+1,cell_y1+1,part_species)+sum_local*gpx*gpy

          !Calculate the particle weights on the grid
          ct(cell_x1-1,cell_y1-1,part_species)=ct(cell_x1-1,cell_y1-1,part_species)+part_weight*gmx*gmy
          ct(cell_x1,cell_y1-1,part_species)=ct(cell_x1,cell_y1-1,part_species)+part_weight*g0x*gmy
          ct(cell_x1+1,cell_y1-1,part_species)=ct(cell_x1+1,cell_y1-1,part_species)+part_weight*gpx*gmy

          ct(cell_x1-1,cell_y1,part_species)=ct(cell_x1-1,cell_y1,part_species)+part_weight*gmx*g0y
          ct(cell_x1,cell_y1,part_species)=ct(cell_x1,cell_y1,part_species)+part_weight*g0x*g0y
          ct(cell_x1+1,cell_y1,part_species)=ct(cell_x1+1,cell_y1,part_species)+part_weight*gpx*g0y

          ct(cell_x1-1,cell_y1+1,part_species)=ct(cell_x1-1,cell_y1+1,part_species)+part_weight*gmx*gpy
          ct(cell_x1,cell_y1+1,part_species)=ct(cell_x1,cell_y1+1,part_species)+part_weight*g0x*gpy
          ct(cell_x1+1,cell_y1+1,part_species)=ct(cell_x1+1,cell_y1+1,part_species )+part_weight*gpx*gpy

          !These are now the weighting factors correct for field weighting
          gmx=0.5_num * (0.5_num + cell_frac_x)**2
          g0x=0.75_num - cell_frac_x**2
          gpx=0.5_num * (0.5_num - cell_frac_x)**2

          gmy=0.5_num * (0.5_num + cell_frac_y)**2
          g0y=0.75_num - cell_frac_y**2
          gpy=0.5_num * (0.5_num - cell_frac_y)**2

          !Particle weighting factors in 2D (2D analogue of 4.140 page 38 of manual)
          !These wieght particle properties onto grid
          !This is used later to calculate J
          Xi0x(-1)=0.5_num * (1.5_num - ABS(cell_frac_x-1.0_num))**2
          Xi0x(+0)=0.75_num - ABS(cell_frac_x)**2
          Xi0x(+1)=0.5_num * (1.5_num - ABS(cell_frac_x + 1.0_num))**2
          Xi0y(-1)=0.5_num * (1.5_num - ABS(cell_frac_y-1.0_num))**2
          Xi0y(+0)=0.75_num - ABS(cell_frac_y)**2
          Xi0y(+1)=0.5_num * (1.5_num - ABS(cell_frac_y + 1.0_num))**2

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

          !Grid weighting factors in 3D (3D analogue of equation 4.77 page 25 of manual)
          !These weight grid properties onto particles
          hmx=0.5_num * (0.5_num + cell_frac_x)**2
          h0x=0.75_num - cell_frac_x**2
          hpx=0.5_num * (0.5_num - cell_frac_x)**2

          hmy=0.5_num * (0.5_num + cell_frac_y)**2
          h0y=0.75_num - cell_frac_y**2
          hpy=0.5_num * (0.5_num - cell_frac_y)**2

          !These are the electric an magnetic fields interpolated to the
          !Particle position. They have been checked and are correct.
          !Actually checking this is messy.
          ex_part=gmy * (hmx*ex(cell_x2-1,cell_y1-1) + h0x*ex(cell_x2,cell_y1-1) + hpx*ex(cell_x2+1,cell_y1-1))&
               +g0y * (hmx*ex(cell_x2-1,cell_y1) + h0x*ex(cell_x2,cell_y1) + hpx*ex(cell_x2+1,cell_y1))&
               +gpy * (hmx*ex(cell_x2-1,cell_y1+1) + h0x*ex(cell_x2,cell_y1+1) + hpx*ex(cell_x2+1,cell_y1+1))
          ey_part=hmy * (gmx*ey(cell_x1-1,cell_y2-1) + g0x*ey(cell_x1,cell_y2-1) + gpx*ey(cell_x1+1,cell_y2-1))&
               +h0y * (gmx*ey(cell_x1-1,cell_y2) + g0x*ey(cell_x1,cell_y2) + gpx*ey(cell_x1+1,cell_y2))&
               +hpy * (gmx*ey(cell_x1-1,cell_y2+1) + g0x*ey(cell_x1,cell_y2+1) + gpx*ey(cell_x1+1,cell_y2+1))
          ez_part=gmy * (gmx*ez(cell_x1-1,cell_y1-1) + g0x*ez(cell_x1,cell_y1-1) + gpx*ez(cell_x1+1,cell_y1-1)) &
               +g0y * (gmx*ez(cell_x1-1,cell_y1) + g0x*ez(cell_x1,cell_y1) + gpx*ez(cell_x1+1,cell_y1)) &
               +gpy * (gmx*ez(cell_x1-1,cell_y1+1) + g0x*ez(cell_x1,cell_y1+1) + gpx*ez(cell_x1+1,cell_y1+1))

          bx_part=hmy * (gmx*bx(cell_x1-1,cell_y2-1) + g0x*Bx(cell_x1,cell_y2-1) + gpx*bx(cell_x1+1,cell_y2-1))&
               +h0y * (gmx*bx(cell_x1-1,cell_y2) + g0x*Bx(cell_x1,cell_y2) + gpx*bx(cell_x1+1,cell_y2))&
               +hpy * (gmx*bx(cell_x1-1,cell_y2+1) + g0x*Bx(cell_x1,cell_y2+1) + gpx*bx(cell_x1+1,cell_y2+1))
          by_part=gmy * (hmx*by(cell_x2-1,cell_y1-1) + h0x*by(cell_x2,cell_y1-1) + hpx*by(cell_x2+1,cell_y1-1)) &
               +g0y * (hmx*by(cell_x2-1,cell_y1) + h0x*by(cell_x2,cell_y1) + hpx*by(cell_x2+1,cell_y1)) &
               +gpy * (hmx*by(cell_x2-1,cell_y1+1) + h0x*by(cell_x2,cell_y1+1) + hpx*by(cell_x2+1,cell_y1+1)) 
          bz_part=hmy * (hmx*bz(cell_x2-1,cell_y2-1) + h0x*bz(cell_x2,cell_y2-1) + hpx*bz(cell_x2+1,cell_y2-1))&
               +h0y * (hmx*bz(cell_x2-1,cell_y2) + h0x*bz(cell_x2,cell_y2) + hpx*bz(cell_x2+1,cell_y2))&
               +hpy *(hmx*bz(cell_x2-1,cell_y2+1) + h0x*bz(cell_x2,cell_y2+1) + hpx*bz(cell_x2+1,cell_y2+1))

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

          !Move particles to end of time step at 2nd order accuracy
          part_x = part_x + part_vx * dt/2.0_num
          part_y = part_y + part_vy * dt/2.0_num

          !Particle has now finished move to end of timestep, so copy back into particle array
          Current%Part_pos(1) = part_x + x_start_local
          Current%Part_pos(2) = part_y + y_start_local
          Current%Part_p  (1) = part_px
          Current%Part_p  (2) = part_py
          Current%Part_p  (3) = part_pz

#ifdef PARTICLE_PROBES
          final_part_x = Current%Part_pos(1)
          final_part_y = Current%Part_pos(2)
#endif

          !If the code is compiled with tracer particle support then put in an
          !If statement so that the current is not calculated for this species
#ifdef TRACER_PARTICLES
          IF (.NOT. ParticleSpecies(iSpecies)%Tracer) THEN
#endif

             !Now advance to t+1.5dt to calculate current. This is detailed in the manual
             !Between pages 37 and 41. The version coded up looks completely different to that
             !In the manual, but is equivalent
             !Use t+1.5 dt so that can update J to t+dt at 2nd order
             part_x = part_x + part_vx * dt/2.0_num
             part_y = part_y + part_vy * dt/2.0_num

             cell_x_r = part_x / dx
             cell_x3  = NINT(cell_x_r)
             cell_frac_x = REAL(cell_x3,num) - cell_x_r
             cell_x3=cell_x3+1

             cell_y_r = part_y / dy
             cell_y3  = NINT(cell_y_r)
             cell_frac_y = REAL(cell_y3,num) - cell_y_r
             cell_y3=cell_y3+1

             IF (cell_y3-cell_y1-1 .LT. -2) PRINT *,"ERROR ",cell_y1,cell_y2,rank

             Xi1x(cell_x3 - cell_x1 - 1) = 0.5_num * (1.5_num - ABS(cell_frac_x - 1.0_num))**2
             Xi1x(cell_x3 - cell_x1 + 0) = 0.75_num - ABS(cell_frac_x)**2
             Xi1x(cell_x3 - cell_x1 + 1) = 0.5_num * (1.5_num - ABS(cell_frac_x + 1.0_num))**2

             Xi1y(cell_y3 - cell_y1 - 1) = 0.5_num * (1.5_num - ABS(cell_frac_y - 1.0_num))**2
             Xi1y(cell_y3 - cell_y1 + 0) = 0.75_num - ABS(cell_frac_y)**2
             Xi1y(cell_y3 - cell_y1 + 1) = 0.5_num * (1.5_num - ABS(cell_frac_y + 1.0_num))**2

             !Now change Xi1* to be Xi1*-Xi0*. This makes the representation of the current update much simpler
             Xi1x = Xi1x - Xi0x
             Xi1y = Xi1y - Xi0y


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

             !Set these to zero due to diffential inside loop
             jxh=0.0_num
             jyh=0.0_num
             jzh=0.0_num

             DO iy=ymin,ymax
                DO ix=xmin,xmax
                   wx = Xi1x(ix) * (Xi0y(iy) + 0.5_num * Xi1y(iy))
                   wy = Xi1y(iy) * (Xi0x(ix) + 0.5_num * Xi1x(ix))
                   wz = Xi0x(ix) * Xi0y(iy) &
                        +0.5_num*xi1x(ix)*xi0y(iy)&
                        +0.5_num*xi0x(ix)*xi1y(iy)&
                        +1.0_num/3.0_num * xi1x(ix) * xi1y(iy)

                   !This is the bit that actually solves d(rho)/dt=-div(J)
                   jxh(ix,iy)=jxh(ix-1,iy) - Part_q * wx * 1.0_num/dt * part_weight/dy 
                   jyh(ix,iy)=jyh(ix,iy-1) - Part_q * wy * 1.0_num/dt * part_weight/dx
                   jzh(ix,iy)=Part_q * Part_vz * wz  * part_weight/(dx*dy)

                   Jx(cell_x1+ix,cell_y1+iy)=Jx(cell_x1+ix,cell_y1+iy)&
                        +jxh(ix,iy)
                   Jy(cell_x1+ix,cell_y1+iy)=Jy(cell_x1+ix,cell_y1+iy)&
                        +jyh(ix,iy)
                   Jz(cell_x1+ix,cell_y1+iy)=Jz(cell_x1+ix,cell_y1+iy)&
                        +jzh(ix,iy)

                ENDDO
             ENDDO
#ifdef TRACER_PARTICLES
          ENDIF
#endif
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

                   ! Probe lines are defined by two points. particles crossing the line are recorded
                   ! assumes transit from left -> right. Put the top at the bottom and bottom at the top
                   ! to catch right-> left moving particles.
                   probe_x1 = current_probe%vertex_bottom(1)
                   probe_y1 = current_probe%vertex_bottom(2)
                   probe_x2 = current_probe%vertex_top(1)
                   probe_y2 = current_probe%vertex_top(2)

                   probe_a = (probe_y1 * probe_x2 - probe_x1 * probe_y2) / (probe_y1 - probe_y2)
                   probe_b = (probe_x1 - probe_x2) / (probe_y1 - probe_y2)

                   ! direction?
                   IF(probe_y2 .LT. probe_y1) THEN
                      probe_temp = init_part_x
                      init_part_x = final_part_x
                      final_part_x = probe_temp
                      probe_temp = init_part_y
                      init_part_y = final_part_y
                      final_part_y = probe_temp
                   ENDIF

                   IF((final_part_y .LT. MAX(probe_y1,probe_y2)) .AND. (final_part_y .GT. min(probe_y1,probe_y2))) THEN
                      IF(init_part_x .LE. probe_a + probe_b * init_part_y) THEN
                         IF(final_part_x .GT. probe_a + probe_b * final_part_y) THEN

                            ! this particle is wanted so copy it to the list associated with this probe
                            ALLOCATE(particle_copy)
                            particle_copy = current
                            CALL add_Particle_To_PartList(current_probe%sampled_particles,particle_copy)
                            NULLIFY(particle_copy)

                         ENDIF
                      ENDIF
                   ENDIF
                ENDIF
             ENDIF
             current_probe => current_probe%next
          ENDDO

#endif
          Current=>Next
       ENDDO
    ENDDO

    !Domain is decomposed. Just add currents at edges
    CALL Processor_Summation_BCS(Jx)
    CALL Field_BC(Jx)
    CALL Processor_Summation_BCS(Jy)
    CALL Field_BC(Jy)
    CALL Processor_Summation_BCS(Jz)
    CALL Field_BC(Jz)


    DO iSpecies=1,nspecies
       CALL Processor_Summation_BCS(ekbar_sum(:,:,iSpecies))
       CALL Field_BC(ekbar_sum(:,:,iSpecies))
       CALL Processor_Summation_BCS(ct(:,:,iSpecies))
       CALL Field_BC(ct(:,:,iSpecies))
    ENDDO


    !Calculate the mean kinetic energy for each species in space
    ekbar=0.0_num
    DO ispecies=1,nspecies
       DO iy=1,ny
          DO ix=1,nx
             mean=ekbar_sum(ix,iy,ispecies)/MAX(ct(ix,iy,ispecies),NONE_ZERO)
             ekbar(ix,iy,ispecies)=mean
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(Xi0x)
    DEALLOCATE(Xi1x)
    DEALLOCATE(Xi0y)
    DEALLOCATE(Xi1y)
    DEALLOCATE(jxh)
    DEALLOCATE(jyh)
    DEALLOCATE(jzh)

    CALL Particle_bcs

    !    JX=0.0_num
    !    JY=0.0_num
    !    Jz=0.0_num

  END SUBROUTINE push_particles

END MODULE particles














