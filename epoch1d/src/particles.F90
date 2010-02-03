MODULE particles
  USE shared_data
  USE boundary
  USE partlist
  USE shape_functions

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
    INTEGER :: cell_x1,cell_x2,cell_x3


    !Xi (space factor see page 38 in manual)
    REAL(num),ALLOCATABLE,DIMENSION(:) :: Xi0x
    REAL(num),ALLOCATABLE,DIMENSION(:) :: Xi1x
    !J from a given particle, can be spread over up to 3 cells in 
    !Each direction due to parabolic weighting. We allocate 4 or 5
    !Cells because the position of the particle at t=t+1.5dt is not
    !known until later. This part of the algorithm could probably be
    !Improved, but at the moment, this is just a straight copy of
    !The core of the PSC algorithm
    REAL(num),ALLOCATABLE,DIMENSION(:) :: jxh,jyh,jzh

    !Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x,part_px,part_py,part_pz,part_q,part_m
    REAL(num) :: root,part_vx,part_vy,part_vz,part_weight
    INTEGER :: part_species

    !Used for particle probes (to see of probe conditions are satisfied)
#ifdef PARTICLE_PROBES
    REAL(num) :: init_part_x, final_part_x
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: probe_temp,probe_energy
#endif

    !Contains the floating point version of the cell number (never actually used)
    REAL(num) :: cell_x_r

    !The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x

    !Weighting factors as Eqn 4.77 page 25 of manual
    !Eqn 4.77 would be written as
    !F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    !Defined at the particle position
    REAL(num),DIMENSION(-2:2) :: gx


    !Weighting factors as Eqn 4.77 page 25 of manual
    !Eqn 4.77 would be written as
    !F(j-1) * hmx + F(j) * h0x + F(j+1) * hpx
    !Defined at the particle position - 0.5 grid cell in each direction
    !This is to deal with the grid stagger
    REAL(num),DIMENSION(-2:2) :: hx

    !Fields at particle location
    REAL(num) :: Ex_part,Ey_part,Ez_part,Bx_part,By_part,Bz_part

    !P+ and P- from Page27 of manual
    REAL(num) :: pxp,pxm,pyp,pym,pzp,pzm

    !Charge to mass ratio modified by normalisation
    REAL(num) :: cmratio

    !Tau variables from Page27 of manual
    REAL(num) :: tau,taux,tauy,tauz 

    !Used by J update
    INTEGER :: xmin,xmax
    REAL(num) :: wx,wy,wz

    !Temporary variables
    REAL(num) :: mean
    INTEGER :: iSpecies

    TYPE(Particle),POINTER :: Current,Next

    ALLOCATE(Xi0x(-3:3))
    ALLOCATE(Xi1x(-3:3))

    ALLOCATE(jxh(-4:3))
    ALLOCATE(jyh(-3:3))
    ALLOCATE(jzh(-3:3))

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

          !Copy the particle properties out for speed
          part_x  = Current%Part_pos - x_start_local
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
          init_part_x = Current%Part_pos
#endif

          !Calculate v(t+0.5dt) from p(t)
          !See PSC manual page (25-27)
          root=1.0_num/SQRT(part_m**2 + (part_px**2 + part_py**2 + part_pz**2)/c**2)
          part_vx = part_px * root
          part_vy = part_py * root
          part_vz = part_pz * root

          !Move particles to half timestep position to first order
          part_x=part_x + part_vx*dt/2.0_num


          !Work out number of grid cells in the particle is
          !Not in general an integer
          cell_x_r = part_x/dx
          !Round cell position to nearest cell
          cell_x1=NINT(cell_x_r)
          !Calculate fraction of cell between nearest cell boundary and particle
          cell_frac_x = REAL(cell_x1,num) - cell_x_r
          cell_x1=cell_x1+1

          !These are now the weighting factors correct for field weighting
			 CALL GridToParticle(cell_frac_x,gx)

          !Particle weighting factors in 1D
          !These wieght particle properties onto grid
          !This is used later to calculate J
			CALL ParticleToGrid(cell_frac_x,Xi0x(-2:2))

          !Now redo shifted by half a cell due to grid stagger.
          !Use shifted version for Ex in X, Ey in Y, Ez in Z
          !And in Y&Z for Bx, X&Z for By, X&Y for Bz
          cell_x_r = part_x/dx - 0.5_num
          cell_x2  = NINT(cell_x_r)
          cell_frac_x = REAL(cell_x2,num) - cell_x_r
          cell_x2=cell_x2+1

          !Grid weighting factors in 3D (3D analogue of equation 4.77 page 25 of manual)
          !These weight grid properties onto particles
			 CALL GridToParticle(cell_frac_x,hx)

			 ex_part=0.0_num
			 ey_part=0.0_num
			 ez_part=0.0_num
			 bx_part=0.0_num
			 by_part=0.0_num
			 bz_part=0.0_num

			 DO ix=-sf_order,sf_order
				 ex_part=ex_part + hx(ix)*ex(cell_x2+ix)
				 ey_part=ey_part + gx(ix)*ey(cell_x1+ix)
				 ez_part=ez_part + gx(ix)*ez(cell_x1+ix)
		
				 bx_part=bx_part + gx(ix)*bx(cell_x1+ix)
				 by_part=by_part + hx(ix)*by(cell_x2+ix)
				 bz_part=bz_part + hx(ix)*bz(cell_x2+ix)
			 ENDDO

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

          !Particle has now finished move to end of timestep, so copy back into particle array
          Current%Part_pos = part_x + x_start_local
          Current%Part_p  (1) = part_px
          Current%Part_p  (2) = part_py
          Current%Part_p  (3) = part_pz

#ifdef PARTICLE_PROBES
          final_part_x = Current%Part_pos
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

             cell_x_r = part_x / dx
             cell_x3  = NINT(cell_x_r)
             cell_frac_x = REAL(cell_x3,num) - cell_x_r
             cell_x3=cell_x3+1

 				 CALL ParticleToGrid(cell_frac_x,Xi1x(cell_x3-cell_x1-2:cell_x3-cell_x1+2))

             !Now change Xi1* to be Xi1*-Xi0*. This makes the representation of the current update much simpler
             Xi1x = Xi1x - Xi0x

             !Remember that due to CFL condition particle can never cross more than one gridcell
             !In one timestep

             IF (cell_x3 == cell_x1) THEN !Particle is still in same cell at t+1.5dt as at t+0.5dt
                xmin = -sf_order
                xmax = +sf_order
             ELSE IF (cell_x3 == cell_x1 - 1) THEN !Particle has moved one cell to left
                xmin = -sf_order-1
                xmax = +sf_order
             ELSE IF (cell_x3 == cell_x1 + 1) THEN !Particle has moved one cell to right
                xmin=-sf_order
                xmax=+sf_order+1
             ENDIF

             !Set these to zero due to diffential inside loop
             jxh=0.0_num
             jyh=0.0_num
             jzh=0.0_num

             DO ix=xmin,xmax
                wx = Xi1x(ix)
                wy = Xi0x(ix) + 0.5_num * Xi1x(ix)
                wz = Xi0x(ix) + 0.5_num * Xi1x(ix)

                !This is the bit that actually solves d(rho)/dt=-div(J)
                jxh(ix)=jxh(ix-1) - Part_q * wx * 1.0_num/dt * Part_Weight 
                jyh(ix)=Part_q * Part_vy * wy  * Part_Weight/dx
                jzh(ix)=Part_q * Part_vz * wz  * Part_Weight/dx


                Jx(cell_x1+ix)=Jx(cell_x1+ix)&
                     +jxh(ix)
                Jy(cell_x1+ix)=Jy(cell_x1+ix)&
                     +jyh(ix)
                Jz(cell_x1+ix)=Jz(cell_x1+ix)&
                     +jzh(ix)
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

                   IF (current_probe%LeftToRight) THEN
                      IF (init_part_x .LT. Current_Probe%Probe_Point .AND. final_part_x .GT. Current_Probe%Probe_Point) THEN
                         ! this particle is wanted so copy it to the list associated with this probe
                         ALLOCATE(particle_copy)
                         particle_copy = current
                         CALL add_Particle_To_PartList(current_probe%sampled_particles,particle_copy)
                         NULLIFY(particle_copy)
                      ENDIF
                   ELSE
                      IF (init_part_x .GT. Current_Probe%Probe_Point .AND. final_part_x .LT. Current_Probe%Probe_Point) THEN
                         ! this particle is wanted so copy it to the list associated with this probe
                         ALLOCATE(particle_copy)
                         particle_copy = current
                         CALL add_Particle_To_PartList(current_probe%sampled_particles,particle_copy)
                         NULLIFY(particle_copy)
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
       CALL Processor_Summation_BCS(ekbar_sum(:,iSpecies))
       CALL Field_BC(ekbar_sum(:,iSpecies))
       CALL Processor_Summation_BCS(ct(:,iSpecies))
       CALL Field_BC(ct(:,iSpecies))
    ENDDO


    !Calculate the mean kinetic energy for each species in space
    ekbar=0.0_num
    DO ispecies=1,nspecies
       DO ix=1,nx
          mean=ekbar_sum(ix,ispecies)/MAX(ct(ix,ispecies),NONE_ZERO)
          ekbar(ix,ispecies)=mean
       ENDDO
    ENDDO

    DEALLOCATE(Xi0x)
    DEALLOCATE(Xi1x)
    DEALLOCATE(jxh)
    DEALLOCATE(jyh)
    DEALLOCATE(jzh)

    CALL Particle_bcs

    !    JX=0.0_num
    !    JY=0.0_num
    !    Jz=0.0_num

  END SUBROUTINE push_particles

END MODULE particles














