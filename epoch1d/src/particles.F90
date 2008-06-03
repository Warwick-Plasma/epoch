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

    !J_temp is used in the MPI
    REAL(num),ALLOCATABLE,DIMENSION(:) :: J_temp

    !Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x,part_y,part_px,part_py,part_pz,part_q,part_m
    REAL(num) :: root,part_vx,part_vy,part_vz,part_weight

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
    REAL(num) :: Ex_part,Ey_part,Ez_part,Bx_part,By_part,Bz_part

    !P+ and P- from Page27 of manual
    REAL(num) :: pxp,pxm,pyp,pym,pzp,pzm

    !Charge to mass ratio modified by normalisation
    REAL(num) :: cmratio

    !Tau variables from Page27 of manual
    REAL(num) :: tau,taux,tauy,tauz 

    !Used by J update
    INTEGER :: xmin,xmax,ymin,ymax
    REAL(num) :: wx,wy,wz,cell_y_r0

    !Kinetic energy calculation
    REAL(num) :: ek_particle,temp
    REAL(num),DIMENSION(1:3) :: ek_particle_dir
    REAL(num),DIMENSION(:),ALLOCATABLE :: temp_1d
    REAL(num),DIMENSION(:,:),ALLOCATABLE :: temp_2d

    REAL(num) :: x0,mean
    INTEGER :: iSpecies

    TYPE (PARTICLE),POINTER :: Cur
    ALLOCATE(Xi0x(-2:2))
    ALLOCATE(Xi1x(-2:2))

    ALLOCATE(jxh(-3:2))
    ALLOCATE(jyh(-2:2))
    ALLOCATE(jzh(-2:2))

    Jx=0.0_num
    Jy=0.0_num
    Jz=0.0_num

    jxh=0.0_num
    jyh=0.0_num
    jzh=0.0_num

    ekbar_sum=0.0_num
    ct=0.0_num

    ipart=0
    Cur=>MainRoot%Head

    en_kinetic=0.0_num
    en_kinetic_species=0.0_num
    en_kinetic_total=0.0_num
    en_kinetic_species_total=0.0_num

    DO WHILE (ASSOCIATED(Cur))

       !Set the weighting functions to zero for each new particle
       Xi0x=0.0
       Xi1x=0.0

       !Copy the particle properties out for speed
       part_x  = cur%Part_pos - x_start_local
       part_px = cur%Part_P(1)
       part_py = cur%Part_P(2)
       part_pz = cur%Part_P(3)
       !Use a lookup table for charge and mass to save memory
       !No reason not to do this (I think), check properly later
#ifdef PER_PARTICLE_CHARGEMASS
       part_q=cur%Charge
       part_m=cur%Mass
#else
       part_q  = species(cur%Part_species,1)
       part_m  = species(cur%Part_species,2)
#endif


       !Particle weighting function
#ifdef PER_PARTICLE_WEIGHT
       part_weight=Cur%Weight
#else
       part_weight=weight
#endif

       !This shouldn't happen, but try to handle bad particles gracefully
       IF (part_weight .EQ. 0 .OR. Cur%Part_Species .EQ. 0) THEN
          ipart=ipart+1
          Cur=>Cur%Next
          CYCLE
       ENDIF
       !Calculate v(t+0.5dt) from p(t)
       !See PSC manual page (25-27)
       root=1.0_num/SQRT(part_m**2 + (part_px**2 + part_py**2 + part_pz**2)/c**2)
       part_vx = part_px * root
       part_vy = part_py * root
       part_vz = part_pz * root


       !PRINT *,part_vx

       !Move particles to half timestep position to first order
       part_x=part_x + part_vx*dt/2.0_num


       !Work out number of grid cells in the particle is
       !Not in general an integer
       cell_x_r = part_x/dx
       !Round cell position to nearest cell
       cell_x1=NINT(cell_x_r)
       !Calculate fraction of cell between nearest cell boundary and particle
       cell_frac_x = REAL(cell_x1,num) - cell_x_r
       !Add one since domain runs from 1->nx
       cell_x1=cell_x1+1

       !Grid weighting factors in 3D (3D analogue of equation 4.77 page 25 of manual)
       !These weight grid properties onto particles
       gmx=0.5_num * (0.5_num + cell_frac_x)**2
       g0x=0.75_num - cell_frac_x**2
       gpx=0.5_num * (0.5_num - cell_frac_x)**2

       !      PRINT *,cell_x1

       IF (cell_x1 .GT. nx+2 .OR. cell_x1 .LT. 0) PRINT *,"ERROR in push",rank, Cur%Part_Pos, part_x, part_px,ipart

       ek_particle=SQRT(((part_px*part_weight)**2+(part_py*part_weight)**2+(part_pz*part_weight)**2)*c**2 + (part_m*part_weight)**2*c**4) - (part_m*part_weight)*c**2
       ek_particle_dir(1)=SQRT(((part_px*part_weight)**2)*c**2 + (part_m*part_weight)**2*c**4) - (part_m*part_weight)*c**2
       ek_particle_dir(2)=SQRT(((part_py*part_weight)**2)*c**2 + (part_m*part_weight)**2*c**4) - (part_m*part_weight)*c**2
       ek_particle_dir(3)=SQRT(((part_pz*part_weight)**2)*c**2 + (part_m*part_weight)**2*c**4) - (part_m*part_weight)*c**2
       !ek_particle=part_weight * part_m/2.0_num * (part_vx**2+part_vy**2+part_vz**2)
       en_kinetic=en_kinetic+ek_particle_dir
       en_kinetic_total=en_kinetic_total+ek_particle
       en_kinetic_species(Cur%Part_species,:)=en_kinetic_species(Cur%Part_species,:)+ek_particle_dir
       en_kinetic_species_total(Cur%Part_species)=en_kinetic_species_total(Cur%Part_species)+ek_particle

       ekbar_sum(cell_x1-1,cur%part_species)=ekbar_sum(cell_x1-1,Cur%part_species) + gmx * ek_particle
       ekbar_sum(cell_x1,cur%part_species)  =ekbar_sum(cell_x1,Cur%part_species)   + g0x * ek_particle
       ekbar_sum(cell_x1+1,cur%part_species)=ekbar_sum(cell_x1+1,Cur%part_species) + gpx * ek_particle

       ct(cell_x1-1,cur%part_species)=ct(cell_x1-1,Cur%part_species) + gmx
       ct(cell_x1,cur%part_species)  =ct(cell_x1,Cur%part_species)   + g0x
       ct(cell_x1+1,cur%part_species)=ct(cell_x1+1,Cur%part_species) + gpx

       !Particle weighting factors in 1D (1D analogue of 4.140 page 38 of manual)
       !These wieght particle properties onto grid
       !This is used later to calculate J
       Xi0x(-1)=0.5_num * (1.5_num - ABS(cell_frac_x-1.0_num))**2
       Xi0x(+0)=0.75_num - ABS(cell_frac_x)**2
       Xi0x(+1)=0.5_num * (1.5_num - ABS(cell_frac_x + 1.0_num))**2


       !Now redo shifted by half a cell due to grid stagger.
       !Use shifted version for Ex in X, Ey in Y, Ez in Z
       !And in Y&Z for Bx, X&Z for By, X&Y for Bz
       cell_x_r = part_x/dx - 0.5_num
       cell_x2  = NINT(cell_x_r)
       cell_frac_x = REAL(cell_x2,num) - cell_x_r
       cell_x2  = cell_x2+1

       !Grid weighting factors in 3D (3D analogue of equation 4.77 page 25 of manual)
       !These weight grid properties onto particles
       hmx=0.5_num * (0.5_num + cell_frac_x)**2
       h0x=0.75_num - cell_frac_x**2
       hpx=0.5_num * (0.5_num - cell_frac_x)**2


       !These are the electric an magnetic fields interpolated to the
       !Particle position. They have been checked and are correct.
       !Actually checking this is messy.
       ex_part=hmx*ex(cell_x2-1) + h0x*ex(cell_x2) + hpx*ex(cell_x2+1)
       ey_part=gmx*ey(cell_x1-1) + g0x*ey(cell_x1) + gpx*ey(cell_x1+1)
       ez_part=gmx*ez(cell_x1-1) + g0x*ez(cell_x1) + gpx*ez(cell_x1+1)

       bx_part=gmx*bx(cell_x1-1) + g0x*Bx(cell_x1) + gpx*bx(cell_x1+1)
       by_part=hmx*by(cell_x2-1) + h0x*by(cell_x2) + hpx*by(cell_x2+1)
       bz_part=hmx*bz(cell_x2-1) + h0x*bz(cell_x2) + hpx*bz(cell_x2+1)

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

       !Particle has now finished move to end of timestep, so copy back into particle array
       cur%Part_pos = part_x + x_start_local
       cur%Part_p  (1) = part_px
       cur%Part_p  (2) = part_py
       cur%Part_p  (3) = part_pz

       !Original code calculates densities of electrons, ions and neutrals here
       !This has been removed to reduce memory footprint

       !Now advance to t+1.5dt to calculate current. This is detailed in the manual
       !Between pages 37 and 41. The version coded up looks completely different to that
       !In the manual, but is equivalent
       !Use t+1.5 dt so that can update J to t+dt at 2nd order
       part_x = part_x + part_vx * dt/2.0_num

       cell_x_r = part_x / dx
       cell_x3  = NINT(cell_x_r)
       cell_frac_x = REAL(cell_x3,num) - cell_x_r
       cell_x3 = cell_x3+1

       Xi1x(cell_x3 - cell_x1 - 1) = 0.5_num * (1.5_num - ABS(cell_frac_x - 1.0_num))**2
       Xi1x(cell_x3 - cell_x1 + 0) = 0.75_num - ABS(cell_frac_x)**2
       Xi1x(cell_x3 - cell_x1 + 1) = 0.5_num * (1.5_num - ABS(cell_frac_x + 1.0_num))**2

       !Now change Xi1* to be Xi1*-Xi0*. This makes the representation of the current update much simpler
       Xi1x = Xi1x - Xi0x


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
       Cur=>Cur%Next
       ipart=ipart+1
    ENDDO

    !PRINT *,"ekbar_sum",Ekbar_sum(nx:,1),MAXVAL(ekbar_sum(:,1))

    !Now have J{x,y,z} on each node, sum them up and scatter to each
    IF (domain == DO_FULL) THEN

       !Domain is shared, used MPI reduction operations to add currents from each processor
       CALL Field_Reduction(Jx)
       CALL Periodic_Summation_BCS(Jx)
       CALL Field_Reduction(Jy)
       CALL Periodic_Summation_BCS(Jy)
       CALL Field_Reduction(Jz)
       CALL Periodic_Summation_BCS(Jz)

       DO iSpecies=1,nspecies
          CALL Field_Reduction(ekbar_sum(:,iSpecies))
          CALL Periodic_Summation_BCS(ekbar_sum(:,iSpecies))
          CALL Field_Reduction(ct(:,iSpecies))
          CALL Periodic_Summation_BCS(ct(:,iSpecies))
       ENDDO

    ELSE
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
    ENDIF


    !Calculate the mean kinetic energy for each species in space
    ek_bar=0.0_num
    DO ispecies=1,nspecies
       DO ix=1,nx
          mean=ekbar_sum(ix,ispecies)/MAX(ct(ix,ispecies),NONE_ZERO)
          ek_bar(ix,ispecies)=mean
       ENDDO
    ENDDO

    !Add the total kinetic energy for the diagnostics
    CALL MPI_REDUCE(en_kinetic_total,temp,1,mpireal,MPI_SUM,0,comm,errcode)
    en_kinetic_total=temp

    ALLOCATE(temp_1d(1:3))
    CALL MPI_REDUCE(en_kinetic,temp_1d,3,mpireal,MPI_SUM,0,comm,errcode)
    en_kinetic=temp_1d
    DEALLOCATE(temp_1d)

    ALLOCATE(temp_1d(1:nspecies))
    CALL MPI_REDUCE(en_kinetic_species_total,temp_1d,nspecies,mpireal,MPI_SUM,0,comm,errcode)   
    en_kinetic_species_total=temp_1d
    DEALLOCATE(temp_1d)

    ALLOCATE(temp_2d(1:nspecies,1:3))
    CALL MPI_REDUCE(en_kinetic_species,temp_2d,nspecies*3,mpireal,MPI_SUM,0,comm,errcode)
    en_kinetic_species=temp_2d
    DEALLOCATE(temp_2d)

    DEALLOCATE(Xi0x)
    DEALLOCATE(Xi1x)
    DEALLOCATE(jxh)
    DEALLOCATE(jyh)
    DEALLOCATE(jzh)

    CALL Particle_bcs


  END SUBROUTINE push_particles
END MODULE particles
