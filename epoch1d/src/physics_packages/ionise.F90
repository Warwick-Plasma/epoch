MODULE ionise

  USE shared_data
  USE partlist
  USE calc_df
  USE helper
  IMPLICIT NONE

CONTAINS

  SUBROUTINE ionise_particles
#ifdef PART_IONISE
    INTEGER :: iSpecies
    TYPE(Particle),POINTER :: Current, Next, New_Part
    REAL(num) :: part_x, part_x2, cell_x_r, cell_frac_x
    REAL(num),DIMENSION(-1:1) :: hx,gx
    REAL(num) :: ex_part,ey_part,ez_part,e_part
    REAL(num) :: number_density_part,ndp_low,ndp_high
    INTEGER :: cell_x1, cell_x2, ix, iy, next_species
    REAL(num),DIMENSION(:),ALLOCATABLE :: Number_Density,ND_Low,ND_High
    REAL(num) :: lambda_db,e_photon,t_eff,saha_rhs,ion_frac,rand
    INTEGER :: idum

    idum=-1445
    rand= random(idum)

    ALLOCATE(Number_Density(-2:nx+3))
    ALLOCATE(ND_Low(-2:nx+3),ND_High(-2:nx+3))
    DO iSpecies=1,nSpecies
       IF (.NOT. ParticleSpecies(iSpecies)%ionise) CYCLE
       CALL calc_number_density(ND_Low,iSpecies)
       CALL calc_number_density(ND_High,ParticleSpecies(iSpecies)%Ionise_to_Species)
       Number_Density=ND_Low+ND_High
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO WHILE(ASSOCIATED(Current))
          Next=>Current%Next
          part_x  = Current%Part_pos - x_start_local

          !Work out number of grid cells in the particle is
          !Not in general an integer
          cell_x_r = part_x/dx
          !Round cell position to nearest cell
          cell_x1=NINT(cell_x_r)
          !Calculate fraction of cell between nearest cell boundary and particle
          cell_frac_x = REAL(cell_x1,num) - cell_x_r
          cell_x1=cell_x1+1

          !These are now the weighting factors correct for field weighting
          gx(-1)=0.5_num * (0.5_num + cell_frac_x)**2
          gx( 0)=0.75_num - cell_frac_x**2
          gx( 1)=0.5_num * (0.5_num - cell_frac_x)**2

          !Now redo shifted by half a cell due to grid stagger.
          !Use shifted version for Ex in X, Ey in Y, Ez in Z
          !And in Y&Z for Bx, X&Z for By, X&Y for Bz
          cell_x_r = part_x/dx - 0.5_num
          cell_x2  = NINT(cell_x_r)
          cell_frac_x = REAL(cell_x2,num) - cell_x_r
          cell_x2=cell_x2+1

          !Grid weighting factors in 3D (3D analogue of equation 4.77 page 25 of manual)
          !These weight grid properties onto particles
          hx(-1)=0.5_num * (0.5_num + cell_frac_x)**2
          hx( 0)=0.75_num - cell_frac_x**2
          hx( 1)=0.5_num * (0.5_num - cell_frac_x)**2

          ex_part=0.0_num
          ey_part=0.0_num
          ez_part=0.0_num
          number_density_part=0.0_num
          ndp_Low=0.0_num
          ndp_High=0.0_num
          DO ix=-1,1
                ex_part=ex_part+hx(ix)*ex(cell_x2+ix)
                ey_part=ey_part+gx(ix)*ex(cell_x1+ix)
                ez_part=ez_part+gx(ix)*ex(cell_x1+ix)
                number_density_part=gx(ix)*&
                     number_density(cell_x1+ix)
                ndp_Low=gx(ix)*ND_Low(cell_x1+ix)
                ndp_High=gx(ix)*ND_High(cell_x1+ix)
          ENDDO
          e_part=SQRT(ex_part**2+ey_part**2+ez_part**2)

          !This is a first attempt at using the 1 level Saha equation to calculate an
          !Ionisation fraction. This isn't really a very good model!

          e_photon=0.5_num * epsilon0 * (e_part)**2 * dx
          t_eff=2.0_num/3.0_num*e_photon/(kb*number_density_part*dx)
          if (t_eff .GT. 1.0e-6_num) THEN
             lambda_db=SQRT(h_planck**2/(2.0_num*pi*m0*kb*t_eff))
             saha_rhs=2.0_num/lambda_db**3 * EXP(-ParticleSpecies(iSpecies)%ionisation_energy/(kb*t_eff))
             ion_frac=0.5_num * (-saha_rhs + SQRT(saha_rhs**2+4.0_num*number_density_part*saha_rhs))
             ion_frac=ion_frac/number_density_part
          ELSE
             ion_frac=0.0_num
          ENDIF
          IF (ion_frac .GT. 1.0_num) ion_frac=1.0_num

          rand=random(idum)
          !After all that, we now know the target ionisation fraction, so subtract the current fraction and ionise
          IF (rand .LT. (ion_frac-ndp_High/MAX(ndp_Low,none_zero))) THEN
!!$
!!$          IF (e_part .GT. ParticleSpecies(iSpecies)%Critical_Field) THEN
             CALL Remove_Particle_From_PartList(ParticleSpecies(iSpecies)%AttachedList,Current)
             next_species=ParticleSpecies(iSpecies)%ionise_to_species
             CALL Add_Particle_To_PartList(ParticleSpecies(next_species)%AttachedList,Current)
             next_species=ParticleSpecies(iSpecies)%Release_Species
             IF (next_species .GT. 0) THEN
                ALLOCATE(New_part)
#ifdef PER_PARTICLE_WEIGHT
                New_Part%weight=Current%weight
#endif
                New_Part%Part_Pos=Current%Part_Pos
                New_Part%Part_P=Current%Part_P
#ifdef PER_PARTICLE_CHARGEMASS
                New_Part%Charge=ParticleSpecies(next_species)%Charge
                New_Part%Mass=ParticleSpecies(next_species)%Mass
#endif
#ifdef PART_DEBUG
                New_Part%Processor=rank
                New_Part%Processor_at_t_0=rank
#endif
                CALL Add_Particle_To_PartList(ParticleSpecies(next_species)%AttachedList,New_Part)
             ENDIF
          ENDIF

          Current=>Next
       ENDDO

    ENDDO

    DEALLOCATE(number_density)
    DEALLOCATE(nd_Low,nd_High)

!!$    PRINT *,m_if
#endif
  END SUBROUTINE ionise_particles


END MODULE ionise
