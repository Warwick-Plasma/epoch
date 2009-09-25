MODULE split_particle
  USE shared_data
  USE partlist
  USE helper
  USE boundary

  IMPLICIT NONE

SAVE

INTEGER(KIND=8) :: npart_per_cell_min=5

#ifdef PARTICLE_CELL_DIVISION
CONTAINS

  SUBROUTINE Reorder_Particles_to_grid

    INTEGER :: iSpecies,cellx,celly
    TYPE(Particle),POINTER :: Current, Next
    INTEGER(KIND=8) :: LocalCount

    DO iSpecies=1,nSpecies
       LocalCount=ParticleSpecies(iSpecies)%AttachedList%Count
       CALL MPI_ALLREDUCE(LocalCount,ParticleSpecies(iSpecies)%GlobalCount,1,mpireal,MPI_SUM,comm,errcode)
       ALLOCATE(ParticleSpecies(iSpecies)%SecondaryList(0:nx+1,0:ny+1))
       DO iy=0,ny+1
          DO ix=0,nx+1
             CALL Create_Empty_PartList(ParticleSpecies(iSpecies)%SecondaryList(ix,iy))
          ENDDO
       ENDDO
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO WHILE(ASSOCIATED(Current))
          Next=>Current%Next
          cellx=INT((Current%Part_pos(1)-x_start_local)/dx)+1
          celly=INT((Current%Part_pos(2)-y_start_local)/dy)+1
          CALL Remove_Particle_From_PartList(ParticleSpecies(iSpecies)%AttachedList,Current)
          CALL Add_Particle_To_PartList(ParticleSpecies(iSpecies)%SecondaryList(cellx,celly),Current)
          Current=>Next
       ENDDO
    ENDDO

  END SUBROUTINE Reorder_Particles_to_grid

  SUBROUTINE Reattach_Particles_to_mainlist

    INTEGER :: iSpecies

    DO iSpecies=1,nSpecies
       DO iy=0,ny+1
          DO ix=0,nx+1
             CALL Append_PartList(ParticleSpecies(iSpecies)%AttachedList,ParticleSpecies(iSpecies)%SecondaryList(ix,iy))
          ENDDO
       ENDDO
       DEALLOCATE(ParticleSpecies(iSpecies)%SecondaryList)
    ENDDO

    CALL Particle_bcs

  END SUBROUTINE Reattach_Particles_to_mainlist

#ifdef PER_PARTICLE_WEIGHT
  SUBROUTINE split_particles

    INTEGER :: iSpecies
    INTEGER(KIND=8) :: count
    TYPE(Particle),POINTER :: Current,New
    INTEGER :: clock,idum
    REAL(num) :: jitter_x,jitter_y

    !Reseed random number generator
    CALL SYSTEM_CLOCK(clock)
    idum=-(clock+rank)

    DO iSpecies=1,nSpecies
       IF (.NOT. ParticleSpecies(iSpecies)%Split) CYCLE
       IF (ParticleSpecies(iSpecies)%nPart_Max .GT. 0 .AND. ParticleSpecies(iSpecies)%GlobalCount .GE. ParticleSpecies(iSpecies)%nPart_Max) CYCLE
       DO iy=0,ny+1
          DO ix=0,nx+1
             Count=ParticleSpecies(iSpecies)%SecondaryList(ix,iy)%Count
             IF (Count .GT. 0 .AND. Count .LE. npart_per_cell_min) THEN
                Current=>ParticleSpecies(iSpecies)%SecondaryList(ix,iy)%Head
                DO WHILE(ASSOCIATED(Current) .AND. Count .LE. npart_per_cell_min .AND. Current%Weight .GE. 1.0_num)
                   Count=ParticleSpecies(iSpecies)%SecondaryList(ix,iy)%Count
                   jitter_x=random(idum)*dx/2.0_num - dx/4.0_num
                   jitter_y=random(idum)*dy/2.0_num - dy/4.0_num
                   Current%Weight=Current%Weight/2.0_num
                   ALLOCATE(New)
                   New=Current
                   New%Part_Pos(1)=Current%Part_Pos(1)+jitter_x
                   New%Part_Pos(2)=Current%Part_Pos(2)+jitter_y
                   CALL Add_Particle_To_PartList(ParticleSpecies(iSpecies)%AttachedList,New)
#ifdef PART_DEBUG
                   !If running with particle debugging, specify that this particle has been split
                   New%Processor_At_T0=-1
#endif
                   NULLIFY(New)

                   Current%Part_Pos(1)=Current%Part_Pos(1)-jitter_x
                   Current%Part_Pos(2)=Current%Part_Pos(2)-jitter_y
                   Current=>Current%Next
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE split_particles
#endif

#endif

END MODULE split_particle
