<<<<<<< .mine
MODULE Probes
  USE shared_data
  USE partlist
  USE iocontrol
  USE output_particle
  USE output
  USE balance
  USE mpi_subtype_control

  SAVE
  TYPE(ParticleList),POINTER,PRIVATE :: CurrentList

CONTAINS
#ifndef PARTICLE_PROBES
  SUBROUTINE Probe_dummy
  END SUBROUTINE Probe_dummy
#else
  SUBROUTINE Init_Probe(Probe)

    TYPE(Particle_Probe),POINTER :: Probe

    NULLIFY(Probe%Next)
    NULLIFY(Probe%probe_species)
    CALL Create_Empty_PartList(Probe%sampled_particles)

  END SUBROUTINE Init_Probe

  SUBROUTINE Attach_Probe(Probe)

    TYPE(Particle_Probe),POINTER :: Probe
    TYPE(Particle_Probe),POINTER :: Current

    Current=>Probe%Probe_Species%AttachedProbes
    IF (.NOT. ASSOCIATED(Current)) THEN
       Probe%Probe_Species%AttachedProbes=>Probe
       RETURN
    ENDIF
    DO WHILE(ASSOCIATED(Current%Next))
       Current=>Current%Next
    ENDDO
    !Now at the last element in the list
    Current%Next=>Probe

  END SUBROUTINE Attach_Probe

  SUBROUTINE Write_Probes(code)
    INTEGER,INTENT(IN) :: code

    TYPE(Particle_Probe),POINTER :: Current_Probe
    CHARACTER(LEN=ENTRYLENGTH) :: probe_name,temp_name
    INTEGER :: iSpecies,subtype_probe_particle_var
    INTEGER(KIND=8) :: npart_probe_local,npart_probe_global,npart_probe_per_it,npart_probe_per_it_local

    DO ispecies=1,nSpecies
       Current_probe=>ParticleSpecies(iSpecies)%AttachedProbes
       DO WHILE(ASSOCIATED(Current_probe))
          !If don't dump this probe currently then just cycle
          IF (IAND(Current_Probe%Dump,code) .EQ. 0) THEN
             Current_Probe=>Current_Probe%Next
             CYCLE
          ENDIF

          CurrentList=>Current_probe%Sampled_Particles

          npart_probe_per_it_local = npart_per_it
          npart_probe_local = current_probe%sampled_particles%count

          IF(npart_probe_local .GT. 0) npart_probe_per_it_local = MIN(npart_probe_local, npart_probe_per_it_local)

          CALL MPI_ALLREDUCE(npart_probe_local,npart_probe_global,1,MPI_INTEGER8,MPI_SUM,comm,errcode)
          CALL MPI_ALLREDUCE(npart_probe_per_it_local,npart_probe_per_it,1,MPI_INTEGER8,MPI_MIN,comm,errcode)
          IF(npart_probe_global .GT. 0) THEN
             subtype_probe_particle_var=Create_Particle_Subtype(npart_probe_local)

             probe_name =  TRIM(ADJUSTL(current_probe%name))

             !Dump particle Positions
             CALL cfd_Write_nD_Particle_Grid_With_Iterator_All(TRIM(probe_name),&
                  "Probe_Grid",Iterate_probe_Particles,2,npart_probe_local,npart_probe_global,npart_probe_per_it,&
                  PARTICLE_CARTESIAN,subtype_probe_particle_var)

             !Dump Px
             WRITE(Temp_Name,'(a,"_Px")') TRIM(probe_name)
             CALL cfd_Write_nD_Particle_Variable_With_Iterator_All(TRIM(temp_name),&
                  TRIM(probe_name),iterate_probe_px,npart_probe_global,npart_probe_per_it,TRIM(probe_name),&
                  "Probe_Grid",subtype_probe_particle_var)

             !Dump Py
             WRITE(Temp_Name,'(a,"_Py")') TRIM(probe_name)
             CALL cfd_Write_nD_Particle_Variable_With_Iterator_All(TRIM(temp_name),&
                  TRIM(probe_name),iterate_probe_py,npart_probe_global,npart_probe_per_it,TRIM(probe_name),&
                  "Probe_Grid",subtype_probe_particle_var)

             !Dump Pz
             WRITE(Temp_Name,'(a,"_Pz")') TRIM(probe_name)
             CALL cfd_Write_nD_Particle_Variable_With_Iterator_All(TRIM(temp_name),TRIM(probe_name)&
                  ,iterate_probe_pz,npart_probe_global,npart_probe_per_it,TRIM(probe_name),"Probe_Grid",&
                  subtype_probe_particle_var)

             !Dump particle weight function
             WRITE(Temp_Name,'(a,"_weight")') TRIM(probe_name)
#ifdef PER_PARTICLE_WEIGHT
             CALL cfd_Write_nD_Particle_Variable_With_Iterator_All(TRIM(temp_name),TRIM(probe_name),iterate_probe_weight,&
                  npart_probe_global,npart_probe_per_it,TRIM(probe_name),"Probe_Grid",subtype_probe_particle_var)
#else
             CALL cfd_Write_Real_Constant(TRIM(temp_name),TRIM(probe_name),weight,0)
#endif

             CALL Destroy_PartList(current_probe%sampled_particles)
             CALL MPI_TYPE_FREE(subtype_probe_particle_var,errcode)
          ENDIF
          current_probe=>current_probe%next

       ENDDO

       NULLIFY(current_probe)
    ENDDO
  END SUBROUTINE Write_Probes


  !Iterator for particle positions
  SUBROUTINE Iterate_Probe_Particles(data,n_points,direction,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL,INTENT(IN) :: start
    INTEGER,INTENT(IN) :: direction
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount


    IF (start)  THEN
       cur => currentList%head
    ENDIF
    partcount=0

    DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       Data(partcount)=Cur%Part_Pos(direction)-Window_Shift(direction)
       Cur=>Cur%Next
    ENDDO

    n_points=partcount

  END SUBROUTINE iterate_probe_particles

  !Iterator for particle momenta
  SUBROUTINE iterate_probe_px(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    IF (start)  THEN
       cur => currentList%head
    ENDIF
    partcount=0

    DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       data(partcount) = Cur%Part_p(1) 
       Cur=>Cur%Next
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_probe_px

  SUBROUTINE iterate_probe_py(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    IF (start)  THEN
       cur => currentList%head
    ENDIF
    partcount=0

    DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       data(partcount) = Cur%Part_p(2) 
       Cur=>Cur%Next
    ENDDO

    n_points=partcount

  END SUBROUTINE iterate_probe_py

  SUBROUTINE iterate_probe_pz(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    IF (start)  THEN
       cur => currentList%head
    ENDIF
    partcount=0

    DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       data(partcount) = Cur%Part_p(3) 
       Cur=>Cur%Next
    ENDDO

    n_points=partcount

  END SUBROUTINE iterate_probe_pz
#ifdef PER_PARTICLE_WEIGHT
  SUBROUTINE iterate_probe_weight(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    IF (start)  THEN
       cur => currentList%head
    ENDIF
    partcount=0

    DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       data(partcount) = Cur%weight
       Cur=>Cur%Next
    ENDDO

    n_points=partcount

  END SUBROUTINE iterate_probe_weight
#endif
#endif

END MODULE Probes
=======
MODULE Probes
  USE shared_data
  USE partlist
  USE iocontrol
  USE output_particle
  USE output
  USE balance

  SAVE
  TYPE(ParticleList),POINTER,PRIVATE :: CurrentList

CONTAINS
#ifndef PARTICLE_PROBES
  SUBROUTINE Probe_dummy
  END SUBROUTINE Probe_dummy
#else
  SUBROUTINE Init_Probe(Probe)

    TYPE(Particle_Probe),POINTER :: Probe

    NULLIFY(Probe%Next)
    NULLIFY(Probe%probe_species)
    CALL Create_Empty_PartList(Probe%sampled_particles)

  END SUBROUTINE Init_Probe

  SUBROUTINE Attach_Probe(Probe)

    TYPE(Particle_Probe),POINTER :: Probe
    TYPE(Particle_Probe),POINTER :: Current

    Current=>Probe%Probe_Species%AttachedProbes
    IF (.NOT. ASSOCIATED(Current)) THEN
       Probe%Probe_Species%AttachedProbes=>Probe
       RETURN
    ENDIF
    DO WHILE(ASSOCIATED(Current%Next))
       Current=>Current%Next
    ENDDO
    !Now at the last element in the list
    Current%Next=>Probe

  END SUBROUTINE Attach_Probe

  SUBROUTINE Write_Probes(code)
    INTEGER,INTENT(IN) :: code

    TYPE(Particle_Probe),POINTER :: Current_Probe
    CHARACTER(LEN=ENTRYLENGTH) :: probe_name,temp_name
    INTEGER :: iSpecies,subtype_probe_particle_var
    INTEGER(KIND=8) :: npart_probe_local,npart_probe_global,npart_probe_per_it,npart_probe_per_it_local

    DO ispecies=1,nSpecies
       Current_probe=>ParticleSpecies(iSpecies)%AttachedProbes
       DO WHILE(ASSOCIATED(Current_probe))
          !If don't dump this probe currently then just cycle
          IF (IAND(Current_Probe%Dump,code) .EQ. 0) THEN
             Current_Probe=>Current_Probe%Next
             CYCLE
          ENDIF

          CurrentList=>Current_probe%Sampled_Particles

          npart_probe_per_it_local = npart_per_it
          npart_probe_local = current_probe%sampled_particles%count

          IF(npart_probe_local .GT. 0) npart_probe_per_it_local = MIN(npart_probe_local, npart_probe_per_it_local)

          CALL MPI_ALLREDUCE(npart_probe_local,npart_probe_global,1,MPI_INTEGER8,MPI_SUM,comm,errcode)
          CALL MPI_ALLREDUCE(npart_probe_per_it_local,npart_probe_per_it,1,MPI_INTEGER8,MPI_MIN,comm,errcode)
          IF(npart_probe_global .GT. 0) THEN
             subtype_probe_particle_var=Create_Particle_Subtype(npart_probe_local)

             probe_name =  TRIM(ADJUSTL(current_probe%name))

             !Dump particle Positions
             CALL cfd_Write_nD_Particle_Grid_With_Iterator_All(TRIM(probe_name),&
                  "Probe_Grid",Iterate_probe_Particles,2,npart_probe_local,npart_probe_global,npart_probe_per_it,&
                  PARTICLE_CARTESIAN,subtype_probe_particle_var)

             !Dump Px
             WRITE(Temp_Name,'(a,"_Px")') TRIM(probe_name)
             CALL cfd_Write_nD_Particle_Variable_With_Iterator_All(TRIM(temp_name),&
                  TRIM(probe_name),iterate_probe_px,npart_probe_global,npart_probe_per_it,TRIM(probe_name),&
                  "Probe_Grid",subtype_probe_particle_var)

             !Dump Py
             WRITE(Temp_Name,'(a,"_Py")') TRIM(probe_name)
             CALL cfd_Write_nD_Particle_Variable_With_Iterator_All(TRIM(temp_name),&
                  TRIM(probe_name),iterate_probe_py,npart_probe_global,npart_probe_per_it,TRIM(probe_name),&
                  "Probe_Grid",subtype_probe_particle_var)

             !Dump Pz
             WRITE(Temp_Name,'(a,"_Pz")') TRIM(probe_name)
             CALL cfd_Write_nD_Particle_Variable_With_Iterator_All(TRIM(temp_name),TRIM(probe_name)&
                  ,iterate_probe_pz,npart_probe_global,npart_probe_per_it,TRIM(probe_name),"Probe_Grid",&
                  subtype_probe_particle_var)

             !Dump particle weight function
             WRITE(Temp_Name,'(a,"_weight")') TRIM(probe_name)
#ifdef PER_PARTICLE_WEIGHT
             CALL cfd_Write_nD_Particle_Variable_With_Iterator_All(TRIM(temp_name),TRIM(probe_name),iterate_probe_weight,&
                  npart_probe_global,npart_probe_per_it,TRIM(probe_name),"Probe_Grid",subtype_probe_particle_var)
#else
             CALL cfd_Write_Real_Constant(TRIM(temp_name),TRIM(probe_name),weight,0)
#endif

             CALL Destroy_PartList(current_probe%sampled_particles)
             CALL MPI_TYPE_FREE(subtype_probe_particle_var,errcode)
          ENDIF
          current_probe=>current_probe%next

       ENDDO

       NULLIFY(current_probe)
    ENDDO
  END SUBROUTINE Write_Probes


  !Iterator for particle positions
  SUBROUTINE Iterate_Probe_Particles(data,n_points,direction,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL,INTENT(IN) :: start
    INTEGER,INTENT(IN) :: direction
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount


    IF (start)  THEN
       cur => currentList%head
    ENDIF
    partcount=0

    DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       Data(partcount)=Cur%Part_Pos(direction)-Window_Shift(direction)
       Cur=>Cur%Next
    ENDDO

    n_points=partcount

  END SUBROUTINE iterate_probe_particles

  !Iterator for particle momenta
  SUBROUTINE iterate_probe_px(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    IF (start)  THEN
       cur => currentList%head
    ENDIF
    partcount=0

    DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       data(partcount) = Cur%Part_p(1) 
       Cur=>Cur%Next
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_probe_px

  SUBROUTINE iterate_probe_py(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    IF (start)  THEN
       cur => currentList%head
    ENDIF
    partcount=0

    DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       data(partcount) = Cur%Part_p(2) 
       Cur=>Cur%Next
    ENDDO

    n_points=partcount

  END SUBROUTINE iterate_probe_py

  SUBROUTINE iterate_probe_pz(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    IF (start)  THEN
       cur => currentList%head
    ENDIF
    partcount=0

    DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       data(partcount) = Cur%Part_p(3) 
       Cur=>Cur%Next
    ENDDO

    n_points=partcount

  END SUBROUTINE iterate_probe_pz
#ifdef PER_PARTICLE_WEIGHT
  SUBROUTINE iterate_probe_weight(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    IF (start)  THEN
       cur => currentList%head
    ENDIF
    partcount=0

    DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       data(partcount) = Cur%weight
       Cur=>Cur%Next
    ENDDO

    n_points=partcount

  END SUBROUTINE iterate_probe_weight
#endif
#endif

END MODULE Probes
>>>>>>> .r31
