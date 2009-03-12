MODULE iterators

  USE Shared_Data
  USE particlepointeradvance
  IMPLICIT NONE

  TYPE :: SettingBlock
     LOGICAL :: Restart
  END TYPE SettingBlock

  SAVE

  TYPE(SettingBlock) :: Iterator_Settings

CONTAINS

  !Iterator for particle positions
  SUBROUTINE iterate_particles(data,n_points,direction,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL,INTENT(IN) :: start
    INTEGER,INTENT(IN) :: direction
    TYPE(Particle),POINTER,SAVE :: Cur
    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily
    INTEGER(8) :: partcount


    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             Data(partcount)=Cur%Part_Pos(direction)-Window_Shift(direction)
             Cur=>Cur%Next
          ENDDO
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          !If this species isn't to be dumped and this isn't a restart dump then
          IF ((.NOT. CurrentFamily%Dump) .AND. (.NOT. Iterator_Settings%Restart)) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_particles

  !Iterator for particle charge
  SUBROUTINE iterate_charge(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root
    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
#ifdef PER_PARTICLE_CHARGEMASS
             data(partcount) = Cur%Charge
#else
             data(partcount) = CurrentFamily%Charge
#endif
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_charge

#ifdef PER_PARTICLE_WEIGHT

  !Iterator for particle weight
  !Only present if you are using the PER_PARTICLE_WEIGHT
  !Precompiler option
  SUBROUTINE iterate_weight(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             data(partcount)=Cur%Weight
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF ((.NOT. CurrentFamily%Dump) .AND. (.NOT. Iterator_Settings%Restart)) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_weight
#endif

  !Iterator for particle mass
  SUBROUTINE iterate_mass(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
#ifdef PER_PARTICLE_CHARGEMASS
             data(partcount) = Cur%Mass
#else
             data(partcount) = CurrentFamily%Mass
#endif
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_mass

#ifdef PART_DEBUG
  !Iterator for particle processor
  SUBROUTINE iterate_processor(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             data(partcount)=REAL(cur%Processor,num)
             IF (cur%Processor .GE. nproc) PRINT *,"Bad Processor"
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_processor

  SUBROUTINE iterate_processor0(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             data(partcount)=REAL(cur%Processor_at_t0,num)
             IF (cur%Processor .GE. nproc) PRINT *,"Bad Processor"
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_processor0
#endif

  !Iterator for particle processor
  SUBROUTINE iterate_species(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             data(partcount)=REAL(CurrentFamily%ID,num)
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF ((.NOT. CurrentFamily%Dump) .AND. (.NOT. Iterator_Settings%Restart)) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_species

  !Iterator for particle velocities
  SUBROUTINE iterate_vx(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root,part_m

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
#ifdef PER_PARTICLE_CHARGEMASS
             part_m = Cur%Mass
#else
             part_m = CurrentFamily%Mass
#endif
#ifdef NEWTONIAN
             root=part_m
#else
             root=SQRT(part_m**2 + (Cur%part_p(1)**2 + Cur%part_p(2)**2 + Cur%part_p(3)**2)/c**2)
#endif
             IF (root .NE. 0.0_num) root=1.0_num/root
             data(partcount) = Cur%Part_p(1) * root
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_vx

  SUBROUTINE iterate_vy(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root,part_m

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
#ifdef PER_PARTICLE_CHARGEMASS
             part_m = Cur%Mass
#else
             part_m = CurrentFamily%Mass
#endif
#ifdef NEWTONIAN
             root=part_m
#else
             root=SQRT(part_m**2 + (Cur%part_p(1)**2 + Cur%part_p(2)**2 + Cur%part_p(3)**2)/c**2)
#endif
             IF (root .NE. 0.0_num) root=1.0_num/root
             data(partcount) = Cur%Part_p(2) * root
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_vy

  SUBROUTINE iterate_vz(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root,part_m
    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
#ifdef PER_PARTICLE_CHARGEMASS
             part_m = Cur%Mass
#else
             part_m = CurrentFamily%Mass
#endif
#ifdef NEWTONIAN
             root=part_m
#else
             root=SQRT(part_m**2 + (Cur%part_p(1)**2 + Cur%part_p(2)**2 + Cur%part_p(3)**2)/c**2)
#endif
             IF (root .NE. 0.0_num) root=1.0_num/root
             data(partcount) = Cur%Part_p(3) * root
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_vz

  !Iterator for particle momenta
  SUBROUTINE iterate_px(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             data(partcount) = Cur%Part_p(1)
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF ((.NOT. CurrentFamily%Dump) .AND. (.NOT. Iterator_Settings%Restart)) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_px

  SUBROUTINE iterate_py(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             data(partcount) = Cur%Part_p(2)
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF ((.NOT. CurrentFamily%Dump) .AND. (.NOT. Iterator_Settings%Restart)) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_py

  SUBROUTINE iterate_pz(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             data(partcount) = Cur%Part_p(3)
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF ((.NOT. CurrentFamily%Dump) .AND. (.NOT. Iterator_Settings%Restart)) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_pz

END MODULE iterators
