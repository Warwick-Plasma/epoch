MODULE probes

#ifdef PARTICLE_PROBES
  USE sdf
  USE mpi_subtype_control
  USE partlist

  IMPLICIT NONE

  SAVE
  TYPE(particle_list), POINTER, PRIVATE :: current_list

CONTAINS

  SUBROUTINE init_probe(probe)

    TYPE(particle_probe), POINTER :: probe

    probe%point = 0.0_num
    probe%normal = 0.0_num
    probe%ek_min = -HUGE(1.0_num)
    probe%ek_max =  HUGE(1.0_num)
    probe%name = blank
    probe%dumpmask = c_io_always
    NULLIFY(probe%next)
    ALLOCATE(probe%use_species(n_species))
    probe%use_species = .FALSE.
    CALL create_empty_partlist(probe%sampled_particles)

  END SUBROUTINE init_probe



  SUBROUTINE attach_probe(probe)

    TYPE(particle_probe), POINTER :: probe
    TYPE(particle_probe), POINTER :: current
    INTEGER :: i

    DO i = 1, n_species
      IF (.NOT. probe%use_species(i)) CYCLE

      current => species_list(i)%attached_probes
      IF (.NOT. ASSOCIATED(current)) THEN
        species_list(i)%attached_probes => probe
        CYCLE
      ENDIF
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      ENDDO
      ! Now at the last element in the list
      current%next => probe
    ENDDO

  END SUBROUTINE attach_probe



  SUBROUTINE write_probes(sdf_handle, code)

    TYPE(sdf_file_handle) :: sdf_handle
    INTEGER, INTENT(IN) :: code

    TYPE(particle_probe), POINTER :: current_probe
    CHARACTER(LEN=string_length) :: probe_name, temp_name
    INTEGER :: ispecies, i
    INTEGER(8) :: npart_probe_global, part_probe_offset
    INTEGER(8), DIMENSION(:), ALLOCATABLE :: npart_probe_per_proc

    ALLOCATE(npart_probe_per_proc(nproc))

    DO ispecies = 1, n_species
      current_probe=>species_list(ispecies)%attached_probes
      DO WHILE(ASSOCIATED(current_probe))
        ! If don't dump this probe currently then just cycle
        IF (IAND(current_probe%dumpmask, code) .EQ. 0) THEN
          current_probe=>current_probe%next
          CYCLE
        ENDIF

        current_list=>current_probe%sampled_particles

        CALL MPI_ALLGATHER(current_probe%sampled_particles%count, 1, &
            MPI_INTEGER8, npart_probe_per_proc, 1, MPI_INTEGER8, comm, errcode)

        npart_probe_global = 0
        DO i = 1, nproc
          IF (rank .EQ. i-1) part_probe_offset = npart_probe_global
          npart_probe_global = npart_probe_global + npart_probe_per_proc(i)
        ENDDO

        IF (npart_probe_global .GT. 0) THEN
          probe_name =  TRIM(ADJUSTL(current_probe%name))

          ! dump particle Positions
          CALL sdf_write_point_mesh(sdf_handle, TRIM(probe_name), &
              'Grid/Probe/' // TRIM(probe_name), npart_probe_global, &
              c_dimension_1d, iterate_probe_particles, part_probe_offset)

          ! dump Px
          WRITE(temp_name, '(a, "/Px")') TRIM(probe_name)
          CALL sdf_write_point_variable(sdf_handle, TRIM(temp_name), &
              TRIM(temp_name), 'Pa', npart_probe_global, TRIM(probe_name), &
              iterate_probe_px, part_probe_offset)

          ! dump Py
          WRITE(temp_name, '(a, "/Py")') TRIM(probe_name)
          CALL sdf_write_point_variable(sdf_handle, TRIM(temp_name), &
              TRIM(temp_name), 'Pa', npart_probe_global, TRIM(probe_name), &
              iterate_probe_py, part_probe_offset)

          ! dump Pz
          WRITE(temp_name, '(a, "/Pz")') TRIM(probe_name)
          CALL sdf_write_point_variable(sdf_handle, TRIM(temp_name), &
              TRIM(temp_name), 'Pa', npart_probe_global, TRIM(probe_name), &
              iterate_probe_pz, part_probe_offset)

          ! dump particle weight function
          WRITE(temp_name, '(a, "/weight")') TRIM(probe_name)
#ifdef PER_PARTICLE_WEIGHT
          CALL sdf_write_point_variable(sdf_handle, TRIM(temp_name), &
              TRIM(temp_name), 'kg', npart_probe_global, TRIM(probe_name), &
              iterate_probe_weight, part_probe_offset)
#else
          CALL sdf_write_srl(sdf_handle, TRIM(temp_name), TRIM(probe_name), &
              species_list(ispecies)%weight)
#endif

          CALL destroy_partlist(current_probe%sampled_particles)
        ENDIF
        current_probe=>current_probe%next

      ENDDO

      NULLIFY(current_probe)
    ENDDO

    DEALLOCATE(npart_probe_per_proc)

  END SUBROUTINE write_probes



  ! iterator for particle positions
  FUNCTION iterate_probe_particles(array, n_points, start, direction)

    REAL(num) :: iterate_probe_particles
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN) :: direction
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER :: part_count

    IF (start)  THEN
      cur=> current_list%head
    ENDIF
    part_count = 0

    DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
      part_count = part_count+1
      array(part_count) = cur%part_pos - window_shift
      cur=>cur%next
    ENDDO

    n_points = part_count

    iterate_probe_particles = 0

  END FUNCTION iterate_probe_particles



  ! iterator for particle momenta
  FUNCTION iterate_probe_px(array, n_points, start)

    REAL(num) :: iterate_probe_px
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER :: part_count

    IF (start)  THEN
      cur=> current_list%head
    ENDIF
    part_count = 0

    DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
      part_count = part_count+1
      array(part_count) = cur%part_p(1)
      cur=>cur%next
    ENDDO
    n_points = part_count

    iterate_probe_px = 0

  END FUNCTION iterate_probe_px



  FUNCTION iterate_probe_py(array, n_points, start)

    REAL(num) :: iterate_probe_py
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER :: part_count

    IF (start)  THEN
      cur=> current_list%head
    ENDIF
    part_count = 0

    DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
      part_count = part_count+1
      array(part_count) = cur%part_p(2)
      cur=>cur%next
    ENDDO

    n_points = part_count

    iterate_probe_py = 0

  END FUNCTION iterate_probe_py



  FUNCTION iterate_probe_pz(array, n_points, start)

    REAL(num) :: iterate_probe_pz
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER :: part_count

    IF (start)  THEN
      cur=> current_list%head
    ENDIF
    part_count = 0

    DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
      part_count = part_count+1
      array(part_count) = cur%part_p(3)
      cur=>cur%next
    ENDDO

    n_points = part_count

    iterate_probe_pz = 0

  END FUNCTION iterate_probe_pz



#ifdef PER_PARTICLE_WEIGHT
  FUNCTION iterate_probe_weight(array, n_points, start)

    REAL(num) :: iterate_probe_weight
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER :: part_count

    IF (start)  THEN
      cur=> current_list%head
    ENDIF
    part_count = 0

    DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
      part_count = part_count+1
      array(part_count) = cur%weight
      cur=>cur%next
    ENDDO

    n_points = part_count

    iterate_probe_weight = 0

  END FUNCTION iterate_probe_weight
#endif
#endif

END MODULE probes
