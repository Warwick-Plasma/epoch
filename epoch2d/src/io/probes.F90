MODULE probes

#ifdef PARTICLE_PROBES
  USE cfd
  USE mpi_subtype_control
  USE partlist

  IMPLICIT NONE

  SAVE
  TYPE(particle_list), POINTER, PRIVATE :: current_list

CONTAINS

  SUBROUTINE init_probe(probe)

    TYPE(particle_probe), POINTER :: probe

    NULLIFY(probe%next)
    NULLIFY(probe%probe_species)
    CALL create_empty_partlist(probe%sampled_particles)

  END SUBROUTINE init_probe



  SUBROUTINE attach_probe(probe)

    TYPE(particle_probe), POINTER :: probe
    TYPE(particle_probe), POINTER :: current

    current=>probe%probe_species%attached_probes
    IF (.NOT. ASSOCIATED(current)) THEN
      probe%probe_species%attached_probes=>probe
      RETURN
    ENDIF
    DO WHILE(ASSOCIATED(current%next))
      current=>current%next
    ENDDO
    ! Now at the last element in the list
    current%next=>probe

  END SUBROUTINE attach_probe



  SUBROUTINE write_probes(code)

    INTEGER, INTENT(IN) :: code

    TYPE(particle_probe), POINTER :: current_probe
    CHARACTER(LEN=string_length) :: probe_name, temp_name
    INTEGER :: ispecies
    INTEGER(8) :: npart_probe_local, npart_probe_global
    INTEGER(8) :: npart_probe_per_it, npart_probe_per_it_local
    INTEGER(KIND=MPI_OFFSET_KIND), DIMENSION(1) :: file_lengths, file_offsets

    DO ispecies = 1, n_species
      current_probe=>particle_species(ispecies)%attached_probes
      DO WHILE(ASSOCIATED(current_probe))
        ! If don't dump this probe currently then just cycle
        IF (IAND(current_probe%dump, code) .EQ. 0) THEN
          current_probe=>current_probe%next
          CYCLE
        ENDIF

        current_list=>current_probe%sampled_particles

        npart_probe_per_it_local = npart_per_it
        npart_probe_local = current_probe%sampled_particles%count

        IF (npart_probe_local .GT. 0) &
            npart_probe_per_it_local = &
                MIN(npart_probe_local, npart_probe_per_it_local)

        CALL MPI_ALLREDUCE(npart_probe_local, npart_probe_global, 1, &
            MPI_INTEGER8, MPI_SUM, comm, errcode)
        CALL MPI_ALLREDUCE(npart_probe_per_it_local, npart_probe_per_it, 1, &
            MPI_INTEGER8, MPI_MIN, comm, errcode)

        IF (npart_probe_global .GT. 0) THEN
          file_lengths(1) = npart_probe_local
          file_offsets(1) = create_particle_offset(npart_probe_local)

          probe_name =  TRIM(ADJUSTL(current_probe%name))

          ! dump particle Positions
          CALL cfd_write_nd_particle_grid_with_iterator_all(&
              TRIM(probe_name), "Probe_Grid", iterate_probe_particles, &
              c_dimension_2d, npart_probe_local, npart_probe_global, &
              npart_probe_per_it, c_particle_cartesian, &
              file_lengths, file_offsets)

          ! dump Px
          WRITE(temp_name, '(a, "_Px")') TRIM(probe_name)
          CALL cfd_write_nd_particle_variable_with_iterator_all(&
              TRIM(temp_name), TRIM(probe_name), iterate_probe_px, &
              npart_probe_global, npart_probe_per_it, TRIM(probe_name), &
              "Probe_Grid", file_lengths, file_offsets)

          ! dump Py
          WRITE(temp_name, '(a, "_Py")') TRIM(probe_name)
          CALL cfd_write_nd_particle_variable_with_iterator_all(&
              TRIM(temp_name), TRIM(probe_name), iterate_probe_py, &
              npart_probe_global, npart_probe_per_it, TRIM(probe_name), &
              "Probe_Grid", file_lengths, file_offsets)

          ! dump Pz
          WRITE(temp_name, '(a, "_Pz")') TRIM(probe_name)
          CALL cfd_write_nd_particle_variable_with_iterator_all(&
              TRIM(temp_name), TRIM(probe_name), iterate_probe_pz, &
              npart_probe_global, npart_probe_per_it, TRIM(probe_name), &
              "Probe_Grid", file_lengths, file_offsets)

          ! dump particle weight function
          WRITE(temp_name, '(a, "_weight")') TRIM(probe_name)
#ifdef PER_PARTICLE_WEIGHT
          CALL cfd_write_nd_particle_variable_with_iterator_all(&
              TRIM(temp_name), TRIM(probe_name), iterate_probe_weight, &
              npart_probe_global, npart_probe_per_it, TRIM(probe_name), &
              "Probe_Grid", file_lengths, file_offsets)
#else
          CALL cfd_write_real_constant(TRIM(temp_name), TRIM(probe_name), &
              weight, 0)
#endif

          CALL destroy_partlist(current_probe%sampled_particles)
        ENDIF
        current_probe=>current_probe%next

      ENDDO

      NULLIFY(current_probe)
    ENDDO

  END SUBROUTINE write_probes



  ! iterator for particle positions
  SUBROUTINE iterate_probe_particles(data, n_points, direction, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN) :: direction
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    IF (start)  THEN
      cur=> current_list%head
    ENDIF
    part_count = 0

    DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
      part_count = part_count+1
      data(part_count) = cur%part_pos(direction)-window_shift(direction)
      cur=>cur%next
    ENDDO

    n_points = part_count

  END SUBROUTINE iterate_probe_particles



  ! iterator for particle momenta
  SUBROUTINE iterate_probe_px(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    IF (start)  THEN
      cur=> current_list%head
    ENDIF
    part_count = 0

    DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
      part_count = part_count+1
      data(part_count) = cur%part_p(1)
      cur=>cur%next
    ENDDO
    n_points = part_count

  END SUBROUTINE iterate_probe_px



  SUBROUTINE iterate_probe_py(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    IF (start)  THEN
      cur=> current_list%head
    ENDIF
    part_count = 0

    DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
      part_count = part_count+1
      data(part_count) = cur%part_p(2)
      cur=>cur%next
    ENDDO

    n_points = part_count

  END SUBROUTINE iterate_probe_py



  SUBROUTINE iterate_probe_pz(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    IF (start)  THEN
      cur=> current_list%head
    ENDIF
    part_count = 0

    DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
      part_count = part_count+1
      data(part_count) = cur%part_p(3)
      cur=>cur%next
    ENDDO

    n_points = part_count

  END SUBROUTINE iterate_probe_pz



#ifdef PER_PARTICLE_WEIGHT
  SUBROUTINE iterate_probe_weight(data, n_points, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(8), INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER(8) :: part_count

    IF (start)  THEN
      cur=> current_list%head
    ENDIF
    part_count = 0

    DO WHILE (ASSOCIATED(cur) .AND. (part_count .LT. n_points))
      part_count = part_count+1
      data(part_count) = cur%weight
      cur=>cur%next
    ENDDO

    n_points = part_count

  END SUBROUTINE iterate_probe_weight
#endif
#endif

END MODULE probes
