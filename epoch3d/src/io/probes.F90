! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE probes

#ifdef NO_PARTICLE_PROBES

CONTAINS

  SUBROUTINE deallocate_probes
  END SUBROUTINE deallocate_probes

#else
  USE partlist
  USE sdf

  IMPLICIT NONE

  SAVE
  TYPE(particle_list), POINTER, PRIVATE :: current_list
  TYPE(particle_species), POINTER, PRIVATE :: current_species

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
    CALL create_empty_partlist(probe%sampled_particles, holds_copies=.TRUE.)

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
      END IF
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      END DO
      ! Now at the last element in the list
      current%next => probe
    END DO

  END SUBROUTINE attach_probe



  SUBROUTINE deallocate_probes

    TYPE(particle_probe), POINTER :: current_probe, next
    INTEGER :: ispecies, j

    DO ispecies = 1, n_species
      current_probe => species_list(ispecies)%attached_probes
      DO j = ispecies + 1, n_species
        IF (ASSOCIATED(current_probe, species_list(j)%attached_probes)) &
            NULLIFY(species_list(j)%attached_probes)
      END DO
      next => current_probe
      DO WHILE(ASSOCIATED(next))
        current_probe => next
        next => current_probe%next
        IF (ASSOCIATED(current_probe%use_species)) &
            DEALLOCATE(current_probe%use_species)
        CALL destroy_partlist(current_probe%sampled_particles)
        DEALLOCATE(current_probe)
      END DO
    END DO

  END SUBROUTINE deallocate_probes



  SUBROUTINE write_probes(sdf_handle, code, mask)

    TYPE(sdf_file_handle) :: sdf_handle
    INTEGER, INTENT(IN) :: code, mask

    TYPE(particle_probe), POINTER :: current_probe
    CHARACTER(LEN=string_length) :: probe_name, temp_name
    INTEGER :: ispecies, i
    INTEGER(i8) :: npart_probe_global, part_probe_offset
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: npart_probe_per_proc
    LOGICAL :: convert

    ALLOCATE(npart_probe_per_proc(nproc))

    DO ispecies = 1, n_species
      current_species => species_list(ispecies)
      current_probe => species_list(ispecies)%attached_probes
      DO WHILE(ASSOCIATED(current_probe))
        ! If don't dump this probe currently then just cycle
        IF (IAND(current_probe%dumpmask, code) == 0) THEN
          current_probe => current_probe%next
          CYCLE
        END IF

        current_list => current_probe%sampled_particles

        CALL MPI_ALLGATHER(current_probe%sampled_particles%count, 1, &
            MPI_INTEGER8, npart_probe_per_proc, 1, MPI_INTEGER8, comm, errcode)

        npart_probe_global = 0
        DO i = 1, nproc
          IF (rank == i-1) part_probe_offset = npart_probe_global
          npart_probe_global = npart_probe_global + npart_probe_per_proc(i)
        END DO

        IF (npart_probe_global > 0) THEN
          convert = (IAND(IOR(mask,current_probe%dumpmask), &
                          c_io_dump_single) /= 0)

          probe_name =  TRIM(ADJUSTL(current_probe%name))

          ! dump particle Positions
          CALL sdf_write_point_mesh(sdf_handle, TRIM(probe_name), &
              'Grid/Probe/' // TRIM(probe_name), TRIM(probe_name), &
              npart_probe_global, c_dimension_3d, it_probe_position, &
              part_probe_offset, convert)

          ! dump Px
          WRITE(temp_name, '(a, ''/Px'')') TRIM(probe_name)
          CALL sdf_write_point_variable(sdf_handle, TRIM(temp_name), &
              TRIM(temp_name), TRIM(probe_name), 'kg.m/s', npart_probe_global, &
              TRIM(probe_name), it_probe_real, c_dump_part_px, &
              part_probe_offset, convert)

          ! dump Py
          WRITE(temp_name, '(a, ''/Py'')') TRIM(probe_name)
          CALL sdf_write_point_variable(sdf_handle, TRIM(temp_name), &
              TRIM(temp_name), TRIM(probe_name), 'kg.m/s', npart_probe_global, &
              TRIM(probe_name), it_probe_real, c_dump_part_py, &
              part_probe_offset, convert)

          ! dump Pz
          WRITE(temp_name, '(a, ''/Pz'')') TRIM(probe_name)
          CALL sdf_write_point_variable(sdf_handle, TRIM(temp_name), &
              TRIM(temp_name), TRIM(probe_name), 'kg.m/s', npart_probe_global, &
              TRIM(probe_name), it_probe_real, c_dump_part_pz, &
              part_probe_offset, convert)

          ! dump ID
#ifdef PARTICLE_ID
          WRITE(temp_name, '(a, ''/ID'')') TRIM(probe_name)
          CALL sdf_write_point_variable(sdf_handle, TRIM(temp_name), &
              TRIM(temp_name), TRIM(probe_name), '', npart_probe_global, &
              TRIM(probe_name), it_probe_integer8, c_dump_part_id, &
              part_probe_offset, convert)
#elif PARTICLE_ID4
          WRITE(temp_name, '(a, ''/ID'')') TRIM(probe_name)
          CALL sdf_write_point_variable(sdf_handle, TRIM(temp_name), &
              TRIM(temp_name), TRIM(probe_name), '', npart_probe_global, &
              TRIM(probe_name), it_probe_integer4, c_dump_part_id, &
              part_probe_offset, convert)
#endif

          ! dump time_at_probe
#ifdef PROBE_TIME
          WRITE(temp_name, '(a, ''/time_at_probe'')') TRIM(probe_name)
          CALL sdf_write_point_variable(sdf_handle, TRIM(temp_name), &
              TRIM(temp_name), TRIM(probe_name), '', npart_probe_global, &
              TRIM(probe_name), it_probe_real, c_dump_probe_time, &
              part_probe_offset, convert)
#endif

          ! dump particle weight function
          WRITE(temp_name, '(a, ''/weight'')') TRIM(probe_name)
#ifndef PER_SPECIES_WEIGHT
          CALL sdf_write_point_variable(sdf_handle, TRIM(temp_name), &
              TRIM(temp_name), TRIM(probe_name), '', npart_probe_global, &
              TRIM(probe_name), it_probe_real, c_dump_part_weight, &
              part_probe_offset, convert)
#else
          CALL sdf_write_srl(sdf_handle, TRIM(temp_name), TRIM(probe_name), &
              species_list(ispecies)%weight)
#endif

          CALL destroy_partlist(current_probe%sampled_particles)
        END IF
        current_probe => current_probe%next

      END DO

      NULLIFY(current_probe)
    END DO

    DEALLOCATE(npart_probe_per_proc)

  END SUBROUTINE write_probes



  ! iterator for particle positions
  FUNCTION it_probe_position(array, npoint_it, start, direction, param)

    REAL(num) :: it_probe_position
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: npoint_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN) :: direction
    INTEGER, INTENT(IN), OPTIONAL :: param
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER :: part_count
    REAL(num) :: window_shift

    IF (start)  THEN
      cur => current_list%head
    END IF
    part_count = 0

    IF (use_offset_grid .AND. direction == c_dir_x) THEN
      window_shift = window_offset
    ELSE
      window_shift = 0.0_num
    END IF

    DO WHILE (ASSOCIATED(cur) .AND. (part_count < npoint_it))
      part_count = part_count + 1
      array(part_count) = cur%part_pos(direction) - window_shift
      cur => cur%next
    END DO

    npoint_it = part_count

    it_probe_position = 0

  END FUNCTION it_probe_position



  FUNCTION it_probe_real(array, npoint_it, start, param)

    REAL(num) :: it_probe_real
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: npoint_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER :: part_count, ndim

    IF (start)  THEN
      cur => current_list%head
    END IF
    part_count = 0

    SELECT CASE (param)
#ifndef PER_SPECIES_WEIGHT
    CASE (c_dump_part_weight)
      DO WHILE (ASSOCIATED(cur) .AND. (part_count < npoint_it))
        part_count = part_count+1
        array(part_count) = cur%weight
        cur => cur%next
      END DO
#endif

#ifdef PROBE_TIME
    CASE (c_dump_probe_time)
      DO WHILE (ASSOCIATED(cur) .AND. (part_count < npoint_it))
        part_count = part_count+1
        array(part_count) = cur%probe_time
        cur => cur%next
      END DO
#endif

    CASE (c_dump_part_px)
      ndim = 1
      DO WHILE (ASSOCIATED(cur) .AND. (part_count < npoint_it))
        part_count = part_count + 1
        array(part_count) = cur%part_p(ndim)
        cur => cur%next
      END DO

    CASE (c_dump_part_py)
      ndim = 2
      DO WHILE (ASSOCIATED(cur) .AND. (part_count < npoint_it))
        part_count = part_count + 1
        array(part_count) = cur%part_p(ndim)
        cur => cur%next
      END DO

    CASE (c_dump_part_pz)
      ndim = 3
      DO WHILE (ASSOCIATED(cur) .AND. (part_count < npoint_it))
        part_count = part_count + 1
        array(part_count) = cur%part_p(ndim)
        cur => cur%next
      END DO

    END SELECT

    npoint_it = part_count

    it_probe_real = 0

  END FUNCTION it_probe_real



#ifdef PARTICLE_ID
  FUNCTION it_probe_integer8(array, npoint_it, start, param)

    INTEGER(i8) :: it_probe_integer8
    INTEGER(i8), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: npoint_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER :: part_count, ndim

    IF (start)  THEN
      cur => current_list%head
      CALL generate_particle_ids(current_list)
    END IF
    part_count = 0

    SELECT CASE (param)

    CASE (c_dump_part_id)
      ndim = 1
      DO WHILE (ASSOCIATED(cur) .AND. (part_count < npoint_it))
        part_count = part_count + 1
        array(part_count) = cur%id
        cur => cur%next
      END DO

    END SELECT

    npoint_it = part_count

    it_probe_integer8 = 0

  END FUNCTION it_probe_integer8
#endif



#ifdef PARTICLE_ID4
  FUNCTION it_probe_integer4(array, npoint_it, start, param)

    INTEGER(i4) :: it_probe_integer4
    INTEGER(i4), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: npoint_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    TYPE(particle), POINTER, SAVE :: cur
    INTEGER :: part_count, ndim

    IF (start)  THEN
      cur => current_list%head
      CALL generate_particle_ids(current_list)
    END IF
    part_count = 0

    SELECT CASE (param)

    CASE (c_dump_part_id)
      ndim = 1
      DO WHILE (ASSOCIATED(cur) .AND. (part_count < npoint_it))
        part_count = part_count + 1
        array(part_count) = cur%id
        cur => cur%next
      END DO

    END SELECT

    npoint_it = part_count

    it_probe_integer4 = 0

  END FUNCTION it_probe_integer4
#endif
#endif

END MODULE probes
