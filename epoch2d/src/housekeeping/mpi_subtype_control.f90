! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2010 Chris Brady <C.S.Brady@warwick.ac.uk>
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

MODULE mpi_subtype_control

  !----------------------------------------------------------------------------
  ! This module contains the subroutines which create the subtypes used in
  ! IO
  !----------------------------------------------------------------------------

  USE shared_data

  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! get_total_local_particles - Returns the number of particles on this
  ! processor.
  !----------------------------------------------------------------------------

  FUNCTION get_total_local_particles()

    ! This subroutine describes the total number of particles on the current
    ! processor. It simply sums over every particle species

    INTEGER(i8) :: get_total_local_particles
    INTEGER :: ispecies

    get_total_local_particles = 0
    DO ispecies = 1, n_species
      get_total_local_particles = get_total_local_particles &
          + species_list(ispecies)%attached_list%count
    END DO

  END FUNCTION get_total_local_particles



  !----------------------------------------------------------------------------
  ! CreateSubtypes - Creates the subtypes used by the main output routines
  ! Run just before output takes place
  !----------------------------------------------------------------------------

  SUBROUTINE create_subtypes

    INTEGER :: i

    ! This subroutines creates the MPI types which represent the data for the
    ! field and particles data. It is used when writing data

    ! Actually create the subtypes
    subtype_field = create_current_field_subtype()
    subarray_field = create_current_field_subarray(ng)
    subarray_field_big = create_current_field_subarray(jng)

    subtype_field_r4 = create_current_field_subtype(MPI_REAL4)
    subarray_field_r4 = create_current_field_subarray(ng, MPI_REAL4)
    subarray_field_big_r4 = create_current_field_subarray(jng, MPI_REAL4)

    DO i = 1, n_subsets
      CALL create_subset_subtypes(i)
    END DO

  END SUBROUTINE create_subtypes



  !----------------------------------------------------------------------------
  ! Creates the subtypes used by the main output routines for subsets
  ! Run every time a time-varying subset is considered for output
  !----------------------------------------------------------------------------

  SUBROUTINE create_subset_subtypes(subset_index)

    INTEGER, INTENT(IN) :: subset_index
    INTEGER :: j, rd, n_min, n_max, mpitype, npd, npdm, n0
    INTEGER, DIMENSION(c_ndims) :: n_global, sn_local, sn_global, starts
    INTEGER, DIMENSION(2,c_ndims) :: ranges
    LOGICAL :: proc_outside_range
    TYPE(subset), POINTER :: sub

    sub => subset_list(subset_index)

    n_global = (/nx_global, ny_global/)
    starts = 0
    n0 = 1

    ranges = cell_global_ranges(global_ranges(sub))
    sn_global = (/ ranges(2,1) - ranges(1,1), ranges(2,2) - ranges(1,2) /)
    ranges = cell_local_ranges(global_ranges(sub))

    ! These calculations rely on the original domain size, so will be wrong
    ! for skipped sets as yet
    proc_outside_range = .FALSE.
    IF (ranges(2,1) - ranges(1,1) <= c_tiny) proc_outside_range = .TRUE.
    IF (ranges(2,2) - ranges(1,2) <= c_tiny) proc_outside_range = .TRUE.
    sn_local =  (/ ranges(2,1) - ranges(1,1), ranges(2,2) - ranges(1,2) /)
    starts = cell_starts(ranges, global_ranges(sub))

    ranges = cell_section_ranges(ranges)

    IF (sub%skip) THEN
      DO j = 1, c_ndims
        rd = sub%skip_dir(j)
        n_min = n_global_min(j)
        n_max = n_global_max(j)
        npd = (n_max - n0) / rd + 1
        npdm = (n_min - 1 - n0) / rd + 1
        IF (n_min < 2) npdm = 0
        sn_global(j) = (n_global(j) - n0) / rd + 1
        sn_local(j) = npd - npdm
        sub%n_start(j) = n0 + npdm * rd - n_min
        starts(j) = npdm
      END DO

      ! Just exit if the subset hasn't changed extents since last time
      IF (equal(sub%n_global, sn_global)) THEN
        IF (equal(sub%n_local, sn_local)) THEN
          IF (equal(sub%starts, starts)) THEN
            RETURN
          END IF
        END IF
      END IF

      sub%n_global(:) = sn_global(:)
      sub%n_local(:) = sn_local(:)
      sub%starts(:) = starts(:)

      IF (sub%subtype /= MPI_DATATYPE_NULL) THEN
        CALL MPI_TYPE_FREE(sub%subtype, errcode)
        CALL MPI_TYPE_FREE(sub%subtype_r4, errcode)
        CALL MPI_TYPE_FREE(sub%subarray, errcode)
        CALL MPI_TYPE_FREE(sub%subarray_r4, errcode)
      END IF

      mpitype = MPI_DATATYPE_NULL
      CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sub%n_global, sub%n_local, &
          starts, MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
      CALL MPI_TYPE_COMMIT(mpitype, errcode)

      sub%subtype = mpitype

      mpitype = MPI_DATATYPE_NULL
      CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sub%n_global, sub%n_local, &
          starts, MPI_ORDER_FORTRAN, MPI_REAL4, mpitype, errcode)
      CALL MPI_TYPE_COMMIT(mpitype, errcode)

      sub%subtype_r4 = mpitype

      starts = 0
      mpitype = MPI_DATATYPE_NULL
      CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sub%n_local, sub%n_local, &
          starts, MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
      CALL MPI_TYPE_COMMIT(mpitype, errcode)

      sub%subarray = mpitype

      mpitype = MPI_DATATYPE_NULL
      CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sub%n_local, sub%n_local, &
          starts, MPI_ORDER_FORTRAN, MPI_REAL4, mpitype, errcode)
      CALL MPI_TYPE_COMMIT(mpitype, errcode)

      sub%subarray_r4 = mpitype
    ELSE
      ! Just exit if the subset hasn't changed extents since last time
      IF (equal(sub%n_global, sn_global)) THEN
        IF (equal(sub%n_local, sn_local)) THEN
          IF (equal(sub%starts, starts)) THEN
            RETURN
          END IF
        END IF
      END IF

      sub%n_global(:) = sn_global(:)
      sub%n_local(:) = sn_local(:)
      sub%starts(:) = starts(:)

      IF (sub%subtype /= MPI_DATATYPE_NULL) THEN
        CALL MPI_TYPE_FREE(sub%subtype, errcode)
        CALL MPI_TYPE_FREE(sub%subtype_r4, errcode)
        CALL MPI_TYPE_FREE(sub%subarray, errcode)
        CALL MPI_TYPE_FREE(sub%subarray_r4, errcode)
      END IF

      mpitype = MPI_DATATYPE_NULL
      IF (proc_outside_range) THEN
        CALL MPI_TYPE_CONTIGUOUS(0, mpireal, mpitype, errcode)
      ELSE
        CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sub%n_global, sub%n_local, &
            starts, MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
      END IF
      CALL MPI_TYPE_COMMIT(mpitype, errcode)

      sub%subtype = mpitype
      mpitype = MPI_DATATYPE_NULL
      IF (proc_outside_range) THEN
        CALL MPI_TYPE_CONTIGUOUS(0, MPI_REAL4, mpitype, errcode)
      ELSE
        CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sub%n_global, sub%n_local, &
            starts, MPI_ORDER_FORTRAN, MPI_REAL4, mpitype, errcode)
      END IF
      CALL MPI_TYPE_COMMIT(mpitype, errcode)

      sub%subtype_r4 = mpitype

      starts = 0
      mpitype = MPI_DATATYPE_NULL
      IF (proc_outside_range) THEN
        CALL MPI_TYPE_CONTIGUOUS(0, mpireal, mpitype, errcode)
      ELSE
        CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sub%n_local, sub%n_local, &
            starts, MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
      END IF
      CALL MPI_TYPE_COMMIT(mpitype, errcode)

      sub%subarray = mpitype

      mpitype = MPI_DATATYPE_NULL
      IF (proc_outside_range) THEN
        CALL MPI_TYPE_CONTIGUOUS(0, MPI_REAL4, mpitype, errcode)
      ELSE
        CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sub%n_local, sub%n_local, &
            starts, MPI_ORDER_FORTRAN, MPI_REAL4, mpitype, errcode)
      END IF
      CALL MPI_TYPE_COMMIT(mpitype, errcode)

      sub%subarray_r4 = mpitype
    END IF

  END SUBROUTINE create_subset_subtypes



  !----------------------------------------------------------------------------
  ! Frees the subtypes created by create_subtypes
  !----------------------------------------------------------------------------

  SUBROUTINE free_subtypes

    CALL MPI_TYPE_FREE(subtype_field, errcode)
    CALL MPI_TYPE_FREE(subarray_field, errcode)
    CALL MPI_TYPE_FREE(subarray_field_big, errcode)

    CALL MPI_TYPE_FREE(subtype_field_r4, errcode)
    CALL MPI_TYPE_FREE(subarray_field_r4, errcode)
    CALL MPI_TYPE_FREE(subarray_field_big_r4, errcode)

  END SUBROUTINE free_subtypes



  !----------------------------------------------------------------------------
  ! create_current_field_subtype - Creates the subtype corresponding to the
  ! current load balanced geometry
  !----------------------------------------------------------------------------

  FUNCTION create_current_field_subtype(basetype_in)

    INTEGER :: create_current_field_subtype
    INTEGER, OPTIONAL, INTENT(IN) :: basetype_in
    INTEGER :: basetype

    IF (PRESENT(basetype_in)) THEN
      basetype = basetype_in
    ELSE
      basetype = mpireal
    END IF

    create_current_field_subtype = &
        create_field_subtype(basetype, nx, ny, nx_global_min, ny_global_min)

  END FUNCTION create_current_field_subtype



  !----------------------------------------------------------------------------
  ! create_current_field_subarray - Creates the subarray corresponding to the
  ! current load balanced geometry
  !----------------------------------------------------------------------------

  FUNCTION create_current_field_subarray(ng, basetype_in)

    INTEGER :: create_current_field_subarray
    INTEGER, INTENT(IN) :: ng
    INTEGER, OPTIONAL, INTENT(IN) :: basetype_in
    INTEGER :: basetype

    IF (PRESENT(basetype_in)) THEN
      basetype = basetype_in
    ELSE
      basetype = mpireal
    END IF

    create_current_field_subarray = create_field_subarray(basetype, ng, nx, ny)

  END FUNCTION create_current_field_subarray



  !----------------------------------------------------------------------------
  ! create_subtypes_for_load - Creates subtypes when the code loads initial
  ! conditions from a file
  !----------------------------------------------------------------------------

  SUBROUTINE create_subtypes_for_load(species_subtypes, species_subtypes_i4, &
      species_subtypes_i8)

    ! This subroutines creates the MPI types which represent the data for the
    ! field and particles data. It is used when reading data.

    INTEGER, ALLOCATABLE :: species_subtypes(:)
    INTEGER, ALLOCATABLE :: species_subtypes_i4(:), species_subtypes_i8(:)
    INTEGER :: i

    subtype_field = create_current_field_subtype()
    subarray_field = create_current_field_subarray(ng)
    subarray_field_big = create_current_field_subarray(jng)

    subtype_field_r4 = create_current_field_subtype(MPI_REAL4)
    subarray_field_r4 = create_current_field_subarray(ng, MPI_REAL4)
    subarray_field_big_r4 = create_current_field_subarray(jng, MPI_REAL4)

    ALLOCATE(species_subtypes(n_species))
    ALLOCATE(species_subtypes_i4(n_species))
    ALLOCATE(species_subtypes_i8(n_species))
    DO i = 1,n_species
      CALL create_particle_subtypes(species_list(i)%attached_list%count, &
          species_subtypes(i), species_subtypes_i4(i), species_subtypes_i8(i))
    END DO

  END SUBROUTINE create_subtypes_for_load



  !----------------------------------------------------------------------------
  ! free_subtypes_for_load - Frees subtypes created by create_subtypes_for_load
  !----------------------------------------------------------------------------

  SUBROUTINE free_subtypes_for_load(species_subtypes, species_subtypes_i4, &
      species_subtypes_i8)

    INTEGER, ALLOCATABLE :: species_subtypes(:)
    INTEGER, ALLOCATABLE :: species_subtypes_i4(:), species_subtypes_i8(:)
    INTEGER :: i

    CALL MPI_TYPE_FREE(subtype_field, errcode)
    CALL MPI_TYPE_FREE(subarray_field, errcode)
    CALL MPI_TYPE_FREE(subarray_field_big, errcode)

    CALL MPI_TYPE_FREE(subtype_field_r4, errcode)
    CALL MPI_TYPE_FREE(subarray_field_r4, errcode)
    CALL MPI_TYPE_FREE(subarray_field_big_r4, errcode)
    DO i = 1,n_species
      CALL MPI_TYPE_FREE(species_subtypes(i), errcode)
      CALL MPI_TYPE_FREE(species_subtypes_i4(i), errcode)
      CALL MPI_TYPE_FREE(species_subtypes_i8(i), errcode)
    END DO
    DEALLOCATE(species_subtypes)
    DEALLOCATE(species_subtypes_i4)
    DEALLOCATE(species_subtypes_i8)

  END SUBROUTINE free_subtypes_for_load



  !----------------------------------------------------------------------------
  ! create_particle_subtype - Creates a subtype representing the local
  ! particles
  !----------------------------------------------------------------------------

  SUBROUTINE create_particle_subtypes(npart_in, subtype, subtype_i4, subtype_i8)

    INTEGER(i8), INTENT(IN) :: npart_in
    INTEGER, INTENT(OUT) :: subtype, subtype_i4, subtype_i8
    INTEGER(i8), DIMENSION(1) :: npart_local
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: npart_each_rank
    INTEGER(KIND=MPI_ADDRESS_KIND) :: particles_to_skip, total_particles
    INTEGER :: i, mpitype, basetype, typesize, intertype

    npart_local = npart_in

    ALLOCATE(npart_each_rank(nproc))

    ! Create the subarray for the particles in this problem: subtype decribes
    ! where this process's data fits into the global picture.
    CALL MPI_ALLGATHER(npart_local, 1, MPI_INTEGER8, &
        npart_each_rank, 1, MPI_INTEGER8, comm, errcode)

    particles_to_skip = 0
    DO i = 1, rank
      particles_to_skip = particles_to_skip &
          + INT(npart_each_rank(i), MPI_ADDRESS_KIND)
    END DO

    total_particles = particles_to_skip
    DO i = rank+1, nproc
      total_particles = total_particles &
          + INT(npart_each_rank(i), MPI_ADDRESS_KIND)
    END DO

    DEALLOCATE(npart_each_rank)

    mpitype = 0
    basetype = mpireal
    CALL MPI_TYPE_SIZE(basetype, typesize, errcode)
    CALL MPI_TYPE_CREATE_HINDEXED(1, [INT(total_particles,i4)], &
        [particles_to_skip * INT(typesize, MPI_ADDRESS_KIND)], basetype, &
        intertype, errcode)
    CALL MPI_TYPE_CREATE_RESIZED(intertype, 0_MPI_ADDRESS_KIND, &
        total_particles * INT(typesize, MPI_ADDRESS_KIND), subtype, errcode)
    CALL MPI_TYPE_COMMIT(subtype, errcode)

    mpitype = 0
    basetype = MPI_INTEGER4
    CALL MPI_TYPE_SIZE(basetype, typesize, errcode)
    CALL MPI_TYPE_CREATE_HINDEXED(1, [INT(total_particles,i4)], &
        [particles_to_skip * INT(typesize, MPI_ADDRESS_KIND)], basetype, &
        intertype, errcode)
    CALL MPI_TYPE_CREATE_RESIZED(intertype, 0_MPI_ADDRESS_KIND, &
        total_particles * INT(typesize, MPI_ADDRESS_KIND), subtype_i4, errcode)
    CALL MPI_TYPE_COMMIT(subtype_i4, errcode)

    mpitype = 0
    basetype = MPI_INTEGER8
    CALL MPI_TYPE_SIZE(basetype, typesize, errcode)
    CALL MPI_TYPE_CREATE_HINDEXED(1, [INT(total_particles,i4)], &
        [particles_to_skip * INT(typesize, MPI_ADDRESS_KIND)], basetype, &
        intertype, errcode)
    CALL MPI_TYPE_CREATE_RESIZED(intertype, 0_MPI_ADDRESS_KIND, &
        total_particles * INT(typesize, MPI_ADDRESS_KIND), subtype_i8, errcode)
    CALL MPI_TYPE_COMMIT(subtype_i8, errcode)

  END SUBROUTINE create_particle_subtypes



  !----------------------------------------------------------------------------
  ! create_field_subtype - Creates a subtype representing the local processor
  ! for any arbitrary arrangement of an array covering the entire spatial
  ! domain. Only used directly during load balancing
  !----------------------------------------------------------------------------

  FUNCTION create_field_subtype(basetype, nx_local, ny_local, &
      cell_start_x_local, cell_start_y_local)

    INTEGER, INTENT(IN) :: basetype
    INTEGER, INTENT(IN) :: nx_local
    INTEGER, INTENT(IN) :: ny_local
    INTEGER, INTENT(IN) :: cell_start_x_local
    INTEGER, INTENT(IN) :: cell_start_y_local
    INTEGER :: create_field_subtype
    INTEGER, DIMENSION(c_ndims) :: n_local, n_global, start

    n_local = (/nx_local, ny_local/)
    n_global = (/nx_global, ny_global/)
    start = (/cell_start_x_local, cell_start_y_local/)

    create_field_subtype = &
        create_2d_array_subtype(basetype, n_local, n_global, start)

  END FUNCTION create_field_subtype



  !----------------------------------------------------------------------------
  ! create_1d_array_subtype - Creates a subtype representing the local fraction
  ! of a completely arbitrary 1D array. Does not assume anything about the
  ! domain at all.
  !----------------------------------------------------------------------------

  FUNCTION create_1d_array_subtype(basetype, n_local, n_global, start) &
      RESULT(vec1d_sub)

    INTEGER, INTENT(IN) :: basetype
    INTEGER, DIMENSION(1), INTENT(IN) :: n_local
    INTEGER, DIMENSION(1), INTENT(IN) :: n_global
    INTEGER, DIMENSION(1), INTENT(IN) :: start
    INTEGER :: vec1d_sub

    vec1d_sub = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(1, n_global, n_local, start - 1, &
        MPI_ORDER_FORTRAN, basetype, vec1d_sub, errcode)
    CALL MPI_TYPE_COMMIT(vec1d_sub, errcode)

  END FUNCTION create_1d_array_subtype



  !----------------------------------------------------------------------------
  ! create_2d_array_subtype - Creates a subtype representing the local fraction
  ! of a completely arbitrary 2D array. Does not assume anything about the
  ! domain at all.
  !----------------------------------------------------------------------------

  FUNCTION create_2d_array_subtype(basetype, n_local, n_global, start) &
      RESULT(vec2d_sub)

    INTEGER, INTENT(IN) :: basetype
    INTEGER, DIMENSION(2), INTENT(IN) :: n_local
    INTEGER, DIMENSION(2), INTENT(IN) :: n_global
    INTEGER, DIMENSION(2), INTENT(IN) :: start
    INTEGER :: vec2d_sub

    vec2d_sub = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(2, n_global, n_local, start - 1, &
        MPI_ORDER_FORTRAN, basetype, vec2d_sub, errcode)
    CALL MPI_TYPE_COMMIT(vec2d_sub, errcode)

  END FUNCTION create_2d_array_subtype



  !----------------------------------------------------------------------------
  ! create_3d_array_subtype - Creates a subtype representing the local fraction
  ! of a completely arbitrary 3D array. Does not assume anything about the
  ! domain at all.
  !----------------------------------------------------------------------------

  FUNCTION create_3d_array_subtype(basetype, n_local, n_global, start) &
      RESULT(vec3d_sub)

    INTEGER, INTENT(IN) :: basetype
    INTEGER, DIMENSION(3), INTENT(IN) :: n_local
    INTEGER, DIMENSION(3), INTENT(IN) :: n_global
    INTEGER, DIMENSION(3), INTENT(IN) :: start
    INTEGER :: vec3d_sub

    vec3d_sub = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_SUBARRAY(3, n_global, n_local, start - 1, &
        MPI_ORDER_FORTRAN, basetype, vec3d_sub, errcode)
    CALL MPI_TYPE_COMMIT(vec3d_sub, errcode)

  END FUNCTION create_3d_array_subtype



  FUNCTION create_field_subarray(basetype, ng, n1, n2, n3)

    INTEGER, INTENT(IN) :: basetype, ng, n1
    INTEGER, INTENT(IN), OPTIONAL :: n2, n3
    INTEGER, DIMENSION(3) :: n_local, n_global, start
    INTEGER :: i, ndim, create_field_subarray

    n_local(1) = n1
    ndim = 1
    IF (PRESENT(n2)) THEN
      n_local(2) = n2
      ndim = 2
    END IF
    IF (PRESENT(n3)) THEN
      n_local(3) = n3
      ndim = 3
    END IF

    DO i = 1, ndim
      start(i) = 1 + ng
      n_global(i) = n_local(i) + 2 * ng
    END DO

    IF (PRESENT(n3)) THEN
      create_field_subarray = &
          create_3d_array_subtype(basetype, n_local, n_global, start)
    ELSE IF (PRESENT(n2)) THEN
      create_field_subarray = &
          create_2d_array_subtype(basetype, n_local, n_global, start)
    ELSE
      create_field_subarray = &
          create_1d_array_subtype(basetype, n_local, n_global, start)
    END IF

  END FUNCTION create_field_subarray



  !----------------------------------------------------------------------------
  ! Clips range of current subset to domain size and cell edge
  !----------------------------------------------------------------------------

  FUNCTION global_ranges(current_subset)

    REAL(NUM), DIMENSION(2,c_ndims) :: global_ranges
    TYPE(subset), INTENT(IN), POINTER :: current_subset
    REAL(num) :: dir_min, dir_max, dir_d
    ! fudge factor allows overshoot of the specified domain extent by about 5%
    REAL(num), PARAMETER :: fudge = 0.019_num
    INTEGER :: idim, n

    global_ranges(1,:) = -HUGE(1.0_num)
    global_ranges(2,:) = HUGE(1.0_num)

    DO idim = 1, c_ndims
      IF (idim == 1) THEN
        dir_d = dx
        dir_min = x_min
        dir_max = x_max
        n = c_subset_x_min
        IF (current_subset%use_restriction(n)) &
            global_ranges(1,idim) = current_subset%restriction(n)
        n = c_subset_x_max
        IF (current_subset%use_restriction(n)) &
            global_ranges(2,idim) = current_subset%restriction(n)
      ELSE
        dir_d = dy
        dir_min = y_min
        dir_max = y_max
        IF (current_subset%use_restriction(n)) &
            global_ranges(1,idim) = current_subset%restriction(n)
        n = c_subset_y_max
        IF (current_subset%use_restriction(n)) &
            global_ranges(2,idim) = current_subset%restriction(n)
      END IF

      IF (global_ranges(2,idim) < global_ranges(1,idim)) THEN
        global_ranges = 0
        RETURN
      END IF

      ! Correct to match cell edges
      global_ranges(1,idim) = dir_min &
          + FLOOR((global_ranges(1,idim) - dir_min) / dir_d + fudge) * dir_d
      global_ranges(2,idim) = dir_min &
          + CEILING((global_ranges(2,idim) - dir_min) / dir_d - fudge) * dir_d

      ! Correct to domain size
      global_ranges(1,idim) = MAX(global_ranges(1,idim), dir_min) &
          + 0.5_num * dir_d
      global_ranges(2,idim) = MIN(global_ranges(2,idim), dir_max) &
          + 0.5_num * dir_d
    END DO

  END FUNCTION global_ranges



  FUNCTION cell_global_ranges(ranges)

    INTEGER, DIMENSION(2,c_ndims) :: cell_global_ranges
    REAL(NUM), DIMENSION(2,c_ndims) :: ranges
    REAL(NUM) :: dir_d, lower_posn
    INTEGER :: idim

    DO idim = 1, c_ndims
      IF (idim == 1) THEN
        dir_d = dx
        lower_posn = x_grid_min
      ELSE IF (idim == 2) THEN
        dir_d = dy
        lower_posn = y_grid_min
      END IF
      cell_global_ranges(1,idim) = &
          NINT((ranges(1,idim) - lower_posn) / dir_d) + 1
      cell_global_ranges(2,idim) = &
          NINT((ranges(2,idim) - lower_posn) / dir_d) + 1
    END DO

  END FUNCTION cell_global_ranges



  !----------------------------------------------------------------------------
  ! Location of current processors section of global array
  !----------------------------------------------------------------------------

  FUNCTION cell_local_ranges(ranges)

    INTEGER, DIMENSION(2,c_ndims) :: cell_local_ranges
    REAL(NUM), DIMENSION(2,c_ndims) :: ranges
    REAL(NUM) :: dir_d, lower_posn
    INTEGER :: idim

    ranges(1,1) = MAX(ranges(1,1), x_min_local)
    ranges(2,1) = MIN(ranges(2,1), x_max_local)
    ranges(1,2) = MAX(ranges(1,2), y_min_local)
    ranges(2,2) = MIN(ranges(2,2), y_max_local)

    DO idim = 1, c_ndims
      IF (idim == 1) THEN
        dir_d = dx
        lower_posn = x_grid_min
      ELSE IF (idim == 2) THEN
        dir_d = dy
        lower_posn = y_grid_min
      END IF
      cell_local_ranges(1,idim) = &
          NINT((ranges(1,idim) - lower_posn) / dir_d) + 1
      cell_local_ranges(2,idim) = &
          NINT((ranges(2,idim) - lower_posn) / dir_d) + 1
    END DO

  END FUNCTION cell_local_ranges



  !----------------------------------------------------------------------------
  ! Convert from local ranges into global array, to section of local array
  ! to use
  !----------------------------------------------------------------------------

  FUNCTION cell_section_ranges(ranges)

    INTEGER, DIMENSION(2,c_ndims) :: cell_section_ranges
    INTEGER, DIMENSION(2,c_ndims) :: ranges
    INTEGER :: idim, min_val

    DO idim = 1, c_ndims
      IF (idim == 1) THEN
        min_val = nx_global_min
      ELSE IF (idim == 2) THEN
        min_val = ny_global_min
      END IF
      cell_section_ranges(1,idim) = (ranges(1,idim) - min_val)
      cell_section_ranges(2,idim) = (ranges(2,idim) - min_val)
    END DO

  END FUNCTION cell_section_ranges



  FUNCTION cell_starts(ranges, global_ranges)

    INTEGER, DIMENSION(c_ndims) :: cell_starts
    INTEGER, DIMENSION(2, c_ndims) :: ranges
    REAL(NUM), DIMENSION(2,c_ndims) :: global_ranges
    INTEGER :: range_global_min

    range_global_min = NINT((global_ranges(1,1) - x_grid_min) / dx)

    ! -1 because ranges is cell indexed and global_ranges isn't
    cell_starts(1) = ranges(1,1) - range_global_min - 1

    range_global_min = NINT((global_ranges(1,2) - y_grid_min) / dy)

    cell_starts(2) = ranges(1,2) - range_global_min - 1

  END FUNCTION cell_starts



  ! Check if two 1-dimensional arrays are equal
  PURE LOGICAL FUNCTION equal(array1, array2)

    INTEGER, DIMENSION(:), INTENT(IN) :: array1, array2
    INTEGER :: i

    equal = SIZE(array1) == SIZE(array2)

    IF (equal) THEN
      DO i = 1, SIZE(array1)
        equal = array1(i) == array2(i)
        IF (.NOT. equal) RETURN
      END DO
    END IF

  END FUNCTION equal

END MODULE mpi_subtype_control
