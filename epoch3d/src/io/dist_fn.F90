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

MODULE dist_fn

  USE mpi_subtype_control
  USE particles, ONLY: f0
  USE sdf

  IMPLICIT NONE

CONTAINS

  SUBROUTINE attach_dist_fn(iblock)

    TYPE(distribution_function_block), POINTER :: iblock
    TYPE(distribution_function_block), POINTER :: current

    current => dist_fns
    IF (.NOT. ASSOCIATED(current)) THEN
      ! This is the first distribution function to add
      dist_fns => iblock
      RETURN
    END IF
    DO WHILE(ASSOCIATED(current%next))
      current => current%next
    END DO
    current%next => iblock

  END SUBROUTINE attach_dist_fn



  SUBROUTINE init_dist_fn(iblock)

    TYPE(distribution_function_block), POINTER :: iblock

    iblock%name = blank
    iblock%ndims = -1
    iblock%dumpmask = c_io_always
    iblock%directions = 0
    iblock%ranges = 1.0_num
    iblock%resolution = 1
    iblock%restrictions = 0.0_num
    iblock%use_restrictions = .FALSE.
    NULLIFY(iblock%next)
    ALLOCATE(iblock%use_species(n_species))
    iblock%use_species = .FALSE.
    iblock%output_deltaf = .FALSE.

  END SUBROUTINE init_dist_fn



  SUBROUTINE deallocate_dist_fns

    TYPE(distribution_function_block), POINTER :: current, next
    INTEGER :: stat

    current => dist_fns
    DO WHILE(ASSOCIATED(current))
      next => current%next
      DEALLOCATE(current%use_species, STAT=stat)
      DEALLOCATE(current)
      current => next
    END DO

  END SUBROUTINE deallocate_dist_fns



  SUBROUTINE write_dist_fns(sdf_handle, code, mask)

    TYPE(sdf_file_handle) :: sdf_handle
    INTEGER, INTENT(IN) :: code, mask

    INTEGER :: ispecies, errcode
    TYPE(distribution_function_block), POINTER :: current
    LOGICAL :: convert

    ! Write the distribution functions
    current => dist_fns
    DO WHILE(ASSOCIATED(current))
      IF (IAND(current%dumpmask, code) /= 0) THEN
        DO ispecies = 1, n_species
          IF (.NOT. current%use_species(ispecies)) CYCLE

          convert = (IAND(IOR(mask,current%dumpmask), c_io_dump_single) /= 0)

          CALL general_dist_fn(sdf_handle, current%name, current%directions, &
              current%ranges, current%resolution, ispecies, &
              current%restrictions, current%use_restrictions, current%ndims, &
              current%output_deltaf, convert, errcode)

          ! If there was an error writing the dist_fn then ignore it in future
          IF (errcode /= 0) current%dumpmask = c_io_never
        END DO
      END IF
      current => current%next
    END DO

  END SUBROUTINE write_dist_fns



  SUBROUTINE general_dist_fn(sdf_handle, name, direction, ranges_in, &
      resolution_in, species, restrictions, use_restrictions, curdims, &
      output_deltaf, convert, errcode)

    TYPE(sdf_file_handle) :: sdf_handle
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, DIMENSION(c_df_maxdims), INTENT(IN) :: direction
    REAL(num), DIMENSION(2,c_df_maxdims), INTENT(IN) :: ranges_in
    INTEGER, DIMENSION(c_df_maxdims), INTENT(IN) :: resolution_in
    INTEGER, INTENT(IN) :: species
    REAL(num), DIMENSION(2,c_df_maxdirs), INTENT(IN) :: restrictions
    LOGICAL, DIMENSION(c_df_maxdirs), INTENT(IN) :: use_restrictions
    INTEGER, INTENT(IN) :: curdims
    LOGICAL, INTENT(IN) :: output_deltaf
    LOGICAL, INTENT(IN) :: convert
    INTEGER, INTENT(OUT) :: errcode

    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: array, array_tmp
    REAL(num), DIMENSION(:), ALLOCATABLE :: grid1, grid2, grid3
    LOGICAL, DIMENSION(c_df_maxdims) :: parallel
    REAL(num), DIMENSION(c_df_maxdims) :: dgrid
    REAL(num) :: current_data, temp_data, theta, p, px, py, pz
    INTEGER :: idim, idir
    LOGICAL, DIMENSION(c_df_maxdims) :: calc_range
    LOGICAL :: calc_ranges

    LOGICAL :: use_x, use_y, use_z, need_reduce
    LOGICAL :: use_xy_angle, use_yz_angle, use_zx_angle
    INTEGER, DIMENSION(c_df_maxdims) :: start_local, global_resolution
    INTEGER, DIMENSION(c_df_maxdims) :: range_global_min
    INTEGER :: color, comm_new
    INTEGER :: new_type, array_type

    REAL(num), DIMENSION(2,c_df_maxdims) :: ranges
    INTEGER, DIMENSION(c_df_maxdims) :: resolution
    INTEGER, DIMENSION(c_df_maxdims) :: cell
    REAL(num) :: part_weight, part_mc, part_mc2, part_u2
    REAL(num) :: gamma_rel, gamma_rel_m1, start
    REAL(num) :: xy_max, yz_max, zx_max
    REAL(num), PARAMETER :: pi2 = 2.0_num * pi
    INTEGER :: rank_local
    INTEGER :: mpireal_new

    TYPE(particle), POINTER :: current, next
    CHARACTER(LEN=string_length) :: var_name
    CHARACTER(LEN=8), DIMENSION(c_df_maxdirs) :: labels, units
    REAL(num), DIMENSION(c_df_maxdirs) :: particle_data
    REAL(num) :: tmpvar
    LOGICAL :: proc_outside_range

    proc_outside_range = .FALSE.
    errcode = 0
    ! Update species count if necessary
    IF (io_list(species)%count_update_step < step) THEN
      CALL MPI_ALLREDUCE(io_list(species)%attached_list%count, &
          io_list(species)%count, 1, MPI_INTEGER8, MPI_SUM, &
          comm, errcode)
      io_list(species)%count_update_step = step
    END IF

    IF (io_list(species)%count < 1) RETURN

    use_x = .FALSE.
    use_y = .FALSE.
    use_z = .FALSE.
    color = 0
    ranges = ranges_in
    resolution = resolution_in
    global_resolution = resolution
    range_global_min = 0
    parallel = .FALSE.
    start_local = 1
    calc_range = .FALSE.
    calc_ranges = .FALSE.
    use_xy_angle = .FALSE.
    use_yz_angle = .FALSE.
    use_zx_angle = .FALSE.

    current_data = 0.0_num
#ifndef PER_PARTICLE_CHARGE_MASS
    part_mc  = io_list(species)%mass * c
    part_mc2 = part_mc * c
#endif
#ifdef PER_SPECIES_WEIGHT
    part_weight = io_list(species)%weight
#endif

    DO idim = 1, curdims
      IF (direction(idim) == c_dir_x) THEN
        use_x = .TRUE.
        IF (ABS(ranges(1,idim) - ranges(2,idim)) <= c_tiny) THEN
          ! If empty range, use whole domain
          ranges(1,idim) = x_grid_min_local - 0.5_num * dx
          ranges(2,idim) = x_grid_max_local + 0.5_num * dx
          global_resolution(idim) = nx_global
          start_local(idim) = nx_global_min
        ELSE
          ! Else use the range including the requested range, but ending
          ! on a cell boundary
          ranges(1,idim) = MAX(ranges(1,idim), x_min)
          ranges(2,idim) = MIN(ranges(2,idim), x_max)
          ranges(1,idim) = x_min_local &
              + FLOOR((ranges(1,idim) - x_min_local) / dx) * dx
          ranges(2,idim) = x_min_local &
              + CEILING((ranges(2,idim) - x_min_local) / dx) * dx
          global_resolution(idim) = NINT((ranges(2,idim) - ranges(1,idim)) / dx)
          range_global_min(idim) = NINT((ranges(1,idim) - x_min) / dx)

          ranges(1,idim) = MAX(ranges(1,idim), x_grid_min_local - 0.5_num * dx)
          ranges(2,idim) = MIN(ranges(2,idim), x_grid_max_local + 0.5_num * dx)

          start_local(idim) = nx_global_min &
              + NINT((ranges(1,idim) - x_min_local) / dx) &
              - range_global_min(idim)
        END IF

        ! resolution is the number of pts
        ! ranges guaranteed to include integer number of grid cells
        resolution(idim) = NINT((ranges(2,idim) - ranges(1,idim)) / dx)
        IF (resolution(idim) <= 0) THEN
          proc_outside_range = .TRUE.
          resolution(idim) = 0
        END IF

        dgrid(idim) = dx
        labels(idim) = 'X'
        units(idim)  = 'm'
        parallel(idim) = .TRUE.
        CYCLE

      ELSE IF (direction(idim) == c_dir_y) THEN
        use_y = .TRUE.
        IF (ABS(ranges(1,idim) - ranges(2,idim)) <= c_tiny) THEN
          ! If empty range, use whole domain
          ranges(1,idim) = y_grid_min_local - 0.5_num * dy
          ranges(2,idim) = y_grid_max_local + 0.5_num * dy
          global_resolution(idim) = ny_global
          start_local(idim) = ny_global_min
        ELSE
          ! Else use the range including the requested range, but ending
          ! on a cell boundary
          ranges(1,idim) = MAX(ranges(1,idim), y_min)
          ranges(2,idim) = MIN(ranges(2,idim), y_max)
          ranges(1,idim) = y_min_local &
              + FLOOR((ranges(1,idim) - y_min_local) / dy) * dy
          ranges(2,idim) = y_min_local &
              + CEILING((ranges(2,idim) - y_min_local) / dy) * dy
          global_resolution(idim) = NINT((ranges(2,idim) - ranges(1,idim)) / dy)
          range_global_min(idim) = NINT((ranges(1,idim) - y_min) / dy)

          ranges(1,idim) = MAX(ranges(1,idim), y_grid_min_local - 0.5_num * dy)
          ranges(2,idim) = MIN(ranges(2,idim), y_grid_max_local + 0.5_num * dy)

          start_local(idim) = ny_global_min &
              + NINT((ranges(1,idim) - y_min_local) / dy) &
              - range_global_min(idim)
        END IF

        ! resolution is the number of pts
        ! ranges guaranteed to include integer number of grid cells
        resolution(idim) = NINT((ranges(2,idim) - ranges(1,idim)) / dy)
        IF (resolution(idim) <= 0) THEN
          proc_outside_range = .TRUE.
          resolution(idim) = 0
        END IF

        dgrid(idim) = dy
        labels(idim) = 'Y'
        units(idim)  = 'm'
        parallel(idim) = .TRUE.
        CYCLE

      ELSE IF (direction(idim) == c_dir_z) THEN
        use_z = .TRUE.
        IF (ABS(ranges(1,idim) - ranges(2,idim)) <= c_tiny) THEN
          ! If empty range, use whole domain
          ranges(1,idim) = z_grid_min_local - 0.5_num * dz
          ranges(2,idim) = z_grid_max_local + 0.5_num * dz
          global_resolution(idim) = nz_global
          start_local(idim) = nz_global_min
        ELSE
          ! Else use the range including the requested range, but ending
          ! on a cell boundary
          ranges(1,idim) = MAX(ranges(1,idim), z_min)
          ranges(2,idim) = MIN(ranges(2,idim), z_max)
          ranges(1,idim) = z_min_local &
              + FLOOR((ranges(1,idim) - z_min_local) / dz) * dz
          ranges(2,idim) = z_min_local &
              + CEILING((ranges(2,idim) - z_min_local) / dz) * dz
          global_resolution(idim) = NINT((ranges(2,idim) - ranges(1,idim)) / dz)
          range_global_min(idim) = NINT((ranges(1,idim) - z_min) / dz)

          ranges(1,idim) = MAX(ranges(1,idim), z_grid_min_local - 0.5_num * dz)
          ranges(2,idim) = MIN(ranges(2,idim), z_grid_max_local + 0.5_num * dz)

          start_local(idim) = nz_global_min &
              + NINT((ranges(1,idim) - z_min_local) / dz) &
              - range_global_min(idim)
        END IF

        ! resolution is the number of pts
        ! ranges guaranteed to include integer number of grid cells
        resolution(idim) = NINT((ranges(2,idim) - ranges(1,idim)) / dz)
        IF (resolution(idim) <= 0) THEN
          proc_outside_range = .TRUE.
          resolution(idim) = 0
        END IF

        dgrid(idim) = dz
        labels(idim) = 'Z'
        units(idim)  = 'm'
        parallel(idim) = .TRUE.
        CYCLE

      END IF

      ! If we're here then this must be a momentum space direction
      ! So determine which momentum space directions are needed
      IF (ABS(ranges(1,idim) - ranges(2,idim)) <= c_tiny) THEN
        calc_range(idim) = .TRUE.
        calc_ranges = .TRUE.
      END IF

      IF (direction(idim) == c_dir_px) THEN
        labels(idim) = 'Px'
        units(idim)  = 'kg.m/s'

      ELSE IF (direction(idim) == c_dir_py) THEN
        labels(idim) = 'Py'
        units(idim)  = 'kg.m/s'

      ELSE IF (direction(idim) == c_dir_pz) THEN
        labels(idim) = 'Pz'
        units(idim)  = 'kg.m/s'

      ELSE IF (direction(idim) == c_dir_mod_p) THEN
        labels(idim) = '|P|'
        units(idim)  = 'kg.m/s'

      ELSE IF (direction(idim) == c_dir_en) THEN
        labels(idim) = 'energy'
        units(idim)  = 'J'

      ELSE IF (direction(idim) == c_dir_gamma_m1) THEN
        labels(idim) = 'gamma-1'
        units(idim)  = ''

      ELSE IF (direction(idim) == c_dir_xy_angle) THEN
        use_xy_angle = .TRUE.
        labels(idim) = 'xy_angle'
        units(idim)  = 'radians'

      ELSE IF (direction(idim) == c_dir_yz_angle) THEN
        use_yz_angle = .TRUE.
        labels(idim) = 'yz_angle'
        units(idim)  = 'radians'

      ELSE IF (direction(idim) == c_dir_zx_angle) THEN
        use_zx_angle = .TRUE.
        labels(idim) = 'zx_angle'
        units(idim)  = 'radians'

      ELSE
        IF (rank == 0) THEN
          WRITE(*,*) '*** WARNING ***'
          WRITE(*,*) 'Unable to write dist_fn. Ignoring.'
        END IF
        errcode = 1
        RETURN

      END IF
    END DO

    ! Calculate range for directions where needed
    IF (calc_ranges) THEN
      DO idim = 1, curdims
        IF (calc_range(idim)) THEN
          ranges(1,idim) =  HUGE(1.0_num) * 0.4_num
          ranges(2,idim) = -HUGE(1.0_num) * 0.4_num
        END IF
      END DO
      next => io_list(species)%attached_list%head

      out1: DO WHILE(ASSOCIATED(next))
        current => next
        next => current%next

#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc  = current%mass * c
        part_mc2 = part_mc * c
#endif
        part_u2 = SUM((current%part_p / part_mc)**2)
        gamma_rel = SQRT(part_u2 + 1.0_num)
        gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)
        px = current%part_p(1)
        py = current%part_p(2)
        pz = current%part_p(3)

        particle_data(1:c_ndims) = current%part_pos
        IF (io_list(species)%species_type == c_species_id_photon) THEN
          particle_data(c_dir_gamma_m1) = 0.0_num
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
          particle_data(c_dir_px) = px
          particle_data(c_dir_py) = py
          particle_data(c_dir_pz) = pz
          particle_data(c_dir_en) = current%particle_energy
          particle_data(c_dir_mod_p) = SQRT(px**2 + py**2 + pz**2)
#else
          particle_data(c_dir_px:c_dir_pz) = 0.0_num
          particle_data(c_dir_en) = 0.0_num
          particle_data(c_dir_mod_p) = 0.0_num
#endif
          ! Can't define gamma for photon so one is as good as anything
          particle_data(c_dir_gamma_m1) = 1.0_num
        ELSE
          particle_data(c_dir_px) = px
          particle_data(c_dir_py) = py
          particle_data(c_dir_pz) = pz
          particle_data(c_dir_en) = gamma_rel_m1 * part_mc2
          particle_data(c_dir_gamma_m1) = gamma_rel_m1
          particle_data(c_dir_mod_p) = SQRT(px**2 + py**2 + pz**2)
        END IF

        IF (use_xy_angle) THEN
          p = SQRT(px**2 + py**2)
          IF (p < c_tiny) CYCLE

          theta = ASIN(py / p)
          IF (px < 0.0_num) theta = SIGN(1.0_num,py) * pi - theta
          particle_data(c_dir_xy_angle) = theta
        END IF

        IF (use_yz_angle) THEN
          p = SQRT(py**2 + pz**2)
          IF (p < c_tiny) CYCLE

          theta = ASIN(pz / p)
          IF (py < 0.0_num) theta = SIGN(1.0_num,pz) * pi - theta
          particle_data(c_dir_yz_angle) = theta
        END IF

        IF (use_zx_angle) THEN
          p = SQRT(pz**2 + px**2)
          IF (p < c_tiny) CYCLE

          theta = ASIN(px / p)
          IF (pz < 0.0_num) theta = SIGN(1.0_num,px) * pi - theta
          particle_data(c_dir_zx_angle) = theta
        END IF

        DO idim = 1, curdims
          IF (calc_range(idim)) THEN
            DO idir = 1, c_df_maxdirs
              IF (use_restrictions(idir) &
                  .AND. (particle_data(idir) < restrictions(1,idir) &
                  .OR. particle_data(idir) > restrictions(2,idir))) &
                      CYCLE out1
            END DO

            current_data = particle_data(direction(idim))
            IF (current_data < ranges(1,idim)) ranges(1,idim) = current_data
            IF (current_data > ranges(2,idim)) ranges(2,idim) = current_data
          END IF
        END DO
      END DO out1

      DO idim = 1, curdims
        IF (.NOT. parallel(idim)) THEN
          ! If not parallel then this is a momentum dimension
          CALL MPI_ALLREDUCE(ranges(1,idim), temp_data, 1, mpireal, MPI_MIN, &
              comm, errcode)
          ranges(1,idim) = temp_data
          CALL MPI_ALLREDUCE(ranges(2,idim), temp_data, 1, mpireal, MPI_MAX, &
              comm, errcode)
          ranges(2,idim) = temp_data
        END IF
      END DO
    END IF

    DO idim = 1, curdims
      ! Fix so that if distribution function is zero then it picks an arbitrary
      ! scale in that direction
      IF (ABS(ranges(1,idim) - ranges(2,idim)) <= c_tiny) THEN
        ranges(1,idim) = -1.0_num
        ranges(2,idim) = 1.0_num
      END IF

      IF (direction(idim) == c_dir_xy_angle) THEN
        xy_max = ranges(2,idim)
      ELSE IF (direction(idim) == c_dir_yz_angle) THEN
        yz_max = ranges(2,idim)
      ELSE IF (direction(idim) == c_dir_zx_angle) THEN
        zx_max = ranges(2,idim)
      END IF

      ! Calculate grid spacing
      IF (.NOT. parallel(idim)) dgrid(idim) = &
          (ranges(2,idim) - ranges(1,idim)) / REAL(resolution(idim), num)
    END DO

    IF (.NOT. proc_outside_range) THEN
      ALLOCATE(array(resolution(1), resolution(2), resolution(3)))
    ELSE
      ! Dummy array
      ALLOCATE(array(1,1,1))
    END IF
    array = 0.0_num

    next => io_list(species)%attached_list%head

    out2: DO WHILE(ASSOCIATED(next))
      current => next
      next => current%next

#ifdef PER_PARTICLE_CHARGE_MASS
      part_mc  = current%mass * c
      part_mc2 = part_mc * c
#endif
#ifndef PER_SPECIES_WEIGHT
      part_weight = current%weight
#endif
#ifdef DELTAF_METHOD
      IF (output_deltaf) THEN
         part_weight = current%weight &
             - current%pvol * f0(species, part_mc / c, current%part_p)
      END IF
#endif
      part_u2 = SUM((current%part_p / part_mc)**2)
      gamma_rel = SQRT(part_u2 + 1.0_num)
      gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)
      px = current%part_p(1)
      py = current%part_p(2)
      pz = current%part_p(3)

      particle_data(1:c_ndims) = current%part_pos
      IF (io_list(species)%species_type == c_species_id_photon) THEN
        particle_data(c_dir_gamma_m1) = 0.0_num
#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
        particle_data(c_dir_px) = px
        particle_data(c_dir_py) = py
        particle_data(c_dir_pz) = pz
        particle_data(c_dir_en) = current%particle_energy
        particle_data(c_dir_mod_p) = SQRT(px**2 + py**2 + pz**2)
#else
        particle_data(c_dir_px:c_dir_pz) = 0.0_num
        particle_data(c_dir_en) = 0.0_num
        particle_data(c_dir_mod_p) = 0.0_num
#endif
        ! Can't define gamma for photon so one is as good as anything
        particle_data(c_dir_gamma_m1) = 1.0_num
      ELSE
        particle_data(c_dir_px) = px
        particle_data(c_dir_py) = py
        particle_data(c_dir_pz) = pz
        particle_data(c_dir_en) = gamma_rel_m1 * part_mc2
        particle_data(c_dir_gamma_m1) = gamma_rel_m1
        particle_data(c_dir_mod_p) = SQRT(px**2 + py**2 + pz**2)
      END IF

      IF (use_xy_angle) THEN
        p = SQRT(px**2 + py**2)
        IF (p < c_tiny) CYCLE

        theta = ASIN(py / p)
        IF (px < 0.0_num) theta = SIGN(1.0_num,py) * pi - theta
        IF ((theta + pi2) < xy_max) theta = theta + pi2
        particle_data(c_dir_xy_angle) = theta
      END IF

      IF (use_yz_angle) THEN
        p = SQRT(py**2 + pz**2)
        IF (p < c_tiny) CYCLE

        theta = ASIN(pz / p)
        IF (py < 0.0_num) theta = SIGN(1.0_num,pz) * pi - theta
        IF ((theta + pi2) < yz_max) theta = theta + pi2
        particle_data(c_dir_yz_angle) = theta
      END IF

      IF (use_zx_angle) THEN
        p = SQRT(pz**2 + px**2)
        IF (p < c_tiny) CYCLE

        theta = ASIN(px / p)
        IF (pz < 0.0_num) theta = SIGN(1.0_num,px) * pi - theta
        IF ((theta + pi2) < zx_max) theta = theta + pi2
        particle_data(c_dir_zx_angle) = theta
      END IF

      DO idir = 1, c_df_maxdirs
        IF (use_restrictions(idir) &
            .AND. (particle_data(idir) < restrictions(1,idir) &
            .OR. particle_data(idir) > restrictions(2,idir))) &
                CYCLE out2
      END DO

      cell = 1
      DO idim = 1, curdims
        current_data = particle_data(direction(idim))
        tmpvar = (current_data - ranges(1,idim)) / dgrid(idim)
        IF (ABS(tmpvar) > REAL(HUGE(1),num)) &
            CYCLE out2
        cell(idim) = FLOOR(tmpvar) + 1
        IF (cell(idim) < 1 .OR. cell(idim) > resolution(idim)) &
            CYCLE out2
      END DO
      IF (.NOT. proc_outside_range) array(cell(1), cell(2), cell(3)) = &
          array(cell(1), cell(2), cell(3)) + part_weight ! * real_space_area
    END DO out2

    need_reduce = .TRUE.
    IF (use_x .AND. use_y .AND. use_z) need_reduce = .FALSE.

    IF (need_reduce) THEN
      ! If using x direction need to reduce across y, z
      IF (use_x) color = color + x_coords
      ! If using y direction need to reduce across x, z
      IF (use_y) color = color + nprocx * y_coords
      ! If using z direction need to reduce across x, y
      IF (use_z) color = color + nprocx * nprocy * z_coords

      CALL MPI_COMM_SPLIT(comm, color, rank, comm_new, errcode)
      IF (.NOT. proc_outside_range) THEN
        ALLOCATE(array_tmp(resolution(1), resolution(2), resolution(3)))
      ELSE
        ALLOCATE(array_tmp(1,1,1))
      END IF
      array_tmp = 0.0_num
      CALL MPI_REDUCE(array, array_tmp, &
          resolution(1)*resolution(2)*resolution(3), mpireal, MPI_SUM, &
          0, comm_new, errcode)
      CALL MPI_COMM_RANK(comm_new, rank_local, errcode)
      array = array_tmp
      DEALLOCATE(array_tmp)
      CALL MPI_COMM_FREE(comm_new, errcode)
    ELSE
      rank_local = 0
    END IF

    ! Create grids
    ALLOCATE(grid1(global_resolution(1)))
    start = ranges(1,1) + 0.5_num * dgrid(1)
    DO idir = 1, global_resolution(1)
      grid1(idir) = start + (idir - 1) * dgrid(1)
    END DO

    IF (curdims >= 2) THEN
      ALLOCATE(grid2(global_resolution(2)))
      start = ranges(1,2) + 0.5_num * dgrid(2)
      DO idir = 1, global_resolution(2)
        grid2(idir) = start + (idir - 1) * dgrid(2)
      END DO
    END IF

    IF (curdims >= 3) THEN
      ALLOCATE(grid3(global_resolution(3)))
      start = ranges(1,3) + 0.5_num * dgrid(3)
      DO idir = 1, global_resolution(3)
        grid3(idir) = start + (idir - 1) * dgrid(3)
      END DO
    END IF

    var_name = TRIM(name) // '/' // TRIM(io_list(species)%name)

    IF (use_offset_grid) THEN
      IF (curdims == 1) THEN
        CALL sdf_write_srl_plain_mesh(sdf_handle, &
            'grid_full/' // TRIM(var_name), 'Grid_Full/' // TRIM(var_name), &
            grid1, convert, labels, units)
      ELSE IF (curdims == 2) THEN
        CALL sdf_write_srl_plain_mesh(sdf_handle, &
            'grid_full/' // TRIM(var_name), 'Grid_Full/' // TRIM(var_name), &
            grid1, grid2, convert, labels, units)
      ELSE IF (curdims == 3) THEN
        CALL sdf_write_srl_plain_mesh(sdf_handle, &
            'grid_full/' // TRIM(var_name), 'Grid_Full/' // TRIM(var_name), &
            grid1, grid2, grid3, convert, labels, units)
      END IF
      IF (parallel(1)) grid1 = grid1 - ranges(1,1)
      IF (parallel(2)) grid2 = grid2 - ranges(1,2)
      IF (parallel(3)) grid3 = grid3 - ranges(1,3)
    END IF

    IF (convert) THEN
      mpireal_new = MPI_REAL4
    ELSE
      mpireal_new = mpireal
    END IF

    IF (curdims == 1) THEN
      CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid/' // TRIM(var_name), &
          'Grid/' // TRIM(var_name), grid1, convert, labels, units)
      DEALLOCATE(grid1)
      new_type = create_1d_array_subtype(mpireal_new, resolution, &
          global_resolution, start_local)
    ELSE IF (curdims == 2) THEN
      CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid/' // TRIM(var_name), &
          'Grid/' // TRIM(var_name), grid1, grid2, convert, labels, units)
      DEALLOCATE(grid1, grid2)
      new_type = create_2d_array_subtype(mpireal_new, resolution, &
          global_resolution, start_local)
    ELSE IF (curdims == 3) THEN
      CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid/' // TRIM(var_name), &
          'Grid/' // TRIM(var_name), grid1, grid2, grid3, convert, labels, &
          units)
      DEALLOCATE(grid1, grid2, grid3)
      new_type = create_3d_array_subtype(mpireal_new, resolution, &
          global_resolution, start_local)
    END IF

    IF (rank_local == 0) THEN
      CALL MPI_TYPE_CONTIGUOUS(resolution(1) * resolution(2) * resolution(3), &
          mpireal_new, array_type, errcode)
      CALL MPI_TYPE_COMMIT(array_type, errcode)
    ELSE
      CALL MPI_TYPE_FREE(new_type, errcode)
      CALL MPI_TYPE_CONTIGUOUS(0, mpireal, new_type, errcode)
      CALL MPI_TYPE_COMMIT(new_type, errcode)
      CALL MPI_TYPE_CONTIGUOUS(0, mpireal, array_type, errcode)
      CALL MPI_TYPE_COMMIT(array_type, errcode)
    END IF

    IF (curdims == 1) THEN
      CALL sdf_write_plain_variable(sdf_handle, TRIM(var_name), &
          'dist_fn/' // TRIM(var_name), 'npart/cell', global_resolution(1:1), &
          c_stagger_vertex, 'grid/' // TRIM(var_name), &
          array(1:resolution(1),1,1), new_type, &
          array_type, convert)
    ELSE IF (curdims == 2) THEN
      CALL sdf_write_plain_variable(sdf_handle, TRIM(var_name), &
          'dist_fn/' // TRIM(var_name), 'npart/cell', global_resolution(1:2), &
          c_stagger_vertex, 'grid/' // TRIM(var_name), &
          array(1:resolution(1),1:resolution(2),1), new_type, &
          array_type, convert)
    ELSE IF (curdims == 3) THEN
      CALL sdf_write_plain_variable(sdf_handle, TRIM(var_name), &
          'dist_fn/' // TRIM(var_name), 'npart/cell', global_resolution, &
          c_stagger_vertex, 'grid/' // TRIM(var_name), array, new_type, &
          array_type, convert)
    END IF

    CALL MPI_TYPE_FREE(new_type, errcode)
    CALL MPI_TYPE_FREE(array_type, errcode)

    DEALLOCATE(array)

  END SUBROUTINE general_dist_fn

END MODULE dist_fn
