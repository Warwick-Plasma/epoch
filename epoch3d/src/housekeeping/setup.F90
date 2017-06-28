! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
! Copyright (C) 2012      Martin Ramsay <M.G.Ramsay@warwick.ac.uk>
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

MODULE setup

  USE fields
  USE mpi_subtype_control
  USE version_data
  USE welcome
  USE split_particle
  USE shunt
  USE laser
  USE window
  USE timer
  USE helper
  USE sdf

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: after_control, minimal_init, restart_data
  PUBLIC :: open_files, close_files, flush_stat_file
  PUBLIC :: setup_species, after_deck_last, set_dt
  PUBLIC :: read_cpu_split

  TYPE(particle), POINTER, SAVE :: iterator_list
#ifndef NO_IO
  CHARACTER(LEN=11+data_dir_max_length), SAVE :: stat_file
#endif

CONTAINS

  SUBROUTINE minimal_init

    INTEGER, PARAMETER :: r4  = SELECTED_REAL_KIND(r=30)
    INTEGER, PARAMETER :: r8  = SELECTED_REAL_KIND(r=300)
    INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(r=3000)

    IF (num == r4) THEN
      ! Should use MPI_SIZEOF() but this breaks on scalimpi
      realsize = 4
      mpireal = MPI_REAL4
    ELSE IF (num == r8) THEN
      realsize = 8
      mpireal = MPI_REAL8
    ELSE
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Cannot determine size of real'
      ENDIF
      CALL abort_code(c_err_terminate)
      STOP
    ENDIF

    dt_plasma_frequency = 0.0_num
    dt_multiplier = 0.95_num
    stdout_frequency = 0
    cpml_thickness = 6
    cpml_kappa_max = 20.0_num
    cpml_a_max = 0.15_num
    cpml_sigma_max = 0.7_num
    cpml_x_min_offset = 0
    cpml_x_max_offset = 0
    cpml_y_min_offset = 0
    cpml_y_max_offset = 0
    cpml_z_min_offset = 0
    cpml_z_max_offset = 0

    window_shift = 0.0_num
    npart_global = -1
    smooth_currents = .FALSE.
    use_balance = .FALSE.
    use_random_seed = .FALSE.
    use_offset_grid = .FALSE.
    use_particle_lists = .FALSE.
    use_multiphoton = .TRUE.
    use_bsi = .TRUE.
    need_random_state = .FALSE.
    force_first_to_be_restartable = .FALSE.
    force_final_to_be_restartable = .FALSE.
    full_dump_every = -1
    restart_dump_every = -1
    nsteps = -1
    t_end = HUGE(1.0_num)
    particles_max_id = 0

    NULLIFY(laser_x_min)
    NULLIFY(laser_x_max)
    NULLIFY(laser_y_max)
    NULLIFY(laser_y_min)
    NULLIFY(laser_z_max)
    NULLIFY(laser_z_min)

    NULLIFY(dist_fns)
    NULLIFY(io_block_list)

    run_date = get_unix_time()

    CALL set_field_order(2)

    CALL timer_init

    ! This array is true if a field component is staggered in the
    ! given direction.
    stagger = .FALSE.
    stagger(c_dir_x,c_stagger_ex) = .TRUE.
    stagger(c_dir_y,c_stagger_ey) = .TRUE.
    stagger(c_dir_z,c_stagger_ez) = .TRUE.

    stagger(c_dir_x:c_dir_z,c_stagger_bx) = .TRUE.
    stagger(c_dir_x:c_dir_z,c_stagger_by) = .TRUE.
    stagger(c_dir_x:c_dir_z,c_stagger_bz) = .TRUE.
    stagger(c_dir_x,c_stagger_bx) = .FALSE.
    stagger(c_dir_y,c_stagger_by) = .FALSE.
    stagger(c_dir_z,c_stagger_bz) = .FALSE.

    CALL set_tokenizer_stagger(c_stagger_centre)

    ALLOCATE(x(1), y(1), z(1))
    x = 0.0_num
    y = 0.0_num
    z = 0.0_num

    ALLOCATE(xb(1), yb(1), zb(1))
    xb = 0.0_num
    yb = 0.0_num
    zb = 0.0_num

    CALL eval_stack_init

  END SUBROUTINE minimal_init



  SUBROUTINE after_control

    CALL setup_grid
    CALL set_initial_values

  END SUBROUTINE after_control



  SUBROUTINE setup_grid

    INTEGER :: iproc, ix, iy, iz
    REAL(num) :: xb_min, yb_min, zb_min

    length_x = x_max - x_min
    dx = length_x / REAL(nx_global-2*cpml_thickness, num)
    x_grid_min = x_min - dx * cpml_thickness
    x_grid_max = x_max + dx * cpml_thickness

    length_y = y_max - y_min
    dy = length_y / REAL(ny_global-2*cpml_thickness, num)
    y_grid_min = y_min - dy * cpml_thickness
    y_grid_max = y_max + dy * cpml_thickness

    length_z = z_max - z_min
    dz = length_z / REAL(nz_global-2*cpml_thickness, num)
    z_grid_min = z_min - dz * cpml_thickness
    z_grid_max = z_max + dz * cpml_thickness

    ! Shift grid to cell centres.
    ! At some point the grid may be redefined to be node centred.

    xb_min = x_grid_min
    yb_min = y_grid_min
    zb_min = z_grid_min
    x_grid_min = x_grid_min + dx / 2.0_num
    x_grid_max = x_grid_max - dx / 2.0_num
    y_grid_min = y_grid_min + dy / 2.0_num
    y_grid_max = y_grid_max - dy / 2.0_num
    z_grid_min = z_grid_min + dz / 2.0_num
    z_grid_max = z_grid_max - dz / 2.0_num

    ! Setup global grid
    DO ix = -2, nx_global + 3
      x_global(ix) = x_grid_min + (ix - 1) * dx
      xb_global(ix) = xb_min + (ix - 1) * dx
      xb_offset_global(ix) = xb_global(ix)
    ENDDO
    DO iy = -2, ny_global + 3
      y_global(iy) = y_grid_min + (iy - 1) * dy
      yb_global(iy) = yb_min + (iy - 1) * dy
      yb_offset_global(iy) = yb_global(iy)
    ENDDO
    DO iz = -2, nz_global + 3
      z_global(iz) = z_grid_min + (iz - 1) * dz
      zb_global(iz) = zb_min + (iz - 1) * dz
      zb_offset_global(iz) = zb_global(iz)
    ENDDO

    DO iproc = 0, nprocx-1
      x_grid_mins(iproc) = x_global(cell_x_min(iproc+1))
      x_grid_maxs(iproc) = x_global(cell_x_max(iproc+1))
    ENDDO
    DO iproc = 0, nprocy-1
      y_grid_mins(iproc) = y_global(cell_y_min(iproc+1))
      y_grid_maxs(iproc) = y_global(cell_y_max(iproc+1))
    ENDDO
    DO iproc = 0, nprocz-1
      z_grid_mins(iproc) = z_global(cell_z_min(iproc+1))
      z_grid_maxs(iproc) = z_global(cell_z_max(iproc+1))
    ENDDO

    x_grid_min_local = x_grid_mins(x_coords)
    x_grid_max_local = x_grid_maxs(x_coords)
    y_grid_min_local = y_grid_mins(y_coords)
    y_grid_max_local = y_grid_maxs(y_coords)
    z_grid_min_local = z_grid_mins(z_coords)
    z_grid_max_local = z_grid_maxs(z_coords)

    x_min_local = x_grid_min_local + (cpml_x_min_offset - 0.5_num) * dx
    x_max_local = x_grid_max_local - (cpml_x_max_offset - 0.5_num) * dx
    y_min_local = y_grid_min_local + (cpml_y_min_offset - 0.5_num) * dy
    y_max_local = y_grid_max_local - (cpml_y_max_offset - 0.5_num) * dy
    z_min_local = z_grid_min_local + (cpml_z_min_offset - 0.5_num) * dz
    z_max_local = z_grid_max_local - (cpml_z_max_offset - 0.5_num) * dz

    ! Setup local grid
    x(-2:nx+3) = x_global(nx_global_min-3:nx_global_max+3)
    y(-2:ny+3) = y_global(ny_global_min-3:ny_global_max+3)
    z(-2:nz+3) = z_global(nz_global_min-3:nz_global_max+3)

    xb(-2:nx+3) = xb_global(nx_global_min-3:nx_global_max+3)
    yb(-2:ny+3) = yb_global(ny_global_min-3:ny_global_max+3)
    zb(-2:nz+3) = zb_global(nz_global_min-3:nz_global_max+3)

  END SUBROUTINE setup_grid



  SUBROUTINE after_deck_last

    INTEGER :: i

    CALL setup_data_averaging
    CALL setup_split_particles
    CALL setup_field_boundaries

    IF (cpml_boundaries) THEN
      CALL allocate_cpml_fields
      CALL set_cpml_helpers(nx, nx_global_min, nx_global_max, &
          ny, ny_global_min, ny_global_max, nz, nz_global_min, nz_global_max)
    ELSE
      cpml_thickness = 0
      cpml_kappa_max = 1.0_num
      cpml_a_max = 0.0_num
      cpml_sigma_max = 0.0_num
      DO i = 1, n_io_blocks
        io_block_list(i)%dumpmask(c_dump_cpml_psi_eyx) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_ezx) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_byx) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_bzx) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_exy) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_ezy) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_bxy) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_bzy) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_exz) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_eyz) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_bxz) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_byz) = c_io_never
      ENDDO
    ENDIF

  END SUBROUTINE after_deck_last



  SUBROUTINE setup_data_averaging()

    INTEGER :: io, ib, nspec_local, nstep_average, mask
    REAL(num) :: dt_average
    TYPE(averaged_data_block), POINTER :: avg

    IF (.NOT. any_average) RETURN

    DO io = 1, n_io_blocks
      IF (io_block_list(io)%any_average) THEN
        dt_min_average = t_end
      ELSE
        dt_min_average = -1.0_num
      ENDIF
    ENDDO

    DO io = 1, num_vars_to_dump
      ib = averaged_var_block(io)
      IF (ib /= 0) THEN
        dt_average  = io_block_list(ib)%dt_average
        nstep_average = io_block_list(ib)%nstep_average

        IF (nstep_average > 0 .AND. dt_average > 0) THEN
          io_block_list(ib)%dt_min_average = &
              MIN(io_block_list(ib)%dt_min_average, &
              dt_average / REAL(nstep_average, num))
        ENDIF

        mask = io_block_list(ib)%dumpmask(io)
        nspec_local = 0
        IF (IAND(mask, c_io_no_sum) == 0) &
            nspec_local = 1
        IF (IAND(mask, c_io_species) /= 0) &
            nspec_local = nspec_local + n_species

        IF (nspec_local <= 0) CYCLE

        avg => io_block_list(ib)%averaged_data(io)
        IF (avg%dump_single) THEN
          ALLOCATE(avg%r4array(-2:nx+3,-2:ny+3,-2:nz+3,nspec_local))
          avg%r4array = 0.0_num
        ELSE
          ALLOCATE(avg%array(-2:nx+3,-2:ny+3,-2:nz+3,nspec_local))
          avg%array = 0.0_num
        ENDIF

        avg%species_sum = 0
        avg%n_species = 0
        IF (IAND(mask, c_io_no_sum) == 0) &
            avg%species_sum = 1
        IF (IAND(mask, c_io_species) /= 0) &
            avg%n_species = n_species
        avg%real_time = 0.0_num
        avg%started = .FALSE.
      ENDIF
    ENDDO

  END SUBROUTINE setup_data_averaging



  SUBROUTINE setup_species

    INTEGER :: ispecies, n

    ALLOCATE(species_list(n_species))
    ALLOCATE(io_list_data(n_species))
    io_list => species_list

    DO ispecies = 1, n_species
      species_list(ispecies)%name = blank
      species_list(ispecies)%mass = -1.0_num
      species_list(ispecies)%charge = 0.0_num
      species_list(ispecies)%weight = 1.0_num
      species_list(ispecies)%dumpmask = c_io_always
      species_list(ispecies)%count = -1
      species_list(ispecies)%id = 0
      species_list(ispecies)%npart_per_cell = -1
      species_list(ispecies)%split = .FALSE.
      species_list(ispecies)%npart_max = 0
      species_list(ispecies)%global_count = 0
      species_list(ispecies)%count_update_step = 0
      species_list(ispecies)%species_type = c_species_id_generic
      species_list(ispecies)%immobile = .FALSE.
      NULLIFY(species_list(ispecies)%next)
      NULLIFY(species_list(ispecies)%prev)
      NULLIFY(species_list(ispecies)%ext_temp_x_min)
      NULLIFY(species_list(ispecies)%ext_temp_x_max)
      NULLIFY(species_list(ispecies)%ext_temp_y_min)
      NULLIFY(species_list(ispecies)%ext_temp_y_max)
      NULLIFY(species_list(ispecies)%ext_temp_z_min)
      NULLIFY(species_list(ispecies)%ext_temp_z_max)
      NULLIFY(species_list(ispecies)%secondary_list)
    ENDDO

    DO ispecies = 1, n_species
      CALL initialise_stack(species_list(ispecies)%density_function)
      CALL set_stack_zero  (species_list(ispecies)%density_function)
      DO n = 1, 3
        CALL initialise_stack(species_list(ispecies)%temperature_function(n))
        CALL set_stack_zero  (species_list(ispecies)%temperature_function(n))
        CALL initialise_stack(species_list(ispecies)%drift_function(n))
        CALL set_stack_zero  (species_list(ispecies)%drift_function(n))
      ENDDO
      species_list(ispecies)%electron = .FALSE.
      species_list(ispecies)%ionise = .FALSE.
      species_list(ispecies)%ionise_to_species = -1
      species_list(ispecies)%release_species = -1
      species_list(ispecies)%n = 0
      species_list(ispecies)%l = 0
      species_list(ispecies)%ionisation_energy = HUGE(0.0_num)
      species_list(ispecies)%migrate%this_species = .FALSE.
      species_list(ispecies)%migrate%fluid = .FALSE.
      species_list(ispecies)%migrate%done = .FALSE.
      species_list(ispecies)%migrate%promoteable = .FALSE.
      species_list(ispecies)%migrate%demoteable = .FALSE.
      species_list(ispecies)%migrate%promote_to_species = 0
      species_list(ispecies)%migrate%demote_to_species = 0
      species_list(ispecies)%migrate%promotion_energy_factor = 1.0_num
      species_list(ispecies)%migrate%demotion_energy_factor = 1.0_num
      species_list(ispecies)%migrate%promotion_density = HUGE(1.0_num)
      species_list(ispecies)%migrate%demotion_density = 0.0_num
#ifndef NO_TRACER_PARTICLES
      species_list(ispecies)%tracer = .FALSE.
#endif
#ifndef NO_PARTICLE_PROBES
      NULLIFY(species_list(ispecies)%attached_probes)
#endif
    ENDDO

  END SUBROUTINE setup_species



  SUBROUTINE setup_field_boundaries

    INTEGER :: nx0, nx1, ny0, ny1, nz0, nz1

    ALLOCATE(ex_x_min(-2:ny+3,-2:nz+3), ex_x_max(-2:ny+3,-2:nz+3))
    ALLOCATE(ey_x_min(-2:ny+3,-2:nz+3), ey_x_max(-2:ny+3,-2:nz+3))
    ALLOCATE(ez_x_min(-2:ny+3,-2:nz+3), ez_x_max(-2:ny+3,-2:nz+3))
    ALLOCATE(bx_x_min(-2:ny+3,-2:nz+3), bx_x_max(-2:ny+3,-2:nz+3))
    ALLOCATE(by_x_min(-2:ny+3,-2:nz+3), by_x_max(-2:ny+3,-2:nz+3))
    ALLOCATE(bz_x_min(-2:ny+3,-2:nz+3), bz_x_max(-2:ny+3,-2:nz+3))

    ALLOCATE(ex_y_min(-2:nx+3,-2:nz+3), ex_y_max(-2:nx+3,-2:nz+3))
    ALLOCATE(ey_y_min(-2:nx+3,-2:nz+3), ey_y_max(-2:nx+3,-2:nz+3))
    ALLOCATE(ez_y_min(-2:nx+3,-2:nz+3), ez_y_max(-2:nx+3,-2:nz+3))
    ALLOCATE(bx_y_min(-2:nx+3,-2:nz+3), bx_y_max(-2:nx+3,-2:nz+3))
    ALLOCATE(by_y_min(-2:nx+3,-2:nz+3), by_y_max(-2:nx+3,-2:nz+3))
    ALLOCATE(bz_y_min(-2:nx+3,-2:nz+3), bz_y_max(-2:nx+3,-2:nz+3))

    ALLOCATE(ex_z_min(-2:nx+3,-2:ny+3), ex_z_max(-2:nx+3,-2:ny+3))
    ALLOCATE(ey_z_min(-2:nx+3,-2:ny+3), ey_z_max(-2:nx+3,-2:ny+3))
    ALLOCATE(ez_z_min(-2:nx+3,-2:ny+3), ez_z_max(-2:nx+3,-2:ny+3))
    ALLOCATE(bx_z_min(-2:nx+3,-2:ny+3), bx_z_max(-2:nx+3,-2:ny+3))
    ALLOCATE(by_z_min(-2:nx+3,-2:ny+3), by_z_max(-2:nx+3,-2:ny+3))
    ALLOCATE(bz_z_min(-2:nx+3,-2:ny+3), bz_z_max(-2:nx+3,-2:ny+3))

    nx0 = 1
    nx1 = nx
    IF (bc_field(c_bd_x_min) == c_bc_cpml_laser) nx0 = cpml_x_min_laser_idx-1
    IF (bc_field(c_bd_x_max) == c_bc_cpml_laser) nx1 = cpml_x_max_laser_idx+1

    ex_x_min = 0.5_num * (ex(nx0,:,:) + ex(nx0-1,:,:))
    ey_x_min = ey(nx0,:,:)
    ez_x_min = ez(nx0,:,:)
    ex_x_max = 0.5_num * (ex(nx1,:,:) + ex(nx1-1,:,:))
    ey_x_max = ey(nx1,:,:)
    ez_x_max = ez(nx1,:,:)

    bx_x_min = bx(nx0,:,:)
    by_x_min = 0.5_num * (by(nx0,:,:) + by(nx0-1,:,:))
    bz_x_min = 0.5_num * (bz(nx0,:,:) + bz(nx0-1,:,:))
    bx_x_max = bx(nx1,:,:)
    by_x_max = 0.5_num * (by(nx1,:,:) + by(nx1-1,:,:))
    bz_x_max = 0.5_num * (bz(nx1,:,:) + bz(nx1-1,:,:))

    ny0 = 1
    ny1 = ny
    IF (bc_field(c_bd_y_min) == c_bc_cpml_laser) ny0 = cpml_y_min_laser_idx-1
    IF (bc_field(c_bd_y_max) == c_bc_cpml_laser) ny1 = cpml_y_max_laser_idx+1

    ex_y_min = ex(:,ny0,:)
    ey_y_min = 0.5_num * (ey(:,ny0,:) + ey(:,ny0-1,:))
    ez_y_min = ez(:,ny0,:)
    ex_y_max = ex(:,ny1,:)
    ey_y_max = 0.5_num * (ey(:,ny1,:) + ey(:,ny1-1,:))
    ez_y_max = ez(:,ny1,:)

    bx_y_min = 0.5_num * (bx(:,ny0,:) + bx(:,ny0-1,:))
    by_y_min = by(:,ny0,:)
    bz_y_min = 0.5_num * (bz(:,ny0,:) + bz(:,ny0-1,:))
    bx_y_max = 0.5_num * (bx(:,ny1,:) + bx(:,ny1-1,:))
    by_y_max = by(:,ny1,:)
    bz_y_max = 0.5_num * (bz(:,ny1,:) + bz(:,ny1-1,:))

    nz0 = 1
    nz1 = nz
    IF (bc_field(c_bd_z_min) == c_bc_cpml_laser) nz0 = cpml_z_min_laser_idx-1
    IF (bc_field(c_bd_z_max) == c_bc_cpml_laser) nz1 = cpml_z_max_laser_idx+1

    bx_z_min = 0.5_num * (bx(:,:,nz0) + bx(:,:,nz0-1))

    ex_z_min = ex(:,:,nz0)
    ey_z_min = ey(:,:,nz0)
    ez_z_min = 0.5_num * (ez(:,:,nz0) + ez(:,:,nz0-1))
    ex_z_max = ex(:,:,nz1)
    ey_z_max = ey(:,:,nz1)
    ez_z_max = 0.5_num * (ez(:,:,nz1) + ez(:,:,nz1-1))

    bx_z_min = 0.5_num * (bx(:,:,nz0) + bx(:,:,nz0-1))
    by_z_min = 0.5_num * (by(:,:,nz0) + by(:,:,nz0-1))
    bz_z_min = bz(:,:,nz0)
    bx_z_max = 0.5_num * (bx(:,:,nz1) + bx(:,:,nz1-1))
    by_z_max = 0.5_num * (by(:,:,nz1) + by(:,:,nz1-1))
    bz_z_max = bz(:,:,nz1)

  END SUBROUTINE setup_field_boundaries



  SUBROUTINE open_files

#ifdef NO_IO
    RETURN
#else
    INTEGER :: errcode
    LOGICAL :: exists

    IF (rank == 0) THEN
      WRITE(stat_file, '(a, ''/epoch3d.dat'')') TRIM(data_dir)
      IF (ic_from_restart) THEN
        INQUIRE(file=stat_file, exist=exists)
        IF (exists) THEN
          OPEN(unit=stat_unit, status='OLD', position='APPEND', &
              file=stat_file, iostat=errcode)
        ELSE
          OPEN(unit=stat_unit, status='NEW', file=stat_file, iostat=errcode)
        ENDIF
      ELSE
        OPEN(unit=stat_unit, status='REPLACE', file=stat_file, iostat=errcode)
      ENDIF
      IF (errcode /= 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Cannot create "epoch3d.dat" output file. The most common ' &
            // 'cause of this problem '
        PRINT*, 'is that the ouput directory does not exist'
        CALL abort_code(c_err_io_error)
        STOP
      ENDIF
      IF (ic_from_restart) THEN
        WRITE(stat_unit,*)
        WRITE(stat_unit,*) 'Restarting from ', TRIM(restart_filename)
        WRITE(stat_unit,*) ascii_header
      ELSE
        WRITE(stat_unit,*) ascii_header
        WRITE(stat_unit,*)
      ENDIF
    ENDIF
#endif

  END SUBROUTINE open_files



  SUBROUTINE flush_stat_file

#ifdef NO_IO
    RETURN
#else
    INTEGER :: errcode

    IF (rank == 0) THEN
      CLOSE(unit=stat_unit)
      OPEN(unit=stat_unit, status='OLD', position='APPEND', &
          file=stat_file, iostat=errcode)
    ENDIF
#endif

  END SUBROUTINE flush_stat_file



  SUBROUTINE close_files

#ifdef NO_IO
    RETURN
#else
    IF (rank == 0) CLOSE(unit=stat_unit)
#endif

  END SUBROUTINE close_files



  SUBROUTINE set_initial_values

    INTEGER :: seed

    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    ! Set up random number seed
    seed = 7842432
    IF (use_random_seed) CALL SYSTEM_CLOCK(seed)
    seed = seed + rank

    CALL random_init(seed)

  END SUBROUTINE set_initial_values



  SUBROUTINE set_plasma_frequency_dt

    INTEGER :: ispecies, ix, iy, iz
    REAL(num) :: min_dt, omega2, omega, k_max, fac1, fac2, clipped_dens

    IF (ic_from_restart) RETURN

    min_dt = 1000000.0_num
    k_max = 2.0_num * pi / MIN(dx, dy, dz)

    ! Identify the plasma frequency (Bohm-Gross dispersion relation)
    ! Note that this doesn't get strongly relativistic plasmas right
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type /= c_species_id_photon) THEN
        fac1 = q0**2 / species_list(ispecies)%mass / epsilon0
        fac2 = 3.0_num * k_max**2 * kb / species_list(ispecies)%mass
        IF (initial_conditions(ispecies)%density_max > 0) THEN
          DO iz = 1, nz
          DO iy = 1, ny
          DO ix = 1, nx
            clipped_dens = MIN(initial_conditions(ispecies)%density(ix,iy,iz), &
                initial_conditions(ispecies)%density_max)
            omega2 = fac1 * clipped_dens &
                + fac2 * MAXVAL(initial_conditions(ispecies)%temp(ix,iy,iz,:))
            IF (omega2 <= c_tiny) CYCLE
            omega = SQRT(omega2)
            IF (2.0_num * pi / omega < min_dt) min_dt = 2.0_num * pi / omega
          ENDDO ! ix
          ENDDO ! iy
          ENDDO ! iz
        ELSE
          DO iz = 1, nz
          DO iy = 1, ny
          DO ix = 1, nx
            omega2 = fac1 * initial_conditions(ispecies)%density(ix,iy,iz) &
                + fac2 * MAXVAL(initial_conditions(ispecies)%temp(ix,iy,iz,:))
            IF (omega2 <= c_tiny) CYCLE
            omega = SQRT(omega2)
            IF (2.0_num * pi / omega < min_dt) min_dt = 2.0_num * pi / omega
          ENDDO ! ix
          ENDDO ! iy
          ENDDO ! iz
        ENDIF
      ENDIF
    ENDDO

    CALL MPI_ALLREDUCE(min_dt, dt_plasma_frequency, 1, mpireal, MPI_MIN, &
        comm, errcode)
    ! Must resolve plasma frequency
    dt_plasma_frequency = dt_plasma_frequency / 2.0_num

  END SUBROUTINE set_plasma_frequency_dt



  SUBROUTINE set_dt        ! sets CFL limited step

    INTEGER :: io

    CALL set_plasma_frequency_dt
    CALL set_laser_dt

    dt = cfl * dx * dy * dz / SQRT((dx*dy)**2 + (dy*dz)**2 + (dz*dx)**2) / c
    IF (dt_plasma_frequency > c_tiny) dt = MIN(dt, dt_plasma_frequency)
    IF (dt_laser > c_tiny) dt = MIN(dt, dt_laser)
    dt = dt_multiplier * dt

    IF (.NOT. any_average) RETURN

    DO io = 1, n_io_blocks
      IF (.NOT. io_block_list(io)%any_average) CYCLE

      io_block_list(io)%average_time = MAX(io_block_list(io)%dt_average, &
          dt * io_block_list(io)%nstep_average)

      IF (io_block_list(io)%dt_min_average > 0 &
          .AND. io_block_list(io)%dt_min_average < dt) THEN
        IF (rank == 0) THEN
          PRINT*,'*** WARNING ***'
          PRINT*,'Time step is too small to satisfy "nstep_average"'
          PRINT*,'Averaging will occur over fewer time steps than specified'
          PRINT*,'Set "dt_multiplier" less than ', &
              dt_multiplier * io_block_list(io)%dt_min_average / dt, &
              ' to fix this'
        ENDIF
        io_block_list(io)%dt_min_average = -1
      ENDIF
    ENDDO

  END SUBROUTINE set_dt



  SUBROUTINE find_species_by_blockid(specname, species_number)

    CHARACTER(LEN=*), INTENT(IN) :: specname
    INTEGER, INTENT(OUT) :: species_number
    INTEGER :: ispecies, i1, i2

    CALL strip_species_blockid(specname, i1, i2)

    species_number = 0
    DO ispecies = 1,n_species
      IF (str_cmp(specname(i1:i2), species_list(ispecies)%name)) THEN
        species_number = ispecies
        EXIT
      ENDIF
    ENDDO

  END SUBROUTINE find_species_by_blockid



  SUBROUTINE find_species_by_id(species_id, species_number)

    CHARACTER(LEN=*), INTENT(IN) :: species_id
    INTEGER, INTENT(OUT) :: species_number
    INTEGER :: ispecies

    species_number = 0
    DO ispecies = 1,n_species
      IF (str_cmp(species_id, species_list(ispecies)%name)) THEN
        species_number = ispecies
        EXIT
      ENDIF
    ENDDO

  END SUBROUTINE find_species_by_id



  SUBROUTINE copy_string(str_in, str_out)

    CHARACTER(LEN=*), INTENT(IN) :: str_in
    CHARACTER(LEN=*), INTENT(INOUT) :: str_out
    INTEGER :: len1, len2, olen, i

    len1 = LEN_TRIM(str_in)
    len2 = LEN(str_out)
    olen = MIN(len1,len2)
    IF (olen > 0) THEN
      str_out(1:olen) = str_in(1:olen)
      DO i = olen+1,len2
        str_out(i:i) = ' '
      ENDDO
    ELSE
      DO i = 1,len2
        str_out(i:i) = ' '
      ENDDO
    ENDIF

  END SUBROUTINE copy_string



  SUBROUTINE find_species_by_id_or_blockid(species_id, block_id, species_number)

    CHARACTER(LEN=*), INTENT(INOUT) :: species_id
    CHARACTER(LEN=*), INTENT(IN) :: block_id
    INTEGER, INTENT(OUT) :: species_number
    INTEGER :: i1, i2

    IF (str_cmp(species_id, '__unknown__') &
        .OR. ICHAR(species_id(1:1)) == 0) THEN
      CALL strip_species_blockid(block_id, i1, i2)
      IF (i2 <= LEN_TRIM(species_id)) THEN
        CALL copy_string(block_id(i1:i2), species_id)
      ELSE
        CALL copy_string('__unknown__', species_id)
      ENDIF
    ENDIF

    CALL find_species_by_id(species_id, species_number)

  END SUBROUTINE find_species_by_id_or_blockid



  SUBROUTINE strip_species_blockid(name, i1, i2)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: i1, i2
    INTEGER :: ii

    i1 = 2
    i2 = LEN_TRIM(name)
    DO ii = 1, i2
      IF (name(ii:ii) == '/') RETURN
      i1 = i1 + 1
    ENDDO
    i1 = 1

  END SUBROUTINE strip_species_blockid



  SUBROUTINE restart_data(step)

    INTEGER, INTENT(OUT) :: step
    CHARACTER(LEN=c_id_length) :: code_name, block_id, mesh_id, str1, str2, str3
    CHARACTER(LEN=c_id_length) :: species_id
    CHARACTER(LEN=c_max_string_length) :: name, len_string
    INTEGER :: blocktype, datatype, code_io_version, string_len, ispecies
    INTEGER :: i, is, iblock, nblocks, ndims, geometry
    INTEGER(i8) :: npart, npart_local
    INTEGER(i8), ALLOCATABLE :: nparts(:), npart_locals(:), npart_proc(:)
    INTEGER, DIMENSION(4) :: dims
    INTEGER, ALLOCATABLE :: random_states_per_proc(:)
    REAL(num), DIMENSION(2*c_ndims) :: extents
    LOGICAL :: restart_flag, got_full
    LOGICAL, ALLOCATABLE :: species_found(:)
    TYPE(sdf_file_handle) :: sdf_handle
    TYPE(particle_species), POINTER :: species
    INTEGER, POINTER :: species_subtypes(:)
    INTEGER, POINTER :: species_subtypes_i4(:), species_subtypes_i8(:)

    got_full = .FALSE.
    npart_global = 0
    step = -1

    IF (rank == 0) THEN
      PRINT*,'Attempting to restart from file: ',TRIM(full_restart_filename)
    ENDIF

    CALL sdf_open(sdf_handle, full_restart_filename, comm, c_sdf_read)

    CALL sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
        string_len, restart_flag)

    IF (.NOT. restart_flag) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file ', TRIM(full_restart_filename), &
            ' is not a restart dump. Unable to continue.'
      ENDIF
      CALL abort_code(c_err_io_error)
      STOP
    ENDIF

    IF (.NOT.str_cmp(code_name, 'Epoch3d')) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF restart file was not generated by Epoch3d. Unable to ', &
            'continue.'
      ENDIF
      CALL abort_code(c_err_io_error)
      STOP
    ENDIF

    IF (string_len > c_max_string_length) THEN
      IF (rank == 0) THEN
        CALL integer_as_string(string_len, len_string)
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file string lengths are too large to read.'
        PRINT*, 'Please increase the size of "c_max_string_length" in ', &
            'shared_data.F90 to ','be at least ',TRIM(len_string)
      ENDIF
      CALL abort_code(c_err_io_error)
      STOP
    ENDIF

    ! Reset io_block parameters
    DO i = 1, n_io_blocks
      IF (io_block_list(i)%dump_first_after_restart) THEN
        io_block_list(i)%dump_first = .TRUE.
      ELSE
        io_block_list(i)%dump_first = .FALSE.
      ENDIF
      IF (io_block_list(i)%dt_snapshot > 0.0_num) THEN
        io_block_list(i)%time_prev = time
      ELSE
        io_block_list(i)%time_prev = 0.0_num
      ENDIF
      IF (io_block_list(i)%nstep_snapshot > 0) THEN
        io_block_list(i)%nstep_prev = step
      ELSE
        io_block_list(i)%nstep_prev = 0
      ENDIF
      IF (ASSOCIATED(io_block_list(i)%dump_at_nsteps)) THEN
        DO is = 1, SIZE(io_block_list(i)%dump_at_nsteps)
          IF (step >= io_block_list(i)%dump_at_nsteps(is)) THEN
            io_block_list(i)%dump_at_nsteps(is) = HUGE(1)
          ENDIF
        ENDDO
      ENDIF
      IF (ASSOCIATED(io_block_list(i)%dump_at_times)) THEN
        DO is = 1, SIZE(io_block_list(i)%dump_at_times)
          IF (time >= io_block_list(i)%dump_at_times(is)) THEN
            io_block_list(i)%dump_at_times(is) = HUGE(1.0_num)
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    nblocks = sdf_read_nblocks(sdf_handle)
    jobid = sdf_read_jobid(sdf_handle)

    IF (rank == 0) THEN
      PRINT*, 'Loading snapshot for time', time
      CALL create_ascii_header
#ifndef NO_IO
      WRITE(stat_unit,*) ascii_header
      WRITE(stat_unit,*)
#endif
    ENDIF

    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    dt_from_restart = 0.0_num

    IF (rank == 0) PRINT*, 'Input file contains', nblocks, 'blocks'

    CALL sdf_read_blocklist(sdf_handle)

    ALLOCATE(nparts(n_species))
    ALLOCATE(npart_locals(n_species))
    ALLOCATE(species_found(n_species))
    nparts = 0
    npart_locals = 0
    species_found = .FALSE.

    ! Scan file for particle species and allocate storage
    ! Both particle grids and CPU split blocks are interrogated. CPU split
    ! blocks take precedence if found and "use_exact_restart" is set.
    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      IF (blocktype == c_blocktype_cpu_split) THEN
        IF (.NOT.use_exact_restart .OR. datatype /= c_datatype_integer8) CYCLE
      ELSE IF (blocktype /= c_blocktype_point_mesh) THEN
        CYCLE
      ENDIF

      CALL sdf_read_point_mesh_info(sdf_handle, npart, geometry, species_id)
      CALL find_species_by_id_or_blockid(species_id, block_id, ispecies)
      IF (ispecies == 0) THEN
        IF (rank == 0) THEN
          IF (.NOT. str_cmp(species_id(1:6), 'subset')) THEN
            PRINT*, '*** WARNING ***'
            PRINT*, 'Particle species "', TRIM(species_id), &
                '" from restart dump ', 'not found in input deck. Ignoring.'
          ENDIF
        ENDIF
        CYCLE
      ENDIF

      IF (species_found(ispecies)) CYCLE

      species => species_list(ispecies)

      IF (blocktype == c_blocktype_cpu_split) THEN
        CALL sdf_read_cpu_split_info(sdf_handle, dims, geometry)

        ALLOCATE(npart_proc(dims(1)))
        CALL sdf_read_srl_cpu_split(sdf_handle, npart_proc)

        npart = 0
        DO i = 1,dims(1)
          npart = npart + npart_proc(i)
        ENDDO
        npart_local = npart_proc(rank+1)
        DEALLOCATE(npart_proc)

        npart_locals(ispecies) = npart_local
        nparts(ispecies) = npart
        species_found(ispecies) = .TRUE.
      ELSE IF (blocktype == c_blocktype_point_mesh) THEN
        CALL sdf_read_point_mesh_info(sdf_handle, npart)

        nparts(ispecies) = npart

        npart_local = npart / nproc
        IF (npart_local * nproc /= npart) THEN
          IF (rank < npart - npart_local * nproc) &
              npart_local = npart_local + 1
        ENDIF

        npart_locals(ispecies) = npart_local
      ENDIF
    ENDDO

    DEALLOCATE(species_found)

    ! Do the species allocation
    DO ispecies = 1, n_species
      species => species_list(ispecies)
      npart_local = npart_locals(ispecies)

      CALL create_allocated_partlist(species%attached_list, npart_local)

      npart_global = npart_global + nparts(ispecies)
      species%count = nparts(ispecies)
    ENDDO

    DEALLOCATE(nparts, npart_locals)

    CALL create_subtypes_for_load(species_subtypes, species_subtypes_i4, &
        species_subtypes_i8)

    CALL sdf_seek_start(sdf_handle)

    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      SELECT CASE(blocktype)
      CASE(c_blocktype_array)
        IF (use_exact_restart .AND. need_random_state &
            .AND. str_cmp(block_id, 'random_states')) THEN
          ALLOCATE(random_states_per_proc(4*nproc))
          CALL sdf_read_srl(sdf_handle, random_states_per_proc)
          CALL set_random_state(random_states_per_proc(4*rank+1:4*(rank+1)))
          DEALLOCATE(random_states_per_proc)
        ELSE IF (str_cmp(block_id, 'file_numbers')) THEN
          CALL sdf_read_array_info(sdf_handle, dims)
          IF (ndims /= 1 .OR. dims(1) /= SIZE(file_numbers)) THEN
            IF (rank == 0) THEN
              PRINT*, '*** WARNING ***'
              PRINT*, 'Output file numbers do not agree. Ignoring.'
            ENDIF
          ELSE
            CALL sdf_read_srl(sdf_handle, file_numbers)
          ENDIF
        ENDIF
      CASE(c_blocktype_constant)
        IF (str_cmp(block_id, 'dt_plasma_frequency')) THEN
          CALL sdf_read_srl(sdf_handle, dt_plasma_frequency)
        ELSE IF (str_cmp(block_id, 'dt')) THEN
          CALL sdf_read_srl(sdf_handle, dt_from_restart)
        ELSE IF (str_cmp(block_id, 'window_shift_fraction')) THEN
          CALL sdf_read_srl(sdf_handle, window_shift_fraction)
        ELSE IF (block_id(1:7) == 'weight/') THEN
          CALL find_species_by_blockid(block_id, ispecies)
          IF (ispecies == 0) CYClE
          CALL sdf_read_srl(sdf_handle, species_list(ispecies)%weight)
        ELSE IF (block_id(1:5) == 'nppc/') THEN
          CALL find_species_by_blockid(block_id, ispecies)
          IF (ispecies == 0) CYClE
          CALL sdf_read_srl(sdf_handle, species_list(ispecies)%npart_per_cell)
        ELSE IF (block_id(1:10) == 'time_prev/') THEN
          DO i = 1, n_io_blocks
            IF (str_cmp(block_id(11:), io_block_list(i)%name)) THEN
              CALL sdf_read_srl(sdf_handle, io_block_list(i)%time_prev)
              EXIT
            ENDIF
          ENDDO
        ELSE IF (block_id(1:11) == 'nstep_prev/') THEN
          DO i = 1, n_io_blocks
            IF (str_cmp(block_id(12:), io_block_list(i)%name)) THEN
              CALL sdf_read_srl(sdf_handle, io_block_list(i)%nstep_prev)
              EXIT
            ENDIF
          ENDDO
        ENDIF
      CASE(c_blocktype_plain_mesh)
        IF (.NOT.got_full) THEN
          IF (str_cmp(block_id, 'grid') &
              .OR. str_cmp(block_id, 'grid_full')) THEN
            CALL sdf_read_plain_mesh_info(sdf_handle, geometry, dims, extents)
            x_min = extents(1)
            x_max = extents(c_ndims+1)
            y_min = extents(2)
            y_max = extents(c_ndims+2)
            z_min = extents(3)
            z_max = extents(c_ndims+3)
            IF (str_cmp(block_id, 'grid_full')) got_full = .TRUE.
          ENDIF
        ENDIF
      CASE(c_blocktype_point_mesh)
        CALL sdf_read_point_mesh_info(sdf_handle, npart, geometry, species_id)

        CALL find_species_by_id_or_blockid(species_id, block_id, ispecies)
        IF (ispecies == 0) CYClE
        species => species_list(ispecies)

        npart_local = species%attached_list%count
        iterator_list => species%attached_list%head

        CALL sdf_read_point_mesh(sdf_handle, npart_local, &
            species_subtypes(ispecies), it_part)

      CASE(c_blocktype_plain_variable)
        CALL sdf_read_plain_variable_info(sdf_handle, dims, str1, mesh_id)

        IF (.NOT.str_cmp(mesh_id, 'grid')) CYCLE

        IF (dims(1) /= nx_global .OR. dims(2) /= ny_global &
            .OR. dims(3) /= nz_global) THEN
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Number of gridpoints in restart dump does not match', &
                ' the input deck.'
            CALL integer_as_string(nx_global, str1)
            CALL integer_as_string(ny_global, str2)
            CALL integer_as_string(nz_global, str3)
            PRINT*, 'Input deck grid: ', TRIM(str1), ',', TRIM(str2), &
                ',', TRIM(str3)
            CALL integer_as_string(dims(1), str1)
            CALL integer_as_string(dims(2), str2)
            CALL integer_as_string(dims(3), str3)
            PRINT*, 'Restart dump grid: ', TRIM(str1), ',', TRIM(str2), &
                ',', TRIM(str3)
          ENDIF
          CALL abort_code(c_err_bad_setup)
          STOP
        ENDIF

        IF (str_cmp(block_id, 'ex')) THEN
          CALL sdf_read_plain_variable(sdf_handle, ex, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'ey')) THEN
          CALL sdf_read_plain_variable(sdf_handle, ey, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'ez')) THEN
          CALL sdf_read_plain_variable(sdf_handle, ez, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'bx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, bx, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'by')) THEN
          CALL sdf_read_plain_variable(sdf_handle, by, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'bz')) THEN
          CALL sdf_read_plain_variable(sdf_handle, bz, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'jx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, jx, &
              subtype_field, subarray_field_big)

        ELSE IF (str_cmp(block_id, 'jy')) THEN
          CALL sdf_read_plain_variable(sdf_handle, jy, &
              subtype_field, subarray_field_big)

        ELSE IF (str_cmp(block_id, 'jz')) THEN
          CALL sdf_read_plain_variable(sdf_handle, jz, &
              subtype_field, subarray_field_big)

        ELSE IF (str_cmp(block_id, 'cpml_psi_eyx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_eyx, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_ezx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_ezx, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_byx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_byx, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_bzx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_bzx, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_exy')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_exy, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_ezy')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_ezy, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_bxy')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_bxy, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_bzy')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_bzy, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_exz')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_exz, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_eyz')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_eyz, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_bxz')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_bxz, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_byz')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_byz, &
              subtype_field, subarray_field)

        ENDIF

      CASE(c_blocktype_point_variable)
        CALL sdf_read_point_variable_info(sdf_handle, npart, mesh_id, &
            str1, species_id)

        CALL find_species_by_id_or_blockid(species_id, mesh_id, ispecies)
        IF (ispecies == 0) CYClE
        species => species_list(ispecies)

        IF (npart /= species%count) THEN
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Malformed restart dump. Number of particle variables', &
                ' does not match grid.'
          ENDIF
          CALL abort_code(c_err_io_error)
          STOP
        ENDIF

        iterator_list => species%attached_list%head
        npart_local = species%attached_list%count

        IF (block_id(1:3) == 'px/') THEN
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_px)

        ELSE IF (block_id(1:3) == 'py/') THEN
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_py)

        ELSE IF (block_id(1:3) == 'pz/') THEN
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_pz)

        ELSE IF (block_id(1:3) == 'id/') THEN
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
          ! Particle IDs may either be 4 or 8-byte integers, depending on the
          ! PARTICLE_ID[4] flag used when writing the file. We must read in
          ! the data using the precision written to file and then convert to
          ! the currently used precision.
          IF (datatype == c_datatype_integer8) THEN
            CALL sdf_read_point_variable(sdf_handle, npart_local, &
                species_subtypes_i8(ispecies), it_id8)
          ELSE
            CALL sdf_read_point_variable(sdf_handle, npart_local, &
                species_subtypes_i4(ispecies), it_id4)
          ENDIF
#else
          IF (rank == 0) THEN
            PRINT*, '*** WARNING ***'
            PRINT*, 'Discarding particle IDs.'
            PRINT*, 'To use, please recompile with the -DPARTICLE_ID option.'
          ENDIF
#endif

        ELSE IF (block_id(1:7) == 'weight/') THEN
#ifndef PER_SPECIES_WEIGHT
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_weight)
#else
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with per particle weight.'
            PRINT*, 'Please recompile without the -DPER_SPECIES_WEIGHT option.'
          ENDIF
          CALL abort_code(c_err_pp_options_missing)
          STOP
#endif

        ELSE IF (block_id(1:14) == 'optical depth/') THEN
#ifdef PHOTONS
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_optical_depth)
#else
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with optical depths.'
            PRINT*, 'Please recompile with the -DPHOTONS option.'
          ENDIF
          CALL abort_code(c_err_pp_options_missing)
          STOP
#endif

        ELSE IF (block_id(1:11) == 'qed energy/') THEN
#ifdef PHOTONS
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_qed_energy)
#else
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with QED energies.'
            PRINT*, 'Please recompile with the -DPHOTONS option.'
          ENDIF
          CALL abort_code(c_err_pp_options_missing)
          STOP
#endif

        ELSE IF (block_id(1:14) == 'trident depth/') THEN
#ifdef TRIDENT_PHOTONS
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_optical_depth_trident)
#else
          IF (rank == 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with Trident optical depths.'
            PRINT*, 'Please recompile with the -DTRIDENT_PHOTONS option.'
          ENDIF
          CALL abort_code(c_err_pp_options_missing)
          STOP
#endif
        ENDIF
      END SELECT
    ENDDO

    CALL sdf_close(sdf_handle)
    CALL free_subtypes_for_load(species_subtypes, species_subtypes_i4, &
        species_subtypes_i8)

    CALL setup_grid
    CALL set_thermal_bcs

    IF (rank == 0) PRINT*, 'Load from restart dump OK'

  END SUBROUTINE restart_data



  SUBROUTINE read_cpu_split

    CHARACTER(LEN=c_id_length) :: code_name, block_id
    CHARACTER(LEN=c_max_string_length) :: name
    INTEGER :: step, code_io_version, string_len, nblocks, ndims
    INTEGER :: blocktype, datatype, geometry, iblock, npx, npy, npz
    INTEGER, DIMENSION(4) :: dims
    LOGICAL :: restart_flag
    TYPE(sdf_file_handle) :: sdf_handle

    CALL sdf_open(sdf_handle, full_restart_filename, comm, c_sdf_read)

    CALL sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
        string_len, restart_flag)

    IF (.NOT. restart_flag) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file ', TRIM(full_restart_filename), &
            ' is not a restart dump. Unable to continue.'
      ENDIF
      CALL abort_code(c_err_io_error)
      STOP
    ENDIF

    IF (.NOT.str_cmp(code_name, 'Epoch3d')) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF restart file was not generated by Epoch3d. Unable to ', &
            'continue.'
      ENDIF
      CALL abort_code(c_err_io_error)
      STOP
    ENDIF

    nblocks = sdf_read_nblocks(sdf_handle)

    CALL sdf_read_blocklist(sdf_handle)

    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      IF (blocktype == c_blocktype_cpu_split &
          .AND. datatype /= c_datatype_integer8) THEN
        CALL sdf_read_cpu_split_info(sdf_handle, dims, geometry)
        npx = dims(1) + 1
        npy = dims(2) + 1
        npz = dims(3) + 1
        IF (npx * npy * npz == nproc) THEN
          nprocx = npx
          nprocy = npy
          nprocz = npz
          ALLOCATE(old_x_max(nprocx))
          ALLOCATE(old_y_max(nprocy))
          ALLOCATE(old_z_max(nprocz))
          CALL sdf_read_srl_cpu_split(sdf_handle, old_x_max, old_y_max, &
              old_z_max)
        ELSE
          IF (rank == 0) THEN
            PRINT*, '*** WARNING ***'
            PRINT'('' SDF restart file was generated using'', &
                & i4,'' CPUs.'')', npx * npy * npz
            PRINT*, 'Ignoring "use_exact_restart" flag.'
          ENDIF
          use_exact_restart = .FALSE.
        ENDIF
        EXIT
      ENDIF
    ENDDO

    CALL sdf_close(sdf_handle)

  END SUBROUTINE read_cpu_split



  FUNCTION it_part(array, npart_this_it, start, direction, param)

    REAL(num) :: it_part
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN) :: direction
    INTEGER, INTENT(IN), OPTIONAL :: param

    INTEGER(i8) :: ipart
    TYPE(particle), POINTER, SAVE :: cur

    IF (start) cur => iterator_list

    DO ipart = 1, npart_this_it
      cur%part_pos(direction) = array(ipart)
      cur => cur%next
    ENDDO

    it_part = 0

  END FUNCTION it_part



  FUNCTION it_px(array, npart_this_it, start, param)

    REAL(num) :: it_px
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%part_p(1) = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_px = 0

  END FUNCTION it_px



  FUNCTION it_py(array, npart_this_it, start, param)

    REAL(num) :: it_py
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%part_p(2) = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_py = 0

  END FUNCTION it_py



  FUNCTION it_pz(array, npart_this_it, start, param)

    REAL(num) :: it_pz
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%part_p(3) = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_pz = 0

  END FUNCTION it_pz



#ifndef PER_SPECIES_WEIGHT
  FUNCTION it_weight(array, npart_this_it, start, param)

    REAL(num) :: it_weight
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%weight = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_weight = 0

  END FUNCTION it_weight
#endif



#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
  FUNCTION it_id4(array, npart_this_it, start, param)

    USE constants
    INTEGER(i4) :: it_id4
    INTEGER(i4), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
#if PARTICLE_ID
      iterator_list%id = INT(array(ipart),i8)
#else
      iterator_list%id = array(ipart)
#endif
      iterator_list => iterator_list%next
    ENDDO

    it_id4 = 0

  END FUNCTION it_id4



  FUNCTION it_id8(array, npart_this_it, start, param)

    USE constants
    INTEGER(i8) :: it_id8
    INTEGER(i8), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
#if PARTICLE_ID
      iterator_list%id = array(ipart)
#else
      iterator_list%id = INT(array(ipart),i4)
#endif
      iterator_list => iterator_list%next
    ENDDO

    it_id8 = 0

  END FUNCTION it_id8
#endif



#ifdef PHOTONS
  FUNCTION it_optical_depth(array, npart_this_it, start, param)

    REAL(num) :: it_optical_depth
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%optical_depth = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_optical_depth = 0

  END FUNCTION it_optical_depth



  FUNCTION it_qed_energy(array, npart_this_it, start, param)

    REAL(num) :: it_qed_energy
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%particle_energy = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_qed_energy = 0

  END FUNCTION it_qed_energy



#ifdef TRIDENT_PHOTONS
  FUNCTION it_optical_depth_trident(array, npart_this_it, start, param)

    REAL(num) :: it_optical_depth_trident
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%optical_depth_tri = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_optical_depth_trident = 0

  END FUNCTION it_optical_depth_trident
#endif
#endif

END MODULE setup
