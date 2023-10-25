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

MODULE file_injectors

  USE partlist
  USE utilities

  IMPLICIT NONE

CONTAINS

  SUBROUTINE open_injector_files(injector)

    ! Called in deck_injector_block after we have read the injector variables,
    ! and only for injectors with the inject_from_file flag. The file units are
    ! chosen such that each variable for each injector has a unique file unit.

    TYPE(injector_block), POINTER :: injector
    INTEGER :: inj_base_unit
    LOGICAL :: file_exists

    inj_base_unit = custom_base_unit + (injector%custom_id-1)*custom_var_num

    ! Injection times (mandatory)
    INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%t_data),&
        EXIST=file_exists)
    IF (file_exists) THEN
      OPEN(UNIT=inj_base_unit + unit_t, &
          FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%t_data))
      IF (rank == 0) THEN
        PRINT*, ''
        PRINT*, 'Successfully opened time file: ', TRIM(data_dir) // '/' // &
            TRIM(injector_filenames%t_data)
      END IF
    ELSE
      IF (rank == 0) THEN
        PRINT*, ''
        PRINT*, '*** ERROR ***'
        PRINT*, 'Unable to locate time file: ', &
            TRIM(data_dir) // '/' // TRIM(injector_filenames%t_data)
      END IF
    END IF

#ifndef PER_SPECIES_WEIGHT
    ! Particle weights (mandatory)
    INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%w_data),&
        EXIST=file_exists)
    IF (file_exists) THEN
      OPEN(UNIT=inj_base_unit + unit_w, &
          FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%w_data))
      IF (rank == 0) THEN
        PRINT*, ''
        PRINT*, 'Successfully opened weight file: ', TRIM(data_dir) // '/' // &
            TRIM(injector_filenames%w_data)
      END IF
    ELSE
      IF (rank == 0) THEN
        PRINT*, ''
        PRINT*, '*** ERROR ***'
        PRINT*, 'Unable to locate weight file: ', &
            TRIM(data_dir) // '/' // TRIM(injector_filenames%w_data)
      END IF
    END IF
#endif

    ! Position (positions that aren't in boundary dimension are mandatory)
    IF (injector%x_data_given) THEN
      INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%x_data),&
          EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(UNIT=inj_base_unit + unit_x, &
            FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%x_data))
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, 'Successfully opened x file: ', TRIM(data_dir) // '/' // &
              TRIM(injector_filenames%x_data)
        END IF
      ELSE
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, '*** ERROR ***'
          PRINT*, 'Unable to locate x file: ', &
              TRIM(data_dir) // '/' // TRIM(injector_filenames%x_data)
        END IF
      END IF
    END IF

    IF (injector%y_data_given) THEN
      INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%y_data),&
          EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(UNIT=inj_base_unit + unit_y, &
            FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%y_data))
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, 'Successfully opened y file: ', TRIM(data_dir) // '/' // &
              TRIM(injector_filenames%y_data)
        END IF
      ELSE
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, '*** ERROR ***'
          PRINT*, 'Unable to locate y file: ', &
              TRIM(data_dir) // '/' // TRIM(injector_filenames%y_data)
        END IF
      END IF
    END IF

    IF (injector%z_data_given) THEN
      INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%z_data),&
          EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(UNIT=inj_base_unit + unit_z, &
            FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%z_data))
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, 'Successfully opened z file: ', TRIM(data_dir) // '/' // &
              TRIM(injector_filenames%z_data)
        END IF
      ELSE
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, '*** ERROR ***'
          PRINT*, 'Unable to locate z file: ', &
              TRIM(data_dir) // '/' // TRIM(injector_filenames%z_data)
        END IF
      END IF
    END IF

    ! Momentum, optional - if missing, will be set to 0
    IF (injector%px_data_given) THEN
      INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%px_data),&
          EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(UNIT=inj_base_unit + unit_px, &
            FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%px_data))
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, 'Successfully opened px file: ', TRIM(data_dir) // '/' // &
              TRIM(injector_filenames%px_data)
        END IF
      ELSE
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, '*** ERROR ***'
          PRINT*, 'Unable to locate px file: ', &
              TRIM(data_dir) // '/' // TRIM(injector_filenames%px_data)
        END IF
      END IF
    END IF

    IF (injector%py_data_given) THEN
      INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%py_data),&
          EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(UNIT=inj_base_unit + unit_py, &
            FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%py_data))
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, 'Successfully opened py file: ', TRIM(data_dir) // '/' // &
              TRIM(injector_filenames%py_data)
        END IF
      ELSE
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, '*** ERROR ***'
          PRINT*, 'Unable to locate py file: ', &
              TRIM(data_dir) // '/' // TRIM(injector_filenames%py_data)
        END IF
      END IF
    END IF

    IF (injector%pz_data_given) THEN
      INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%pz_data),&
          EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(UNIT=inj_base_unit + unit_pz, &
            FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%pz_data))
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, 'Successfully opened pz file: ', TRIM(data_dir) // '/' // &
              TRIM(injector_filenames%pz_data)
        END IF
      ELSE
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, '*** ERROR ***'
          PRINT*, 'Unable to locate pz file: ', &
              TRIM(data_dir) // '/' // TRIM(injector_filenames%pz_data)
        END IF
      END IF
    END IF

#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
    ! Particle ID, optional - if missing, particle ID will be left empty
    IF (injector%id_data_given) THEN
      INQUIRE(FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%id_data),&
          EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(UNIT=inj_base_unit + unit_id, &
            FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%id_data))
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, 'Successfully opened id file: ', TRIM(data_dir) // '/' // &
              TRIM(injector_filenames%id_data)
        END IF
      ELSE
        IF (rank == 0) THEN
          PRINT*, ''
          PRINT*, '*** ERROR ***'
          PRINT*, 'Unable to locate id file: ', &
              TRIM(data_dir) // '/' // TRIM(injector_filenames%id_data)
        END IF
      END IF
    END IF
#endif

  END SUBROUTINE open_injector_files



  SUBROUTINE run_file_injection(injector)

    ! This subroutine reads particles from files (opened after the full injector
    ! block has been read, in a call from deck_injector_block.f90). We read from
    ! the injection time file (t_data) until we hit a particle injected after
    ! the NEXT timestep (time + dt). We position injected particles such that
    ! they pass the injection boundary at the time specified in t_data (assuming
    ! constant velocity), only injecting with ranks at the appropriate position

    TYPE(injector_block), POINTER :: injector
    REAL(num) :: mass, inv_m2c2
    REAL(num) :: x_in, y_in, z_in, px_in, py_in, pz_in
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: w_in
#endif
#ifdef PARTICLE_ID4
    INTEGER :: id_in
#elif PARTICLE_ID
    INTEGER(i8) :: id_in
#endif
    INTEGER :: boundary
    REAL(num) :: next_time, time_to_bdy
    REAL(num) :: vx, vy, vz, gamma, inv_gamma_mass, iabs_p
    REAL(num) :: x_start, y_start, z_start
    TYPE(particle), POINTER :: new
    TYPE(particle_list) :: plist
    LOGICAL :: no_particles_added, skip_processor

    ! Set to true if any of the associated files have reached the end of file
    IF (injector%file_finished) RETURN

    mass = species_list(injector%species)%mass
    IF (mass > c_tiny) THEN
      inv_m2c2 = 1.0_num/(mass*c)**2
    END IF

    no_particles_added = .TRUE.

    ! Add particles until we reach an injection time greater than the next
    ! timestep (particles must pass the injection boundary in the following
    ! particle push), or until there are no more particles to add
    DO
      ! We always start with injector%next_time known. Global time is a half
      ! timestep ahead of particle time when this is called
      IF (.NOT. injector%next_time < time + 0.5_num*dt ) EXIT

      ! If on x boundary read the y and z position, and vice versa
      boundary = injector%boundary
      IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
        CALL read_injector_real(unit_y, y_start, injector)
        CALL read_injector_real(unit_z, z_start, injector)
      ELSE IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
        CALL read_injector_real(unit_x, x_start, injector)
        CALL read_injector_real(unit_z, z_start, injector)
      ELSE IF (boundary == c_bd_z_min .OR. boundary == c_bd_z_max) THEN
        CALL read_injector_real(unit_x, x_start, injector)
        CALL read_injector_real(unit_y, y_start, injector)
      END IF

#ifndef PER_SPECIES_WEIGHT
      ! Read weight data
      CALL read_injector_real(unit_w, w_in, injector)
#endif

      ! Read momentum data
      IF (injector%px_data_given) THEN
        CALL read_injector_real(unit_px, px_in, injector)
      ELSE
        px_in = 0.0_num
      END IF

      IF (injector%py_data_given) THEN
        CALL read_injector_real(unit_py, py_in, injector)
      ELSE
        py_in = 0.0_num
      END IF

      IF (injector%pz_data_given) THEN
        CALL read_injector_real(unit_pz, pz_in, injector)
      ELSE
        pz_in = 0.0_num
      END IF

#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
      IF (injector%id_data_given) THEN
#ifdef PARTICLE_ID4
        CALL read_injector_int4(unit_id, id_in, injector)
#elif PARTICLE_ID
        CALL read_injector_int(unit_id, id_in, injector)
#endif
      ELSE
        id_in = 0
      END IF
#endif

      ! Ensure we still have values for each variable here
      IF (injector%file_finished) EXIT

      ! Identify processors which aren't adding this particle to the simulation
      boundary = injector%boundary
      skip_processor = .FALSE.
      IF (boundary == c_bd_x_min) THEN
        IF (.NOT. x_min_boundary) skip_processor = .TRUE.
      ELSE IF (boundary == c_bd_x_max) THEN
        IF (.NOT. x_max_boundary) skip_processor = .TRUE.
      ELSE IF (boundary == c_bd_y_min) THEN
        IF (.NOT. y_min_boundary) skip_processor = .TRUE.
      ELSE IF (boundary == c_bd_y_max) THEN
        IF (.NOT. y_max_boundary) skip_processor = .TRUE.
      ELSE IF (boundary == c_bd_z_min) THEN
        IF (.NOT. z_min_boundary) skip_processor = .TRUE.
      ELSE IF (boundary == c_bd_z_max) THEN
        IF (.NOT. z_max_boundary) skip_processor = .TRUE.
      END IF

      ! Skip the following phases if this rank isn't on the right boundary
      IF (skip_processor) THEN
        ! Calculate time for the next particle to be injected
        CALL read_injector_real(unit_t, next_time, injector)
        ! If there are no more particles to add, exit the loop, otherwise save
        ! the next injection time
        IF (injector%file_finished) EXIT
        injector%next_time = next_time
        CYCLE
      END IF

      ! If code has been restarted from an output dump, the position in our file
      ! will be lost. All particles with injection time < simulation time were
      ! injected in the previous run, and shouldn't be added again
      !
      ! Note that particles injected in the step before the restart will have
      ! been deleted, and must be re-added
      IF (injector%next_time < time - 0.5_num*dt) THEN
        CALL read_injector_real(unit_t, next_time, injector)
        IF (.NOT. injector%file_finished) injector%next_time = next_time
        CYCLE
      END IF

      ! Only ranks on the same boundary as the particle can reach here
      ! Calculate particle velocity
      IF (mass > c_tiny) THEN
        gamma = SQRT(1.0_num + (px_in**2 + py_in**2 + pz_in**2)*inv_m2c2)
        inv_gamma_mass = 1.0_num/(gamma*mass)
        vx = px_in*inv_gamma_mass
        vy = py_in*inv_gamma_mass
        vz = pz_in*inv_gamma_mass
      ELSE 
        iabs_p = 1.0_num / SQRT(px_in**2 + py_in**2 + pz_in**2)
        vx = px_in * iabs_p * c 
        vy = py_in * iabs_p * c
        vz = pz_in * iabs_p * c  
      END IF

      ! Calculate position of injection such that paritlces reach the boundary
      ! at next_time. Note that global time is a half timestep ahead of the time
      ! our particles are at
      time_to_bdy = (injector%next_time - (time-0.5_num*dt))
      IF (boundary == c_bd_x_min) THEN
        x_in = x_min - time_to_bdy * vx
        y_in = y_start - time_to_bdy * vy
        z_in = z_start - time_to_bdy * vz
      ELSE IF (boundary == c_bd_x_max) THEN
        x_in = x_max - time_to_bdy * vx
        y_in = y_start - time_to_bdy * vy
        z_in = z_start - time_to_bdy * vz
      ELSE IF (boundary == c_bd_y_min) THEN
        x_in = x_start - time_to_bdy * vx
        y_in = y_min - time_to_bdy * vy
        z_in = z_start - time_to_bdy * vz
      ELSE IF (boundary == c_bd_y_max) THEN
        x_in = x_start - time_to_bdy * vx
        y_in = y_max - time_to_bdy * vy
        z_in = z_start - time_to_bdy * vz
      ELSE IF (boundary == c_bd_z_min) THEN
        x_in = x_start - time_to_bdy * vx
        y_in = y_start - time_to_bdy * vy
        z_in = z_min - time_to_bdy * vz
      ELSE IF (boundary == c_bd_z_max) THEN
        x_in = x_start - time_to_bdy * vx
        y_in = y_start - time_to_bdy * vy
        z_in = z_max - time_to_bdy * vz
      END IF

      ! Now we have start position, find the processor which will be adding the
      ! particle
      IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
        ! Skip all processors which are at the wrong y position
        IF (y_in <= y_min_local .OR. y_in > y_max_local) THEN
          skip_processor = .TRUE.
        END IF
        ! Skip all processors which are at the wrong z position
        IF (z_in <= z_min_local .OR. z_in > z_max_local) THEN
          skip_processor = .TRUE.
        END IF

      ELSE IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
        ! Skip all processors which are at the wrong x position
        IF (x_in <= x_min_local .OR. x_in > x_max_local) THEN
          skip_processor = .TRUE.
        END IF
        ! Skip all processors which are at the wrong z position
        IF (z_in <= z_min_local .OR. z_in > z_max_local) THEN
          skip_processor = .TRUE.
        END IF

      ELSE IF (boundary == c_bd_z_min .OR. boundary == c_bd_z_max) THEN
        ! Skip all processors which are at the wrong x position
        IF (x_in <= x_min_local .OR. x_in > x_max_local) THEN
          skip_processor = .TRUE.
        END IF
        ! Skip all processors which are at the wrong y position
        IF (y_in <= y_min_local .OR. y_in > y_max_local) THEN
          skip_processor = .TRUE.
        END IF

      END IF

      ! Skip the particle adding phase if this rank isn't adding the particle
      IF (skip_processor) THEN
        ! Calculate time for the next particle to be injected
        CALL read_injector_real(unit_t, next_time, injector)
        ! If there are no more particles to add, exit the loop, otherwise save
        ! the next injection time
        IF (injector%file_finished) EXIT
        injector%next_time = next_time
        CYCLE
      END IF

      ! Create the particle and assign properties
      CALL create_particle(new)
      new%part_pos(1) = x_in
      new%part_pos(2) = y_in
      new%part_pos(3) = z_in
      new%part_p(1) = px_in
      new%part_p(2) = py_in
      new%part_p(3) = pz_in
#ifdef PER_PARTICLE_CHARGE_MASS
      new%charge = species_list(injector%species)%charge
      new%mass = mass
#endif
#ifndef PER_SPECIES_WEIGHT
      new%weight = w_in
#endif
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
      new%id = id_in
#endif

      ! Add the particle to our dummy particle list plist
      IF (no_particles_added) THEN
        ! Create an empty particle list the first time a particle is added
        CALL create_empty_partlist(plist)
        no_particles_added = .FALSE.
      END IF
      CALL add_particle_to_partlist(plist, new)

      ! Calculate the next time
      CALL read_injector_real(unit_t, next_time, injector)
      IF (injector%file_finished) THEN
        EXIT
      ELSE
        injector%next_time = next_time
      END IF
    END DO

    ! Append particles to the main particle list for this species on this
    ! processor if we have added particles
    IF (.NOT. no_particles_added) THEN
      CALL append_partlist(species_list(injector%species)%attached_list, plist)
    END IF

  END SUBROUTINE run_file_injection



  SUBROUTINE read_injector_real(unit_code, value, injector)

    ! Opens the file specified by the unit_code, reads a real variable and
    ! changes value to match. Use the unit codes defined in shared_data.F90 -
    ! e.g. unit_t for the t_data file. Also checks for the end of file condition

    TYPE(injector_block), POINTER :: injector
    REAL(num) :: value
    INTEGER :: unit_code, eof

    ! Read file, checking a value has been assigned to our variable, and that we
    ! have not reached the end of the file
    READ(custom_base_unit + (injector%custom_id-1)*custom_var_num + unit_code, &
        *, IOSTAT = eof) value

    ! File read successfully
    IF (eof == 0) RETURN

    ! Triggered if there are no more values to read
    IF (eof < 0) THEN
      injector%file_finished = .TRUE.

    ! Triggered if illegal value entered
    ELSE
      IF(rank == 0) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Illegal value found in read_injector_real'
        PRINT*,'Injected particle will not behave as expected'
        CALL abort_code(c_err_bad_value)
      END IF
    END IF

  END SUBROUTINE read_injector_real



#ifdef PARTICLE_ID
  SUBROUTINE read_injector_int(unit_code, value, injector)

    ! Opens the file specified by the unit_code, reads an integer variable and
    ! changes value to match. Use the unit codes defined in shared_data.F90 -
    ! e.g. unit_t for the t_data file. Also checks for the end of file condition

    TYPE(injector_block), POINTER :: injector
    INTEGER(i8) :: value
    INTEGER :: unit_code, eof

    ! Read file, checking a value has been assigned to our variable, and that we
    ! have not reached the end of the file
    READ(custom_base_unit + (injector%custom_id-1)*custom_var_num + unit_code, &
        *, IOSTAT = eof) value

    ! File read successfully
    IF (eof == 0) RETURN

    ! Triggered if there are no more values to read
    IF (eof < 0) THEN
      injector%file_finished = .TRUE.

    ! Triggered if illegal value entered
    ELSE
      IF(rank == 0) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Illegal value found in read_injector_int'
        PRINT*,'Injected particle will not behave as expected'
        CALL abort_code(c_err_bad_value)
      END IF
    END IF

  END SUBROUTINE read_injector_int
#endif



#ifdef PARTICLE_ID4
  SUBROUTINE read_injector_int4(unit_code, value, injector)

    ! Opens the file specified by the unit_code, reads an integer variable and
    ! changes value to match. Use the unit codes defined in shared_data.F90 -
    ! e.g. unit_t for the t_data file. Also checks for the end of file condition

    TYPE(injector_block), POINTER :: injector
    INTEGER :: value
    INTEGER :: unit_code, eof

    ! Read file, checking a value has been assigned to our variable, and that we
    ! have not reached the end of the file
    READ(custom_base_unit + (injector%custom_id-1)*custom_var_num + unit_code, &
        *, IOSTAT = eof) value

    ! File read successfully
    IF (eof == 0) RETURN

    ! Triggered if there are no more values to read
    IF (eof < 0) THEN
      injector%file_finished = .TRUE.

    ! Triggered if illegal value entered
    ELSE
      IF(rank == 0) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Illegal value found in read_injector_int4'
        PRINT*,'Injected particle will not behave as expected'
        CALL abort_code(c_err_bad_value)
      END IF
    END IF

  END SUBROUTINE read_injector_int4
#endif

END MODULE file_injectors
