! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
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
!
!------------------------------------------------------------------------------
!
! hybrid.F90
!
! This is the main interface to the hybrid routines. When running in hybrid
! mode, we ignore the main PIC loop and control is diverted to the hybrid PIC
! loop. This is also the script which initialises the hybrid arrays, and checks
! solids have been correctly defined.

MODULE hybrid

#ifdef HYBRID
  USE balance
  USE boundary
  USE current_smooth
  USE diagnostics
  USE injectors
  USE particles
  USE particle_migration
  USE random_generator
#ifdef PHOTONS
  USE photons
#endif
#ifdef BREMSSTRAHLUNG
  USE bremsstrahlung
#endif
  USE hy_elastic
  USE hy_fields
  USE hy_heating
  USE hy_ionisation_loss
  USE hy_resistivity
#endif

  IMPLICIT NONE

CONTAINS

#ifdef HYBRID

  SUBROUTINE run_hybrid_pic(push, halt, force_dump)

    LOGICAL :: push, halt, force_dump

    ! A copy of the EPOCH PIC loop, modified to run in the hybrid scheme
    DO
      ! Timing information.
      ! Functions/subroutines found in housekeeping/timer.f90
      IF (timer_collect) THEN
        CALL timer_stop(c_timer_step)
        CALL timer_reset
        timer_first(c_timer_step) = timer_walltime
      END IF

      ! Radiation scripts
#ifdef PHOTONS
      ! Non-linear Compton scatter/synchrotron radiation calculation (photons
      ! can be generated)
      IF (use_qed .AND. time > qed_start_time .AND. push) THEN
        CALL qed_update_optical_depth()
      END IF
#endif

      ! Evaluate fields a half timestep ahead of the particles
      IF (use_hybrid_fields) CALL run_hybrid_fields

      ! Logical flag set to true when particles can start to move
      push = (time >= particle_push_start_time)

      ! The following scripts will only be executed if particles can move
      IF (push) THEN

        ! Inject particles into the simulation
        CALL run_injectors

        ! .FALSE. this time to use load balancing threshold
        IF (use_balance) CALL balance_workload(.FALSE.)

        ! Move particles, leapfrogging over E and B
        CALL push_particles

        ! Pass current to neighbouring ranks (housekeeping/current_smooth.F90),
        ! and write new currents to neighbouring ghost cells
        CALL current_finish
        CALL field_bc(jx, ng)
        CALL field_bc(jy, ng)
        CALL field_bc(jz, ng)

        ! Calculate scattering from elastic collisions
        IF (use_hybrid_scatter) CALL elastic_scatter

        ! Obtain heat capacity to calculate the temperature change in
        ! hybrid_collisions and ohmic_heating
        CALL get_heat_capacity

        ! Calculates ionisational energy loss, and updates grid temperature
        IF (use_hybrid_collisions) CALL run_ionisation_loss

        ! Updates grid temperature due to Ohmic heating
        IF (use_ohmic) CALL ohmic_heating
        CALL clear_heat_capacity

        ! Now that temperature has been fully updated, re-evaluate resistivity
        IF (use_hy_ionisation) CALL update_ionisation
        IF (use_hy_cou_log) CALL update_coulomb_logarithm
        IF (use_ion_temp) CALL thermal_equilibration
        CALL update_resistivity

        ! Migrate particle species if they pass the migration criteria
        IF (use_particle_migration) CALL migrate_particles(step)

        ! See housekeeping/partlist.F90
        CALL update_particle_count

      END IF

      CALL check_for_stop_condition(halt, force_dump)
      IF (halt) EXIT
      step = step + 1
      time = time + dt / 2.0_num
      CALL output_routines(step)
      ! Check we have not passed the end condition
      IF ((step >= nsteps .AND. nsteps >= 0) .OR. (time >= t_end) &
          .OR. halt) EXIT
      time = time + dt / 2.0_num

      ! Iterate B and E, such that they are evaluated at the same time as the
      ! particles (note: main PIC loop also does this after the output dump)
      IF (use_hybrid_fields) CALL run_hybrid_fields

    END DO

  END SUBROUTINE run_hybrid_pic



  SUBROUTINE initialise_hybrid

    ! This subroutine initialises the hybrid arrays, and sets the values of some
    ! constants, to speed up the code

    REAL(num) :: resistivity_init, max_ne, z_real, iex
    INTEGER :: ix, iy, iz, i_sol
    INTEGER :: max_id
    INTEGER :: io, iu
    LOGICAL :: sendbuf(2), recvbuf(2)
    INTEGER :: err

    IF (use_hybrid) THEN

      ! Ensure solids have been called with all necessary variables specified
      CALL check_solids

      ! Initialise variables for extra physics
      IF (use_hybrid_collisions) CALL setup_hy_ionisation_loss
      IF (use_ohmic .OR. use_hybrid_collisions) CALL setup_heating

      ! Preset useful constants for solids
      ALLOCATE(hy_sum_ne(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      hy_sum_ne = 0
      DO i_sol = 1, solid_count
        z_real = REAL(solid_array(i_sol)%z, num)
        iex = solid_array(i_sol)%iex
        hy_sum_ne = hy_sum_ne + solid_array(i_sol)%el_density

        solid_array(i_sol)%el_density = z_real * solid_array(i_sol)%ion_density
        solid_array(i_sol)%z_prime = z_real**(-4.0_num/3.0_num) *  kelvin_to_ev
        solid_array(i_sol)%iex_term = 2.0_num  / (iex / m0c2)**2
        solid_array(i_sol)%dedx_c = 1.0_num + 2.0_num &
            * LOG(iex / (h_bar * q0) * SQRT(epsilon0 * m0))
        solid_array(i_sol)%theta_fac = z_real * q0**4  / (2.0_num * pi &
            * epsilon0**2)
        solid_array(i_sol)%ln_s = 4.0_num * epsilon0 * h_planck &
            / (z_real**(1.0_num/3.0_num) * m0 * q0**2)
      END DO

      ! Allocate additional arrays for running in hybrid mode. These require
      ! extra remapping scripts in balance.F90 (for domain change in
      ! load-balance)
      ALLOCATE(resistivity(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(resistivity_model(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))

      ! Assign a resistivity model to each cell based on present solids
      ! (default to vacuum)
      resistivity_model = c_resist_vacuum
      DO iz = 1-ng, nz+ng
        DO iy = 1-ng, ny+ng
          DO ix = 1-ng, nx+ng

            ! Use the same resistivity model as the solid with the highest
            ! local electron number density
            max_ne = 0.0_num
            DO i_sol = 1, solid_count
              IF (solid_array(i_sol)%el_density(ix,iy,iz) > max_ne) THEN
                max_ne = solid_array(i_sol)%el_density(ix,iy,iz)
                resistivity_model(ix,iy,iz) = solid_array(i_sol)%res_model
              END IF
            END DO

            ! Do we need ionisation and Coulomb logarithm scripts?
            ! If using Lee-More resistivity
            IF (resistivity_model(ix,iy,iz) == c_resist_rlm) THEN
              use_hy_cou_log = .TRUE.
              use_hy_ionisation = .TRUE.
            END IF

          END DO
        END DO
      END DO

      ! If one processor needs to calculate a Coulomb logarithm and deal with
      ! background ionisation, then all processors must
      sendbuf(1) = use_hy_cou_log
      sendbuf(2) = use_hy_ionisation
      CALL MPI_ALLREDUCE(sendbuf, recvbuf, 2, MPI_LOGICAL, MPI_LOR, comm, err)
      use_hy_cou_log = recvbuf(1)
      use_hy_ionisation = recvbuf(2)

      ! Allocate additional arrays for extra physics processes
      IF (use_hy_ionisation) THEN
        ALLOCATE(ion_charge(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
        ALLOCATE(ion_z_avg(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
        ALLOCATE(ion_reduced_density(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
        ion_z_avg = 0.0_num
        ion_reduced_density = 0.0_num
      END IF
      IF (use_hy_cou_log) THEN
        ALLOCATE(ion_cou_log(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      END IF
      IF (use_ion_temp) THEN
        ALLOCATE(ion_A(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
        ion_A = 0.0_num
      END IF
      IF (use_ion_temp .OR. use_hy_ionisation) THEN
        ALLOCATE(ion_ni(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
        ion_ni = 0.0_num
      END IF

      ! Initialise global background variables (some routines use an "average"
      ! solid background). Note reduced density is in g/cm³
      DO i_sol = 1, solid_count
        IF (ALLOCATED(ion_ni)) ion_ni = ion_ni + solid_array(i_sol)%ion_density
        IF (ALLOCATED(ion_z_avg)) ion_z_avg = ion_z_avg + solid_array(i_sol)%z &
            * solid_array(i_sol)%ion_density
        IF (ALLOCATED(ion_a)) ion_a = ion_a + solid_array(i_sol)%mass_no &
            * solid_array(i_sol)%ion_density
        IF (ALLOCATED(ion_reduced_density)) ion_reduced_density = &
            ion_reduced_density + 0.001_num * amu &
            * solid_array(i_sol)%ion_density / solid_array(i_sol)%z
      END DO
      IF (ALLOCATED(ion_z_avg)) ion_z_avg = ion_z_avg / ion_ni
      IF (ALLOCATED(ion_a)) ion_a = ion_a / ion_ni

      ! Initialise resistivity
      IF (use_hy_ionisation) CALL update_ionisation
      IF (use_hy_cou_log) CALL update_coulomb_logarithm
      CALL update_resistivity

      IF (rank == 0) PRINT*, 'Code is running in hybrid mode'

    ELSE
      ! Do not try to output hybrid variables if we aren't running in hybrid
      ! mode
      IF (rank == 0) THEN
        PRINT*, ''
        PRINT*, 'Code is not running in hybrid mode'
        PRINT*, 'Any attempts to output hybrid-only variables (like ', &
           'resistivity) will be ignored'
        PRINT*, 'Switch off the -DHYBRID pre-processor flag to stop this ', &
            'message from printing'
        PRINT*, ''
      END IF
    END IF

  END SUBROUTINE initialise_hybrid



  SUBROUTINE deallocate_hybrid

    ! Deallocates all hybrid global arrays, and solid arrays

    INTEGER :: i_sol

    ! Solids
    DO i_sol = 1, solid_count
      DEALLOCATE(solid_array(i_sol)%ion_density, solid_array(i_sol)% el_density)
    END DO
    DEALLOCATE(solid_array)

    ! Global arrays
    DEALLOCATE(hy_te, resistivity, resistivity_model)
    DEALLOCATE(jbx, jby, jbz)

    ! Ionisation/resistivity optional arrays
    IF (use_hy_cou_log) DEALLOCATE(ion_cou_log)
    IF (use_ion_temp) DEALLOCATE(ion_a)
    IF (use_ion_temp .OR. use_hy_ionisation) DEALLOCATE(ion_ni)
    IF (ALLOCATED(hy_ti)) DEALLOCATE(hy_ti)
    IF (use_hy_ionisation) &
        DEALLOCATE(ion_charge, ion_z_avg, ion_reduced_density)

  END SUBROUTINE deallocate_hybrid



  SUBROUTINE check_solids

    INTEGER :: i_sol, io, iu

    ! Do we have a solid background?
    IF (solid_count == 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No solids specified! Please set all solid ion ', &
              'number densities to positive values.'
          WRITE(io,*) 'Code will terminate.'
        END DO
      END IF
      errcode = c_err_bad_value + c_err_terminate
    END IF

    ! Check the solid parameters have been entered correctly
    DO i_sol = 1, solid_count

      ! Check ion number density is positive
      IF (MINVAL(solid_array(i_sol)%ion_density) < 0.0_num) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Please ensure the solid ion number densities are ', &
                'positive values in all cells.'
            WRITE(io,*) 'Code will terminate.'
          END DO
        END IF
        errcode = c_err_bad_value + c_err_terminate
      END IF

      ! Check the atomic number is positive
      IF (solid_array(i_sol)%Z < 1) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Please set all solid atomic numbers to positive ', &
                'values greater than 0.'
            WRITE(io,*) 'Code will terminate.'
          END DO
        END IF
        errcode = c_err_bad_value + c_err_terminate
      END IF

      ! Check the mean excitation energy is positive (only used in collisions)
      IF (solid_array(i_sol)%iex < 0.0_num &
          .AND. use_hybrid_collisions) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Please set all solid mean excitation energies to ', &
                'positive values for hybrid collisions'
            WRITE(io,*) 'Code will terminate.'
          END DO
        END IF
        errcode = c_err_bad_value + c_err_terminate
      END IF

    END DO

  END SUBROUTINE check_solids



#else
  ! If hybrid precompiler option is off, let the routines called by other
  ! modules do nothing
  SUBROUTINE run_hybrid_pic(push, halt, force_dump)

    LOGICAL :: push, halt, force_dump

  END SUBROUTINE run_hybrid_pic



  SUBROUTINE deallocate_hybrid

  END SUBROUTINE deallocate_hybrid

#endif

END MODULE hybrid
