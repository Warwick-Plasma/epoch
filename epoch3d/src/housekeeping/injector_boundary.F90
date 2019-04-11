MODULE injector_boundary

  USE shared_data
  USE injectors

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setup_injector_boundaries, create_boundary_injectors, &
      finish_setup_injector_boundaries

  ! Temporary injector list, while setting up
  TYPE(injector_block), POINTER :: auto_injectors
  TYPE(injector_block), POINTER :: prev_inj

CONTAINS

  SUBROUTINE setup_injector_boundaries

    INTEGER :: i
    TYPE(particle_species), POINTER :: current

    DO i = 1, n_species
      current => species_list(i)
      IF (ANY(current%bc_particle == c_bc_continue)) THEN
        ! Create and setup injector to continue species
        CALL create_boundary_injectors(current%bc_particle, i)
      END IF

    END DO

    CALL reset_auto_list

  END SUBROUTINE setup_injector_boundaries



  ! Point prev_inj back to the start of the auto list
  SUBROUTINE reset_auto_list

    prev_inj => auto_injectors

  END SUBROUTINE reset_auto_list


  ! Attach an auto injector to the list
  SUBROUTINE attach_auto_injector(working_injector)

    TYPE(injector_block), INTENT(IN), POINTER :: working_injector

    NULLIFY(working_injector%next)

    IF (ASSOCIATED(prev_inj)) THEN
      prev_inj%next => working_injector
    ELSE
      auto_injectors => working_injector
    END IF

    prev_inj => working_injector

  END SUBROUTINE attach_auto_injector



  ! Create an injector for a return-boundary species

  SUBROUTINE create_boundary_injectors(bcs, ispecies)

    INTEGER, INTENT(IN) :: ispecies
    INTEGER, DIMENSION(c_ndims*2), INTENT(IN) :: bcs
    TYPE(injector_block), POINTER :: working_injector
    INTEGER :: i

    DO i = 1, 2
      IF (bcs(i) /= c_bc_continue) CYCLE
      use_injectors = .TRUE.
      need_random_state = .TRUE.

      ALLOCATE(working_injector)
      CALL init_injector(i, working_injector)

      CALL copy_stack(species_list(ispecies)%drift_function(1), &
          working_injector%drift_function(1))
      CALL copy_stack(species_list(ispecies)%drift_function(2), &
          working_injector%drift_function(2))
      CALL copy_stack(species_list(ispecies)%drift_function(3), &
          working_injector%drift_function(3))
      CALL copy_stack(species_list(ispecies)%density_function, &
          working_injector%density_function)
      CALL copy_stack(species_list(ispecies)%temperature_function(1), &
          working_injector%temperature_function(1))
      CALL copy_stack(species_list(ispecies)%temperature_function(2), &
          working_injector%temperature_function(2))
      CALL copy_stack(species_list(ispecies)%temperature_function(3), &
          working_injector%temperature_function(3))

      working_injector%species = ispecies
      working_injector%npart_per_cell = 1 ! Temporary, fixed after load

      CALL attach_auto_injector(working_injector)
    END DO

  END SUBROUTINE create_boundary_injectors



  ! Fix up the ppc now we've loaded the particles
  ! attach the injectors to the main injector lists
  ! and fix up boundary values

  SUBROUTINE finish_setup_injector_boundaries

    TYPE(injector_block), POINTER :: working_injector, next
    INTEGER, DIMENSION(c_ndims*2) :: bcs
    INTEGER :: bnd, species

    working_injector => auto_injectors

    DO WHILE(ASSOCIATED(working_injector))
      next => working_injector%next

      species = working_injector%species
      bcs = species_list(species)%bc_particle
      bnd = working_injector%boundary

      IF (bcs(bnd) == c_bc_continue) THEN
        ! Fix nppc
        working_injector%npart_per_cell = &
            FLOOR(species_list(species)%npart_per_cell)
        ! Set boundary to open for future
        species_list(species)%bc_particle(bnd) = c_bc_open

        ! Attach to main injector list
        CALL attach_injector(working_injector)

      END IF

      working_injector => next

    END DO
    NULLIFY(auto_injectors)

  END SUBROUTINE finish_setup_injector_boundaries


END MODULE
