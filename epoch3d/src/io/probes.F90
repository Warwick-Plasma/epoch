MODULE probes

#ifdef PARTICLE_PROBES
  USE mpi
  USE sdf
  USE mpi_subtype_control
  USE partlist

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
    probe%radial = .FALSE.
    NULLIFY(probe%next)
    ALLOCATE(probe%use_species(n_species))
    probe%use_species = .FALSE.

  END SUBROUTINE init_probe



  SUBROUTINE attach_probe(probe)

    TYPE(particle_probe), POINTER :: probe
    TYPE(particle_probe), POINTER :: current
    INTEGER :: i

    ALLOCATE(probe%sampled_particles(1:n_species))

    DO i = 1, n_species
      probe%sampled_particles(i)=species_list(i)
      CALL create_empty_partlist(probe%sampled_particles(i)%attached_list)
    ENDDO
    current=>attached_probes
    IF (.NOT. ASSOCIATED(attached_probes)) THEN
      attached_probes=>current
    ELSE
      DO WHILE (ASSOCIATED(current%next))
        current=>current%next
      ENDDO
    ENDIF
    current%next=>probe
    n_probes=n_probes+1

  END SUBROUTINE attach_probe



  SUBROUTINE create_probe_subsets()

    TYPE(particle_probe), POINTER :: current
    TYPE(subset), POINTER :: current_subset
    TYPE(subset), DIMENSION(:), POINTER :: temp_subsets
    INTEGER :: iprobe

    IF (n_probes .EQ. 0) RETURN
    ALLOCATE(temp_subsets(n_subsets))
    temp_subsets = subset_list
    DEALLOCATE(subset_list)
    ALLOCATE(subset_list(n_subsets + n_probes))
    subset_list(1:n_subsets)=temp_subsets
    DEALLOCATE(temp_subsets)

    current=>attached_probes
    iprobe=0
    DO WHILE (ASSOCIATED(current))
      iprobe=iprobe+1
      current_subset=>subset_list(n_subsets+iprobe)
      current_subset%connected_probe=>current
      current_subset%name = current%name
      current_subset%dumpmask = dumpmask(c_dump_probes)
    ENDDO

    n_subsets = n_subsets + n_probes

  END SUBROUTINE create_probe_subsets

#endif

END MODULE probes
