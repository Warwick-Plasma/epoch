MODULE prefetch

  USE shared_data

CONTAINS

  !Uses Intel 12 prefetch call to move particle data into cache.
  SUBROUTINE prefetch_particle(p)

    TYPE(particle), INTENT(INOUT) :: p

#ifdef PREFETCH
    CALL mm_prefetch(p%part_p(1)) !, 1)
    CALL mm_prefetch(p%weight) !, 1)
#endif

  END SUBROUTINE prefetch_particle

END MODULE prefetch
