PRO ParticlePhasePlane, positions, variable, direction=dir, _extra=extra

  IF (N_ELEMENTS(dir) EQ 0) THEN dir = 0
  PLOT, positions.particlepositions[*,dir], variable, psym=3, _extra=extra

END
