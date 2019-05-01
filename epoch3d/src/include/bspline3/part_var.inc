  part_var = &
        gz(-2) * (gy(-2) * (gx(-2) * grid_var(cell_x1-2,cell_y1-2,cell_z1-2) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1-2,cell_z1-2) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1-2,cell_z1-2) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1-2,cell_z1-2) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1-2,cell_z1-2)) &
      +           gy(-1) * (gx(-2) * grid_var(cell_x1-2,cell_y1-1,cell_z1-2) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1-1,cell_z1-2) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1-1,cell_z1-2) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1-1,cell_z1-2) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1-1,cell_z1-2)) &
      +           gy( 0) * (gx(-2) * grid_var(cell_x1-2,cell_y1  ,cell_z1-2) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1  ,cell_z1-2) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1  ,cell_z1-2) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1  ,cell_z1-2) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1  ,cell_z1-2)) &
      +           gy( 1) * (gx(-2) * grid_var(cell_x1-2,cell_y1+1,cell_z1-2) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1+1,cell_z1-2) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1+1,cell_z1-2) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1+1,cell_z1-2) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1+1,cell_z1-2)) &
      +           gy( 2) * (gx(-2) * grid_var(cell_x1-2,cell_y1+2,cell_z1-2) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1+2,cell_z1-2) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1+2,cell_z1-2) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1+2,cell_z1-2) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1+2,cell_z1-2)))
  part_var = part_var &
      + gz(-1) * (gy(-2) * (gx(-2) * grid_var(cell_x1-2,cell_y1-2,cell_z1-1) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1-2,cell_z1-1) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1-2,cell_z1-1) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1-2,cell_z1-1) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1-2,cell_z1-1)) &
      +           gy(-1) * (gx(-2) * grid_var(cell_x1-2,cell_y1-1,cell_z1-1) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1-1,cell_z1-1) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1-1,cell_z1-1) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1-1,cell_z1-1) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1-1,cell_z1-1)) &
      +           gy( 0) * (gx(-2) * grid_var(cell_x1-2,cell_y1  ,cell_z1-1) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1  ,cell_z1-1) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1  ,cell_z1-1) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1  ,cell_z1-1) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1  ,cell_z1-1)) &
      +           gy( 1) * (gx(-2) * grid_var(cell_x1-2,cell_y1+1,cell_z1-1) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1+1,cell_z1-1) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1+1,cell_z1-1) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1+1,cell_z1-1) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1+1,cell_z1-1)) &
      +           gy( 2) * (gx(-2) * grid_var(cell_x1-2,cell_y1+2,cell_z1-1) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1+2,cell_z1-1) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1+2,cell_z1-1) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1+2,cell_z1-1) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1+2,cell_z1-1)))
  part_var = part_var &
      + gz( 0) * (gy(-2) * (gx(-2) * grid_var(cell_x1-2,cell_y1-2,cell_z1  ) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1-2,cell_z1  ) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1-2,cell_z1  ) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1-2,cell_z1  ) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1-2,cell_z1  )) &
      +           gy(-1) * (gx(-2) * grid_var(cell_x1-2,cell_y1-1,cell_z1  ) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1-1,cell_z1  ) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1-1,cell_z1  ) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1-1,cell_z1  ) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1-1,cell_z1  )) &
      +           gy( 0) * (gx(-2) * grid_var(cell_x1-2,cell_y1  ,cell_z1  ) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1  ,cell_z1  ) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1  ,cell_z1  ) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1  ,cell_z1  ) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1  ,cell_z1  )) &
      +           gy( 1) * (gx(-2) * grid_var(cell_x1-2,cell_y1+1,cell_z1  ) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1+1,cell_z1  ) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1+1,cell_z1  ) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1+1,cell_z1  ) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1+1,cell_z1  )) &
      +           gy( 2) * (gx(-2) * grid_var(cell_x1-2,cell_y1+2,cell_z1  ) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1+2,cell_z1  ) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1+2,cell_z1  ) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1+2,cell_z1  ) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1+2,cell_z1  )))
  part_var = part_var &
      + gz( 1) * (gy(-2) * (gx(-2) * grid_var(cell_x1-2,cell_y1-2,cell_z1+1) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1-2,cell_z1+1) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1-2,cell_z1+1) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1-2,cell_z1+1) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1-2,cell_z1+1)) &
      +           gy(-1) * (gx(-2) * grid_var(cell_x1-2,cell_y1-1,cell_z1+1) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1-1,cell_z1+1) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1-1,cell_z1+1) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1-1,cell_z1+1) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1-1,cell_z1+1)) &
      +           gy( 0) * (gx(-2) * grid_var(cell_x1-2,cell_y1  ,cell_z1+1) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1  ,cell_z1+1) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1  ,cell_z1+1) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1  ,cell_z1+1) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1  ,cell_z1+1)) &
      +           gy( 1) * (gx(-2) * grid_var(cell_x1-2,cell_y1+1,cell_z1+1) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1+1,cell_z1+1) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1+1,cell_z1+1) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1+1,cell_z1+1) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1+1,cell_z1+1)) &
      +           gy( 2) * (gx(-2) * grid_var(cell_x1-2,cell_y1+2,cell_z1+1) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1+2,cell_z1+1) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1+2,cell_z1+1) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1+2,cell_z1+1) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1+2,cell_z1+1)))
  part_var = part_var &
      + gz( 2) * (gy(-2) * (gx(-2) * grid_var(cell_x1-2,cell_y1-2,cell_z1+2) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1-2,cell_z1+2) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1-2,cell_z1+2) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1-2,cell_z1+2) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1-2,cell_z1+2)) &
      +           gy(-1) * (gx(-2) * grid_var(cell_x1-2,cell_y1-1,cell_z1+2) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1-1,cell_z1+2) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1-1,cell_z1+2) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1-1,cell_z1+2) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1-1,cell_z1+2)) &
      +           gy( 0) * (gx(-2) * grid_var(cell_x1-2,cell_y1  ,cell_z1+2) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1  ,cell_z1+2) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1  ,cell_z1+2) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1  ,cell_z1+2) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1  ,cell_z1+2)) &
      +           gy( 1) * (gx(-2) * grid_var(cell_x1-2,cell_y1+1,cell_z1+2) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1+1,cell_z1+2) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1+1,cell_z1+2) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1+1,cell_z1+2) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1+1,cell_z1+2)) &
      +           gy( 2) * (gx(-2) * grid_var(cell_x1-2,cell_y1+2,cell_z1+2) &
      +                     gx(-1) * grid_var(cell_x1-1,cell_y1+2,cell_z1+2) &
      +                     gx( 0) * grid_var(cell_x1  ,cell_y1+2,cell_z1+2) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1+2,cell_z1+2) &
      +                     gx( 2) * grid_var(cell_x1+2,cell_y1+2,cell_z1+2)))
