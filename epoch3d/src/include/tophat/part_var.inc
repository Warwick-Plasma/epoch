  part_var = &
        gz( 0) * (gy( 0) * (gx( 0) * grid_var(cell_x1  ,cell_y1  ,cell_z1  ) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1  ,cell_z1  )) &
      +           gy( 1) * (gx( 0) * grid_var(cell_x1  ,cell_y1+1,cell_z1  ) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1+1,cell_z1  ))) &
      + gz( 1) * (gy( 0) * (gx( 0) * grid_var(cell_x1  ,cell_y1  ,cell_z1+1) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1  ,cell_z1+1)) &
      +           gy( 1) * (gx( 0) * grid_var(cell_x1  ,cell_y1+1,cell_z1+1) &
      +                     gx( 1) * grid_var(cell_x1+1,cell_y1+1,cell_z1+1)))
