        ex_part = &
              gz( 0) * (gy( 0) * (hx( 0) * ex(cell_x2  ,cell_y1  ,cell_z1  ) &
            +                     hx( 1) * ex(cell_x2+1,cell_y1  ,cell_z1  )) &
            +           gy( 1) * (hx( 0) * ex(cell_x2  ,cell_y1+1,cell_z1  ) &
            +                     hx( 1) * ex(cell_x2+1,cell_y1+1,cell_z1  ))) &
            + gz( 1) * (gy( 0) * (hx( 0) * ex(cell_x2  ,cell_y1  ,cell_z1+1) &
            +                     hx( 1) * ex(cell_x2+1,cell_y1  ,cell_z1+1)) &
            +           gy( 1) * (hx( 0) * ex(cell_x2  ,cell_y1+1,cell_z1+1) &
            +                     hx( 1) * ex(cell_x2+1,cell_y1+1,cell_z1+1)))

        ey_part = &
              gz( 0) * (hy( 0) * (gx( 0) * ey(cell_x1  ,cell_y2  ,cell_z1  ) &
            +                     gx( 1) * ey(cell_x1+1,cell_y2  ,cell_z1  )) &
            +           hy( 1) * (gx( 0) * ey(cell_x1  ,cell_y2+1,cell_z1  ) &
            +                     gx( 1) * ey(cell_x1+1,cell_y2+1,cell_z1  ))) &
            + gz( 1) * (hy( 0) * (gx( 0) * ey(cell_x1  ,cell_y2  ,cell_z1+1) &
            +                     gx( 1) * ey(cell_x1+1,cell_y2  ,cell_z1+1)) &
            +           hy( 1) * (gx( 0) * ey(cell_x1  ,cell_y2+1,cell_z1+1) &
            +                     gx( 1) * ey(cell_x1+1,cell_y2+1,cell_z1+1)))

        ez_part = &
              hz( 0) * (gy( 0) * (gx( 0) * ez(cell_x1  ,cell_y1  ,cell_z2  ) &
            +                     gx( 1) * ez(cell_x1+1,cell_y1  ,cell_z2  )) &
            +           gy( 1) * (gx( 0) * ez(cell_x1  ,cell_y1+1,cell_z2  ) &
            +                     gx( 1) * ez(cell_x1+1,cell_y1+1,cell_z2  ))) &
            + hz( 1) * (gy( 0) * (gx( 0) * ez(cell_x1  ,cell_y1  ,cell_z2+1) &
            +                     gx( 1) * ez(cell_x1+1,cell_y1  ,cell_z2+1)) &
            +           gy( 1) * (gx( 0) * ez(cell_x1  ,cell_y1+1,cell_z2+1) &
            +                     gx( 1) * ez(cell_x1+1,cell_y1+1,cell_z2+1)))
