        bx_part = &
              hz(-1) * (hy(-1) * (gx(-1) * bx(cell_x1-1,cell_y2-1,cell_z2-1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2-1,cell_z2-1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2-1,cell_z2-1)) &
            +           hy( 0) * (gx(-1) * bx(cell_x1-1,cell_y2  ,cell_z2-1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2  ,cell_z2-1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2  ,cell_z2-1)) &
            +           hy( 1) * (gx(-1) * bx(cell_x1-1,cell_y2+1,cell_z2-1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2+1,cell_z2-1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2+1,cell_z2-1))) &
            + hz( 0) * (hy(-1) * (gx(-1) * bx(cell_x1-1,cell_y2-1,cell_z2  ) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2-1,cell_z2  ) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2-1,cell_z2  )) &
            +           hy( 0) * (gx(-1) * bx(cell_x1-1,cell_y2  ,cell_z2  ) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2  ,cell_z2  ) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2  ,cell_z2  )) &
            +           hy( 1) * (gx(-1) * bx(cell_x1-1,cell_y2+1,cell_z2  ) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2+1,cell_z2  ) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2+1,cell_z2  ))) &
            + hz( 1) * (hy(-1) * (gx(-1) * bx(cell_x1-1,cell_y2-1,cell_z2+1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2-1,cell_z2+1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2-1,cell_z2+1)) &
            +           hy( 0) * (gx(-1) * bx(cell_x1-1,cell_y2  ,cell_z2+1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2  ,cell_z2+1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2  ,cell_z2+1)) &
            +           hy( 1) * (gx(-1) * bx(cell_x1-1,cell_y2+1,cell_z2+1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2+1,cell_z2+1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2+1,cell_z2+1)))

        by_part = &
              hz(-1) * (gy(-1) * (hx(-1) * by(cell_x2-1,cell_y1-1,cell_z2-1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1-1,cell_z2-1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1-1,cell_z2-1)) &
            +           gy( 0) * (hx(-1) * by(cell_x2-1,cell_y1  ,cell_z2-1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1  ,cell_z2-1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1  ,cell_z2-1)) &
            +           gy( 1) * (hx(-1) * by(cell_x2-1,cell_y1+1,cell_z2-1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1+1,cell_z2-1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1+1,cell_z2-1))) &
            + hz( 0) * (gy(-1) * (hx(-1) * by(cell_x2-1,cell_y1-1,cell_z2  ) &
            +                     hx( 0) * by(cell_x2  ,cell_y1-1,cell_z2  ) &
            +                     hx( 1) * by(cell_x2+1,cell_y1-1,cell_z2  )) &
            +           gy( 0) * (hx(-1) * by(cell_x2-1,cell_y1  ,cell_z2  ) &
            +                     hx( 0) * by(cell_x2  ,cell_y1  ,cell_z2  ) &
            +                     hx( 1) * by(cell_x2+1,cell_y1  ,cell_z2  )) &
            +           gy( 1) * (hx(-1) * by(cell_x2-1,cell_y1+1,cell_z2  ) &
            +                     hx( 0) * by(cell_x2  ,cell_y1+1,cell_z2  ) &
            +                     hx( 1) * by(cell_x2+1,cell_y1+1,cell_z2  ))) &
            + hz( 1) * (gy(-1) * (hx(-1) * by(cell_x2-1,cell_y1-1,cell_z2+1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1-1,cell_z2+1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1-1,cell_z2+1)) &
            +           gy( 0) * (hx(-1) * by(cell_x2-1,cell_y1  ,cell_z2+1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1  ,cell_z2+1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1  ,cell_z2+1)) &
            +           gy( 1) * (hx(-1) * by(cell_x2-1,cell_y1+1,cell_z2+1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1+1,cell_z2+1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1+1,cell_z2+1)))

        bz_part = &
              gz(-1) * (hy(-1) * (hx(-1) * bz(cell_x2-1,cell_y2-1,cell_z1-1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2-1,cell_z1-1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2-1,cell_z1-1)) &
            +           hy( 0) * (hx(-1) * bz(cell_x2-1,cell_y2  ,cell_z1-1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2  ,cell_z1-1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2  ,cell_z1-1)) &
            +           hy( 1) * (hx(-1) * bz(cell_x2-1,cell_y2+1,cell_z1-1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2+1,cell_z1-1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2+1,cell_z1-1))) &
            + gz( 0) * (hy(-1) * (hx(-1) * bz(cell_x2-1,cell_y2-1,cell_z1  ) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2-1,cell_z1  ) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2-1,cell_z1  )) &
            +           hy( 0) * (hx(-1) * bz(cell_x2-1,cell_y2  ,cell_z1  ) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2  ,cell_z1  ) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2  ,cell_z1  )) &
            +           hy( 1) * (hx(-1) * bz(cell_x2-1,cell_y2+1,cell_z1  ) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2+1,cell_z1  ) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2+1,cell_z1  ))) &
            + gz( 1) * (hy(-1) * (hx(-1) * bz(cell_x2-1,cell_y2-1,cell_z1+1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2-1,cell_z1+1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2-1,cell_z1+1)) &
            +           hy( 0) * (hx(-1) * bz(cell_x2-1,cell_y2  ,cell_z1+1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2  ,cell_z1+1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2  ,cell_z1+1)) &
            +           hy( 1) * (hx(-1) * bz(cell_x2-1,cell_y2+1,cell_z1+1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2+1,cell_z1+1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2+1,cell_z1+1)))
