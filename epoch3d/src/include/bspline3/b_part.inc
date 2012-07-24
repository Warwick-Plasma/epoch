        bx_part = &
              hz(-2) * (hy(-2) * (gx(-2) * bx(cell_x1-2,cell_y2-2,cell_z2-2) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2-2,cell_z2-2) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2-2,cell_z2-2) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2-2,cell_z2-2) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2-2,cell_z2-2)) &
            +           hy(-1) * (gx(-2) * bx(cell_x1-2,cell_y2-1,cell_z2-2) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2-1,cell_z2-2) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2-1,cell_z2-2) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2-1,cell_z2-2) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2-1,cell_z2-2)) &
            +           hy( 0) * (gx(-2) * bx(cell_x1-2,cell_y2  ,cell_z2-2) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2  ,cell_z2-2) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2  ,cell_z2-2) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2  ,cell_z2-2) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2  ,cell_z2-2)) &
            +           hy( 1) * (gx(-2) * bx(cell_x1-2,cell_y2+1,cell_z2-2) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2+1,cell_z2-2) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2+1,cell_z2-2) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2+1,cell_z2-2) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2+1,cell_z2-2)) &
            +           hy( 2) * (gx(-2) * bx(cell_x1-2,cell_y2+2,cell_z2-2) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2+2,cell_z2-2) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2+2,cell_z2-2) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2+2,cell_z2-2) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2+2,cell_z2-2)))
        bx_part = bx_part &
            + hz(-1) * (hy(-2) * (gx(-2) * bx(cell_x1-2,cell_y2-2,cell_z2-1) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2-2,cell_z2-1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2-2,cell_z2-1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2-2,cell_z2-1) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2-2,cell_z2-1)) &
            +           hy(-1) * (gx(-2) * bx(cell_x1-2,cell_y2-1,cell_z2-1) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2-1,cell_z2-1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2-1,cell_z2-1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2-1,cell_z2-1) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2-1,cell_z2-1)) &
            +           hy( 0) * (gx(-2) * bx(cell_x1-2,cell_y2  ,cell_z2-1) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2  ,cell_z2-1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2  ,cell_z2-1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2  ,cell_z2-1) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2  ,cell_z2-1)) &
            +           hy( 1) * (gx(-2) * bx(cell_x1-2,cell_y2+1,cell_z2-1) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2+1,cell_z2-1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2+1,cell_z2-1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2+1,cell_z2-1) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2+1,cell_z2-1)) &
            +           hy( 2) * (gx(-2) * bx(cell_x1-2,cell_y2+2,cell_z2-1) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2+2,cell_z2-1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2+2,cell_z2-1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2+2,cell_z2-1) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2+2,cell_z2-1)))
        bx_part = bx_part &
            + hz( 0) * (hy(-2) * (gx(-2) * bx(cell_x1-2,cell_y2-2,cell_z2  ) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2-2,cell_z2  ) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2-2,cell_z2  ) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2-2,cell_z2  ) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2-2,cell_z2  )) &
            +           hy(-1) * (gx(-2) * bx(cell_x1-2,cell_y2-1,cell_z2  ) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2-1,cell_z2  ) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2-1,cell_z2  ) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2-1,cell_z2  ) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2-1,cell_z2  )) &
            +           hy( 0) * (gx(-2) * bx(cell_x1-2,cell_y2  ,cell_z2  ) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2  ,cell_z2  ) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2  ,cell_z2  ) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2  ,cell_z2  ) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2  ,cell_z2  )) &
            +           hy( 1) * (gx(-2) * bx(cell_x1-2,cell_y2+1,cell_z2  ) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2+1,cell_z2  ) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2+1,cell_z2  ) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2+1,cell_z2  ) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2+1,cell_z2  )) &
            +           hy( 2) * (gx(-2) * bx(cell_x1-2,cell_y2+2,cell_z2  ) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2+2,cell_z2  ) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2+2,cell_z2  ) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2+2,cell_z2  ) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2+2,cell_z2  )))
        bx_part = bx_part &
            + hz( 1) * (hy(-2) * (gx(-2) * bx(cell_x1-2,cell_y2-2,cell_z2+1) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2-2,cell_z2+1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2-2,cell_z2+1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2-2,cell_z2+1) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2-2,cell_z2+1)) &
            +           hy(-1) * (gx(-2) * bx(cell_x1-2,cell_y2-1,cell_z2+1) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2-1,cell_z2+1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2-1,cell_z2+1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2-1,cell_z2+1) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2-1,cell_z2+1)) &
            +           hy( 0) * (gx(-2) * bx(cell_x1-2,cell_y2  ,cell_z2+1) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2  ,cell_z2+1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2  ,cell_z2+1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2  ,cell_z2+1) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2  ,cell_z2+1)) &
            +           hy( 1) * (gx(-2) * bx(cell_x1-2,cell_y2+1,cell_z2+1) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2+1,cell_z2+1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2+1,cell_z2+1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2+1,cell_z2+1) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2+1,cell_z2+1)) &
            +           hy( 2) * (gx(-2) * bx(cell_x1-2,cell_y2+2,cell_z2+1) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2+2,cell_z2+1) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2+2,cell_z2+1) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2+2,cell_z2+1) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2+2,cell_z2+1)))
        bx_part = bx_part &
            + hz( 2) * (hy(-2) * (gx(-2) * bx(cell_x1-2,cell_y2-2,cell_z2+2) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2-2,cell_z2+2) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2-2,cell_z2+2) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2-2,cell_z2+2) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2-2,cell_z2+2)) &
            +           hy(-1) * (gx(-2) * bx(cell_x1-2,cell_y2-1,cell_z2+2) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2-1,cell_z2+2) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2-1,cell_z2+2) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2-1,cell_z2+2) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2-1,cell_z2+2)) &
            +           hy( 0) * (gx(-2) * bx(cell_x1-2,cell_y2  ,cell_z2+2) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2  ,cell_z2+2) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2  ,cell_z2+2) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2  ,cell_z2+2) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2  ,cell_z2+2)) &
            +           hy( 1) * (gx(-2) * bx(cell_x1-2,cell_y2+1,cell_z2+2) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2+1,cell_z2+2) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2+1,cell_z2+2) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2+1,cell_z2+2) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2+1,cell_z2+2)) &
            +           hy( 2) * (gx(-2) * bx(cell_x1-2,cell_y2+2,cell_z2+2) &
            +                     gx(-1) * bx(cell_x1-1,cell_y2+2,cell_z2+2) &
            +                     gx( 0) * bx(cell_x1  ,cell_y2+2,cell_z2+2) &
            +                     gx( 1) * bx(cell_x1+1,cell_y2+2,cell_z2+2) &
            +                     gx( 2) * bx(cell_x1+2,cell_y2+2,cell_z2+2)))

        by_part = &
              hz(-2) * (gy(-2) * (hx(-2) * by(cell_x2-2,cell_y1-2,cell_z2-2) &
            +                     hx(-1) * by(cell_x2-1,cell_y1-2,cell_z2-2) &
            +                     hx( 0) * by(cell_x2  ,cell_y1-2,cell_z2-2) &
            +                     hx( 1) * by(cell_x2+1,cell_y1-2,cell_z2-2) &
            +                     hx( 2) * by(cell_x2+2,cell_y1-2,cell_z2-2)) &
            +           gy(-1) * (hx(-2) * by(cell_x2-2,cell_y1-1,cell_z2-2) &
            +                     hx(-1) * by(cell_x2-1,cell_y1-1,cell_z2-2) &
            +                     hx( 0) * by(cell_x2  ,cell_y1-1,cell_z2-2) &
            +                     hx( 1) * by(cell_x2+1,cell_y1-1,cell_z2-2) &
            +                     hx( 2) * by(cell_x2+2,cell_y1-1,cell_z2-2)) &
            +           gy( 0) * (hx(-2) * by(cell_x2-2,cell_y1  ,cell_z2-2) &
            +                     hx(-1) * by(cell_x2-1,cell_y1  ,cell_z2-2) &
            +                     hx( 0) * by(cell_x2  ,cell_y1  ,cell_z2-2) &
            +                     hx( 1) * by(cell_x2+1,cell_y1  ,cell_z2-2) &
            +                     hx( 2) * by(cell_x2+2,cell_y1  ,cell_z2-2)) &
            +           gy( 1) * (hx(-2) * by(cell_x2-2,cell_y1+1,cell_z2-2) &
            +                     hx(-1) * by(cell_x2-1,cell_y1+1,cell_z2-2) &
            +                     hx( 0) * by(cell_x2  ,cell_y1+1,cell_z2-2) &
            +                     hx( 1) * by(cell_x2+1,cell_y1+1,cell_z2-2) &
            +                     hx( 2) * by(cell_x2+2,cell_y1+1,cell_z2-2)) &
            +           gy( 2) * (hx(-2) * by(cell_x2-2,cell_y1+2,cell_z2-2) &
            +                     hx(-1) * by(cell_x2-1,cell_y1+2,cell_z2-2) &
            +                     hx( 0) * by(cell_x2  ,cell_y1+2,cell_z2-2) &
            +                     hx( 1) * by(cell_x2+1,cell_y1+2,cell_z2-2) &
            +                     hx( 2) * by(cell_x2+2,cell_y1+2,cell_z2-2)))
        by_part = by_part &
            + hz(-1) * (gy(-2) * (hx(-2) * by(cell_x2-2,cell_y1-2,cell_z2-1) &
            +                     hx(-1) * by(cell_x2-1,cell_y1-2,cell_z2-1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1-2,cell_z2-1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1-2,cell_z2-1) &
            +                     hx( 2) * by(cell_x2+2,cell_y1-2,cell_z2-1)) &
            +           gy(-1) * (hx(-2) * by(cell_x2-2,cell_y1-1,cell_z2-1) &
            +                     hx(-1) * by(cell_x2-1,cell_y1-1,cell_z2-1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1-1,cell_z2-1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1-1,cell_z2-1) &
            +                     hx( 2) * by(cell_x2+2,cell_y1-1,cell_z2-1)) &
            +           gy( 0) * (hx(-2) * by(cell_x2-2,cell_y1  ,cell_z2-1) &
            +                     hx(-1) * by(cell_x2-1,cell_y1  ,cell_z2-1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1  ,cell_z2-1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1  ,cell_z2-1) &
            +                     hx( 2) * by(cell_x2+2,cell_y1  ,cell_z2-1)) &
            +           gy( 1) * (hx(-2) * by(cell_x2-2,cell_y1+1,cell_z2-1) &
            +                     hx(-1) * by(cell_x2-1,cell_y1+1,cell_z2-1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1+1,cell_z2-1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1+1,cell_z2-1) &
            +                     hx( 2) * by(cell_x2+2,cell_y1+1,cell_z2-1)) &
            +           gy( 2) * (hx(-2) * by(cell_x2-2,cell_y1+2,cell_z2-1) &
            +                     hx(-1) * by(cell_x2-1,cell_y1+2,cell_z2-1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1+2,cell_z2-1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1+2,cell_z2-1) &
            +                     hx( 2) * by(cell_x2+2,cell_y1+2,cell_z2-1)))
        by_part = by_part &
            + hz( 0) * (gy(-2) * (hx(-2) * by(cell_x2-2,cell_y1-2,cell_z2  ) &
            +                     hx(-1) * by(cell_x2-1,cell_y1-2,cell_z2  ) &
            +                     hx( 0) * by(cell_x2  ,cell_y1-2,cell_z2  ) &
            +                     hx( 1) * by(cell_x2+1,cell_y1-2,cell_z2  ) &
            +                     hx( 2) * by(cell_x2+2,cell_y1-2,cell_z2  )) &
            +           gy(-1) * (hx(-2) * by(cell_x2-2,cell_y1-1,cell_z2  ) &
            +                     hx(-1) * by(cell_x2-1,cell_y1-1,cell_z2  ) &
            +                     hx( 0) * by(cell_x2  ,cell_y1-1,cell_z2  ) &
            +                     hx( 1) * by(cell_x2+1,cell_y1-1,cell_z2  ) &
            +                     hx( 2) * by(cell_x2+2,cell_y1-1,cell_z2  )) &
            +           gy( 0) * (hx(-2) * by(cell_x2-2,cell_y1  ,cell_z2  ) &
            +                     hx(-1) * by(cell_x2-1,cell_y1  ,cell_z2  ) &
            +                     hx( 0) * by(cell_x2  ,cell_y1  ,cell_z2  ) &
            +                     hx( 1) * by(cell_x2+1,cell_y1  ,cell_z2  ) &
            +                     hx( 2) * by(cell_x2+2,cell_y1  ,cell_z2  )) &
            +           gy( 1) * (hx(-2) * by(cell_x2-2,cell_y1+1,cell_z2  ) &
            +                     hx(-1) * by(cell_x2-1,cell_y1+1,cell_z2  ) &
            +                     hx( 0) * by(cell_x2  ,cell_y1+1,cell_z2  ) &
            +                     hx( 1) * by(cell_x2+1,cell_y1+1,cell_z2  ) &
            +                     hx( 2) * by(cell_x2+2,cell_y1+1,cell_z2  )) &
            +           gy( 2) * (hx(-2) * by(cell_x2-2,cell_y1+2,cell_z2  ) &
            +                     hx(-1) * by(cell_x2-1,cell_y1+2,cell_z2  ) &
            +                     hx( 0) * by(cell_x2  ,cell_y1+2,cell_z2  ) &
            +                     hx( 1) * by(cell_x2+1,cell_y1+2,cell_z2  ) &
            +                     hx( 2) * by(cell_x2+2,cell_y1+2,cell_z2  )))
        by_part = by_part &
            + hz( 1) * (gy(-2) * (hx(-2) * by(cell_x2-2,cell_y1-2,cell_z2+1) &
            +                     hx(-1) * by(cell_x2-1,cell_y1-2,cell_z2+1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1-2,cell_z2+1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1-2,cell_z2+1) &
            +                     hx( 2) * by(cell_x2+2,cell_y1-2,cell_z2+1)) &
            +           gy(-1) * (hx(-2) * by(cell_x2-2,cell_y1-1,cell_z2+1) &
            +                     hx(-1) * by(cell_x2-1,cell_y1-1,cell_z2+1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1-1,cell_z2+1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1-1,cell_z2+1) &
            +                     hx( 2) * by(cell_x2+2,cell_y1-1,cell_z2+1)) &
            +           gy( 0) * (hx(-2) * by(cell_x2-2,cell_y1  ,cell_z2+1) &
            +                     hx(-1) * by(cell_x2-1,cell_y1  ,cell_z2+1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1  ,cell_z2+1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1  ,cell_z2+1) &
            +                     hx( 2) * by(cell_x2+2,cell_y1  ,cell_z2+1)) &
            +           gy( 1) * (hx(-2) * by(cell_x2-2,cell_y1+1,cell_z2+1) &
            +                     hx(-1) * by(cell_x2-1,cell_y1+1,cell_z2+1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1+1,cell_z2+1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1+1,cell_z2+1) &
            +                     hx( 2) * by(cell_x2+2,cell_y1+1,cell_z2+1)) &
            +           gy( 2) * (hx(-2) * by(cell_x2-2,cell_y1+2,cell_z2+1) &
            +                     hx(-1) * by(cell_x2-1,cell_y1+2,cell_z2+1) &
            +                     hx( 0) * by(cell_x2  ,cell_y1+2,cell_z2+1) &
            +                     hx( 1) * by(cell_x2+1,cell_y1+2,cell_z2+1) &
            +                     hx( 2) * by(cell_x2+2,cell_y1+2,cell_z2+1)))
        by_part = by_part &
            + hz( 2) * (gy(-2) * (hx(-2) * by(cell_x2-2,cell_y1-2,cell_z2+2) &
            +                     hx(-1) * by(cell_x2-1,cell_y1-2,cell_z2+2) &
            +                     hx( 0) * by(cell_x2  ,cell_y1-2,cell_z2+2) &
            +                     hx( 1) * by(cell_x2+1,cell_y1-2,cell_z2+2) &
            +                     hx( 2) * by(cell_x2+2,cell_y1-2,cell_z2+2)) &
            +           gy(-1) * (hx(-2) * by(cell_x2-2,cell_y1-1,cell_z2+2) &
            +                     hx(-1) * by(cell_x2-1,cell_y1-1,cell_z2+2) &
            +                     hx( 0) * by(cell_x2  ,cell_y1-1,cell_z2+2) &
            +                     hx( 1) * by(cell_x2+1,cell_y1-1,cell_z2+2) &
            +                     hx( 2) * by(cell_x2+2,cell_y1-1,cell_z2+2)) &
            +           gy( 0) * (hx(-2) * by(cell_x2-2,cell_y1  ,cell_z2+2) &
            +                     hx(-1) * by(cell_x2-1,cell_y1  ,cell_z2+2) &
            +                     hx( 0) * by(cell_x2  ,cell_y1  ,cell_z2+2) &
            +                     hx( 1) * by(cell_x2+1,cell_y1  ,cell_z2+2) &
            +                     hx( 2) * by(cell_x2+2,cell_y1  ,cell_z2+2)) &
            +           gy( 1) * (hx(-2) * by(cell_x2-2,cell_y1+1,cell_z2+2) &
            +                     hx(-1) * by(cell_x2-1,cell_y1+1,cell_z2+2) &
            +                     hx( 0) * by(cell_x2  ,cell_y1+1,cell_z2+2) &
            +                     hx( 1) * by(cell_x2+1,cell_y1+1,cell_z2+2) &
            +                     hx( 2) * by(cell_x2+2,cell_y1+1,cell_z2+2)) &
            +           gy( 2) * (hx(-2) * by(cell_x2-2,cell_y1+2,cell_z2+2) &
            +                     hx(-1) * by(cell_x2-1,cell_y1+2,cell_z2+2) &
            +                     hx( 0) * by(cell_x2  ,cell_y1+2,cell_z2+2) &
            +                     hx( 1) * by(cell_x2+1,cell_y1+2,cell_z2+2) &
            +                     hx( 2) * by(cell_x2+2,cell_y1+2,cell_z2+2)))

        bz_part = &
              gz(-2) * (hy(-2) * (hx(-2) * bz(cell_x2-2,cell_y2-2,cell_z1-2) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2-2,cell_z1-2) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2-2,cell_z1-2) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2-2,cell_z1-2) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2-2,cell_z1-2)) &
            +           hy(-1) * (hx(-2) * bz(cell_x2-2,cell_y2-1,cell_z1-2) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2-1,cell_z1-2) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2-1,cell_z1-2) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2-1,cell_z1-2) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2-1,cell_z1-2)) &
            +           hy( 0) * (hx(-2) * bz(cell_x2-2,cell_y2  ,cell_z1-2) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2  ,cell_z1-2) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2  ,cell_z1-2) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2  ,cell_z1-2) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2  ,cell_z1-2)) &
            +           hy( 1) * (hx(-2) * bz(cell_x2-2,cell_y2+1,cell_z1-2) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2+1,cell_z1-2) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2+1,cell_z1-2) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2+1,cell_z1-2) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2+1,cell_z1-2)) &
            +           hy( 2) * (hx(-2) * bz(cell_x2-2,cell_y2+2,cell_z1-2) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2+2,cell_z1-2) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2+2,cell_z1-2) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2+2,cell_z1-2) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2+2,cell_z1-2)))
        bz_part = bz_part &
            + gz(-1) * (hy(-2) * (hx(-2) * bz(cell_x2-2,cell_y2-2,cell_z1-1) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2-2,cell_z1-1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2-2,cell_z1-1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2-2,cell_z1-1) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2-2,cell_z1-1)) &
            +           hy(-1) * (hx(-2) * bz(cell_x2-2,cell_y2-1,cell_z1-1) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2-1,cell_z1-1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2-1,cell_z1-1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2-1,cell_z1-1) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2-1,cell_z1-1)) &
            +           hy( 0) * (hx(-2) * bz(cell_x2-2,cell_y2  ,cell_z1-1) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2  ,cell_z1-1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2  ,cell_z1-1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2  ,cell_z1-1) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2  ,cell_z1-1)) &
            +           hy( 1) * (hx(-2) * bz(cell_x2-2,cell_y2+1,cell_z1-1) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2+1,cell_z1-1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2+1,cell_z1-1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2+1,cell_z1-1) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2+1,cell_z1-1)) &
            +           hy( 2) * (hx(-2) * bz(cell_x2-2,cell_y2+2,cell_z1-1) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2+2,cell_z1-1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2+2,cell_z1-1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2+2,cell_z1-1) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2+2,cell_z1-1)))
        bz_part = bz_part &
            + gz( 0) * (hy(-2) * (hx(-2) * bz(cell_x2-2,cell_y2-2,cell_z1  ) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2-2,cell_z1  ) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2-2,cell_z1  ) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2-2,cell_z1  ) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2-2,cell_z1  )) &
            +           hy(-1) * (hx(-2) * bz(cell_x2-2,cell_y2-1,cell_z1  ) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2-1,cell_z1  ) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2-1,cell_z1  ) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2-1,cell_z1  ) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2-1,cell_z1  )) &
            +           hy( 0) * (hx(-2) * bz(cell_x2-2,cell_y2  ,cell_z1  ) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2  ,cell_z1  ) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2  ,cell_z1  ) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2  ,cell_z1  ) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2  ,cell_z1  )) &
            +           hy( 1) * (hx(-2) * bz(cell_x2-2,cell_y2+1,cell_z1  ) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2+1,cell_z1  ) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2+1,cell_z1  ) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2+1,cell_z1  ) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2+1,cell_z1  )) &
            +           hy( 2) * (hx(-2) * bz(cell_x2-2,cell_y2+2,cell_z1  ) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2+2,cell_z1  ) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2+2,cell_z1  ) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2+2,cell_z1  ) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2+2,cell_z1  )))
        bz_part = bz_part &
            + gz( 1) * (hy(-2) * (hx(-2) * bz(cell_x2-2,cell_y2-2,cell_z1+1) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2-2,cell_z1+1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2-2,cell_z1+1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2-2,cell_z1+1) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2-2,cell_z1+1)) &
            +           hy(-1) * (hx(-2) * bz(cell_x2-2,cell_y2-1,cell_z1+1) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2-1,cell_z1+1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2-1,cell_z1+1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2-1,cell_z1+1) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2-1,cell_z1+1)) &
            +           hy( 0) * (hx(-2) * bz(cell_x2-2,cell_y2  ,cell_z1+1) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2  ,cell_z1+1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2  ,cell_z1+1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2  ,cell_z1+1) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2  ,cell_z1+1)) &
            +           hy( 1) * (hx(-2) * bz(cell_x2-2,cell_y2+1,cell_z1+1) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2+1,cell_z1+1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2+1,cell_z1+1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2+1,cell_z1+1) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2+1,cell_z1+1)) &
            +           hy( 2) * (hx(-2) * bz(cell_x2-2,cell_y2+2,cell_z1+1) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2+2,cell_z1+1) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2+2,cell_z1+1) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2+2,cell_z1+1) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2+2,cell_z1+1)))
        bz_part = bz_part &
            + gz( 2) * (hy(-2) * (hx(-2) * bz(cell_x2-2,cell_y2-2,cell_z1+2) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2-2,cell_z1+2) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2-2,cell_z1+2) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2-2,cell_z1+2) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2-2,cell_z1+2)) &
            +           hy(-1) * (hx(-2) * bz(cell_x2-2,cell_y2-1,cell_z1+2) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2-1,cell_z1+2) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2-1,cell_z1+2) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2-1,cell_z1+2) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2-1,cell_z1+2)) &
            +           hy( 0) * (hx(-2) * bz(cell_x2-2,cell_y2  ,cell_z1+2) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2  ,cell_z1+2) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2  ,cell_z1+2) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2  ,cell_z1+2) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2  ,cell_z1+2)) &
            +           hy( 1) * (hx(-2) * bz(cell_x2-2,cell_y2+1,cell_z1+2) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2+1,cell_z1+2) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2+1,cell_z1+2) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2+1,cell_z1+2) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2+1,cell_z1+2)) &
            +           hy( 2) * (hx(-2) * bz(cell_x2-2,cell_y2+2,cell_z1+2) &
            +                     hx(-1) * bz(cell_x2-1,cell_y2+2,cell_z1+2) &
            +                     hx( 0) * bz(cell_x2  ,cell_y2+2,cell_z1+2) &
            +                     hx( 1) * bz(cell_x2+1,cell_y2+2,cell_z1+2) &
            +                     hx( 2) * bz(cell_x2+2,cell_y2+2,cell_z1+2)))
