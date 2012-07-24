        ex_part = &
              hx( 0) * ex(cell_x2  ) &
            + hx( 1) * ex(cell_x2+1)

        ey_part = &
              gx( 0) * ey(cell_x1  ) &
            + gx( 1) * ey(cell_x1+1)

        ez_part = &
              gx( 0) * ez(cell_x1  ) &
            + gx( 1) * ez(cell_x1+1)
