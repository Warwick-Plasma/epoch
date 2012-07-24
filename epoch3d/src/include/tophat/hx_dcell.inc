! IMPORTANT NOTE
! These weight need to be multiplied by 1/2
        hx(dcellx  ) = 0.5_num + cell_frac_x
        hx(dcellx+1) = 0.5_num - cell_frac_x

        hy(dcelly  ) = 0.5_num + cell_frac_y
        hy(dcelly+1) = 0.5_num - cell_frac_y

        hz(dcellz  ) = 0.5_num + cell_frac_z
        hz(dcellz+1) = 0.5_num - cell_frac_z
