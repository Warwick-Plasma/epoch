        cf2 = cell_frac_x**2
        hx(dcellx-1) = 0.25_num + cf2 + cell_frac_x
        hx(dcellx  ) = 1.5_num - 2.0_num * cf2
        hx(dcellx+1) = 0.25_num + cf2 - cell_frac_x

        cf2 = cell_frac_y**2
        hy(dcelly-1) = 0.25_num + cf2 + cell_frac_y
        hy(dcelly  ) = 1.5_num - 2.0_num * cf2
        hy(dcelly+1) = 0.25_num + cf2 - cell_frac_y

        cf2 = cell_frac_z**2
        hz(dcellz-1) = 0.25_num + cf2 + cell_frac_z
        hz(dcellz  ) = 1.5_num - 2.0_num * cf2
        hz(dcellz+1) = 0.25_num + cf2 - cell_frac_z
