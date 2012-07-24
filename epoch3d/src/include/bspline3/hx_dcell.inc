! IMPORTANT NOTE
! These weight need to be multiplied by 1/24
        cf2 = cell_frac_x**2
        hx(dcellx-2) = (0.5_num + cell_frac_x)**4
        hx(dcellx-1) = 4.75_num + 11.0_num * cell_frac_x &
            + 4.0_num * cf2 * (1.5_num - cell_frac_x - cf2)
        hx(dcellx  ) = 14.375_num + 6.0_num * cf2 * (cf2 - 2.5_num)
        hx(dcellx+1) = 4.75_num - 11.0_num * cell_frac_x &
            + 4.0_num * cf2 * (1.5_num + cell_frac_x - cf2)
        hx(dcellx+2) = (0.5_num - cell_frac_x)**4

        cf2 = cell_frac_y**2
        hy(dcelly-2) = (0.5_num + cell_frac_y)**4
        hy(dcelly-1) = 4.75_num + 11.0_num * cell_frac_y &
            + 4.0_num * cf2 * (1.5_num - cell_frac_y - cf2)
        hy(dcelly  ) = 14.375_num + 6.0_num * cf2 * (cf2 - 2.5_num)
        hy(dcelly+1) = 4.75_num - 11.0_num * cell_frac_y &
            + 4.0_num * cf2 * (1.5_num + cell_frac_y - cf2)
        hy(dcelly+2) = (0.5_num - cell_frac_y)**4

        cf2 = cell_frac_z**2
        hz(dcellz-2) = (0.5_num + cell_frac_z)**4
        hz(dcellz-1) = 4.75_num + 11.0_num * cell_frac_z &
            + 4.0_num * cf2 * (1.5_num - cell_frac_z - cf2)
        hz(dcellz  ) = 14.375_num + 6.0_num * cf2 * (cf2 - 2.5_num)
        hz(dcellz+1) = 4.75_num - 11.0_num * cell_frac_z &
            + 4.0_num * cf2 * (1.5_num + cell_frac_z - cf2)
        hz(dcellz+2) = (0.5_num - cell_frac_z)**4
