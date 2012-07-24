! IMPORTANT NOTE
! These weight need to be multiplied by 1/24
        cf2 = cell_frac_x**2
        gx(-2) = (0.5_num + cell_frac_x)**4
        gx(-1) = 4.75_num + 11.0_num * cell_frac_x &
            + 4.0_num * cf2 * (1.5_num - cell_frac_x - cf2)
        gx( 0) = 14.375_num + 6.0_num * cf2 * (cf2 - 2.5_num)
        gx( 1) = 4.75_num - 11.0_num * cell_frac_x &
            + 4.0_num * cf2 * (1.5_num + cell_frac_x - cf2)
        gx( 2) = (0.5_num - cell_frac_x)**4

        cf2 = cell_frac_y**2
        gy(-2) = (0.5_num + cell_frac_y)**4
        gy(-1) = 4.75_num + 11.0_num * cell_frac_y &
            + 4.0_num * cf2 * (1.5_num - cell_frac_y - cf2)
        gy( 0) = 14.375_num + 6.0_num * cf2 * (cf2 - 2.5_num)
        gy( 1) = 4.75_num - 11.0_num * cell_frac_y &
            + 4.0_num * cf2 * (1.5_num + cell_frac_y - cf2)
        gy( 2) = (0.5_num - cell_frac_y)**4

        cf2 = cell_frac_z**2
        gz(-2) = (0.5_num + cell_frac_z)**4
        gz(-1) = 4.75_num + 11.0_num * cell_frac_z &
            + 4.0_num * cf2 * (1.5_num - cell_frac_z - cf2)
        gz( 0) = 14.375_num + 6.0_num * cf2 * (cf2 - 2.5_num)
        gz( 1) = 4.75_num - 11.0_num * cell_frac_z &
            + 4.0_num * cf2 * (1.5_num + cell_frac_z - cf2)
        gz( 2) = (0.5_num - cell_frac_z)**4
