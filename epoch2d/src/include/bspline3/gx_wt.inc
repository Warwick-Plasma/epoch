! IMPORTANT NOTE
! These weight need to be multiplied by 1/24
        wt_var1 = wt_facx * (wt_dtx + cell_frac_x)**4
        wt_var2 = wt_facx * (wt_dtx - cell_frac_x)**4
        wt_var3 = SIGN(wt_var1, wt_dtx + cell_frac_x) &
            + SIGN(wt_var2, wt_dtx - cell_frac_x)
        wt_var4 = ((3.0_num - cell_frac_x) * cell_frac_x + 3.0_num - wt_dtx2) &
            * cell_frac_x + 1.0_num + wt_dtx2
        wt_var5 = ((3.0_num + cell_frac_x) * cell_frac_x - 3.0_num + wt_dtx2) &
            * cell_frac_x + 1.0_num + wt_dtx2
        gx(-2) = wt_var3 + wt_var1 - wt_var2
        gx(-1) = 4.0_num * (wt_var4 - wt_var3)
        gx( 0) = 24.0_num + 6.0_num * wt_var3 - 4.0_num * (wt_var4 + wt_var5)
        gx( 1) = 4.0_num * (wt_var5 - wt_var3)
        gx( 2) = wt_var3 + wt_var2 - wt_var1

        wt_var1 = wt_facy * (wt_dty + cell_frac_y)**4
        wt_var2 = wt_facy * (wt_dty - cell_frac_y)**4
        wt_var3 = SIGN(wt_var1, wt_dty + cell_frac_y) &
            + SIGN(wt_var2, wt_dty - cell_frac_y)
        wt_var4 = ((3.0_num - cell_frac_y) * cell_frac_y + 3.0_num - wt_dty2) &
            * cell_frac_y + 1.0_num + wt_dty2
        wt_var5 = ((3.0_num + cell_frac_y) * cell_frac_y - 3.0_num + wt_dty2) &
            * cell_frac_y + 1.0_num + wt_dty2
        gy(-2) = wt_var3 + wt_var1 - wt_var2
        gy(-1) = 4.0_num * (wt_var4 - wt_var3)
        gy( 0) = 24.0_num + 6.0_num * wt_var3 - 4.0_num * (wt_var4 + wt_var5)
        gy( 1) = 4.0_num * (wt_var5 - wt_var3)
        gy( 2) = wt_var3 + wt_var2 - wt_var1
