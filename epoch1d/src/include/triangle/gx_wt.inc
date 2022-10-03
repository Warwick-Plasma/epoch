! IMPORTANT NOTE
! These weight need to be multiplied by 1/2
        wt_var1 = wt_facx * ( &
            ABS(wt_dtx + cell_frac_x) * (wt_dtx + cell_frac_x) &
            + ABS(wt_dtx - cell_frac_x) * (wt_dtx - cell_frac_x))
        gx(-1) = wt_var1 + cell_frac_x
        gx( 0) = 2.0_num - 2.0_num * wt_var1
        gx( 1) = wt_var1 - cell_frac_x
