        cf2 = cell_frac_x**2
        hx(dcellx-1) = 0.25_num + cf2 + cell_frac_x
        hx(dcellx  ) = 1.5_num - 2.0_num * cf2
        hx(dcellx+1) = 0.25_num + cf2 - cell_frac_x
