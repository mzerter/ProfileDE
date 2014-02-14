function chemo0 = make_chemo0()

chemo0.fn      = @chemo0_fun;
chemo0.fn_ode  = @chemo0_fun_ode;
chemo0.dfdx    = @chemo0_dfdx;
chemo0.dfdp    = @chemo0_dfdp;
chemo0.d2fdx2  = @chemo0_d2fdx2;
chemo0.d2fdxdp = @chemo0_d2fdxdp;

end

