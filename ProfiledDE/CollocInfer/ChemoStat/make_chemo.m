function chemo = make_chemo()

chemo.fn      = @chemo_fun;
chemo.fn_ode  = @chemo_fun_ode;
chemo.dfdx    = @chemo_dfdx;
chemo.dfdp    = @chemo_dfdp;
chemo.d2fdx2  = @chemo_d2fdx2;
chemo.d2fdxdp = @chemo_d2fdxdp;

end

