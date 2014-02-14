function fhn = make_fhn()

fhn.fn      = @fhn_fun;
fhn.fn_ode  = @fhn_fun_ode;
fhn.dfdx    = @fhn_dfdx;
fhn.dfdp    = @fhn_dfdp;
fhn.d2fdx2  = @fhn_d2fdx2;
fhn.d2fdxdp = @fhn_d2fdxdp;

end
