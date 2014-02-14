function genlin = make_genlin()

genlin.fun_ode = @genlin_fun_ode;
genlin.fn      = @genlin_fn;
genlin.dfdx    = @genlin_dfdx;
genlin.dfdp    = @genlin_dfdp;
genlin.d2fdx2  = @genlin_d2fdx2;
genlin.d2fdxdp = @genlin_d2fdxdp;

end
