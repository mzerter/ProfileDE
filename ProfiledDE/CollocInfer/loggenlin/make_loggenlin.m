function loggenlin = make_loggenlin()

loggenlin.fun_ode = @loggenlin_fun_ode;
loggenlin.fn      = @loggenlin_fn;
loggenlin.dfdx    = @loggenlin_dfdx;
loggenlin.dfdp    = @loggenlin_dfdp;
loggenlin.d2fdx2  = @loggenlin_d2fdx2;
loggenlin.d2fdxdp = @loggenlin_d2fdxdp;

end
