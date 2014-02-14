function SIR = make_SIR_dynamic()

SIR.ode.fn  = @SIR_dynamic_ode;
SIR.fn      = @SIR_dynamic_fun;
SIR.dfdx    = @SIR_dynamic_dfdx;
SIR.dfdp    = @SIR_dynamic_dfdp;
SIR.d2fdx2  = @SIR_dynamic_d2fdx2;
SIR.d2fdxdp = @SIR_dynamic_d2fdxdp;

end
