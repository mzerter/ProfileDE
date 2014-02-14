function SIR = make_SIR()

SIR.ode.fn  = @SIR_ode;
SIR.fn      = @SIR_fun;
SIR.dfdx    = @SIR_dfdx;
SIR.dfdp    = @SIR_dfdp;
SIR.d2fdx2  = @SIR_d2fdx2;
SIR.d2fdxdp = @SIR_d2fdxdp;

end
