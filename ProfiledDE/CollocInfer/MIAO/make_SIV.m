function SIV = make_SIV()

SIV.ode.fn  = @SIV_ode;
SIV.fn      = @SIV_fun;
SIV.dfdx    = @SIV_dfdx;
SIV.dfdp    = @SIV_dfdp;
SIV.d2fdx2  = @SIV_d2fdx2;
SIV.d2fdxdp = @SIV_d2fdxdp;

end
