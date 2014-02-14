function SEIR = make_SEIR()

SEIR.ode.fn  = @SEIR_ode;
SEIR.fn      = @SEIR_fun;
SEIR.dfdx    = @SEIR_dfdx;
SEIR.dfdp    = @SEIR_dfdp;
SEIR.d2fdx2  = @SEIR_d2fdx2;
SEIR.d2fdxdp = @SEIR_d2fdxdp;

end
