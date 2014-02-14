function Lorenz = make_Lorenz()

Lorenz.ode.fn  = @Lorenz_ode;
Lorenz.fn      = @Lorenz_fun;
Lorenz.dfdx    = @Lorenz_dfdx;
Lorenz.dfdp    = @Lorenz_dfdp;
Lorenz.d2fdx2  = @Lorenz_d2fdx2;
Lorenz.d2fdxdp = @Lorenz_d2fdxdp;

end
