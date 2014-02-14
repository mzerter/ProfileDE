%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite differencing for ODE RHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function findif_ode = make_findif_ode()

findif_ode.fn      = @findif_ode_fn;
findif_ode.dfdx    = @findif_ode_dfdx;
findif_ode.dfdy    = @findif_ode_dfdy;
findif_ode.dfdp    = @findif_ode_dfdp;
findif_ode.d2fdx2  = @findif_ode_d2fdx2;
findif_ode.d2fdxdp = @findif_ode_d2fdxdp;

end


