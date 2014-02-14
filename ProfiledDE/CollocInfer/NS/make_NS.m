function NS = make_NS()
NS.ode.fn  = @NS_ode;
NS.fn      = @NS_fun;
NS.dfdx    = @NS_dfdx;
NS.dfdp    = @NS_dfdp;
NS.d2fdx2  = @NS_d2fdx2;
NS.d2fdxdp = @NS_d2fdxdp;
end
