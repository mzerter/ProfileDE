function findif_loglik = make_findif_loglik()

findif_loglik.fn      = @findif_loglik_fun;
findif_loglik.dfdx    = @findif_loglik_dfdx;
findif_loglik.dfdy    = @findif_loglik_dfdy;
findif_loglik.dfdp    = @findif_loglik_dfdp;
findif_loglik.d2fdx2  = @findif_loglik_d2fdx2;
findif_loglik.d2fdy2  = @findif_loglik_d2fdy2;
findif_loglik.d2fdxdy = @findif_loglik_d2fdxdy;
findif_loglik.d2fdxdp = @findif_loglik_d2fdxdp;
findif_loglik.d2fdydp = @findif_loglik_d2fdydp;

end

