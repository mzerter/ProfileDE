function logstate_lik = make_logstate_lik()

logstate_lik.fn      = @logstate_lik_fun;
logstate_lik.dfdx    = @logstate_lik_dfdx;
logstate_lik.dfdp    = @logstate_lik_dfdp;
logstate_lik.d2fdx2  = @logstate_lik_d2fdx2;
logstate_lik.d2fdxdp = @logstate_lik_d2fdxdp;

end
