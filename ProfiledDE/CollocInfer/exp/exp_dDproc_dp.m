function dpval = exp_dDproc_dp(coefs,bvals,pars,more)

[m,n] = size(bvals);
devals  = exp(bvals(1:m-1,:) * coefs);
ddevals = exp(bvals(2:m,:)   * coefs);
% colnames(devals) = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
g = more.dfdp(ddevals,more.qpts,devals,pars,more.more);
dpval = applysum(g);

end

