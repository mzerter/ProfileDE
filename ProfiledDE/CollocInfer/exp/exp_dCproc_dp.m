function dpval = exp_dCproc_dp(coefs,bvals,pars,more)

devals  = exp(bvals.bvals * coefs);
ddevals = (bvals.dbvals * coefs)*devals;
% colnames(devals) = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
g = more.dfdp(ddevals,more.qpts,devals,pars,more.more);
dpval = applysum(g);

end

