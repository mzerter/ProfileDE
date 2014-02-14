function fnval = exp_Cproc(coefs,bvals,pars,more)

devals  = exp(bvals.bvals * coefs);
ddevals = (bvals.dbvals * coefs)*devals;
% colnames(devals) = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
f = more.fn(ddevals,more.qpts,devals,pars,more.more);
fnval = sum(f);

end

