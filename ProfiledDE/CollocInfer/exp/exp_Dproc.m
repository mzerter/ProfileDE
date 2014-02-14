function fnval = exp_Dproc(coefs,bvals,pars,more)

[m,n] = size(bvals);
devals  = exp(bvals(1:m-1,:) * coefs);
ddevals = exp(bvals(2:m,:)   * coefs);
% colnames(devals) = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
f = more.fn(ddevals,more.qpts,devals,pars,more.more);

fnval = sum(f);

end

