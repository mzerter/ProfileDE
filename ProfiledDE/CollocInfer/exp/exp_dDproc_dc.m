function dcval = exp_dDproc_dc(coefs,bvals,pars,more)

[m,n] = size(bvals);
devals  = exp(bvals(1:m-1,:) * coefs);
ddevals = exp(bvals(2:m,:)   * coefs);
% colnames(devals) = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
g1 = more.dfdx(ddevals,more.qpts,devals,pars,more.more)*devals;
g2 = more.dfdy(ddevals,more.qpts,devals,pars,more.more)*ddevals;
ind = 1:m-1;
g = as.vector( bvals(ind,:)' * g1 + bvals(ind+1,:)' * g2 );
dcval = g;

end

