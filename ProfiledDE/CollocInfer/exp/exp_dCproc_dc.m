function dcval = exp_dCproc_dc(coefs,bvals,pars,more)

[m,n] = size(bvals);
devals  = bvals(1:m-1,:) * coefs;
ddevals = bvals(2:m,:  ) * coefs;
devals  = exp(devals);
ddevals = (devals * coefs)*ddevals;
% colnames(devals) = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
g1 = more.dfdx(ddevals,more.qpts,devals,pars,more.more);g2 = more.dfdy(ddevals,more.qpts,devals,pars,more.more);
dcval = as.vector( bvals.bvals'  * (g1*devals) + ...
                   bvals.dbvals' * (g2*devals) + ...
                   bvals.bvals'  * (g2*ddevals) );

end

