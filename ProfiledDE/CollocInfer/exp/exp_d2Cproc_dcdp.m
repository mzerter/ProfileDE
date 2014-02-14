function dcdpval = exp_d2Cproc_dcdp(coefs,bvals,pars,more)

bmat = bvals.bvals;
dmat = bvals.dbvals;
devals = exp(bmat * coefs);
ddevals = (dmat * coefs)*devals;
% colnames(devals)  = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
H1 = more.d2fdxdp(ddevals,more.qpts,devals,pars,more.more);
H2 = more.d2fdydp(ddevals,more.qpts,devals,pars,more.more);
H = [];
for i = 1:length(pars)
    H = [H,as.vector(bmat' * (devals*H1(:,:,i) + ddevals*H2(:,:,i)) + ...
                     dmat' * (devals*H2(:,:,i)))];
end
dcdpval = H;

end

