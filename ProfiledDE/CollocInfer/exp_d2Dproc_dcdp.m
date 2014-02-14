function dcdpval = exp_d2Dproc_dcdp(coefs,bvals,pars,more)

[m,n] = size(bvals);
devals  = exp(bvals(1:m-1,:) * coefs);
ddevals = exp(bvals(2:m,:)   * coefs);
% colnames(devals) = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames; 
H1 = more.d2fdxdp(ddevals,more.qpts,devals,pars,more.more);
H2 = more.d2fdydp(ddevals,more.qpts,devals,pars,more.more);
H = [];
ind = 1:m-1;
for i = 1:length(pars)
    H = [H,as.vector(bvals(ind,  :)' * ( devals*H1(:,:,i)) + ...
                     bvals(ind+1,:)' * (ddevals*H2(:,:,i)))];
end
dcdpval = H;

end

