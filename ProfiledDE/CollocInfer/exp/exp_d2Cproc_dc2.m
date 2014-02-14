function dc2val = exp_d2Cproc_dc2(coefs,bvals,pars,more)

bmat = bvals.bvals;
dmat = bvals.dbvals;
devals  = exp(bmat * coefs);
ddevals = (dmat * coefs)*devals;
% colnames(devals) = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
g1 = more.dfdx(ddevals,more.qpts,devals,pars,more.more);
g2 = more.dfdy(ddevals,more.qpts,devals,pars,more.more);
H1 = more.d2fdx2(ddevals,more.qpts,devals,pars,more.more);
H2 = more.d2fdxdy(ddevals,more.qpts,devals,pars,more.more);
H3 = more.d2fdy2(ddevals,more.qpts,devals,pars,more.more);
H = cell(size(bmat,2));
for i = 1:size(devals,2)
    devi    = diag(devals(:,i));
    ibvals  = diag(devi) * bmat;
    idbvals = diag(ddevals(:,i)) * bmat + devi * dmat;
    for j = 1:size(devals,2)
        H{i,j} = cell(size(devals));
        devj = diag(devals(:,j));
        jbvals = diag(devj) * bmat;
        jdbvals = diag(ddevals(:,j)) * bmat + devj * dmat;
        H{i,j} = ibvals'  * diag(H1(:,i,j)) * jbvals  + ...
                 ibvals'  * diag(H2(:,i,j)) * jdbvals + ...
                 idbvals' * diag(H2(:,j,i)) * jbvals  + ...
                 idbvals' * diag(H3(:,i,j)) * jdbvals;
    end
    H{i,i} = H{i,i} + ibvals'  * diag(g1(:,i)) * bmat + ...
                        idbvals' * diag(g2(:,i)) * bmat + ...
                        ibvals'  * diag(g2(:,i)) * dmat;   
end

dc2val = blocks2mat(H);

end

