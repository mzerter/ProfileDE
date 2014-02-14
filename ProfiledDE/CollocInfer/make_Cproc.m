function Cproc = make_Cproc()

Cproc.fn      = @Cproc_fn;
Cproc.dfdc    = @dCproc_dfdc;
Cproc.dfdp    = @dCproc_dfdp;
Cproc.d2fdc2  = @d2Cproc_d2fdc2;
Cproc.d2fdcdp = @d2Cproc_d2fdcdp;

end

function fnval = Cproc_fn(coefs,bvals,pars,more)

devals  = bvals.bvals  * coefs;
ddevals = bvals.dbvals * coefs;
% colnames(devals)  = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
fnval = sum(more.fn(ddevals,more.qpts,devals,pars,more.more));

end

function dfdcval = dCproc_dfdc(coefs,bvals,pars,more)

devals  = bvals.bvals  * coefs;
ddevals = bvals.dbvals * coefs;
% colnames(devals)  = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
g1 = more.dfdx(ddevals,more.qpts,devals,pars,more.more);
g2 = more.dfdy(ddevals,more.qpts,devals,pars,more.more);
dfdcval = as.vector( bvals.bvals' * g1 + bvals.dbvals' * g2 );

end

function dfdpval = dCproc_dfdp(coefs,bvals,pars,more)

devals  = bvals.bvals  * coefs;
ddevals = bvals.dbvals * coefs;
% colnames(devals)  = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
g = more.dfdp(ddevals,more.qpts,devals,pars,more.more);
dfdpval = applysum(g);

end

function d2dc2val = d2Cproc_d2fdc2(coefs,bvals,pars,more)

devals  = bvals.bvals  * coefs;
ddevals = bvals.dbvals * coefs;
% colnames(devals)  = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
H1 = more.d2fdx2(ddevals,more.qpts,devals,pars,more.more);
H2 = more.d2fdxdy(ddevals,more.qpts,devals,pars,more.more);
H3 = more.d2fdy2(ddevals,more.qpts,devals,pars,more.more);
m = size(bvals.bvals,2);
n = size(devals,2);
H = cell(m,n);
for i = 1:m
    for j = 1:n
        H{i,j} = bvals.bvals'  * diag(H1(:,i,j)) * bvals.bvals  + ...
                 bvals.bvals'  * diag(H2(:,i,j)) * bvals.dbvals + ...
                 bvals.dbvals' * diag(H2(:,j,i)) * bvals.bvals  + ...
                 bvals.dbvals' * diag(H3(:,i,j)) * bvals.dbvals;
    end
end
d2dc2val = blocks2mat(H);

end

function d2fdcdpval = d2Cproc_d2fdcdp(coefs,bvals,pars,more)

devals  = bvals.bvals  * coefs;
ddevals = bvals.dbvals * coefs;
% colnames(devals)  = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
H1 = more.d2fdxdp(ddevals,more.qpts,devals,pars,more.more);
H2 = more.d2fdydp(ddevals,more.qpts,devals,pars,more.more);
H = [];
for i = 1:length(pars)
    H = [H, bvals.bvals'  * H1(:,:,i) + bvals.dbvals' * H2(:,:,i)];
end
d2fdcdpval = H;

end

