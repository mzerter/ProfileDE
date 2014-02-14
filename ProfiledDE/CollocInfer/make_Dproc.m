%%%%%%%%%%% Likelihood for Discrete-Time Dynamics %%%%%%%%%%%%%%%% 

function Dproc = make_Dproc()

Dproc.fn      = @Dproc_fn;
Dproc.dfdc    = @dDproc_dfdc;
Dproc.dfdp    = @dDproc_dfdp;
Dproc.d2fdc2  = @d2Dproc_d2fdc2;
Dproc.d2fdcdp = @d2Dproc_d2fdcdp;

end 

function funval = Dproc_fn(coefs,bvals,pars,more)

[m,n] = size(bvals);
devals  = bvals(1:m-1,:) * coefs;
ddevals = bvals(2:m,:)   * coefs;
%    colnames(devals) = more.names
%    colnames(ddevals) = more.names
%    names(pars) = more.parnames
f = more.fn(ddevals,more.qpts,devals,pars,more.more);
funval = sum(f);

end

function dfdcval = dDproc_dfdc(coefs,bvals,pars,more)

[m,n] = size(bvals);
devals  = bvals(1:m-1,:) * coefs;
ddevals = bvals(2:m,:)   * coefs;
% colnames(devals) = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
g1 = more.dfdx(ddevals,more.qpts,devals,pars,more.more);
g2 = more.dfdy(ddevals,more.qpts,devals,pars,more.more);
g = reshape(bvals(1:m-1,:)'*g1 + bvals(2:m,:)'*g2,n*size(g2,2),1);
dfdcval = g;

end

function dfdpval = dDproc_dfdp(coefs,bvals,pars,more)

[m,n] = size(bvals);
devals  = bvals(1:m-1,:) * coefs;
ddevals = bvals(2:m,:  ) * coefs;
% colnames(devals) = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
g = more.dfdp(ddevals,more.qpts,devals,pars,more.more);
g = applysum(g);
dfdpval = g;

end

function d2fdc2val = d2Dproc_d2fdc2(coefs,bvals,pars,more)

[m,n] = size(bvals);
devals  = bvals(1:m-1,:) * coefs;
ddevals = bvals(2:m,:  ) * coefs;
% colnames(devals) = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
H1 = more.d2fdx2( ddevals,more.qpts,devals,pars,more.more);
H2 = more.d2fdxdy(ddevals,more.qpts,devals,pars,more.more);
H3 = more.d2fdy2( ddevals,more.qpts,devals,pars,more.more);
H = cell(n);
k = size(devals,2);
for i = 1:k
    for j = 1:k
        H{i,j }= bvals(1:m-1,:)' * diag(H1(:,i,j)) * bvals(1:m-1,:) + ...
                 bvals(1:m-1,:)' * diag(H2(:,i,j)) * bvals(2:m,:  ) + ...
                 bvals(2:m,:  )' * diag(H2(:,j,i)) * bvals(1:m-1,:) + ...
                 bvals(2:m,:  )' * diag(H3(:,i,j)) * bvals(1:m-1,:);
    end
end
H = blocks2mat(H);
d2fdc2val = H;

end

function d2fdcdpval = d2Dproc_d2fdcdp(coefs,bvals,pars,more)

nbvals  = size(bvals,1);
ind     = 1:(nbvals-1);
devals  = bvals(ind,  :) * coefs;
ddevals = bvals(ind+1,:) * coefs;
% colnames(devals) = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
H1 = more.d2fdxdp(ddevals,more.qpts,devals,pars,more.more);H2 = more.d2fdydp(ddevals,more.qpts,devals,pars,more.more);
H = [];
for i = 1:length(pars)
    H = [H,as.vector(bvals(ind,  :)' * H1(:,:,i) + ...
                     bvals(ind+1,:)' * H2(:,:,i))];
end
d2fdcdpval = H;

end
