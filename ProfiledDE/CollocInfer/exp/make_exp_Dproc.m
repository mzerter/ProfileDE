
%%%%%%%%%%% Likelihood for Discrete-Time Dynamics %%%%%%%%%%%%%%%%

function expDproc = make_exp_Dproc()

expDproc.fn      = @exp_Dproc;
expDproc.dfdc    = @exp_dDproc_dc;
expDproc.dfdp    = @exp_dDproc_dp;
expDproc.d2fdc2  = @exp_d2Dproc_dc2;
expDproc.d2fdcdp = @exp_d2Dproc_dcdp;

end

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


function dpval = exp_dDproc_dp(coefs,bvals,pars,more)

[m,n] = size(bvals);
devals  = exp(bvals(1:m-1,:) * coefs);
ddevals = exp(bvals(2:m,:)   * coefs);
% colnames(devals) = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
g = more.dfdp(ddevals,more.qpts,devals,pars,more.more);
dpval = applysum(g);

end

function dc2val = exp_d2Dproc_dc2(coefs,bvals,pars,more)

ind = 1:(size(bvals,1)-1);
devals  = exp(bvals(ind,:) * coefs);
ddevals = exp(bvals(ind+1,:) * coefs);
% colnames(devals) = more.names;
% colnames(ddevals) = more.names;
% names(pars) = more.parnames;
g1 = more.dfdx(ddevals,more.qpts,devals,pars,more.more)*devals;
g2 = more.dfdy(ddevals,more.qpts,devals,pars,more.more)*ddevals;
H1 = more.d2fdx2(ddevals,more.qpts,devals,pars,more.more);
H2 = more.d2fdxdy(ddevals,more.qpts,devals,pars,more.more);
H3 = more.d2fdy2(ddevals,more.qpts,devals,pars,more.more);
H = cell(size(bvals,2),1);
ndev = size(devals,2);
for i = 1:ndev
    devi  = devals(:,i);
    ddevi = ddevals(:,i);
    H{i,j} = cell(size(devals))
    for j = 1:ndev
        devj  = devals(:,j);
        ddevj = ddevals(:,j);
        H{i,j} = ...
       bvals(ind,  :)' * diag(devi *H1(:,i,j)*devj)  * bvals(ind,  :) + ...
       bvals(ind,  :)' * diag(devi *H2(:,i,j)*ddevj) * bvals(ind+1,:) + ...
       bvals(ind+1,:)' * diag(ddevi*H2(:,j,i)*devj)  * bvals(ind,  :) + ...
       bvals(ind+1,:)' * diag(ddevi*H3(:,i,j)*ddevj) * bvals(ind+1,:);
    end
    H{i,j} = H{i,j} + ...
        bvals(ind,  :)' * diag(g1(:,i)) * bvals(ind,  :) + ...
        bvals(ind+1,:)' * diag(g2(:,i)) * bvals(ind+1,:);
end
dc2val = blocks2mat(H);

end

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

