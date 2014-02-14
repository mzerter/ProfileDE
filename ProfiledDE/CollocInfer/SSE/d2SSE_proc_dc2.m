function dc2val = d2SSE_proc_dc2(coefs,bvals,pars,more)

devals  = bvals.bvals * coefs;
ddevals = bvals.dbvals * coefs;
tmp = make_SSElik;
H1  = tmp.d2fdx2(ddevals,more.qpts,devals,pars,more);
H2  = more.dfdx(more.qpts,devals,pars,more.more);
if isfield(more,'whichobs')
    whichobs = more.whichobs;
else
    whichobs = [];
end
if isfield(more,'which')
    whichobs = more.which;
else
    whichobs = [];
end
weights = checkweights(more.weights,whichobs,H1(:,:,1));
m = size(bvals.bvals,2);
n = size(devals,2);
H = cell(n,n);
for i = 1:n
    for j = 1:n
        H{i,j} = bvals.bvals' *diag(H1(:,i,j))              *bvals.bvals - ...
               2*bvals.dbvals'*diag(H2(:,i,j).*weights(:,i))*bvals.bvals - ...
               2*bvals.bvals' *diag(H2(:,j,i).*weights(:,j))*bvals.dbvals;
    end
    H{i,i} = H{i,i} + 2*bvals.dbvals'*diag(weights(:,i))*bvals.dbvals;
end
Hmat = blocks2mat(H);
dc2val = Hmat;

end

