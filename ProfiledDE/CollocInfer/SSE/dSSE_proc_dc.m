function dcval = dSSE_proc_dc(coefs,bvals,pars,more)

if ~isfield(more,'more')
    more.more = [];
end

m       = size(bvals.bvals,2);
n       = size(coefs,2);
devals  = bvals.bvals *coefs;
ddevals = bvals.dbvals*coefs;

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

tmp     = make_SSElik;
g1      = tmp.dfdx(ddevals,more.qpts,devals,pars,more);
weights = checkweights(more.weights,whichobs,g1);
term2   = more.fn(more.qpts,devals,pars,more.more);
g2      = weights.*(ddevals - term2);
dcval   = reshape(bvals.bvals'*g1 + 2*bvals.dbvals'*g2,m,n,1);

end

