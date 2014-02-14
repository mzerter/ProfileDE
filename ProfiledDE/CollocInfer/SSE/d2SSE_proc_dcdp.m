%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dcdpval = d2SSE_proc_dcdp(coefs,bvals,pars,more)

if ~isfield(more,'more')
    more = [];
end
devals  = bvals.bvals  * coefs;
ddevals = bvals.dbvals * coefs;
tmp = make_SSElik;
H1 = tmp.d2fdxdp(ddevals,more.qpts,devals,pars,more);
H2 = 2*more.dfdp(more.qpts,devals,pars,more.more);
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
H = zeros(size(coefs,2)*size(bvals.bvals,2),length(pars));
for i = 1:length(pars)
    term1 = bvals.bvals' *squeeze(H1(:,:,i));
    term2 = bvals.dbvals'*(weights.*squeeze(H2(:,:,i)));
    term = term1 - term2;
    H(:,i) = term(:);
end
dcdpval = H;

end