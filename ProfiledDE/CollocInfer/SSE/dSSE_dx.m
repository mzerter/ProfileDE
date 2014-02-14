function dxval = dSSE_dx(data,times,devals,pars,more)

if ~isfield(more,'more')
    more.more = [];
end

fdevals = more.fn(times,devals,pars,more.more);
dfdx    = more.dfdx(times,devals,pars,more.more);
[l,m,n] = size(dfdx);
if any(size(data) ~= size(fdevals))
    disp(['size of data:     ',num2str(size(data))])
    disp(['size of fdevals:  ',num2str(size(fdevals))])
end
difs = data - fdevals;
difs(isnan(difs)) = 0;

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

weights = checkweights(more.weights,whichobs,difs);

difs = weights.*difs;
g    = zeros(l,m);
for i = 1:n
    g(:,i) = sum((difs.*squeeze(dfdx(:,:,i))),2);
end
dxval = -2*g;

end

