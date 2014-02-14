%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dpval = dSSE_dp(data,times,devals,pars,more)

if ~isfield(more,'more')
    more.more = [];
end

fdevals = more.fn(times,devals,pars,more.more);
difs    = data-fdevals;
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

difs    = weights.*difs;
dfdp = more.dfdp(times,devals,pars,more.more);
[l,m,n] = size(dfdp);
g = zeros(l,m);
for i = 1:n
    temp = sum(difs.*squeeze(dfdp(:,:,i)),2);
    g(:,i)  = temp;
end
dpval = -2*g;

end

