function dydpval = d2SSE_dydp(data,times,devals,pars,more)

if ~isfield(more,'more')
    more.more = [];
end

dfdp = more.dfdp(times,devals,pars,more.more);

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

weights = checkweights(more.weights,more.whichobs, ...
                       squeeze(dfdp(:,:,1)));
weights(isnan(data)) = 0;

for i = 1:size(dfdp,3)
    dfdp(:,:,i) = weights.*squeeze(dfdp(:,:,i));
end
dydpval = -2*dfdp;

end
