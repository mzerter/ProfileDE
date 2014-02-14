%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdyval = d2SSE_dxdy(data,times,devals,pars,more)

if ~isfield(more,'more')
    more.more = [];
end

dfdx = more.dfdx(times,devals,pars,more.more);

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

weights = checkweights(more.weights,whichobs, ...
                       squeeze(dfdx(:,:,1)));
weights(isnan(data)) = 0;

for i = 1:size(dfdx,3)
    dfdx(:,:,i) = weights.*dfdx(:,:,i);
end
dxdyval = -2*permute(dfdx,[1,3,2]);

end

