%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy2val = d2SSE_dy2(data,times,devals,pars,more)

[m,n] = size(data);

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

weights = checkweights(more.weights,more.whichobs,difs);

ind = [repmat((1:m)',1,n), kron([1:n,1:n],ones(m,1))];
r       = zeros(m,n,n);
r(ind)  = weights;
dy2val  = 2*r;

end

