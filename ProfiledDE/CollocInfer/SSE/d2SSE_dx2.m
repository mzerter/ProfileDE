%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dx2val = d2SSE_dx2(data,times,devals,pars,more)

if ~isfield(more,'more')
    more.more = [];
end

fdevals = more.fn(times,devals,pars,more.more);

difs    = data - fdevals;
dfdx    = more.dfdx(times,devals,pars,more.more);
d2fdx2  = more.d2fdx2(times,devals,pars,more.more);
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
[m,n] = size(devals);
H = zeros(m,n,n);
for i = 1:n
    for j = 1:n
        H(:,i,j) = sum(-difs.*squeeze(d2fdx2(:,:,i,j)) +  ...
            weights.*squeeze(dfdx(:,:,j).*dfdx(:,:,i)),2);
    end
end
dx2val = 2*H;

end

