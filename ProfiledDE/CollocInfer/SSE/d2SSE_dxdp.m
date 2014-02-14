%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdpval = d2SSE_dxdp(data,times,devals,pars,more)

if ~isfield(more,'more')
    more.more = [];
end

fdevals = more.fn(times,devals,pars,more.more);
difs    = data - fdevals;
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
dfdx    = more.dfdx(times,devals,pars,more.more);
dfdp    = more.dfdp(times,devals,pars,more.more);
d2fdxdp = more.d2fdxdp(times,devals,pars,more.more);

[nt,nout,nx,np] = size(d2fdxdp);
H = zeros(nt,nx,np);
for i = 1:nx
    dfdxi = squeeze(dfdx(:,:,i));
    for j = 1:np
        dfdpj = squeeze(dfdp(:,:,j));
        d2fdxdpij = squeeze(d2fdxdp(:,:,i,j));
        term1 = -difs.*d2fdxdpij;
        term2 = weights.*dfdxi.*dfdpj;
        H(:,i,j) = sum(term1+term2,2);
    end
end
dxdpval = 2*H;

end

