function diagnostics = make_diagnostics()

    diagnostics.fn      = @diagnostics_fn;
    diagnostics.dfdx    = @diagnostics_dfdx;
    diagnostics.dfdp    = @diagnostics_dfdp;
    diagnostics.d2fdx2  = @diagnostics_d2fdx2;
    diagnostics.d2fdxdp = @diagnostics_d2fdxdp;
    diagnostics.more    = NULL;

end

function fnval = diagnostics_fn(t,x,p,more)

force = zeros(size(x));
force(:,more.which) = ...
    more.psi*reshape(p,size(more.psi,2),length(more.which));
fnval = more.fn(t,x,more.p,more.more) + force;

end


function dfdxval = diagnostics_dfdx(t,x,p,more)

dfdxval = more.dfdx(t,x,more.p,more.more);

end

function force = diagnostics_dfdp(t,x,p,more)

force = zeros(size(x),length(p));
k = size(more.psi,2);
whichp = 0;
for i = 1:length(more.which)
    force(:,more.which(i),whichp+(1:k)) = more.psi;
    whichp = whichp + k;
end

end

function d2fdx2val = diagnostics_d2fdx2(t,x,p,more)

d2fdx2val = more.d2fdx2(t,x,more.p,more.more);

end

function d2fdxdpval = diagnostics_d2fdxdp(t,x,p,more)

[m,n] = size(x);
d2fdxdpval = zeros(m,n,n,length(p));

end

