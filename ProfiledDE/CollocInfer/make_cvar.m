%%% CONSTANT VARIANCE FUNCTIONS

function cvar = make_cvar()

cvar.fn      = @cvar_fn;
cvar.dfdx    = @cvar_dfdx;
cvar.dfdp    = @cvar_dfdp;
cvar.d2fdxdp = @cvar_d2fdxdp;
cvar.d2fdx2  = @cvar_d2fdx2;

end 

function pmat = cvar_fn(t,y,pmat,more)

n     = size(y,2);
[l,m] = size(pmat);
more = checkmore_cvar(more,n);
pmat = 0.*more.mat;
pmat(more.sub(:,1:2)) = pmat(more.sub(:,1:2)) + p(more.sub(:,3));
if l > 1
    pmat = pmat - diag(diag(pmat));
    pmat(more.sub(:,c(2,1))) = pmat(more.sub(:,c(2,1))) + p(more.sub(:,3));
end
pmat = pmat + more.mat;
pmat = reshape(pmat,1,l,m);

end

function dfdxval = cvar_dfdx(t,y,pmat,more)
        
if ~isempty(more.mat) 
    ny = size(more.mat,1); 
else
    ny = size(y,2);
    nx = size(y,2);
end
    dfdxval = zeros(1,ny,ny,nx);

end

function r = cvar_dfdp(t,y,pmat,more)

n = size(y,2);
more = checkmore_cvar(more,n);
r = zeros(1,size(more.mat),length(p));
ind = [ones(2*size(more.sub,1),1), ...
       [more.sub; more.sub(:,2,1,3)]];
r(ind) = 1;

end

function d2fdxdpval = cvar_d2fdxdp(t,y,p,more)

if ~isempty(more.mat)
    ny = size(more.mat);
else
    ny = size(y,2);
    nx = size(y,2);
end
d2fdxdpval = zeros(1,ny,ny,nx,length(p));

end

function d2fdx2val = cvar_d2fdx2(t,y,p,more)

if ~isempty(more.mat)
    ny = size(more.mat,1);
else
    ny = size(y,2);
end
nx = size(y,2);
d2fdx2val = zeros(1,ny,ny,nx,nx);

end

