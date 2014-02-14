%%%% A general multivariate normal function; some of this 
% should be coded as C

function multinorm = make_multinorm()

multinorm.fn      = @multinorm_fn;
multinorm.dfdx    = @multinorm_dfdx;
multinorm.dfdy    = @multinorm_dfdy;
multinorm.dfdp    = @multinorm_dfdp;
multinorm.d2fdx2  = @multinorm_d2fdx2;
multinorm.d2fdy2  = @multinorm_d2fdy2;
multinorm.d2fdxdy = @multinorm_d2fdxdy;
multinorm.d2fdxdp = @multinorm_d2fdxdp;
multinorm.d2fdydp = @multinorm_d2fdydp;
    
end

function fnval = multinorm_fn(y,times,x,pars,more)

%  computes negative log density for multivariate normal distribution
F = more.fn(times,x,pars,more.f.more);
S = more.var.fn(times,x,pars,more.v.more);
n = size(F,1);
d = zeros(n,1);
if  size(S,1) == 1
    %  if variance is constant, inverse of S is computed once
    Sinv = inv(squeeze(S(1,:,:)));
    lognormconst = -0.5*log(det(Sinv));
end
for i = 1:n
    if  size(S,1) > 1
        Sinv = inv(squeeze(S(i,:,:)));
        lognormconst = -0.5*log(det(Sinv));
    end
    resi = y(i,:)-F(i,:);
    d(i) = 0.5.*resi'*Sinv*resi + lognormconst;
end
fnval = d;

end

function dfdyval = multinorm_dfdy(y,times,x,pars,more)

F = more.fn(times,x,pars,more.f.more);
S = more.var.fn(times,x,pars,more.v.more);
d = zeros(size(y));
m = size(S,1);
n = size(y,1);
if m == 1
    Sinv = inv(squeeze(S(1,:,:)));
end
for i = 1:n
    if  m > 1
        Sinv = inv(squeeze(S(i,:,:)));
    end
    d(i,:) = Sinv*(y(i,:)-F(i,:));
end
dfdyval = d;

end

function dfdxval = multinorm_dfdx(y,times,x,pars,more)

F  = more.fn(times,x,pars,more.f.more);
dF = more.dfdx(times,x,pars,more.f.more);
S  = more.var.fn(times,x,pars,more.v.more);
dS = more.var.dfdx(times,x,pars,more.v.more);
d = zeros(size(x));
l = size(S,1);
m = size(dS,2:4);
n = size(dF,2:3);
if  l == 1
    SS = inv(squeeze(S(1,:,:)));
    dSS = array(dS(1,:,:,:),m);
end
for i = 1:size(x,1)
    if  l > 1
        SS = inv(squeeze(S(i,:,:)));
        dSS = dS(i,m);
    end
    wdifs = SS*(y(i,:)-F(i,:));
    d(i,:) = -squeeze(dF(i,m))'*wdifs';
    for j = 1:size(x,2)
        d(i,j) = d(i,j) - 0.5*wdifs'*squeeze(dSS(:,:,j))*wdifs + ...
                          0.5*sum(diag(SS*squeeze(dSS(:,:,j))));
    end
end
dfdxval = d;

end

function dfdpval = multinorm_dfdp(y,times,x,pars,more)

F  = more.fn(times,x,pars,more.f.more);
dF = more.dfdp(times,x,pars,more.f.more);
S  = more.var.fn(times,x,pars,more.v.more);
dS = more.var.dfdp(times,x,pars,more.v.more);
k = size(x,1);
l = size(S,1);
m = size(dS,2:4);
n = size(dF,2:3);
d = zeros(k,length(pars));
if  l == 1
    SS = inv(squeeze(S(1,:,:)));
    dSS = dS(1,m);
end
for i = 1:k
    if  l > 1
        SS = inv(squeeze(S(i,:,:)));
        dSS = dS(i,m);
    end
    wdifs = SS*(y(i,:)-F(i,:));
    d(i,:) = -squeeze(dF(i,n))'*wdifs;
    for j = 1:length(pars)
        d(i,j) = d(i,j) - 0.5*wdifs'*squeeze(dSS(:,:,j)*wdifs + ...
                          0.5*sum(diag(SS*squeeze(dSS(:,:,j)))));
    end
end
dfdpval = d;

end

function d2fdx2val = multinorm_d2fdx2(y,times,x,pars,more)

F   = more.fn(times,x,pars,more.f.more);
dF  = more.dfdx(times,x,pars,more.f.more);
d2F = more.d2fdx2(times,x,pars,more.f.more);
S   = more.var.fn(times,x,pars,more.v.more);
dS  = more.var.dfdx(times,x,pars,more.v.more);
d2S = more.var.d2fdx2(times,x,pars,more.v.more);
[k1,k2] = size(x);
l = size(S,1);
m = size(dS,2:4);
n = size(dF,2:3);
d = zeros(k1,k2,k2);
if  l == 1
    SS = inv(squeeze(S(1,:,:)));
    dSS = array(dS(1,m));
    d2SS = array(d2S(1,size(d2S,2:5)));
end
for i = 1:k1
    if  l > 1
        SS = inv(squeeze(S(i,:,:)));
        dSS = array(dS(i,m));
        d2SS = array(d2S(i,size(d2S,2:5)));
    end
    wdifs = SS * (y(i,:)-F(i,:));
    for j = 1:k2
        dFj  = dF(i,:,j);
        dSSj = squeeze(dSS(:,:,j));
        for k = j:k2
            dFk = dF(i,:,k);
            dSSk = squeeze(dSS(:,:,k));
            d2SSjk = squeeze(d2SS(:,:,j,k));
            %  this statement is probably wrong ... check
            d(i,j,k) =  ...
                -(d2F(i,:,j,k)-dFj*SS*dSSk)-dFk*SS*dSSj*wdifs + ...
                0.5*wdifs'*(dSSj*SS*dSSk+dSSk*SS*dSSj-d2SSjk)*wdifs + ...
                dFj'*SS*dFk - 0.5*sum(diag(SS*(dSSj*SS*dSSk-d2SSjk)));
            d(i,k,j) = d(i,j,k);
        end
    end
end
d2fdx2val = d;

end

function d2fdy2val = multinorm_d2fdy2(y,times,x,pars,more)

S = more.var.fn(times,x,pars,more.v.more);
[k1,k2] = size(x);
l = size(S,1);
d = zeros(size(y),size(y,2));
if  l == 1 
    SS = inv(squeeze(S(1,:,:)));
end
for i = 1:k1
    if  l > 1 
        SS = inv(squeeze(S(i,:,:)));
    end
    d(i,:,:) = 0.5*(SS' + SS);
end
d2fdy2val = d;

end

function d2fdxdyval = multinorm_d2fdxdy(y,times,x,pars,more)

F  = more.fn(times,x,pars,more.f.more);
dF = more.dfdx(times,x,pars,more.f.more);
S  = more.var.fn(times,x,pars,more.v.more);
dS = more.var.dfdx(times,x,pars,more.v.more);
[k1,k2] = size(x);
l = size(S,1);
m = size(dS,2:4);
d = zeros(k1,k2,size(y,2));
if  l == 1
    SS = inv(squeeze(S(1,:,:)));
    dSS = array(dS(1,m));
end
for i = 1:k1
    if  l > 1
        SS = inv(squeeze(S(i,:,:)));
        dSS = array(dS(i,size(dS,2:5)));
    end
    for j = 1:k2
        d(i,j,:) = -SS*squeeze(dSS(:,:,j))*SS*(y(i,:)-F(i,:)) - ...
            SS*dF(i,:,j);
    end
end
d2fdxdyval = d;

end

function d2fdxdpval = multinorm_d2fdxdp(y,times,x,pars,more)

F   = more.fn(times,x,pars,more.f.more);
dFx = more.dfdx(times,x,pars,more.f.more);
dFp = more.dfdp(times,x,pars,more.f.more);
d2F = more.d2fdxdp(times,x,pars,more.f.more);
S   = more.var.fn(times,x,pars,more.v.more);
dSx = more.var.dfdx(times,x,pars,more.v.more);
dSp = more.var.dfdp(times,x,pars,more.v.more);
d2S = more.var.d2fdxdp(times,x,pars,more.v.more);
[k1,k2] = size(x);
l = size(S,1);
m = size(dS,2:4);
n = size(dF,2:3);
d = zeros(k1,k2,length(pars));
if  l == 1
    SS = inv(squeeze(S(1,:,:)));
    dSSx = array(dSx(1,size(dSx,2:4)));
    dSSp = array(dSp(1,size(dSp,2:4)));
    d2SS = array( d2S(1,size(d2S,2:5)));
end
for i = 1:k1
    if  l > 1
        SS = inv(squeeze(S(i,:,:)));
        dSSx = array(dSx(i,size(dSx,2:4)));
        dSSp = array(dSp(i,size(dSp,2:4)));
        d2SS = array(d2S(i,size(d2S,2:5)));
    end
    wdifs = SS*(y(i,:)-F(i,:));
    for j = 1:k2
        dFxj = dFx(i,:,j);
        dSSxj = squeeze(dSSx(:,:,j));
        for k = 1:length(pars)
            dFpk = dFp(i,:,k);         
            dSSpk = squeeze(dSSp(:,:,k));
            d2Fjk = d2F(i,:,j,k);
            d2SSjk = squeeze(d2SS(:,:,j,k));
            d(i,j,k) = ...
                -(d2Fjk-dFx(i,:,j)*SS*dSSpk-dFpk*SS*dSSxj*wdifs + ...
            0.5*wdifs'*(dSSxj*SS*dSSpk+dSSpk*SS*dSSxj-d2SSjk*wdifs + ...
            dFxj'*SS*dFpK - 0.5*sum(diag(SS*(dSSxJ*SS*dSSpk - d2SSjk)))));
        end
    end
end
d2fdxdpval = d;

end

function d2fdydpval = multinorm_d2fdydp(y,times,x,pars,more)

F  = more.fn(times,x,pars,more.f.more);
dF = more.dfdp(times,x,pars,more.f.more);
S  = more.var.fn(times,x,pars,more.v.more);
dS = more.var.dfdp(times,x,pars,more.v.more);
[j1,j2] = size(y);
m = size(dS,2:4);
d = zeros(j1,j2,length(pars));
if  l == 1
    SS = inv(squeeze(S(1,:,:)));
    dSS = array(dS(1,m));
end
for i = 1:j1
    if  l > 1
        SS = inv(squeeze(S(i,:,:)));
        dSS = array(dS(i,m));
    end
    for j = 1:length(pars)
        dSSj = squeeze(dSS(:,:,j));
        d(i,:,j) = -SS*(dSSj*SS*(y(i,:)-F(i,:)) + dF(i,:,j));
    end
end
d2fdydpval = d;

end
