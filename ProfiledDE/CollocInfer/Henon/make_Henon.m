function Henon = make_Henon()

Henon.ode     = @Henon_ode;
Henon.fun     = @Henon_fun;
Henon.dfdx    = @Henon_dfdx;
Henon.dfdp    = @Henon_dfdp;
Henon.d2fdx2  = @Henon_d2fdx2;
Henon.d2fdxdp = @Henon_d2fdxdp;

end

function fnval = Henon_ode(times,y,p,more)
r = y;
r(1) = 1- p(1)*y(1)^2 + y(2);
r(2) = p(2)*y(1);
fnval = r;

end

function fnval = Henon_fun(times,y,p,more)

r = y;
r(:,1) = 1- p(1)*y(:,1)^2 + y(:,2);
r(:,2) = p(2)*y(:,1);
fnval = r;

end

function dfdxval = Henon_dfdx(times,y,p,more)

r = zeros(length(times),size(y,2),size(y,2));
r(:,1,1) = -2*p(1)*y(:,1);
r(:,1,2) = 1;
r(:,2,1) = p(2);
dfdxval = r;

end

function dfdpval = Henon_dfdp(times,y,p,more)

r = zeros(length(times),size(y,2),length(p));
r(:,1,1) = -y(:,1).^2;
r(:,2,2) = y(:,1);
dfdpval = r;

end

function d2fdx2val = Henon_d2fdx2(times,y,p,more)

r = zeros(length(times),rep(size(y,2),3));
r(:,1,1,1) = -2*p(1);
d2fdx2val = r;

end

function d2fdxdpval = Henon_d2fdxdp(times,y,p,more)

r = zeros(length(times),size(y,2),size(y,2),length(p));
r(:,1,1,1) = -2*y(:,1);
r(:,2,1,2) = 1;
d2fdxdpval = r;

end
