function d2fdxdpval = SIVlog_d2fdxdp(t,x,p,more)
%  x-p-derivative for SIV differential equation
%  To be called by CollocInfer functions
%  x(1) = S,    x(2) = I,     x(3) = V
%  p(1) = beta, p(2) = gamma, p(3) = mu
%  p(1) = rho, p(2) = beta, p(3) = delta, p(4) = pi, p(5) = c
p = exp(p);
x = 10.^x;
logten = log(10);
r = zeros(length(t),3,3,5);
r(:,1,1,2) = (-x(:,3))*p(2);
r(:,1,3,2) = (-x(:,1))*p(2);
r(:,2,1,2) = ( x(:,3)./(logten*x(:,2)))*p(2);
r(:,2,2,2) = (-x(:,1).*x(:,3)./(logten*x(:,2).^2))*p(2);
r(:,2,3,2) = ( x(:,1)./(logten*x(:,2)))*p(2);
r(:,3,2,4) = ( 1./x(:,3))*p(4);
r(:,3,3,4) = (-x(:,2)./x(:,3).^2)*p(4);
d2fdxdpval = r;

end
