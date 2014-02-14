function fnval = SIVlog_fun(t,x,p,more)
%  Right side function for SIV differential equation
%  x%  To be called by CollocInfer functions
%  x(1) = S,    x(2) = I,     x(3) = V
%  p(1) = rho, p(2) = beta, p(3) = delta, p(4) = pi, p(5) = c
p = exp(p);
x = 10.^x;
logten = log(10);
r = x;
r(:,1) = (p(1)                 - p(2).*x(:,3))./logten;
r(:,2) = (p(2).*x(:,1).*x(:,3) - p(3).*x(:,2))./(logten*x(:,2));
r(:,3) = (p(4).*x(:,2)         - p(5).*x(:,3))./(logten*x(:,3));
fnval = r;

end

