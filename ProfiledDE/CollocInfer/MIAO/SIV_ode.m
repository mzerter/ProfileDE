function fnval = SIV_ode(t,x,p,more)
%  Right side function for SIV differential equation
%  To be called by CollocInfer functions
%  x(1) = S,    x(2) = I,     x(3) = V
%  p(1) = rho, p(2) = beta, p(3) = delta, p(4) = pi, p(5) = c
r = x;
r(:,1) = p(1).*x(:,1)         - p(2).*x(:,1).*x(:,3);
r(:,2) = p(2).*x(:,1).*x(:,3) - p(3).*x(:,2);
r(:,3) = p(4).*x(:,2)         - p(5).*x(:,3);
fnval = r;

end

