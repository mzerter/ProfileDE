function fnval = SIR_ode(t,x,p)
%  Right side function for SIR differential equation
%  To be called by CollocInfer functions
%  x(1) = S,    x(2) = I,     x(3) = R
%  p(1) = beta, p(2) = gamma, p(3) = mu
r = x;
r(:,1) = -p(1).*x(:,1).*x(:,2) - p(3)*x(:,1) + p(3);
r(:,2) =  p(1).*x(:,1).*x(:,2) - (p(2)+p(3)).*x(:,2);
r(:,3) =  p(2).*x(:,2)         - p(3).*x(:,3);
fnval = r;

end

