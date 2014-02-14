function fnval = SIR_fun(t,x,p,more)
%  Right side function for SIR differential equation
%  x%  To be called by CollocInfer functions
% x(1) = S,    x(2) = I,     x(3) = R
%  p(1) = beta, p(2) = gamma, p(3) = mu
p = exp(p);
n = size(x,2);
if n == 1, x = x'; end
r = x;
r(:,1) = -p(1).*x(:,1).*x(:,2) -       p(3) .*x(:,1) + p(3);
r(:,2) =  p(1).*x(:,1).*x(:,2) - (p(2)+p(3)).*x(:,2);
r(:,3) =  p(2).*x(:,2)         -       p(3) .*x(:,3);
fnval = r;
if n == 1, fnval = fnval'; end

end

