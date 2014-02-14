function fnval = Lorenz_fun(t,x,p,more)
%  Right side function for Lorenz differential equation
%  x(1) = x,    x(2) = y,     x(3) = z
%  p(1) = sigma, p(2) = r, p(3) = b
p = exp(p);
n = size(x,2);
if n == 1, x = x'; end
r = x;
r(:,1) = -p(1).*x(:,1) + p(1).*x(:,2);
r(:,2) =  p(2).*x(:,1) - x(:,2) + x(:,1).*x(:,3);
r(:,3) = -p(3).*x(:,3) + x(:,1).*x(:,2);
fnval = r;
if n == 1, fnval = fnval'; end

end

