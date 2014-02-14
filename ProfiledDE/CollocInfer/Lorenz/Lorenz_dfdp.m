function dfdpval = Lorenz_dfdp(t,x,p,more)
%  p-derivative of right side function for Lorenz equation
%  x(1) = x,    x(2) = y,     x(3) = z
%  p(1) = sigma, p(2) = r, p(3) = b
p = exp(p);
r = zeros(length(t),3,3);
r(:,1,3) =  x(:,2) - x(:,1);
r(:,2,2) =  x(:,1);
r(:,3,3) = -x(:,3);
dfdpval  = r;

end

