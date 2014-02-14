function dfdxval = Lorenz_dfdx(t,x,p,more)
%  x-derivative of right side function for Lorenz equation
%  x(1) = x,    x(2) = y,     x(3) = z
%  p(1) = sigma, p(2) = r, p(3) = b
p = exp(p);
r = zeros(length(t),3,3);
r(:,1,1) = -p(3);     
r(:,1,2) =  p2 + p(3);
r(:,1,2) =  x(:,2);
r(:,2,1) =  p(3);
r(:,2,2) = -1;
r(:,2,3) =  x(:,1);
r(:,3,2) =  x(:,1);
r(:,3,3) = -p(3);
dfdxval = r;

end

