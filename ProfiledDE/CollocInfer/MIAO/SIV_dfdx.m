function dfdxval = SIV_dfdx(t,x,p,more)
%  x-derivative for SIV differential equation
%  To be called by CollocInfer functions
%  x(1) = S,    x(2) = I,     x(3) = V
%  p(1) = rho, p(2) = beta, p(3) = delta, p(4) = pi, p(5) = c
r = zeros(length(t),3,3);
r(:,1,1) =  p(1) - p(2).*x(:,3);
r(:,1,3) = -p(2).*x(:,1);
r(:,2,1) =  p(2).*x(:,3);
r(:,2,2) = -p(3);
r(:,2,3) =  p(2).*x(:,1);
r(:,3,2) =  p(4);
r(:,3,3) = -p(5);
dfdxval = r;

end

