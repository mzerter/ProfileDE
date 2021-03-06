function dfdpval = SIV_dfdp(t,x,p,more)
%  p-derivative of for SIV differential equation
%  To be called by CollocInfer functions
%  x(1) = S,    x(2) = I,     x(3) = V
%  p(1) = rho, p(2) = beta, p(3) = delta, p(4) = pi, p(5) = c
r = zeros(length(t),3,5);
r(:,1,1) =  x(:,1);
r(:,1,2) = -x(:,1).*x(:,3);
r(:,2,2) =  x(:,1).*x(:,3);
r(:,2,3) = -x(:,2);
r(:,3,4) =  x(:,2);
r(:,3,5) = -x(:,3);
dfdpval  = r;

end

