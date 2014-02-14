function dfdpval = SIVlog_dfdp(t,x,p,more)
%  p-derivative of for SIV differential equation
%  To be called by CollocInfer functions
%  x(1) = S,    x(2) = I,     x(3) = V
%  p(1) = rho, p(2) = beta, p(3) = delta, p(4) = pi, p(5) = c
p = exp(p);
x = 10.^x;
logten = log(10);
r = zeros(length(t),3,5);
r(:,1,1) = ( 1/logten)*p(1);
r(:,1,2) = (-x(:,3)./logten)*p(2);
r(:,2,2) = ( x(:,1).*x(:,3)./(logten*x(:,2)))*p(2);
r(:,2,3) = (-1/logten)*p(3);
r(:,3,4) = ( x(:,2)./x(:,3))*p(4);
r(:,3,5) = (-1)*p(5);
dfdpval  = r;

end

