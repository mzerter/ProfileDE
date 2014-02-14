function dfdxval = SIVlog_dfdx(t,x,p,more)
%  x-derivative for SIV differential equation
%  To be called by CollocInfer functions
%  x(1) = S,    x(2) = I,     x(3) = V
%  p(1) = rho, p(2) = beta, p(3) = delta, p(4) = pi, p(5) = c
p = exp(p);
x = 10.^x;
logten = log(10);
r = zeros(length(t),3,3);
r(:,1,3) = -p(2)./logten;
r(:,2,1) = (p(2).*x(:,3))./(logten*x(:,2));
r(:,2,2) = -p(2).*x(:,1).*x(:,3)./(logten*x(:,2).^2);
r(:,2,3) =  p(2).*x(:,1)./(logten*x(:,2));
r(:,3,2) =  p(4)./x(:,3);
r(:,3,3) = -p(4).*x(:,2)./x(:,3).^2;
dfdxval = r;

end

