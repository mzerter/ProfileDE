function dfdxval = SIR_dynamic_dfdx(t,x,p,more)
%  x-derivative of right side function for SIR differential equation
%  To be called by CollocInfer functions
%  x(1) = S,    x(2) = I,     x(3) = R
%  p(1) = beta, p(2) = nu, p(3) = mu
% Value of N passed by more.N

N = more.N;


r = zeros(length(t),3,3);

r(:,1,1) = -p(1).*x(:,2)/N + p(3);
r(:,1,2) = -p(1).*x(:,1)/N;
r(:,2,1) =  p(1).*x(:,2)/N;
r(:,2,2) =  p(1).*x(:,1)/N - (p(2)+p(3));
r(:,3,2) =  p(2);
r(:,3,3) = -p(3);
dfdxval = r;

end

