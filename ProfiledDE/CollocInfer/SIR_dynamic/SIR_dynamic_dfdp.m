function dfdpval = SIR_dynamic_dfdp(t,x,p,more)
%  p-derivative of right side function for SIR differential equation
%  To be called by CollocInfer functions
%  x(1) = S,    x(2) = I,     x(3) = R
%  p(1) = beta, p(2) = nu, p(3) = mu
% Value of N passed by more.N

N = more.N;


r = zeros(length(t),3,3);
r(:,1,1) = -x(:,1).*x(:,2)/N;
r(:,1,3) = (N-x(:,1));
r(:,2,1) =  x(:,1).*x(:,2)/N;
r(:,2,2) = -x(:,2);
r(:,2,3) = -x(:,2);
r(:,3,2) =  x(:,2);
r(:,3,3) = -x(:,3);
dfdpval  = r;

end

