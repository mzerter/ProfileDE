function d2fdxdpval = SIR_dynamic_d2fdxdp(t,x,p,more)
%  x-p-derivative of right side function for SIR differential equation
%  To be called by CollocInfer functions
%  x(1) = S,    x(2) = I,     x(3) = R
%  p(1) = beta, p(2) = nu, p(3) = mu
% Value of N passed by more.N

N = more.N;

r = zeros(length(t),3,3,3);
r(:,1,1,1) = -x(:,2)/N;
r(:,1,1,3) = 1;

r(:,1,2,1) = -x(:,1)/N;

r(:,2,1,1) =  x(:,2)/N;
r(:,2,2,1) =  x(:,1)/N;
r(:,2,2,2) = -1;
r(:,2,2,3) = -1;

r(:,3,2,2) =  1;
r(:,3,3,3) =  -1;

d2fdxdpval = r;


end
