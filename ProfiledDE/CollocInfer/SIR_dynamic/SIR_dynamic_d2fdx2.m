function d2fdx2val = SIR_dynamic_d2fdx2(t,x,p,more)
%  second x-derivative of right side function for SIR differential equation
%  To be called by CollocInfer functions
%  x(1) = S,    x(2) = I,     x(3) = R
%  p(1) = beta, p(2) = nu, p(3) = mu
% Value of N passed by more.N

N = more.N;


r = zeros(length(t),3,3,3);
r(:,1,1,2) = -p(1)/N;
r(:,1,2,1) = -p(1)/N;
r(:,2,1,2) =  p(1)/N;
r(:,2,2,1) =  p(1)/N;
d2fdx2val = r;

end

