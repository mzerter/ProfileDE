function d2fdx2val = SIV_d2fdx2(t,x,p,more)
%  second x-derivative for SIV differential equation
%  To be called by CollocInfer functions
%  x(1) = S,    x(2) = I,     x(3) = V
%  p(1) = rho, p(2) = beta, p(3) = delta, p(4) = pi, p(5) = c
r = zeros(length(t),3,3,3);
r(:,1,1,3) =  -p(2);
r(:,1,3,1) =  -p(2);
r(:,2,1,3) =   p(2);
r(:,2,3,1) =   p(2);
d2fdx2val = r;

end

