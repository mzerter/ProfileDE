function d2fdx2val = Lorenz_d2fdx2(t,x,p,more)
%  second x-derivative of right side function for Lorenz equation
%  x(1) = x,    x(2) = y,     x(3) = z
%  p(1) = sigma, p(2) = r, p(3) = b
p = exp(p);
r = zeros(length(t),3,3,3);
r(:,2,3,1) = 1;
r(:,2,1,3) = 1;
r(:,3,1,2) = 1;
r(:,3,2,1) = 1;

end

