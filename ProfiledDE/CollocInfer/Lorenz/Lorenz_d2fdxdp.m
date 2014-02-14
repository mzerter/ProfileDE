function d2fdxdpval = Lorenz_d2fdxdp(t,x,p,more)
%  x-p-derivative of right side function for Lorenz equation
%  x(1) = x,    x(2) = y,     x(3) = z
%  p(1) = sigma, p(2) = r, p(3) = b
p = exp(p);
r = zeros(length(t),3,3,3);
r(:,1,1,3) = -1;
r(:,1,2,3) =  1;
r(:,2,1,2) =  1;
r(:,3,3,3) = -1;

d2fdxdpval = r;

end
