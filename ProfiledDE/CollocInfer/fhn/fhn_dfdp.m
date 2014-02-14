function dfdpval = fhn_dfdp(times,y,p,more)

[m,n] = size(y);
r = zeros(m,n,length(p));
r(:,1,3) = y(:,1) - y(:,1).^3./3 + y(:,2);
r(:,2,1) = 1/p(3);
r(:,2,2) = -y(:,2)./p(3);
r(:,2,3) = (y(:,1) - p(1)+p(2).*y(:,2))./(p(3)^2);
dfdpval  = r;

end
