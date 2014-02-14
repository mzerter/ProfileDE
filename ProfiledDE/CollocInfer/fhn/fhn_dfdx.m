function dfdxval = fhn_dfdx(times,y,p,more)

[m,n] = size(y);
r = zeros(m,n,2);
r(:,1,1) =  p(3) - p(3).*y(:,1).^2;
r(:,1,2) =  p(3);
r(:,2,1) = -1/p(3);
r(:,2,2) = -p(2)/p(3);
dfdxval  = r;

end
