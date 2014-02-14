function d2fdx2val = fhn_d2fdx2(times,y,p,more)

[m,n] = size(y);
r = zeros(m,n,2,2);
r(:,1,1,1) = -2*p(3).*y(:,1);
d2fdx2val = r;

end
