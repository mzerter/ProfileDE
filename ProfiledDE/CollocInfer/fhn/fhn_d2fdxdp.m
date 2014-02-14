function d2fdxdpval = fhn_d2fdxdp(times,y,p,more)

[m,n] = size(y);
r = zeros(m,n,2,length(p));
r(:,1,1,3) = 1 - y(:,1).^2;
r(:,1,2,3) = 1;
r(:,2,1,3) = 1/p(3).^2;
r(:,2,2,2) = -1/p(3);
r(:,2,2,3) = p(2)/p(3).^2;
d2fdxdpval = r;

end
