function d2fdxdpval = findif_ode_d2fdxdp(times,y,p,more)

x1 = findif_ode_dfdx(times,y,p,more);
[l,m,n] = size(x1);
x = zeros(l,m,n,length(p));

for i = 1:length(p)
    tp = p;
    tp(i) = p(i) + more.eps;
    x(:,:,:,i) = (findif_ode_dfdx(times,y,tp,more)-x1)/more.eps;
end
d2fdxdpval = x;

end
