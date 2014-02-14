function d2fdxdpval = findif_var_d2fdxdp(times,y,p,more)

x1 = findif_var_dfdx(times,y,p,more);
x  = zeros(size(x1),length(p));

for i = 1:length(p)
    tp = p;
    tp(i) = p(i) + more.eps;
    x(:,:,:,:,i) = (findif_var_dfdx(times,y,tp,more)-x1)/more.eps;
end
d2fdxdpval = x;

end

