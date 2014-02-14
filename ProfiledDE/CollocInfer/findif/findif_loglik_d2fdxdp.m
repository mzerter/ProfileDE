function d2fdxdpval = findif_loglik_d2fdxdp(data,times,y,p,more)

x1 = findif_loglik_dfdx(data,times,y,p,more);
[m,n] = size(x1);
x = zeros(m,n,length(p));
for i = 1:length(p)
    tp = p;
    tp(i) = p(i) + more.eps;
    x(:,:,i) = (findif_loglik_dfdx(data,times,y,tp,more)-x1)/more.eps;
end
d2fdxdpval = x;

end

