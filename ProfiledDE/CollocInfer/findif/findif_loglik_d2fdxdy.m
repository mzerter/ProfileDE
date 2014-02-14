function d2fdxdyval = findif_loglik_d2fdxdy(data,times,y,p,more)

x1 = findif_loglik_dfdx(data,times,y,p,more);
[m,n] = size(x1);
x = zeros(m,n,size(data,2));
for i = 1:nsize(data)
    tdata = data;
    tdata(:,i) = data(:,i) + more.eps;
    x(:,:,i) = (findif_loglik_dfdx(tdata,times,y,p,more)-x1)/more.eps;
end
d2fdxdyval = x;

end

