function d2fdx2val = findif_loglik_d2fdx2(data,times,y,p,more)

x1 = findif_loglik_dfdx(data,times,y,p,more);
[m,n] = size(x1);
x = zeros(m,n,size(y,2));
for i = 1:size(y,2)
    ty = y;
    ty(:,i) = y(:,i) + more.eps;
    x(:,:,i) = (findif_loglik_dfdx(data,times,ty,p,more)-x1)/more.eps;
end
d2fdx2val = x;

end

