function d2fdx2val = findif_var_d2fdx2(times,y,p,more)
n = size(y,2);
x1 = findif_var_dfdx(times,y,p,more);
x = zeros(size(x1),n);
for i = 1:n
    ty = y;
    ty(:,i) = y(:,i) + more.eps;
    x(:,:,:,:,i) = (findif_var_dfdx(times,ty,p,more)-x1)/more.eps;
end
d2fdx2val = x;

end

