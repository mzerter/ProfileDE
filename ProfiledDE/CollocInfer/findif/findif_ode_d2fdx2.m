function d2fdx2val = findif_ode_d2fdx2(times,y,p,more)

x1 = findif_ode_dfdx(times,y,p,more);
[l,m,n] = size(x1);
x = zeros(l,m,n,size(y,2));

for i = 1:size(y,2)
    ty = y; 
    ty(:,i) = y(:,i) + more.eps;
    x(:,:,:,i) = (findif_ode_dfdx(times,ty,p,more)-x1)/more.eps;
end
d2fdx2val = x;

end

