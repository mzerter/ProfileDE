function dfdxval = findif_ode_dfdx(times,y,p,more)

x1 = more.fn(times,y,p,more.more);
[m,n] = size(x1);
x = zeros(m,n,size(y,2));

for i = 1:size(y,2)
    ty = y;
    ty(:,i) = y(:,i) + more.eps;
    x(:,:,i) = (more.fn(times,ty,p,more.more)-x1)/more.eps;
end
dfdxval = x;

end

