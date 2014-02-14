function dfdxval = findif_var_dfdx(times,y,p,more)

x1 = more.var.fn(times,y,p,more.more);
x = zeros(size(x1),size(y,2));

for i = 1:size(y,2)
    ty = y;
    ty(:,i) = y(:,i) + more.eps;
    x(:,:,:,i) = (more.var.fn(times,ty,p,more.more)-x1)/more.eps;
end
dfdxval = x;

end

