function dfdpval = findif_var_dfdp(times,y,p,more)

x1 = more.var.fn(times,y,p,more.more);
x = zeros(size(x1),length(p));

for i = 1:length(p)
    tp = p;
    tp(i) = p(i) + more.eps;
    x(:,:,:,i) = (more.var.fn(times,y,tp,more.more)-x1)/more.eps;
end
dfdpval = x;

end

