function dfdpval = findif_ode_dfdp(times,y,p,more)

x1 = more.fn(times,y,p,more.more);
[m,n] = size(x1);
x = zeros(m,n,length(p));

for i = 1:length(p)
    tp = p;
    tp(i) = p(i) + more.eps;
    x(:,:,i) = (more.fn(times,y,tp,more.more)-x1)/more.eps;
end
dfdpval = x;

end

