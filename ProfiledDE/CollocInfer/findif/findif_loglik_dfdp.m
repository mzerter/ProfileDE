function dfdpval = findif_loglik_dfdp(data,times,y,p,more)

x1 = more.fn(data,times,y,p,more.more);
x = zeros(length(x1),length(p));
for i = 1:length(p)
    tp = p;
    tp(i) = p(i) + more.eps;
    x(:,i) = (more.fn(data,times,y,tp,more.more)-x1)/more.eps;
end
dfdpval = x;

end

