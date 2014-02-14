function dfdxval = findif_loglik_dfdx(data,times,y,p,more)  

x1 = more.fn(data,times,y,p,more.more);
x = y;

for i = 1:size(y,2)
    ty = y;
    ty(:,i) = y(:,i) + more.eps;
    x(:,i) = (more.fn(data,times,ty,p,more.more)-x1)/more.eps;
end
dfdxval = x;

end

