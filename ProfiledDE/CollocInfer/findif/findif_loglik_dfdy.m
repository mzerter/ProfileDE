function dfdyval = findif_loglik_dfdy(data,times,y,p,more)  % deriv wrt response

x1 = more.fn(data,times,y,p,more.more);
x = data;

for i = 1:size(data,2)
    tdata = data;
    tdata(:,i) = data(:,i) + more.eps;
    x(:,i) = (more.fn(tdata,times,y,p,more.more)-x1)/more.eps;
end
dfdyval = x;

end

