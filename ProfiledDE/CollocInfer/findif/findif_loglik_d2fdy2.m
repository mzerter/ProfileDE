function d2fdy2 = findif_loglik_d2fdy2(data,times,y,p,more)

x1 = findif_loglik_dfdy(data,times,y,p,more);
x = zeros(size(x1),size(data,2));
for i = 1:size(data,2)
    tdata = data;
    tdata(:,i) = data(:,i) + more.eps;
    x(:,:,i) = (findif_loglik_dfdy(tdata,times,y,p,more)-x1)/more.eps;
end
d2fdy2 = x;

end

