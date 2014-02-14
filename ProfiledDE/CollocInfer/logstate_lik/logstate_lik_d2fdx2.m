function d2fdx2val = logstate_lik_d2fdx2(data,times,y,p,more)

x1 = logstate_lik_dfdx(data,times,y,p,more);
y  = exp(y);
x  = more.d2fdx2(data,times,y,p,more.more);
for i = 1:size(x,2)
    for j = 1:size(x,3)
        x(:,i,j) = x(:,i,j).*y(:,i).*y(:,j);
    end
    x(:,i,i) = x(:,i,i) + x1(:,i);
end
d2fdx2val = x;

end
