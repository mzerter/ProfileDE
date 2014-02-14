function dfdxval = logtrans_dfdx(times,y,p,more)

x1 = logtrans_fun(times,y,p,more);
y = exp(y);
x = more.dfdx(times,y,p,more.more);
[l,m,n] = size(x);
for i = 1:m
    for j = 1:n
        x(:,i,j) = x(:,i,j).*y(:,j)./y(:,i);
    end
    x(:,i,i) = x(:,i,i) - x1(:,i);
end
dfdxval = x;

end
