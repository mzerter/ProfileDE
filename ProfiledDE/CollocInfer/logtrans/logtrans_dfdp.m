function dfdpval = logtrans_dfdp(times,y,p,more)

y = exp(y);
x = more.dfdp(times,y,p,more.more);

[l,m,n] = size(x);
for i = 1:m
    for j = 1:n
        x(:,i,j) = x(:,i,j)./y(:,i);
    end
end

dfdpval = x;

end
