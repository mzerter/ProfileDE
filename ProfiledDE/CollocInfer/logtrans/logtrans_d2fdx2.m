function d2fdx2val = logtrans_d2fdx2(times,y,p,more)

x1 = logtrans_dfdx(times,y,p,more);
y = exp(y);
x = more.d2fdx2(times,y,p,more.more);

[kk,ll,mm,nn] = size(x);
for i = 1:ll
    for j = 1:mm
        for k = 1:nn
            x(:,i,j,k) = x(:,i,j,k).*y(:,j).*y(:,k)./y(:,i);
        end
    end
end

for i = 1:ll
    for j = 1:mm
        x(:,i,i,j) = x(:,i,i,j) - x1(:,i,j);
        x(:,i,j,i) = x(:,i,j,i) - x1(:,i,j);
        x(:,i,j,j) = x(:,i,j,j) + x1(:,i,j);
    end
end

d2fdx2val = x;

end
