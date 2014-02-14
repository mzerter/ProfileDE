function d2fdxdpval = logtrans_d2fdxdp(times,y,p,more)

x1 = logtrans_dfdp(times,y,p,more);
y = exp(y);
x = more.d2fdxdp(times,y,p,more.more);
[kk,ll,mm,nn] = size(x);

for i = 1:ll
    for j = 1:mm
        for k = 1:nn
            x(:,i,j,k) = x(:,i,j,k).*y(:,j)./y(:,i);
        end
    end
end

for i = 1:ll
    for j = 1:nn
        x(:,i,i,j) = x(:,i,i,j) - x1(:,i,j);
    end
end

d2fdxdpval = x;

end
