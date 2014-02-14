function d2fdx2val = exptrans_d2fdx2(times,y,p,more)

d2fdx2val1 = exptrans_dfdx(times,y,p,more);
y = exp(y);
d2fdx2val = more.d2fdx2(times,y,p,more.more);
[k,l,m,n] = size(d2fdx2val);
for i = 1:l
    for j = 1:m
        for k = 1:n
            d2fdx2val(:,i,j,k) = d2fdx2val(:,i,j,k).*y(:,j).*y(:,k);
        end
        d2fdx2val(:,i,j,j) = d2fdx2val(:,i,j,j) + d2fdx2val1(:,i,j);
    end
end

end
