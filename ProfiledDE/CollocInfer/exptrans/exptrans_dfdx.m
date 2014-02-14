function dfdxval = exptrans_dfdx(times,y,p,more)

y  = exp(y);
dfdxval  = more.dfdx(times,y,p,more.more);
[l,m,n] = size(dfdx);
for i = 1:m
    for j = 1:n
        dfdxval(:,i,j) = dfdxval(:,i,j).*y(:,j);
    end
end

end
