function d2fdxdpval = exptrans_d2fdxdp(times,y,p,more)

y = exp(y);
d2fdxdpval = more.d2fdxdp(times,y,p,more.more);
[k,l,m,n] = size(d2fdxdpval);
for i = 1:l
    for j = 1:m
        for k = 1:n
            d2fdxdpval(:,i,j,k) = d2fdxdpval(:,i,j,k).*y(:,j);
        end
    end
end

end

