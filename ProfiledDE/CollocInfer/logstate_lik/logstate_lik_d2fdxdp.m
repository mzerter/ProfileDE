function d2fdxdpval = logstate_lik_d2fdxdp(data,times,y,p,more)

y = exp(y);
x = more.d2fdxdp(data,times,y,p,more.more);
for i = 1:size(x,2)
    for j = 1:size(x,3)
        x(:,i,j) = x(:,i,j).*y(:,i);
    end
end
d2fdxdpval = x;

end
