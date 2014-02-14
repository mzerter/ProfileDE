function d2fdx2val = exp_d2fdx2(t,x,pars,more)

n = size(x,2);
[l,m] = size(x);
g =  zeros(l,m,n,n);
for i = 1:n
    g(:,i,i,i) = exp(x(:,i));
end
d2fdx2val = g;

end

