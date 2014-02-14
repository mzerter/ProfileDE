function dfdxval = id_dfdx(t,x,pars,more)

[m,n] = size(x);
g = zeros(m,n,n);
for i = 1:n
    g(:,i,i) = 1;
end
dfdxval = g;

end

