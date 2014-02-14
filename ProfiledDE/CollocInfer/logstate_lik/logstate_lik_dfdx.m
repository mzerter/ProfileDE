function dfdxval = logstate_lik_dfdx(data,times,y,p,more)

y = exp(y);
x = more.dfdx(data,times,y,p,more.more);
for i = 1:size(x,2)
    x(:,i) = x(:,i).*y(:,i);
end
dfdxval = x;

end
