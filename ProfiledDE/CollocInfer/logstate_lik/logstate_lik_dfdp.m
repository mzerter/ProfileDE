function dfdpval = logstate_lik_dfdp(data,times,y,p,more)

y = exp(y);
dfdpval = more.dfdp(data,times,y,p,more.more);

end
