function fnval = logstate_lik_fun(data,times,y,p,more)

y = exp(y);
fnval = more.fn(data,times,y,p,more.more);

end
