function fnval = findif_loglik_fun(data,times,y,p,more)

fnval = more.fn(data,times,y,p,more.more);

end

