function funval = exptrans_fun(times,y,p,more)

y = exp(y);
funval = more.fn(times,y,p,more.more);

end
