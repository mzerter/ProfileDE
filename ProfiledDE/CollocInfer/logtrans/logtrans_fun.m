function xovery = logtrans_fun(times,y,p,more)
y = exp(y);
x = more.fn(times,y,p,more.more);
xovery = x./y;

end
