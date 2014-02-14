function dfdpval = exptrans_dfdp(times,y,p,more)

y = exp(y);
dfdpval = more.dfdp(times,y,p,more.more);
     
end
