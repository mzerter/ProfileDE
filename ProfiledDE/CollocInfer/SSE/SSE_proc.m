function fnval = SSE_proc(coefs,bvals,pars,more)

devals  = bvals.bvals  * coefs;
ddevals = bvals.dbvals * coefs;
tmp     = make_SSElik;
f       = tmp.fn(ddevals,more.qpts,devals,pars,more);
fnval   = sum(f);

end

