%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dpval = dSSE_proc_dp(coefs,bvals,pars,more)

if ~isfield(more,'more')
    more.more = [];
end
devals  = bvals.bvals * coefs;
ddevals = bvals.dbvals * coefs;
tmp = make_SSElik;
g = tmp.dfdp(ddevals,more.qpts,devals,pars,more);
g = sum(g);
dpval = g;

end

