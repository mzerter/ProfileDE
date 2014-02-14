function logtrans = make_logtrans()

logtrans.fn      = @logtrans_fun;
logtrans.dfdx    = @logtrans_dfdx;
logtrans.dfdp    = @logtrans_dfdp;
logtrans.d2fdx2  = @logtrans_d2fdx2;
logtrans.d2fdxdp = @logtrans_d2fdxdp;
logtrans.extras  = [];

end
