function exptransstruct = make_exptrans()

exptransstruct.fn      = @exptrans_fun;
exptransstruct.dfdx    = @exptrans_dfdx;
exptransstruct.dfdp    = @exptrans_dfdp;
exptransstruct.d2fdx2  = @exptrans_d2fdx2;
exptransstruct.d2fdxdp = @exptrans_d2fdxdp;
exptransstruct.extras  = [];

end
