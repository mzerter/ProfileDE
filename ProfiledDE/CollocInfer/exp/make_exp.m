%%%%%%%%%%% exponential of identity functions   %%%%%%%%%%%%%%%%%%

function expstruct = make_exp()

expstruct.fn      = @exp_fn;
expstruct.dfdx    = @exp_dfdx;
expstruct.dfdp    = @exp_dfdp;
expstruct.d2fdx2  = @exp_d2fdx2;
expstruct.d2fdxdp = @exp_d2fdxdp;

end

