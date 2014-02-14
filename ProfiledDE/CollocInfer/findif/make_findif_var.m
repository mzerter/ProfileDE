%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite differencing for process variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function findif_var = make_findif_var()

findif_var.fn      = @findif_var_fun;
findif_var.dfdx    = @findif_var_dfdx;
findif_var.dfdp    = @findif_var_dfdp;
findif_var.d2fdx2  = @findif_var_d2fdx2;
findif_var.d2fdxdp = @findif_var_d2fdxdp;

end

