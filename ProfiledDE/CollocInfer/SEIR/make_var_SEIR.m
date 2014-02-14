
function var_SEIR = make_var_SEIR()

var_SEIR.var.fn      = @SEIR_var_fun;
var_SEIR.var.dfdx    = @SEIR_var_dfdx;
var_SEIR.var.dfdp    = @SEIR_var_dfdp;
var_SEIR.var.d2fdx2  = @SEIR_var_d2fdx2;
var_SEIR.var.d2fdxdp = @SEIR_var_d2fdxdp;

end

