function expCproc = make_exp_Cproc()

expCproc.fn      = @exp_Cproc;
expCproc.dfdc    = @exp_dCproc_dc;
expCproc.dfdp    = @exp_dCproc_dp;
expCproc.d2fdc2  = @exp_d2Cproc_dc2;
expCproc.d2fdcdp = @exp_d2Cproc_dcdp;

end

