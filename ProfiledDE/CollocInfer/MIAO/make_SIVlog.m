function SIVlog = make_SIVlog()

SIVlog.ode.fn  = @SIVlog_ode;
SIVlog.fn      = @SIVlog_fun;
SIVlog.dfdx    = @SIVlog_dfdx;
SIVlog.dfdp    = @SIVlog_dfdp;
SIVlog.d2fdx2  = @SIVlog_d2fdx2;
SIVlog.d2fdxdp = @SIVlog_d2fdxdp;

end
