function SSEproc = make_SSEproc()

SSEproc.fn      = @SSE_proc;
SSEproc.dfdc    = @dSSE_proc_dc;
SSEproc.dfdp    = @dSSE_proc_dp;
SSEproc.d2fdc2  = @d2SSE_proc_dc2;
SSEproc.d2fdcdp = @d2SSE_proc_dcdp;

end

