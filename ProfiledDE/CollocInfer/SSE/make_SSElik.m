function SSElik = make_SSElik()

SSElik.fn      = @SSE;
SSElik.dfdx    = @dSSE_dx;
SSElik.dfdy    = @dSSE_dy;
SSElik.dfdp    = @dSSE_dp;
SSElik.d2fdx2  = @d2SSE_dx2;
SSElik.d2fdxdy = @d2SSE_dxdy;
SSElik.d2fdy2  = @d2SSE_dy2;
SSElik.d2fdxdp = @d2SSE_dxdp;
SSElik.d2fdydp = @d2SSE_dydp;

end

