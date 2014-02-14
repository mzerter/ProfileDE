% -----------  proc object for loggenlin analysis  ---------------------
% 
% proc                 %  SSE functions for evaluating fit to DE
% proc.more            %  findif   functions for difference-derivatives
% proc.more.more       %  logtrans functions for exponential transformation
% proc.more.more.more  %  chemo functions for right side of DE
% 
% proc = 
% 
%          fn: @SSE_proc
%        dfdc: @dSSE_proc_dc
%        dfdp: @dSSE_proc_dp
%      d2fdc2: @d2SSE_proc_dc2
%     d2fdcdp: @d2SSE_proc_dcdp
%       bvals: [1x1 struct]
%        more: [1x1 struct]
% 
% 
% ans = 
% 
%           fn: @findif_ode_fn
%         dfdx: @findif_ode_dfdx
%         dfdy: @findif_ode_dfdy
%         dfdp: @findif_ode_dfdp
%       d2fdx2: @findif_ode_d2fdx2
%      d2fdxdp: @findif_ode_d2fdxdp
%         qpts: [1x202 double]
%      weights: [5x1 double]
%        names: [5x2 char]
%     parnames: [16x6 char]
%         more: [1x1 struct]
% 
% 
% ans = 
% 
%       fn: @logtrans_fun
%      eps: 1.0000e-008
%     more: [1x1 struct]
% 
% 
% ans = 
% 
%     fn: @chemo_fun
%     
%     
% ----------------  lik object for loggenlin analysis  -------------------
% 
% lik                 %  logstate functions for exponentiating log state
% lik.more            %  SSE functions for evaluating fit
% lik.more.more       %  genlin functions fitting composite observations
% lik.more.more.more  %  mat and sub matrices for genlin functions
% 
% lik = 
% 
%          fn: @logstate_lik_fun
%        dfdx: @logstate_lik_dfdx
%        dfdp: @logstate_lik_dfdp
%      d2fdx2: @logstate_lik_d2fdx2
%     d2fdxdp: @logstate_lik_d2fdxdp
%        more: [1x1 struct]
%       bvals: [101x203 double]
% 
% 
% ans = 
% 
%          fn: @SSE
%        dfdx: @dSSE_dx
%        dfdy: @dSSE_dy
%        dfdp: @dSSE_dp
%      d2fdx2: @d2SSE_dx2
%     d2fdxdy: @d2SSE_dxdy
%      d2fdy2: @d2SSE_dy2
%     d2fdxdp: @d2SSE_dxdp
%     d2fdydp: @d2SSE_dydp
%        more: [1x1 struct]
% 
% 
% ans = 
% 
%     fun_ode: @loggenlin_fun_ode
%          fn: @loggenlin_fn
%        dfdx: @loggenlin_dfdx
%        dfdp: @loggenlin_dfdp
%      d2fdx2: @loggenlin_d2fdx2
%     d2fdxdp: @loggenlin_d2fdxdp
%        more: [1x1 struct]
%     weights: [100 1]
% 
% 
% ans = 
% 
%     mat: [2x5 double]
%     sub: [4x3 double]
% 
