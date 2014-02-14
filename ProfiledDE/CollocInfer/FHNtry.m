
%  add paths to required functions

addpath('../../fdaM')
addpath('fhn')
addpath('SSE')
addpath('id')

%  set up objects required for analysis

times = (0:0.5:20)';
FHN_xnames = ['V','R'];
FHN_knots  = (0:2:20)';
FHN_fn = make_fhn;
FHN_order  = 3;
FHN_nbasis = length(FHN_knots) + FHN_order - 2;
FHN_range  = [0,20];
FHN_basis  = create_bspline_basis(FHN_range, FHN_nbasis, ...
                                  FHN_order, FHN_knots);
FHN_pars0 = [0.2, 0.2, 3.0];
% odeopts = odeset('RelTol',1e-13);
% [x,y] = ode45(FHN_fn.fn_ode, times, [-1,1], odeopts, FHN_pars0);
% FHN_data = y + 0.5.*randn(41,2);
load FHN_data.txt
FHN_fdnames = cell(1,3);                                  
FHN_fdnames{3} = FHN_xnames;
FHN_xfd = smooth_basis(times, FHN_data, FHN_basis, ones(41,1), ...
                       FHN_fdnames);
FHN_coefs0 = getcoef(FHN_xfd);
FHN_lambda = 1e4*[1,1]';

%  set up lik and proc objects for least squares fitting

[lik, proc, coefs] = ...
    LS_setup(FHN_fn, times, FHN_coefs0, FHN_basis, FHN_lambda);

%  compute initial function and derivatives

[f, dfdc, d2fdc2, d2fdcdp] = ...
    SplineCoefs(FHN_coefs0, times, FHN_data, FHN_pars0, lik, proc);

%  optimize initial coefficient values by Newton-Raphson method

options_in.maxit      = 100;
options_in.trace      = 2;

coefs_opt = ...
    SplineEst_NewtRaph(times, FHN_data, FHN_coefs0, FHN_pars0, ...
                       lik, proc, options_in);

%  optimizing using function inneropt

global INNEROPT_COEFS0
INNEROPT_COEFS0 = FHN_coefs0;

options_in.Display    = 'on';
options_in.MaxIter    = 200;
options_in.LargeScale = 'on';
options_in.GradObj    = 'on';
options_in.Hessian    = 'on';

tic;
coefs_opt = inneropt(times, FHN_data, FHN_pars0, lik, proc, options_in);
toc

%  Profile_LS argument correspondences

data      = FHN_data;
fn        = FHN_fn;
pars      = FHN_pars0;
coefs     = FHN_coefs0;
coefs_0   = FHN_coefs0;
basisvals = FHN_basis;
lambda    = FHN_lambda;
fdobj       = [];
more        = [];
weights     = []; 
quadrature  = [];
active      = [];
in_method   = [];
out_method  = [];
options_in  = [];
options_out = [];
eps         = 1e-6;
poslik      = 0; 
posproc     = 0; 
discrete    = 0; 

%  default values within Profile_GaussNewt

active=1:length(FHN_pars0);   

coefs0 = coefs;

options_in = optimset('Large', 'on', 'Display','iter', 'TolFun', 1e-5, ...
                             'GradObj','on', 'Hessian','on',   ...
                             'maxit', 100, 'DiffMaxChange', 1e-4);

in_method = 'SplineEst';
in_method = [];
options_out.trace = 2;

tic;
[newpars, res0, F0, coefs] = ...
    Profile_GausNewt(times, FHN_data, FHN_coefs0, FHN_pars0, lik, proc, ...  
                     in_method, options_in, options_out, active);
toc

out_method = [];

tic;
Profile_sse_struct = ...
    Profile_LS(FHN_fn, times, FHN_data, FHN_coefs0, FHN_pars0, ...
               FHN_basis, FHN_lambda, [], [], [], [], active, ...
               in_method, out_method, options_in, options_out);
toc


               