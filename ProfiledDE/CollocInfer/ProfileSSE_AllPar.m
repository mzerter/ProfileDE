function [res, Dres, ncoefs, dcdp] = ...
             ProfileSSE_AllPar(pars, times, data, coefs0, lik, proc, ...
                               in_method, options_in, dcdp, oldpars)
%  Inner optimization of coefficient values using the Gauss-Newton method 
%  for nonlinear least squares problems
%
%  Arguments:
%    The first six arguments are required.  The remaining arguments each
%    have default values, can be omitted.  All arguments after the first
%    omitted argument must be also omitted.
%  PARS       ...  A vector of parameter values.  These are initial values
%                  that are used to start the outer iterations.  Some
%                  of these may be fixed, and the remainder are optimized.
%                  The indices of the optimized parameters are contained
%                  in argument ACTIVE.
%  TIMES      ...  The times of observation of those variables that are
%                  measured.  These times are assumed to be common to all
%                  measured variables.
%  DATA       ...  Data matrix.  The number of columns is equal to the
%                  number of variables (including unobserved variables).
%                  The number of columns is equal to the number of 
%                  observation points per variable times the number of
%                  replications.
%  COEFS0     ...  A matrix containing coefficients for the expansion.  
%                  It has NVAR columns, and its number of rows is
%                  NBASIS*NVAR.
%  LIK        ...  A struct object set up by function LS_SETUP that
%                  defines the loss function for fitting the data
%  PROC       ...  A struct object set up by function LS_SETUP that
%                  defines the loss function for fitting the 
%                  differential equation.
%  IN_METHOD  ...  A string identifing the optimization function to be
%                  used in the outer optimization loop.  It defaults to [].
%  OPTIONS_IN ...  A struct object containing values controlling the 
%                  inner optimization process.
%  DCDP       ...  Matrix of derivatives of the coefficients with respect
%                  to the parameters.  This will not be available on
%                  the first call to the function, but is optionally
%                  available subsequently, and if present, is used for
%                  a single gradient step prior to the inner optimization.
%  OLDPARS    ...  A vector of previous parameter values.  
%                  This will not be available on
%                  the first call to the function, but is optionally
%                  available subsequently, and if present, is used for
%                  a single gradient step prior to the inner optimization.

%  Last modified 13 January 2011

if nargin < 10, oldpars    = [];  end
if nargin <  9, dcdp       = [];  end
if nargin <  8, options_in = [];  end
if nargin <  7, in_method  = [];  end

global INNEROPT_COEFS0

if ~isfield(lik, 'report')
    lik.report = [];
end

%  Evaluate the objective function at the current value of coefs

f1 = SplineCoefs(coefs0, times, data, pars, lik, proc);

% disp('lik object')
% disp(lik)
% disp(lik.more)
% disp(lik.more.more)
% 
% disp('proc object')
% disp(proc)
% disp(proc.more)
% disp(proc.more.more)

if isnan(f1)
    error('Initial evaluation of inner opt. criterion is NaN');
end

%  If a matrix DCDP containing derivatives of the coefficients with respect
%  to the parameter is available, as well as previous parameter values,
%  try a single gradient step prior to the inner optimization, 

if  ~isempty(dcdp)
    tcoefs = coefs0(:) - dcdp*(pars(:) - oldpars(:));
    f2 = SplineCoefs(tcoefs, times, data, pars, lik, proc);
    %  if this results in a reduced function value, update COEFS
    if f2 < f1
        coefs0 = tcoefs;
        f1 = f2;
    end
end

%  Invoke the inner optimization at the current values of COEFS
%  in order to initialize coefs for the outer optimization loop

ncoefs = inneropt(times, data, pars, lik, proc, ...
                  in_method, options_in);
            
INNEROPT_COEFS0 = ncoefs;

devals = lik.bvals * ncoefs;

%  Partial second derivatives

[f, dfdc, d2fdc2, d2fdcdp] = ...
    SplineCoefs(ncoefs, times, data, pars, lik, proc);

%  Compute the weighted residuals

if isfield(more,'whichobs')
    whichobs = more.whichobs;
else
    whichobs = [];
end
weights  = checkweights(lik.more.weights, whichobs, data);
residual = data - lik.more.fn(times, devals, pars, lik.more.more);
res      = residual.*sqrt(weights);
res      = res(:);
isnaf    = isnan(res);
res(isnaf) = 0;

% Matrix of derivatives of fit values with respect to parameters

dlikdp  = lik.more.dfdp(times, devals, pars, lik.more.more);
[l,m,n] = size(dlikdp);
dlikdp  = reshape(dlikdp,l*m,n);

%  Matrix of derivatives of fit values with respect to path

dlikdx  = lik.more.dfdx(times, devals, pars, lik.more.more);

%  Matrix of derivatives of fit values with respect to coefficients

dlikdc = [];
[l,m,n]  = size(dlikdx);
for i = 1:m
    tH = [];
    for j = 1:n
        tH = [tH, diag(dlikdx(:,i,j)) * lik.bvals];
    end
    dlikdc = [dlikdc; tH];
end

%  Matrix DCDP of total derivatives of coefficients with respect to
%  parameters

dcdp = d2fdc2\d2fdcdp;

%  Matrix DCDP of total derivatives of residuals with respect to
%  parameters

Dres = dlikdc * dcdp + dlikdp;
Dres(isnan(Dres)) = 0;

end
