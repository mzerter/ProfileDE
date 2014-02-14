%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pars_opt, ncoefs, value, gradient] = ...
         outeropt(times, data, coefs, allpars, lik, proc, active, ...
                  in_method, options_in, out_method, options_out)
%  The multivariate data = argment data are observations of one or more
%  functional variables or processes at time points = argument times.
%
%  The data are smoothed using differential operators.  
%
%  The operators typically depend on a number of parameters whose values
%  are = argument pars.  The function optimizes these parameter values.
%
%  The fitting criterion for both the data and the differential equation
%  is user-defined = arguments lik and proc, respectively.
%
%  The function optimizes the coefficients of the basis expansions for the
%  representations of the variables or processes = what is called the
%  inner optimization loop, implemented by function Inneropt.  It optimizes the
%  parameter values = the outer optimization loop.
%
%  The optimization method for the inner loop is defined by argument in_meth.
%  The optimization method for the outer loop is defined by argument out_meth.
%
%  The function returns the optimized parameter values, a functional data object defining 
%  the smoothing functions, the likelihood and process functions, and the
%  residual values.
%
%  Arguments:
%    The first six arguments are required.  The remaining arguments each
%    have default values, can be omitted.  All arguments after the first
%    omitted argument must be also omitted.
%  TIMES      ...  The times of observation of those variables that are
%                  measured.  These times are assumed to be common to all
%                  measured variables.
%  DATA       ...  Data matrix.  The number of columns is equal to the
%                  number of variables (including unobserved variables).
%                  The number of columns is equal to the number of 
%                  observation points per variable times the number of
%                  replications.
%  COEFS      ...  A matrix containing coefficients for the expansion.  
%                  It has NVAR columns, and its number of rows is
%                  NBASIS*NVAR.
%  ALLPARS    ...  A vector of parameter values.  These are initial values
%                  that are used to start the outer iterations.  Some
%                  of these may be fixed, and the remainder are optimized.
%                  The indices of the optimized parameters are contained
%                  in argument ACTIVE.
%  LIK        ...  A struct object set up by function LS_SETUP that
%                  defines the loss function for fitting the data
%  PROC       ...  A struct object set up by function LS_SETUP that
%                  defines the loss function for fitting the 
%                  differential equation.
%  ACTIVE     ...  A vector containing the indices of the subset of
%                  parameters to be optimized.  If empty, it defaults
%                  to all indices.
%                  Defaults to [].
%  IN_METHOD  ...  A string identifing the optimization function to be
%                  used in the outer optimization loop.  It defaults to [].
%                  Defaults to [].
%  OPTIONS_IN ...  A struct object containing values controlling the 
%                  inner optimization process.
%  OUT_METHOD ...  A string identifing the optimization function to be
%                  used in the inner optimization loop.  It defaults to [].
%                  Defaults to [].
%  OPTIONS_OUT...  A struct object containing values controlling the 
%                  outer optimization process.  
%                  If not empty, field names are:
%                      EPS or TolFun:  convergence criterion for function 
%                          value.  
%                          Defaults to 1e-6
%                      maxit or MaxIter:  limit on the number of iterations
%                          Defaults to 100
%                  Defaults to [].

%  Last modified 10 January 2011

if nargin < 11, options_out = [];  end
if nargin < 10, out_method  = [];  end
if nargin <  9, options_in  = [];  end
if nargin <  8, in_method   = [];  end
if nargin <  7, active      = [];  end

global INNEROPT_COEFS0

INNEROPT_COEFS0 = coefs;

if isempty(active)
    active = 1:length(allpars);
end

if isempty(options_out)
    options_out = optimset('GradObj', 'on', 'Hessian', 'off', ...
        'MaxIter', 100, 'display', 'iter', 'TolFun', 1e-5);
else
    options_out.Hessian = 'off';
end

pars0 = allpars(active);
newpars = fminunc(@ProfileErr, pars0, options_out, ...
                  times, data, coefs, allpars, lik, proc, ...
                  active, in_method, options_in);
pars_opt = allpars;
pars_opt(active) = newpars;

[value, gradient, ncoefs] = ...
    ProfileErr_AllPar(pars_opt, times, data, coefs, lik, proc, ...
                  in_method, options_in);

gradient = gradient(active);

end

