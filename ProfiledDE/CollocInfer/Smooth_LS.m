%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fdobj, ncoefs, lik, proc] = ...
    Smooth_LS(fn, times, data,  coefs, pars, basisvals, lambda, ...
              fdobj, more, obsweights, quadrature, ...
              in_method, options_in, diffeps, poslik, posproc, discrete)
%  The multivariate data = argment data are observations of one or more
%  functional variables or processes at time points = argument times.
%
%  The data are smoothed using differential operators.  These operators are
%  defined by the functions = argument fn.  The smooth can optionally be
%  required to be positive, set by argument pos.
%
%  The operators typically depend on a number of parameters whose values
%  are = argument pars.
%
%  The level of smoothing is defined by a scalar or vector of smoothing
%  parameters = argument lambda.
%
%  The fitting criterion for both the data and the differential equation
%  is weighted error sum of squares, and consequently the function first calls the
%  function LS_setup.
%
%  The function optimizes the coefficients of the basis expansions for the
%  representations of the variables or processes = what is called the
%  inner optimization loop, implemented by function Inneropt.  The
%  optimization method is defined by argument in_meth.
%
%  The function returns a functional data object defining the smoothing
%  functions, along with the likelihood of process functions and the
%  residual values.
%
%  call LS_setup to set up likelihood and process functions for
%  least squares estimation using ODE evaluation functions = fn
%
%  Arguments:
%    The first six arguments are required.  The remaining arguments each
%    have default values, can be omitted.  All arguments after the first
%    omitted argument must be also omitted.
%  FN         ...  A struct object containing handles to functions that
%                  evaluate the value of the right side of the 
%                  differential equation as well as the values of a 
%                  number of its partial derivatives.  These functions
%                  must be provided by the user.  However, difference
%                  approximations of derivatives may also be requested.
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
%  PARS       ...  A vector of parameter values.  These are initial values
%                  that are used to start the outer iterations.  Some
%                  of these may be fixed, and the remainder are optimized.
%                  The indices of the optimized parameters are contained
%                  in argument ACTIVE.
%  BASISVALS  ...  Either a basis object to be used to approximate each
%                  solution, or a struct object containing the values
%                  of the basis function at each observation time.
%  LAMBDA     ...  A single constant, or a vector of length NVAR containing
%                  smoothing parameter values for each variable controlling
%                  how closely a variable comes to solving the equation.
%                  Defaults to 0.
%  FDOBJ      ...  An alternative means of contributing a basis system.
%                  When this object is supplied, its basis overides that
%                  supplied in argument BASISVALS.
%                  Defaults to [].
%  MORE       ...  Any additional information required by the functions
%                  referenced in argument FN.
%                  Defaults to [].
%  OBSWEIGHTS ...  Weights to be applied in computing error sums of 
%                  squares over the observed variables.
%                  Defaults to [].
%  QUADRATURE ...  Locations of quadrature points for approximating the
%                  integrals defining the penalty terms.  If not 
%                  supplied, quadrature points are positioned at the 
%                  centers of the inter-knot intervals.
%                  Defaults to [].
%  ACTIVE     ...  A vector containing the indices of the subset of
%                  parameters to be optimized.  If empty, it defaults
%                  to all indices.
%                  Defaults to [].
%  IN_METHOD  ...  A string identifing the optimization function to be
%                  used in the outer optimization loop.  It defaults to [].
%                  Defaults to [].
%  OPTIONS_IN ...  A struct object containing values controlling the 
%                  inner optimization process.
%  DIFFEPS    ...  A small constant value used to compute difference 
%                  approximations to derivatives if this is requested.
%  POSLIK     ...  If nonzero, the functions approximating the data are
%                  constrained to be positive.
%  POSPROC    ...  If nonzero, the functions used in the penalty terms
%                  are constrained to be positive.
%  DISCRETE   ...  If nonzero, a discrete time model is used instead of
%                  the default continuous time model.

%  Last modified 8 January 2011

if nargin < 17, discrete    = 0;     end
if nargin < 16, posproc     = 0;     end
if nargin < 15, poslik      = 0;     end
if nargin < 14, diffeps     = 1e-6;  end
if nargin < 13, options_in  = [];    end
if nargin < 12, in_method   = [];    end
if nargin < 11, quadrature  = [];    end
if nargin < 10, obsweights     = [];    end
if nargin <  9, more        = [];    end
if nargin <  8, fdobj       = [];    end
if nargin <  7, basisvals   = [];    end

[lik, proc, coefs] = ...
    LS_setup(fn, times, coefs, basisvals, lambda, fdobj,...
             more, obsweights, quadrature, ...
             diffeps, poslik, posproc, discrete);

%  Carry out an inner optimization to estimate coefficients for
%  processes given parameter values

ncoefs = inneropt(times, data, pars, lik, proc, options_in, in_method);
nbasis = size(lik.bvals,2);
ncoefs  = reshape(ncoefs,nbasis,numel(ncoefs)/nbasis);

%  define functional data object for optimized coefficients

if ~isempty(basisvals) && isa_basis(basisvals)
    fdobj = fd(ncoefs, basisvals);
else
    if ~isempty(fdobj) && isa_fd(fdobj)
        fdobj = fd(ncoefs, fdobj.basis);
    end
end

end

