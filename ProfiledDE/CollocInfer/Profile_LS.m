function Profile_LS_struct = ...
    Profile_LS(fn, times, data, coefs, allpars, basisobj, lambda, ...
               fdobj, more, obsweights, quadrature, active, ...
               in_method, options_in, out_method, options_out, ...
               diffeps, poslik, posproc, discrete)
%  The multivariate data = argment data are observations of one or more
%  functional variables or processes at time points = argument times.
%
%  The data are smoothed using differential operators.  These operators are
%  defined by the functions = argument fn.  The smooth can optionally be
%  required to be positive, set by argument pos.
%
%  The operators typically depend on a number of parameters whose values
%  are = argument pars.  The function optimizes these parameter values
%  The level of smoothing is defined by a scalar or vector of smoothing
%  parameters = argument lambda.
%
%  The fitting criterion for both the data and the differential equation
%  is weighted error sum of squares, and consequently the function first 
%  calls the function sse_setup.
%
%  The function optimizes the coefficients of the basis expansions for the
%  representations of the variables or processes = what is called the
%  inner optimization loop, implemented by function Inneropt.  It optimizes 
%  the parameter values = the outer optimization loop.
%
%  The optimization method for the inner loop is defined by in_meth.
%  The optimization method for the outer loop is defined by out_meth.
%  The function returns the optimized parameter values, a functional data 
%  object defining the smoothing functions, the likelihood and process 
%  functions, and the residual values.
%
%  Arguments:  In the following description of the arguments, N is the
%    number of observations per measured variable, NREP is the number
%    of replications of the observations (usually 1), NBASIS is the
%    number of basis functions used to represent the approximate solution
%    for each variable, and NVAR is the number of variables (including
%    unmeasured variables.)
%    The first six arguments are required.  The remaining arguments each
%    have default values, can be omitted.  All arguments after the first
%    omitted argument must be also omitted.
%
%  Arguments:
%    The first seven arguments are required.  The remaining arguments each
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
%                  The number of rows is equal to the number of 
%                  observation points per variable times the number of
%                  replications.
%                  Unobserved variables have all NaN's in the respective
%                  column.
%  COEFS      ...  A matrix containing coefficients for the expansion.  
%                  It has NVAR columns, and its number of rows is
%                  NBASIS*NVAR.
%  ALLPARS    ...  A vector of parameter values.  These are initial values
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
%  OPTIONS_OUT...  A struct object containing values controlling the 
%                  outer optimization process.  
%                  If not empty, field names are:
%                      EPS or TolFun:  convergence criterion for function 
%                          value.  
%                          Defaults to 1e-6
%                      maxit or MaxIter:  limit on the number of iterations
%                          Defaults to 100
%                  Defaults to [].
%  DIFFEPS    ...  A small constant value used to compute difference 
%                  approximations to derivatives if this is requested.
%  POSLIK     ...  If nonzero, the functions approximating the data are
%                  constrained to be positive.
%  POSPROC    ...  If nonzero, the functions used in the penalty terms
%                  are constrained to be positive.
%  DISCRETE   ...  If nonzero, a discrete time model is used instead of
%                  the default continuous time model.

%  Last modified 2 November 2013

if nargin < 7
    error('Less than seven arguments supplied.');
end

if nargin < 20 || isempty(discrete), discrete = 0;      end
if nargin < 19 || isempty(posproc),  posproc  = 0;      end
if nargin < 18 || isempty(poslik),   poslik   = 0;      end
if nargin < 17 || isempty(diffeps),  diffeps  = 1e-6;   end
if nargin < 16,   options_out = [];  end
if nargin < 15,   out_method  = [];  end
if nargin < 14,   options_in  = [];  end
if nargin < 13,   in_method   = [];  end
if nargin < 12,   active      = [];  end
if nargin < 11,   quadrature  = [];  end
if nargin < 10,   obsweights  = [];  end
if nargin <  9,   more        = [];  end
if nargin <  8,   fdobj       = [];  end

% ---------------------------------
%  check the first seven arguments
% ---------------------------------

%  check fn ... it must be a struct variable and the fields must be
%  function handles

if ~isstruct(fn)
    error('Argument fn is not a struct object.');
end

fieldnames = cell(1,5);
fieldnames{1} = 'fn';
fieldnames{2} = 'dfdx';
fieldnames{3} = 'dfdp';
fieldnames{4} = 'd2fdx2';
fieldnames{5} = 'd2fdxdp';

tp = isfield(fn, fieldnames);
if ~all(tp)
    error('Required field names in argument fn are not present.');
end

if ~strcmp(class(fn.fn), 'function_handle')
    error('Argument fn.fn is not a function handle.');
end
if ~strcmp(class(fn.dfdx), 'function_handle')
    error('Argument fn.dfdx is not a function handle.');
end
if ~strcmp(class(fn.dfdp), 'function_handle')
    error('Argument fn.dfdp is not a function handle.');
end
if ~strcmp(class(fn.d2fdx2), 'function_handle')
    error('Argument fn.d2fdx2 is not a function handle.');
end
if ~strcmp(class(fn.d2fdxdp), 'function_handle')
    error('Argument fn.d2fdxdp is not a function handle.');
end

%  check times

if ~isnumeric(times)
    error('Argument times is not numeric');
end
[nobs,ntmp] = size(times);
if nobs > 1 && ntmp > 1
    error('Argument times is not a vector');
end
times = times(:);
nobs = length(times);

%  check data for consistency with times

if ~isnumeric(data)
    error('Argument data is not numeric');
end
[ntmp, nvar] = size(data);
if ntmp ~= nobs
    error(['Length of argument times not equal to ', ...
           'number rows in argument data.']);
end

% check coefs

if ~isnumeric(coefs)
    error('Argument coefs is not numeric');
end
[nbasis, ntmp] = size(coefs);
if ntmp ~= nvar
    error(['Number of columns in argument coefs not equal to ', ...
           'number columns in argument data.']);
end

% check allpars

if ~isnumeric(allpars)
    error('Argument allpars is not numeric');
end
npar = length(allpars);

%  check basisobj

if ~(strcmp(class(basisobj), 'basis') || ...
     strcmp(class(basisobj), 'struct'))
    error(['Argument basisobj is neither a basis object nor ', ...
           'a struct object.']);
end

if strcmp(class(basisobj), 'basis')
    ntmp = getnbasis(basisobj);
    if nbasis ~= ntmp
        error(['Number of basis functions is not consistent with ', ...
               'number of columns of argument coef.']);
    end
else
    ntmp = size(basisobj,2);
    if nbasis ~= ntmp
        error(['Number of basis functions is not consistent with ', ...
               'number of columns of argument coef.']);
    end
end

%  check lambda

nvar = size(data,2);
if ~isnumeric(lambda)
    error('Argument lambda is not numeric');
end
ntmp = length(lambda);
if ~(ntmp == nvar || ntmp == 1)
    error(['Length of argument lambda is neither 1 nor ', ...
           'the number of state variables.']);
end

%  ----------------------------------------------------------
%  set convergence tolerance and maximum number of iterations
%  ----------------------------------------------------------

if isempty(options_out) 
    eps   = 1e-6;
    maxit = 500;
else
    if      isfield(options_out, 'eps')
        eps = options_out.eps;
    elseif  isfield(options_out, 'TolFun')
        eps = options_out.TolFun;
    else
        eps = 1e-6;
    end
    if      isfield(options_out, 'maxit')
        maxit = options_out.maxit;
    elseif  isfield(options_out, 'MaxIter')
        maxit = options_out.MaxIter;
    else
        maxit = 500;
    end
end

%  set default value for ACTIVE or check input values

if isempty(active)
    active = 1:npar;
else
    if any(active < 1) || any(active > npar) || any(diff(active) <= 0)
        disp('ACTIVE')
        disp(active)
        error('Argument ACTIVE is not admissible.');
    end
end

%  ----------------------------------------------------------
%  Set up likelihood and process functions for
%  least squares estimation using ODE evaluation functions = fn
%  ----------------------------------------------------------

[lik, proc, coefs] = ...
            LS_setup(fn, times, coefs, basisobj, lambda, ...
                     fdobj, more, obsweights, quadrature, ...
                     diffeps, poslik, posproc, discrete);
                           
%  ----------------------------------------------------------
%  Select optimization method and optimize with respect to parameters
%  The default method is Gauss-Newton optimizing using Matlab function
%  lsnonlin, and the alternative is using home-grown function
%  Profile_GausNewt.  
%  The proceed with optimization.
%  ----------------------------------------------------------

pars0 = allpars(active);

% disp('optimization')

if strcmp(out_method,'ProfileGN')
    %%%  Gauss-Newton method  using function Profile_GausNewt  %%%
    [newpars, res, F, ncoefs] = ...
        Profile_GausNewt(times, data, coefs, pars0, lik, proc, ...
                         in_method, options_in, options_out, active);
                     
else
%     disp('lsqnonlin')
    if isempty(options_out)
        options_out = optimset('Jacobian', 'on', 'MaxIter', maxit, ...
                               'Display', 'iter', ...
                               'TolFun', eps);
    end
%     disp('ProfileSSE')
%     res0 = ProfileSSE(pars0, times, data, coefs, allpars, ... 
%                       lik, proc, active, in_method, options_in);
%     disp(res0)         
%     disp(sum(res0.^2))
              
    newpars = lsqnonlin(@ProfileSSE, pars0, [], [], options_out, ...
                        times, data, coefs, allpars, lik, proc, ...
                        active, in_method, options_in);
    allpars(active) = newpars;
    %  wrap up by computing final residuals, Jacobian matrix, and
    %  coefficients
    [res, Dres, ncoefs] = ...
             ProfileSSE_AllPar(allpars, times, data, coefs, lik, proc, ...
                               in_method, options_in);
end

allpars(active) = newpars;

%  ----------------------------------------------------------
%  Set up output struct object Profile_LS_struct
%  ----------------------------------------------------------

if ~isempty(fdobj)
    fdobj = fd(ncoefs,getbasis(fdobj));
    Profile_LS_struct.fd = fdobj;
else
    Profile_LS_struct.coefs = ncoefs;
end

Profile_LS_struct.pars  = allpars;
Profile_LS_struct.lik   = lik;
Profile_LS_struct.proc  = proc;
Profile_LS_struct.res   = res;
Profile_LS_struct.coefs = ncoefs;

end
