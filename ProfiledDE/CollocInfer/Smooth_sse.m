%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Smooth_sse_struct = ...
    Smooth_sse(fn, data, times, pars, lambda, coefs, basisvals, ...
    fdobj, more, weights, quadrature, ...
    in_method, options_in, eps, pos, discrete)

%  The multivariate data = argment data are observations of one or more
%  functional variables or processes at time points = argument times.

%  The data are smoothed using differential operators.  These operators are
%  defined by the functions = argument fn.  The smooth can optionally be
%  required to be positive, set by argument pos.

%  The operators typically depend on a number of parameters whose values
%  are = argument pars.

%  The level of smoothing is defined by a scalar or vector of smoothing
%  parameters = argument lambda.

%  The fitting criterion for both the data and the differential equation
%  is weighted error sum of squares, and consequently the function first calls the
%  function sse_setup.

%  The function optimizes the coefficients of the basis expansions for the
%  representations of the variables or processes = what is called the
%  inner optimization loop, implemented by function Inneropt.  The
%  optimization method is defined by argument in_meth.

%  The function returns a functional data object defining the smoothing
%  functions, along with the likelihood of process functions and the
%  residual values.

%  call sse_setup to set up likelihood and process functions for
%  least squares estimation using ODE evaluation functions = fn

if nargin < 16, discrete = 0;      end
if nargin < 15, pos      = 0;      end
if nargin < 14, eps      = 1e-6;   end
if nargin < 13, options_in  = [];  end
if nargin < 12, in_method   = [];  end
if nargin < 11, quadrature  = [];  end
if nargin < 10, weights     = [];  end
if nargin <  9, more        = [];  end
if nargin <  8, fdobj       = [];  end
if nargin <  7, basisvals   = [];  end
if nargin <  6, coefs       = [];  end

profile_obj = sse_setup(pars, coefs, fn, basisvals, lambda, fdobj,...
    more, weights, times, quadrature, ...
    eps, pos, discrete);

lik   = profile_obj.lik;
proc  = profile_obj.proc;
coefs = profile_obj.coefs;
dims = size(data);
if length(dims) > 2
    [l,m,n] = size(data);
    newdata = reshape(data,l*m,n);
    dims(1) = length(coefs)/(m*n);
else
    newdata = data;
    dims = size(coefs);
end
%  carry out an inner optimization to estimate coefficients for
%  processes given parameter values
Ires = inneropt(newdata, times, pars, coefs, lik, proc, ...
                options_in, in_method);
ncoefs = Ires.coefs;
Ires   = Ires.res;
ncoefs = reshape(ncoefs,dims);
% if ~isempty(proc.more.names))
%     colnames(ncoefs) = proc.more.names
% end
%  define functional data object for optimized coefficients
if ~isempty(fdobj)
    fdobj = fd(ncoefs, fdobj.basis);
    Smooth_sse_struct.fd = fdobj;
else
    Smooth_sse_struct.coefs = ncoefs;
end
Smooth_sse_struct.lik  = lik;
Smooth_sse_struct.proc = proc;
Smooth_sse_struct.res  = Ires;

end

