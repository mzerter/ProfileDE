function Profile_sse_struct = ...
    Profile_sse(fn, data, times, pars, lambda, coefs, basisvals, ...
                fdobj, more, weights, quadrature, active, ...
                in_method, out_method, options_in, options_out, ...
                eps, pos, discrete)

%  The multivariate data = argment data are observations of one or more
%  functional variables or processes at time points = argument times.

%  The data are smoothed using differential operators.  These operators are
%  defined by the functions = argument fn.  The smooth can optionally be
%  required to be positive, set by argument pos.

%  The operators typically depend on a number of parameters whose values
%  are = argument pars.  The function optimizes these parameter values

%  The level of smoothing is defined by a scalar or vector of smoothing
%  parameters = argument lambda.

%  The fitting criterion for both the data and the differential equation
%  is weighted error sum of squares, and consequently the function first 
%  calls the function sse_setup.

%  The function optimizes the coefficients of the basis expansions for the
%  representations of the variables or processes = what is called the
%  inner optimization loop, implemented by function Inneropt.  It optimizes 
%  the parameter values = the outer optimization loop.

%  The optimization method for the inner loop is defined by in_meth.
%  The optimization method for the outer loop is defined by out_meth.

%  The function returns the optimized parameter values, a functional data 
%  object defining the smoothing functions, the likelihood and process 
%  functions, and the residual values.

if nargin < 19, discrete = 0;      end
if nargin < 18, pos      = 0;      end
if nargin < 17, eps      = 1e-6;   end
if nargin < 16, options_out = [];  end
if nargin < 15, options_in  = [];  end
if nargin < 14, out_method  = [];  end
if nargin < 13, in_method   = [];  end
if nargin < 12, active      = [];  end
if nargin < 11, quadrature  = [];  end
if nargin < 10, weights     = [];  end
if nargin <  9, more        = [];  end
if nargin <  8, fdobj       = [];  end
if nargin <  7, basisvals   = [];  end
if nargin <  6, coefs       = [];  end

if isempty(active)
    active = 1:length(pars);
end

%  Call sse_setup to set up likelihood and process functions for
%  least squares estimation using ODE evaluation functions = fn

[lik, proc, coefs] = ...
            sse_setup(fn, times, pars, lambda, coefs, basisvals, ...
                      fdbj, more, weights, quadrature, eps, pos, discrete);
                           
%  format data

[newdata, dims] = dataformat(data, coefs);

%  clear any old temporary files
% if (file.exists('curcoefs.tmp')) file.remove('curcoefs.tmp')
% if (file.exists('optcoefs.tmp')) file.remove('optcoefs.tmp')
% if (file.exists('counter.tmp'))  file.remove('counter.tmp')

%  initial inner optimization to get coefficients

ncoefs = inneropt(newdata, times, pars, coefs, lik, proc, ...
                  options_in, in_method);
              
%  write coefficients to temporary files
% write.table(ncoefs,file='optcoefs.tmp',col.names=FALSE,row.names=FALSE)
% write.table(ncoefs,file='curcoefs.tmp',col.names=FALSE,row.names=FALSE)

%  Select optimization method and optimize with respect to parameters

if strcmp(out_method,'ProfileGN')
    %%%  Gauss-Newton method  using function Profile_GausNewt  %%%
    [newpars, res] = ...
        Profile_GausNewt(pars, times, data, ncoefs, active, lik, proc, ...
                         in_method, options_in, options_out);
    apars  = newpars(active);  % optimal parameter values
    ncoefs = res.in.res.coefs;
    g      = res.in.res.df;
    resid  = res.in.res.f;
else
    error('no method defined');
end

% if file.exists('curcoefs.tmp'))
%     ncoefs = read.table(file='curcoefs.tmp'))
% else
ncoefs = coefs;
% end

pars(active) = apars;
ncoefs = reshape(ncoefs,dims);

%  remove temporary files
% if file.exists('curcoefs.tmp')
%     file.remove('curcoefs.tmp')
% end
% if file.exists('optcoefs.tmp')
%     file.remove('optcoefs.tmp')
% end
% if file.exists('counter.tmp')
%     file.remove('counter.tmp')
% end

if ~isempty(fdobj)
    fdobj = fd(ncoefs,fdobj.basis);
    Profile_sse_struct.fd = fdobj;
else
    Profile_sse_struct.coefs = ncoefs;
end
Profile_sse_struct.pars = pars;
Profile_sse_struct.lik  = lik;
Profile_sse_struct.proc = proc;
Profile_sse_struct.res  = res;

end
