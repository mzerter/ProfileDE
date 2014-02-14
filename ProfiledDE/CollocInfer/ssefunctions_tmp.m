
%  This file defines three functions:

%  Smooth_sse:   Smoothing data using the inner coefficient optimization
%  Profile_sse:  Optimizing parameter values using outer optimization
%  sse_setup:    Set up error sum of squares functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Smooth_sse_struct = ...
    Smooth_sse(fn, data, times, pars, coefs, lambda, basisvals, fdobj, ...
               more, weights, quadrature, in_method, options_in, ...
               eps, pos, discrete)

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

profile_obj = sse_setup(pars, coefs, fn, basisvals, lambda, fdobj, ...
                        more, weights, times, quadrature, ...
                        eps, pos, discrete);
                    
lik   = profile_obj.lik;
proc  = profile_obj.proc;
coefs = profile_obj.coefs;
%  format data
[newdata, dims] = dataformat(data);
%  carry out an inner optimization to estimate coefficients for
%  processes given parameter values
ncoefs = inneropt(newdata, times, pars, coefs, lik, proc, ...
                  options_in, in_method);
% if ~isempty(proc.more.names)
%     colnames(ncoefs) = proc.more.names
% end
%  define functional data object for optimized coefficients
if ~isempty(fdobj)
    fdobj = fd(ncoefs,fdobj.basis);
    Smooth_sse_struct.fd = fdobj;
else
    Smooth_sse_struct.coefs = ncoefs;
end
Smooth_sse_struct.lik  = lik;
Smooth_sse_struct.proc = proc;
Smooth_sse_struct.res  = Ires;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Profile_sse_struct = ...
    Profile_sse(fn, data, times, pars, lambda, coefs, basisvals, fdobj, ...
                more, weights, quadrature, active, ...
                in_method, out_method, options_in, options_out, ...
                eps, pos, discrete)

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
%  is weighted error sum of squares, and consequently the function first calls the
%  function sse_setup.

%  The function optimizes the coefficients of the basis expansions for the
%  representations of the variables or processes = what is called the
%  inner optimization loop, implemented by function Inneropt.  It optimizes the
%  parameter values = the outer optimization loop.

%  The optimization method for the inner loop is defined by argument in_meth.
%  The optimization method for the outer loop is defined by argument out_meth.

%  The function returns the optimized parameter values, a functional data object defining
%  the smoothing functions, the likelihood and process functions, and the
%  residual values.

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

%  call sse_setup to set up likelihood and process functions for
%  least squares estimation using ODE evaluation functions = fn

profile_obj = sse_setup(pars, coefs, fn, basisvals, lambda, fdobj, ...
                        more, weights, times, quadrature, ...
                        eps, pos, discrete);
lik   = profile_obj.lik;
proc  = profile_obj.proc;
coefs = profile_obj.coefs;
%  format data
[newdata, dims] = dataformat(data);
%  clear any old temporary files
%     if (file.exists('curcoefs.tmp')) file.remove('curcoefs.tmp')
%     if (file.exists('optcoefs.tmp')) file.remove('optcoefs.tmp')
%     if (file.exists('counter.tmp'))  file.remove('counter.tmp')
%  initial inner optimization to get coefficients
ncoefs = inneropt(newdata, times, pars, coefs, lik, proc, ...
                  options_in, in_method);
%  write coefficients to temporary files
%   write.table(ncoefs,file='optcoefs.tmp',col.names=FALSE,row.names=FALSE)
%   write.table(ncoefs,file='curcoefs.tmp',col.names=FALSE,row.names=FALSE)
apars = pars(active);
% aparamnames = names(apars);
%  Select optimization method and optimize with respect to parameters
if strcmp(out_method,'ProfileGN')
     %%%  Gauss-Newton method  using function Profile_GausNewt  %%%  
    res = Profile_GausNewt(pars, times, newdata, ncoefs, lik, proc, ...
                           active, in_method, options_in, options_out);    
    apars  = res.pars(active);  % optimal parameter values
    ncoefs = res.in.res.coefs;
    g      = res.in.res.df;
    resid  = res.in.res.f;
else
    error('no method specified.');
end
% names(apars) = aparamnames
pars(active) = apars;
ncoefs = reshape(ncoefs,dims);
% if ~isempty(proc.more.names)
%     colnames(ncoefs) = proc.more.names
% end
%  remove temporary files
%     if file.exists('curcoefs.tmp')) file.remove('curcoefs.tmp')
%     if file.exists('optcoefs.tmp')) file.remove('optcoefs.tmp')
%     if file.exists('counter.tmp'))  file.remove('counter.tmp')
%  return results

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lik, proc,coefs] = ...
    sse_setup(pars, times, fn, lambda, coefs, basisvals, fdobj, ...
              more, weights, quadrature, eps, pos, discrete)

if nargin < 13, discrete = 0;      end
if nargin < 12, pos      = 0;      end
if nargin < 11, eps      = 1e-6;   end
if nargin < 10, quadrature  = [];  end
if nargin <  9, weights     = [];  end
if nargin <  8, more        = [];  end
if nargin <  7, fdobj       = [];  end
if nargin <  6, basisvals   = [];  end
if nargin <  5, coefs       = [];  end
% If an fd object is provided, it overrides the basis and function values
if ~isempty(fdobj)
    basisvals = fdobj.basis;
    if ~isempty(fdobj.coefs)
        coefs = fdobj.coefs;
    end
end
% colnames = NULL
%  See file SSE_lik_R for function make_SSElik
%  this file defines functions for evaluating the values of the
%  error sum of squares for the data and its derivatives
%  names(lik)
%  (1) 'fn'      'dfdx'    'dfdy'    'dfdp'    'd2fdx2'  'd2fdxdy'
%  (7) 'd2fdy2'  'd2fdxdp' 'd2fdydp' 'more'    'bvals'
lik = make_SSElik();
%  add the more slot to lik that defines the transformation of the process
%  to fit the data.  pos == 0 means no transformation, otherwise an
%  exponential transformation to provide a positive fit.
if pos==0
    lik.more = make_id();
else
    lik.more  = make_exp();
end
%  determine number of repeated measurements nrep
if ndims(coefs) > 2
%     if isempty(colnames)
%         colnames = dimnames(coefs)((3))
%     end
    [nbasis, nrep, nvar] = size(coefs);
    coefs = reshape(coefs,nbasis*nrep,nvar);
else
    nrep = 1;
%     colnames = colnames(coefs)
end
%  see file SSE.proc.R for function make_SSEproc
%  this file defines functions for evaluating the values of the
%  error sum of squares for the ODE and its derivatives
proc = make_SSEproc();
if isstruct(fn)
    procmore = fn;
    procmore.more = more;
elseif strcmp(class(fn), 'fd')
    procmore = make_findif_ode();
    procmore.more.fn = fn;
    procmore.more.more = more;
    procmore.more.eps = eps;
else
    error('fn must be either a list of functions or a function');
end
if pos==0
    proc.more = procmore;
else
    proc.more      = make_logtrans();
    proc.more.more = procmore;
end
% proc.more.names    = colnames
% proc.more.parnames = names(pars)
n = length(times);
if isa_basis
    lik.bvals = kron(diag(ones(nep,1)), eval_basis(times,basisvals));
    if isempty(quadrature) || isempty(quadrature.qpts)
        knots = [basisvals.rangeval(1), basisvals.params, ...
                 basisvals.rangeval(2)];
        qpts = knots(1:length(knots)-1) + diff(knots);
    else
        qpts = quadrature.qpts;
    end
    if discrete==0
        proc.bvals.bvals  = kron(diag(ones(nep,1)), ...
                                 eval_basis(qpts,basisvals));
        proc.bvals.dbvals = kron(diag(ones(nep,1)), ...
                                 eval_basis(qpts,basisvals,1));
        proc.more.weights = (1/length(qpts)).*ones(length(qpts)*nrep,nrep);
        proc.more.qpts = qpts;
    else
        basis = eval_basis(times,basisvals);
        proc.bvals.bvals  = basis(1:(n-1),:);
        proc.bvals.dbvals = basis(2:n,:);
        proc.more.weights = (1/(n-1)).*ones((n-1)*nrep,nrep);
        proc.more.qpts = times(1:(n-1));
    end
else                                   % quadrature is ignored if basisvals
    % is not a basis object
    if discrete && (isnumeric(basisvals) || isempty(basisvals))
        if isempty(basisvals) 
            basisvals = Diagonal(size(coefs,1)); 
        end
        lik.bvals = basisvals;
        proc.bvals.bvals  = basisvals(1:n-1,:);
        proc.bvals.dbvals = basisvals(2:n,:);
        proc.more.weights = (1/(size(basisvals,1)-1)).*ones(n-1,nrep);
        proc.more.qpts = times(1:n-1);
    else
        lik.bvals = basisvals.bvals.obs;
        proc.bvals.bvals  = basisvals.bvals;
        proc.bvals.dbvals =basisvals.dbvals;
        proc.more.qpts = basisvals.qpts;
        if ~isempty(basisvals.qwts)
            proc.more.weights = ...
                reshape(basisvals.qwts,length(basisvals.qwts)*nrep,nrep);
        else
            proc.more.weights = ...
                (1/length(proc.more.qpts)).* ...
                ones(length(proc.more.qpts)*nrep,nrep);
        end
    end
    %  define weights
    if  isempty(weights)
        lik.more.weights = ones(nrep,1);
    else
        lik.more.weights = reshape(weights,size(lik.bvals,1)*nrep,nrep);
    end
    %  define vector of smoothing parameters lambda
    if length(lambda)==1
        lambda = lambda.*ones(nrep,1);
    end    
    proc.more.weights = proc.more.weights * diag(lambda);
    
end

end
