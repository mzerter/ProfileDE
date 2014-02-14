function [lik, proc, coefs] = ...
    LS_setup(fn, times, pars, lambda, coefs, basisvals, fdobj, ...
              more, weights, quadrature, eps, poslik, posproc, discrete)

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
    basisvals = getbasis(fdobj);
    if ~isempty(getcoef(fdobj));
        coefs = getcoef(fdobj);
    end
end

%  define functions for evaluating the values of the
%  error sum of squares for the data and its derivatives
%  names(lik)
%  (1) 'fn'      'dfdx'    'dfdy'    'dfdp'    'd2fdx2'  'd2fdxdy'
% (7) 'd2fdy2'  'd2fdxdp' 'd2fdydp' 'more'    'bvals'

lik = make_SSElik();

%  add the more slot to lik that defines the transformation of the process
%  to fit the data.  pos == 0 means no transformation, otherwise an
%  exponential transformation to provide a positive fit.

if pos==0
    lik.more = make_id;
else
    lik.more = make_exp();
end

%  determine number of repeated measurements nrep

[coefs, nrep] = coefsformat(coefs);

%  Define functions for evaluating the values of the
%  error sum of squares for the ODE and its derivatives

proc = make_SSEproc;

if isstruct(fn)
    procmore = fn;
    procmore.more = more;
elseif strcmp(class(fn), 'fd')
    procmore = make_findif_ode();
    procmore.more.fn   = fn;
    procmore.more.more = more;
    procmore.more.eps  = eps;
else
    error('fn must be either a struct of function handles or a function');
end
if pos==0
    proc.more = procmore;
else
    proc.more      = make_logtrans();
    proc.more.more = procmore;
end

%  set up basis for representing solution

if strcmp(class(basisvals),'basis')
    %  basis object supplied as an argument
    basisobj = basisvals;
    if isempty(times)
        error(['if basisvals is a basis object, ' ...
               'you must specify the observation times']);
    end
    %  basis value array for sampling points
    if nrep > 1
        lik.bvals = kron(sparse(diag(ones(nrep,1))), ...
                         eval_basis(times,basisobj));
    else
        lik.bvals = eval_basis(times,basisobj);
    end
    %  define quadrature points
    if isempty(quadrature) || isempty(quadrature.qpts)
        rangeval = getbasisrange(basisobj);
        params   = getbasispar(basisobj);
        knots = [rangeval(1), params, rangeval(2)];
        qpts = knots(1:length(knots)-1) + diff(knots);
    else
        qpts = quadrature.qpts;
    end
    %  basis value array for quadrature points
    if discrete==0
        if nrep > 1
            proc.bvals.bvals  = kron(sparse(diag(ones(nrep,1))), ...
                                     eval_basis(times,basisvals));
            proc.bvals.dbvals = kron(sparse(diag(ones(nrep,1))), ...
                                     eval_basis(qpts,basisvals));
        else
            proc.bvals.bvals  = eval_basis(times,basisvals);
            proc.bvals.dbvals = eval_basis(qpts,basisvals);
        end
        proc.more.weights = ...
            repmat(1/length(qpts),length(qpts)*nrep,size(coefs,2));
        proc.more.qpts = qpts;
    else
        n = length(times);
        basismat          = eval_basis(times,basisvals);
        proc.bvals.bvals  = basismat(1:(n-1),:);
        proc.bvals.dbvals = basismat(2:n,:);
        proc.more.weights = repmat(1/(n-1),(n-1)*nrep,size(coefs,1));
        proc.more.qpts    = times(1:(n-1));
    end
else
    %  basis values supplied  as argument
    if discrete && (isnumeric(basisvals) || isempty(basisvals))
        if isempty(basisvals)
            basisvals = Diagonal(size(coefs,1));
        end
        lik.bvals = basisvals;
        n = size(basisvals,1);
        proc.bvals.bvals  = basisvals(1:n-1,:);
        proc.bvals.dbvals = basisvals(2:n,:);
        proc.more.weights = (1/(n-1)).*ones(n-1,nrep);
        proc.more.qpts    = times(1:(length(times)-1));
    else
        lik.bvals         = basisvals.bvals.obs;
        proc.bvals.bvals  = basisvals.bvals;
        proc.bvals.dbvals = basisvals.dbvals;
        proc.more.qpts    = basisvals.qpts;
        if ~isempty(basisvals.qwts)
            proc.more.weights = ...
                reshape(basisvals.qwts, ...
                        length(basisvals.qwts)* nrep,nrep);
        else
            proc.more.weights = ...
                repmat(1/length(proc.more.qpts), ...
                       length(proc.more.qpts)*nrep,nrep);
        end
    end
end

%  define weights

if  isempty(weights)
    lik.more.weights = ones(nrep,1);
else
    lik.more.weights = repmat(weights,size(lik.bvals,2)*nrep,nrep);
end

%  define vector of smoothing parameters lambda

if length(lambda)==1
    lambda = lambda.*ones(nrep,1);
end

proc.more.weights = proc.more.weights * diag(lambda);

end
