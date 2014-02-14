
function multinorm = ...
    Profile_multinorm(fn, data, times, pars, coefs, basisvals, ...
                      var, fdobj, more, quadrature, ...
                      in_method, out_method, options_in, options_out, ...
                      eps, active, pos, discrete)

if nargin < 18, discrete = 0;     end                                
if nargin < 17, pos = 0;          end 
if nargin < 16, active = [];      end
if nargin < 15, eps = 1e-6;       end
if nargin < 14, options.out = []; end
if nargin < 13, options.in = [];  end
if nargin < 12, out_method  = []; end
if nargin < 11, in_method  = [];  end
if nargin < 10, quadrature = [];  end
if nargin <  9, more = [];        end
if nargin <  8, fdobj = [];       end
if nargin <  7, var = [1, 0.01];  end
if nargin <  6, basisvals = [];   end
if nargin <  5, coefs=[];         end

if isempty(active)
    active = 1:length(pars);
end

profile_obj = multinorm_setup(pars, coefs, fn, basisvals, var, ...
                              fdobj, more, times, quadrature, ...
                              1e-6, pos, discrete);

lik   = profile_obj.lik;
proc  = profile_obj.proc;
coefs = profile_obj.coefs;

coef_opt = inneropt(data, times, pars, coefs, lik, proc, ...
                    options_in, in_method);

apars = pars(active);
%     aparamnames = names(apars);

[pars, ncoefs, res, counter] = ...
    outeropt(data, times, apars, coef_opt, lik, proc,active, ...
               in_method, out_method, options_in, options_out);

apars = res.pars(active);
%     names(apars) = aparamnames

pars(active) = apars;

if ~isempty(fdobj)
    fdobj = fd(ncoefs,fdobj.basis);
    multinorm.pars = pars;
    multinorm.fd   = fdobj;
    multinorm.lik  = lik;
    multinorm.proc = proc;
    multinorm.res  = res;
else
    multinorm.pars  = pars;
    multinorm.coefs = ncoefs;
    multinorm.lik   = lik;
    multinorm.proc  = proc;
    multinorm.res   = res;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function profile = Smooth_multinorm(fn, data, times, pars, coefs, ...
                                    basisvals, var, fdobj, more, ...
                                    quadrature, options_in, in_method, ...
                                    eps, pos, discrete)
                                
if nargin < 13, discrete = 0;     end                                
if nargin < 13, pos = 0;          end                                
if nargin < 12, eps = 1e-6;       end
if nargin < 11, in_method  = [];  end
if nargin < 10, options.in = [];  end
if nargin <  9, quadrature = [];  end
if nargin <  8, more = [];        end
if nargin <  7, fdobj = [];       end
if nargin <  6, var = [1, 0.01];  end
if nargin <  5, basisvals = [];   end
if nargin <  4, coefs=[];         end

profile_obj = multinorm_setup(pars, coefs, fn, basisvals, var, fdobj, ...
                              more, times, quadrature, eps, pos, discrete);

lik   = profile_obj.lik;
proc  = profile_obj.proc;
coefs = profile_obj.coefs;

dims = size(data);
if length(dims) > 2
    data = reshape(data, dims(1)*dims(2),dims(3));
    dims(1) = length(coefs)/(dims(2)*dims(3));
else
    dims = size(coefs);
end

Ires = inneropt(data, times, pars, coefs, lik, proc, ...
                options_in, in_method);

ncoefs = Ires.coefs;
Ires   = Ires.res;
ncoefs = array(ncoefs,dims);

if ~isempty(fdobj)
    fdobj = fd(ncoefs,fdobj.basis);
    profile.fd   = fdobj;
    profile.lik  = lik;
    profile.proc = proc;
    profile.res  = Ires;
else
    profile.coefs = ncoefs;
    profile.lik   = lik;
    profile.proc  = proc;
    profile.res   = Ires;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lik, proc, coefs] = ...
    multinorm_setup(pars, fn, coefs, basisval, var, fdobj, ...
                    more, times, quadrature, eps, pos, discrete)

if nargin < 12, discrete = 0;     end                                
if nargin < 11, pos = 0;          end                                
if nargin < 10, eps = 1e-6;       end
if nargin <  9, quadrature = [];  end
if nargin <  8, times = [];       end
if nargin <  7, more = [];        end
if nargin <  6, fdobj = [];       end
if nargin <  5, var = [1, 0.01];  end
if nargin <  4, basisvals = [];   end
if nargin <  3, coefs=[];         end

if ~isempty(fdobj)                 
    % If an fd object is provided, it overrides
    % the basis and function values
    basisvals = fdobj.basis;
    if ~isempty(fdobj.coefs)
        coefs = fdobj.coefs;
    end
%     colnames = fdobj.fdnames((3))
% else
%     colnames = NULL
end

lik = make_multinorm();

if pos==0
    lik.more = c(make_id(),make_cvar());
else
    lik.more = c(make_exp(),make_cvar());
end

lik.more.f.more = [];
lik.more.v.more.mat = var(1);
lik.more.v.more.sub = zeros(3,1);

if length(size(coefs)) > 2
%     if isempty(colnames))
%         colnames = dimnames(coefs)((3))
%     end
    [nbasis,nrep,nvar] = size(coefs);
    coefs = reshape(coefs,basis*nrep,nvar);
else
    nrep = 1;
%     colnames = colnames(coefs)
end

if pos==0
    if discrete==0
        proc = make_Cproc();
    else
        proc = make_Dproc();
    end
else
    if discrete==0
        proc = make_expCproc();
    else
        proc = make_expDproc();
    end
end

proc.more = make_multinorm();

if isstruct(fn)
    proc.more.more = c(fn,make_cvar());
    proc.more.more.f.more = NULL;
    proc.more.more.v.more.mat = var(2)*diag(ones(2,1));
    proc.more.more.v.more.sub = zeros(3,1);
    proc.more.more.more = more;
elseif strcmp(class(fn),'fd')
    proc.more.more = [make_findif_ode(),make_cvar()];
    proc.more.more.f.more.eps = eps;
    proc.more.more.f.more.fn = fn;
    proc.more.more.more = more;
else
    error('fn must be either a list of functions or a function');
end

% proc.more.names = colnames
% proc.more.parnames = names(pars)

if is_basis(basisvals)
    if isempty(times)
        error('observation times not specified for basisvals');
    end
    
    lik.bvals = sparse(diag(ones(nrep,1))); 
    
    if isempty(quadrature) || isempty(quadrature.qpts)
        knots = [basisvals.rangeval(1), basisvals.params, ...
                 basisvals.rangeval(2)];
        qpts = knots(1:length(knots)-1) + diff(knots);
    else
        qpts = quadrature.qpts;
    end
    
    if discrete==0
        proc.bvals.bvals  = lik.bvals;
        proc.bvals.dbvals = lik.bvals;
        proc.more.qpts = qpts;
    else
        len = length(times);
        proc.bvals = lik.bvals;
        proc.more.qpts = times(1:(len-1));
    end
else
    if discrete && (is_reshape(basisvals) || isempty(basisvals))
        if isempty(basisvals) 
            basisvals = sparse(diag(ones(size(coef,1),1))); 
        end
        lik.bvals  = basisvals;
        proc.bvals = basisvals;
        proc.more.qpts = times(1:(length(times)-1));
    else
        lik.bvals = basisvals.bvals.obs;
        proc.bvals.bvals  = basisvals.bvals;
        proc.bvals.dbvals = basisvals.dbvals;
        proc.more.qpts = basisvals.qpts;
    end
end

end
