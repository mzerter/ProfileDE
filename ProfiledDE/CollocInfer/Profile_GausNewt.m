function [newpars, res0, F0, coefs] = ...
    Profile_GausNewt(times, data, coefs0, pars, lik, proc, ...  
                     in_method, options_in, options_out, active)
% Gauss-Newton optimization routine for squared error outer objective
%  Arguments:  In the follow description of the arguments, N is the
%    number of observations per measured variable, NREP is the number
%    of replications of the observations (usually 1), NBASIS is the
%    number of basis functions used to represent the approximate solution
%    for each variable, and NVAR is the number of variables (including
%    unmeasured variables.)
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
%  COEFS0     ...  A matrix containing coefficients for the expansion.  
%                  It has NVAR columns, and its number of rows is
%                  NBASIS*NVAR.
%  PARS       ...  A vector of parameter values.  These are initial values
%                  that are used to start the outer iterations.  Some
%                  of these may be fixed, and the remainder are optimized.
%                  The indices of the optimized parameters are contained
%                  in argument ACTIVE.
%  LIK        ...  A struct object set up by function LS_SETUP that
%                  defines the loss function for fitting the data
%  PROC       ...  A struct object set up by function LS_SETUP that
%                  defines the loss function for fitting the 
%                  differential equation.
%  IN_METHOD  ...  A string identifing the optimization function to be
%                  used in the outer optimization loop.  It defaults to [].
%  OUT_METHOD ...  A string identifing the optimization function to be
%                  used in the inner optimization loop.  It defaults to [].
%  OPTIONS_IN ...  A struct object containing values controlling the 
%                  inner optimization process.
%  OPTIONS_OUT...  A struct object containing values controlling the 
%                  outer optimization process.  
%                  If not empty, field names are:
%                  EPS or TolFun:  convergence criterion for function 
%                       value.  
%                       Defaults to 1e-6
%                  maxit or MaxIter:  limit on the number of iterations
%                       Defaults to 100
%  DISCRETE   ...  If nonzero, a discrete time model is used instead of
%                  the default continuous time model.
%  ACTIVE     ...  A vector containing the indices of the subset of
%                  parameters to be optimized.  If empty, it defaults
%                  to all indices.

%  Last modified 8 January 2011

%  Set default argument values

if nargin < 10, active=1:length(pars);  end
if nargin <  9, options_out = [];       end
if nargin <  8, options_in  = [];       end
if nargin <  7, in_method   = [];       end

%  set convergence and number of iteration values

if isempty(options_out)
    reltol = 1e-6;
    maxit  = 100;
    maxtry = 15;
    trace  = 1;
else
    if isfield(options_out, 'reltol')
        reltol = options_out.reltol;
    else
        reltol = 1e-6;
    end
    if isfield(options_out, 'maxit')
        maxit = options_out.maxit;
    else
        maxit = 100;
    end
    if isfield(options_out, 'maxtry')
        maxtry = options_out.maxtry;
    else
        maxtry = 15;
    end
    if isfield(options_out, 'trace')
        trace = options_out.trace;
    else
        trace = 1;
    end
end

[res0, gradient, coefs] = ...
    ProfileSSE_AllPar(pars, times, data, coefs0, lik, proc, ...
                         in_method, options_in);
                     
ind    = 1:length(pars);
ind    = ind(active);                                        
pars0  = pars;
pars1  = pars;
F0     = sum(res0.^2);
F1     = F0;
iter   = 0;
fundif = 1;

%  iterate through Gauss-Newton steps

gradnorm1 = 1;
dcdp = [];
ncoefs = coefs;
while gradnorm1 > reltol && ...
      fundif    > reltol && ...
      iter      < maxit
    iter = iter + 1;
    Dpars = (gradient(:,ind)'*gradient(:,ind)) \ ...
            (gradient(:,ind)'*res0);
    %  half step size until function decreases
    ntry = 0;
    while F1           >= F0     && ...
          Dpars'*Dpars >  reltol && ...
          ntry         <  maxtry
        pars1(ind) = pars0(ind) - 0.5*Dpars';
        [res0, gradient, ncoefs, dcdp] = ...
           ProfileSSE_AllPar(pars1, times, data, coefs, lik, proc, ...
                             in_method, options_in, dcdp, pars0);
        F1    = sum(res0.^2);
        Dpars = Dpars/2;
        ntry  = ntry + 1;
    end
    gradnorm1 = abs(mean(gradient(:,ind)' * res0));
    fundif    = (F0 - F1)/abs(F0);
    pars0     = pars1;
    F0        = F1;
    coefs     = ncoefs(:);
    if trace > 1 
        disp([iter,ntry,F0,pars0]) 
    end
end
newpars = pars0;

nbasis = size(lik.bvals,2);
coefs  = reshape(coefs,nbasis,numel(coefs)/nbasis);
n      = length(times);
res0   = reshape(res0, n,     numel(res0 )/n);

end