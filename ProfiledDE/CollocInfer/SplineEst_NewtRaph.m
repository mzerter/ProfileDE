%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function coefs_opt = ...
    SplineEst_NewtRaph(times, data, coefs0, pars, lik, proc, options_in)
% Newton-Raphson routine for  inner optimization
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
%  OPTIONS_IN ...  A struct object containing values controlling the 
%                  inner optimization process.

%  Last modified 8 January 2011

if nargin <  7, options_in = [];  end

%  set values of iteration control parameters

if isempty(options_in)
    reltol = 1e-6;
    maxit  = 100;
    maxtry = 15;
    trace  = 1;
else
    if isfield(options_in, 'reltol')
        reltol = options_in.reltol;
    else
        reltol = 1e-6;
    end
    if isfield(options_in, 'maxit')
        maxit = options_in.maxit;
    else
        maxit = 100;
    end
    if isfield(options_in, 'maxtry')
        maxtry = options_in.maxtry;
    else
        maxtry = 15;
    end
    if isfield(options_in, 'trace')
        trace = options_in.trace;
    else
        trace = 1;
    end
end

%  Initial evaluation of the objective function and its derivatives

[f0, dfdc, d2fdc2] = SplineCoefs(coefs0, times, data, pars, lik, proc);
gradnorm1 = 1;
fundif    = 1;
iter      = 0;
f1        = f0;
coefs0    = coefs0(:);

while gradnorm1 > reltol && fundif > 0 && iter < maxit
    iter = iter + 1;
    DC = -d2fdc2\dfdc;
    ntry = 0;
    coefs1 = coefs0;    
    while f1 >= f0 && DC'*DC > reltol && ntry < maxtry
        coefs1 = coefs0 + DC;
        f1 = SplineCoefs(coefs1, times, data, pars, lik, proc);
        DC = DC/2;
        ntry = ntry + 1;
    end
    coefs0 = coefs1;
    [f1, dfdc, d2fdc2] = SplineCoefs(coefs0, times, data, pars, lik, proc);
    gradnorm1 = mean(abs(DC));
    fundif = (f0-f1)/abs(f0);   
    f0 = f1;
    if trace > 1 
        disp([iter, ntry, f0, gradnorm1, fundif]) 
    end    
end

coefs_opt = coefs0;

if trace > 0
    disp([f0, gradnorm1, iter])
end

end
