function coefs_opt = inneropt(times, data, pars, lik, proc, ...
                              in_method, options_in)
%  INNEROPT optimizes the coefficients of the basis expansions for the
%  representations of the variables or processes in what is called the
%  inner optimization loop.
%
%  The fitting criterion for both the data and the differential equation
%  is user-defined = arguments lik and proc, respectively.
%
%  The options for the Matlab unconstrained optimization function FMINUNC
%  are defined in the argument OPTIONS_IN, which is in turned defined
%  by a call to Matlab function OPTIMSET.
%
%  The optimization method is defined by argument IN_METHOD.
%
%  The initial values of matrix COEFS0 is defined outside of the function
%  and passed in as the global variable INNEROPT_COEFS0.
%
%  The function returns the optimized coefficients.
%
%  Arguments:
%    The first five arguments are required.  The remaining arguments each
%    have default values, can be omitted.  All arguments after the first
%    omitted argument must be also omitted.
%  DATA       ...  Data matrix.  The number of columns is equal to the
%                  number of variables (including unobserved variables).
%                  The number of columns is equal to the number of 
%                  observation points per variable times the number of
%                  replications.
%  TIMES      ...  The times of observation of those variables that are
%                  measured.  These times are assumed to be common to all
%                  measured variables.
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
%  IN_METHOD  ...  A string identifing the optimization function to be
%                  used in the outer optimization loop.  It defaults to [].

%  Last modified 13 January 2011

global INNEROPT_COEFS0

if isempty(INNEROPT_COEFS0)
    error('Initial coefficients matrix is empty.');
end
%%

if nargin < 7, in_method  = [];  end
if nargin < 6, options_in = [];           end

if strcmp(in_method,'SplineEst')
    
    %  optimization by Newton-Raphson using COLLOCINFER function 
    %  SplineEst_NewtRaph
    
    coefs_opt  = ...
        SplineEst_NewtRaph(times, data, INNEROPT_COEFS0, pars, lik, proc,...
                           options_in);
    
else
    
    %  Optimization using Matlab function FMINUNC to optimize the
    %  function value returned by COLLOCINFER function SPLINECOEFSERR,
    %  which may, optionally, also return the gradient and hessian as
    %  the second and third returned objects, respectively.

    if isempty(options_in)
        options_in = optimset('LargeScale', 'on', 'GradObj', 'on', ...
                              'Hessian', 'on', 'Diagnostics', 'off', ...
                              'Display', 'off');
    end
    coefs_opt  = fminunc(@SplineCoefs, INNEROPT_COEFS0, options_in, ...
                         times, data, pars, lik, proc);
end

%  return COEFS_OPT as a matrix

npars = length(coefs_opt(:));
n     = size(lik.bvals,2);
coefs_opt = reshape(coefs_opt,n,npars/n);

end

