function [lik, proc, coefs] = ...
    LS_setup(fn, times, coefs, basisvals, lambda, fdobj, ...
             more, obsweights, quadrature, ...
             diffeps, poslik, posproc, discrete)
%  LS_SETUP sets up the LIK, PROC and COEF objects for a CollocInfer
%  analysis where squared error measures of fit to the data and equation
%  are used.   In this common situation, the effort required of the user
%  is considerably reduced by using this function.
%
%  Arguments:  In the follow description of the arguments, N is the
%    number of observations per measured variable, NREP is the number
%    of replications of the observations (usually 1), NBASIS is the
%    number of basis functions used to represent the approximate solution
%    for each variable, and NVAR is the number of variables (including
%    unmeasured variables.)
%
%  FN         ...  A struct object containing handles to functions that
%                  evaluate the value of the right side of the 
%                  differential equation as well as the values of a 
%                  number of its partial derivatives.  These functions
%                  must be provided by the user.  However, difference
%                  approximations of derivatives may also be requested.
%  TIMES      ...  The times of observation of those variables that are
%                  measured.  These times are assumed to be common to all
%                  measured variables.
%  LAMBDA     ...  A single constant, or a vector of length NVAR containing
%                  smoothing parameter values for each variable controlling
%                  how closely a variable comes to solving the equation.
%  COEFS      ...  A matrix containing coefficients for the expansion. 
%                  As an input argument is either an NBASIS by NVAR matrix,
%                  in which case NREP = 1, or an NBASIS by NVAR by NREP 
%                  array.
%                  On output it has NVAR columns, and its number of rows is
%                  NBASIS*NREP.
%  BASISVALS  ...  Either a basis object to be used to approximate each
%                  solution, or a struct object containing the values
%                  of the basis function at each observation time.
%                  Defaults to [].
%  FDOBJ      ...  An alternative means of contributing a basis system.
%                  When this object is supplied, its basis overides that
%                  supplied in argument BASISVALS.  But one or the other
%                  must be supplied.
%                  Defaults to [];
%  MORE       ...  Any additional information required by the functions
%                  referenced in argument FN.
%  OBSWEIGHTS ...  Weights to be applied in computing error sums of 
%                  squares over the observed variables.
%  QUADRATURE ...  Locations of quadrature points for approximating the
%                  integrals defining the penalty terms.  If not 
%                  supplied, quadrature points are positioned at the 
%                  centers of the inter-knot intervals.
%  DIFFEPS    ...  A small constant value used to compute difference 
%                  approximations to derivatives if this is requested.
%  POSLIK     ...  If nonzero, the functions approximating the data are
%                  constrained to be positive.
%  POSPROC    ...  If nonzero, the functions used in the penalty terms
%                  are constrained to be positive.
%  DISCRETE   ...  If nonzero, a discrete time model is used instead of
%                  the default continuous time model.

%  Note:  the R version also inputs argument PARS, but only uses it to 
%  get variable names.  This argument is left out here since Matlab
%  cannot address array elements by name.

%  Last modified 1 November 2013

if nargin < 3
    error('Number of arguments is less than three.');
end

%  ------------------------------------------------------------------------
%             Define default values for the arguments
%  ------------------------------------------------------------------------

if nargin < 13, discrete   = 0;     end
if nargin < 12, posproc    = 0;     end
if nargin < 11, poslik     = 0;     end
if nargin < 10, diffeps    = [];    end
if nargin <  9, quadrature = [];    end
if nargin <  8, obsweights = [];    end
if nargin <  7, more       = [];    end
if nargin <  6, fdobj      = [];    end
if nargin <  5, lambda     = [];    end
if nargin <  4, basisvals  = [];    end

if isempty(basisvals) && isempty(fdobj)
    error('Both arguments BASISVALS and FDOBJ are empty.');
end

if isempty(diffeps),  diffeps = 1e-6;  end

% If an fd object is provided, it overrides values of arguments
%      BASISVALUES and COEFS, if supplied.

if ~isempty(fdobj)
    basisvals = getbasis(fdobj);
    if ~isempty(getcoef(fdobj));
        coefs = getcoef(fdobj);
    end
end

%  ------------------------------------------------------------------------
%  Define the dimensions of the analysis:
%  Determine from coefficient array COEF: 
%      number of basis functions       NBASIS
%      number of repeated measurements NREP
%      number of variables             NVAR
%  Determine from TIMES the number of observations per variable N
%  ------------------------------------------------------------------------

if ndims(coefs) > 2
    %  NREP > 1
    [nbasis, nrep, nvar] = size(coefs);
    %  Put coefs array into matrix format:
    %     Number of rows = NBASIS*NREP
    %     Number of cols = NVAR
    coefs = reshape(coefs,nbasis*nrep,nvar);
else
    %  COEFS 2-dimensional.  In this case NREP set to 1
    %  and length of second dimension defined as number of variables 
    [nbasis, nvar] = size(coefs);
    nrep = 1;
end
n = length(times);

%  ------------------------------------------------------------------------
%  Define struct object LIK containing handles for functions for 
%  evaluating the values of the error sum of squares for the data and 
%  their derivatives.
%  Names of the members of the struct object LIK are:
%  'fn'      'dfdx'    'dfdy'    'dfdp'    'd2fdx2'  'd2fdxdy'
%  'd2fdy2'  'd2fdxdp' 'd2fdydp' 'more'    'bvals'
%  ------------------------------------------------------------------------

lik = make_SSElik;

%  Add the member MORE to LIK that defines the transformation of the 
%  process to fit the data.  POS == 0 means no transformation, otherwise an
%  exponential transformation is applied to provide a positive fit.

if poslik==0
    lik.more = make_id;
else
    lik.more = make_exp;
end

%  Members MORE and WHICHOBS are referenced as arguments but 
%  usually not used

lik.more.more     = [];   
lik.more.whichobs = [];

%  ------------------------------------------------------------------------
%  Define struct object PROC containing functions for evaluating the
%  penalty term and its derivatives
%  ------------------------------------------------------------------------

proc = make_SSEproc;

%  Define a struct object PROCMORE containing handles to code for
%  evalution of right side of differential equation

if isstruct(fn)
    %  Argument FN is a struct object with members being function handles 
    %  for functions for evaluating the right side of the
    %  differential equation and its derivatives
    procmore = fn;
    procmore.more = more;  %  contains any additional information req'd
elseif strcmp(class(fn), 'function_handle')
    %  Argument FN is a functional data object for evaluating only 
    %  the value of the right side of the differential equation.
    %  Required derivatives are computed by differencing.
    %  Set up handles to the differencing code
    procmore = make_findif_ode;  
    %  Put right side evaluation code as a member of slot more
    procmore.more.fn   = fn;
    %  Add any additional required information supplied in argument MORE
    procmore.more.more = more;
    %  Add an amount defining differences
    procmore.more.eps  = diffeps;
else
    error('fn must be either a struct of function handles or a function');
end

%  Add  member MORE to PROC containing:
%    PROCMORE if variables are not constrained to have positive values
%    the struct object returned by function MAKE_LOGTRANS with
%    an additional member containing PROCMORE

if posproc==0
    proc.more = procmore;
else
    proc.more      = make_logtrans;
    proc.more.more = procmore;
end
proc.more.whichobs = [];

%  ------------------------------------------------------------------------
%  Define the basis for representing approximate solutions to the 
%  differential equation, and 
%     add member BVALS to the LIK object
%           to contain basis values and first derivative values
%     define quadrature points and add as member MORE to PROC object
%     define basis values and first derivative
%           values and add as member BVALS to PROC object
%  ------------------------------------------------------------------------

if isa_basis(basisvals)
    %  Basis object supplied as an argument.  
    basisobj = basisvals;
    if isempty(times)
        error(['if basisvals is a basis object, ' ...
               'you must specify the observation times']);
    end
    
    %  basis value array for sampling points.  This is always 2-D,
    %  but if there are replications, it is block diagonal with
    %  NREP blocks of basis function values.
    
    if nrep > 1
        %  replications present ... put into block diagonal
        %  format using the kronecker product operator
        lik.bvals = kron(sparse(diag(ones(nrep,1))), ...
                         eval_basis(times,basisobj));
    else
        %  no replications ... straightforward evaluation
        lik.bvals = eval_basis(times,basisobj);
    end
    
    %  define quadrature points
    
    if isempty(quadrature)
        rangeval = getbasisrange(basisobj);
        params   = getbasispar(basisobj);
        knots    = [rangeval(1), params, rangeval(2)];
        quadvals = MakeQuadPoints(knots,5);
        qpts     = quadvals(:,1);
%         qpts = knots(1:length(knots)-1) + diff(knots)/2;
    else
        qpts = quadrature;
    end
    nqpts = length(qpts);
    
    %  basis value array for quadrature points
    
    if discrete==0
        %  continuous time model
        if nrep > 1
            %  replications present ... put into block diagonal
            %  format using the kronecker product operator
            %  values of basis functions
            proc.bvals.bvals  = kron(sparse(diag(ones(nrep,1))), ...
                                     eval_basis(qpts,basisvals));
            %  values of basis function first derivatives
            proc.bvals.dbvals = kron(sparse(diag(ones(nrep,1))), ...
                                     eval_basis(qpts,basisvals,1));
        else
            %  no replications ... straightforward evaluation
            %  values of basis functions
            proc.bvals.bvals  = eval_basis(qpts,basisvals);
            %  values of basis function first derivatives
            proc.bvals.dbvals = eval_basis(qpts,basisvals,1);
        end
        proc.more.weights = repmat(1/nqpts,nqpts*nrep,nvar);
        proc.more.qpts    = qpts;
    else
        %  discrete time model
        basismat          = eval_basis(times,basisvals);
        proc.bvals.bvals  = basismat(1:(n-1),:);
        proc.bvals.dbvals = basismat(2:n,:);
        proc.more.weights = repmat(1/(n-1),(n-1)*nrep,nbasis*nrep);
        proc.more.qpts    = times(1:(n-1));
    end
else
    %  Basis values supplied  as argument
    %  quadrature is ignored
    if discrete && (is.numeric(basisvals) || isempty(basisvals))
        if isempty(basisvals)
            basisvals = Diagonal(size(coefs,1));
        end
        lik.bvals = basisvals;
        n = size(basisvals,1);
        proc.bvals.bvals  = basisvals(1:n-1,:);
        proc.bvals.dbvals = basisvals(2:n,:);
        proc.more.weights = (1/(n-1)).*ones(n-1,nrep);
        proc.more.qpts    = times(1:n-1);
    else
        lik.bvals         = basisvals.bvals.obs;
        proc.bvals.bvals  = basisvals.bvals;
        proc.bvals.dbvals = basisvals.dbvals;
        proc.more.qpts    = basisvals.qpts;
        if ~isempty(basisvals.qwts)
            proc.more.weights = ...
                reshape(basisvals.qwts, ...
                        length(basisvals.qwts)*nrep,nrep);
        else
            nqpts = n;
            proc.more.weights = ones(nqpts*nrep,nrep)/nqpts;
        end
    end
end

%  ------------------------------------------------------------------------
%  Define weights for observations and add them in member MORE
%    of the LIK object
%  ------------------------------------------------------------------------

if  isempty(obsweights)
    obsweights = ones(1,nvar);
end

if numel(obsweights) == nvar
    obsweights = obsweights(:)';
    lik.more.weights = repmat(obsweights,n*nrep*nrep,1);
elseif all(size(obsweights) == [n,nvar])
    lik.more.weights = repmat(obsweights,nrep*nrep,1);
else
    error('Dimensions of observation OBSWEIGHTS is not correct.');
end

%  ------------------------------------------------------------------------
%  Define vector LAMBDA of smoothing parameters and add as a
%  member of member MORE of the PROC object as weights multiplied by 
%  values in LAMBDA
%  ------------------------------------------------------------------------

if length(lambda) > 1
%     disp(size(lambda))
%     disp(size(proc.more.weights))
    proc.more.weights = proc.more.weights*diag(lambda);
else
    if ~isempty(lambda)
        proc.more.weights = lambda.*proc.more.weights;
    end
end

%  ------------------------------------------------------------------------

end
