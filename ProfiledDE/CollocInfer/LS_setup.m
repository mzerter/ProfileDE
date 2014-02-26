function [lik, proc, coefs] = ...
    LS_setup(fn, times, coefs, basisvals, lambda, fdobj, ...
             more, data, obsweights, quadrature, ...
             diffeps, poslik, posproc, discrete, likfn, likmore)
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
%  LAMBDA     ...  A single constant, or a vector of length NVAR containing
%                  smoothing parameter values for each variable controlling
%                  how closely a variable comes to solving the equation.
%  FDOBJ      ...  An alternative means of contributing a basis system.
%                  When this object is supplied, its basis overides that
%                  supplied in argument BASISVALS.  But one or the other
%                  must be supplied.
%                  Defaults to [];
%  MORE       ...  Any additional information required by the functions
%                  referenced in argument FN.
%  DATA       ...  Array of data.  Defaults to [];
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
%  LIKFN      ...  Function handle to a defined function for mapping 
%                  trajectories into observations.  Defaults to make_id.
%  LIKMORE    ...  Additional information for LIKFN

%  Last modified 25 Feburary 2014 (Mathieu)
disp('LS_Setup Version: 25 Feburary 2014 (Mathieu)')

if nargin < 3
    error('Number of arguments is less than three.');
end

%  ------------------------------------------------------------------------
%             Define default values for the arguments
%  ------------------------------------------------------------------------

if nargin < 16, likmore    = [];    end
if nargin < 15 || isempty(likfn),   likfn    = make_id;end
if nargin < 14 || isempty(discrete), discrete = 0;      end
if nargin < 13 || isempty(posproc),  posproc  = 0;      end
if nargin < 12 || isempty(poslik),   poslik   = 0;      end
if nargin < 11 || isempty(diffeps),  diffeps  = 1e-6;   end
if nargin < 10, quadrature = [];    end
if nargin <  9, obsweights = [];    end
if nargin <  8, data       = [];    end
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

n = length(times);

%  set up matrix or array COEFS along with dimensons NBASIS and NREP

if ~isempty(coefs) && strcmp(class(coefs),'double')
    if ndims(coefs) > 2
        %  NREP > 1
        [nbasis,nvar,nrep] = size(coefs);
        %  Put coefs array into matrix format:
        %     Number of rows = NBASIS*NREP
        %     Number of cols = NVAR
        coefs = reshape(coefs,nbasis*nrep,nvar);
    else
        %  COEFS 2-dimensional.  In this case NREP set to 1
        %  and length of second dimension defined as number of variables
        [nbasis,nvar] = size(coefs);
        nrep = 1;
    end
end

%  If COEFS is empty, set up coefs as zero matrix with NREP = 1

if isempty(coefs)
    if ~isempty(basisvals)
        if isa_basis(basisvals)
            nbasis = getnbasis(basisvals);
        else
            if strcmp(class(basisvals),'double') && ndims(basisvals) == 2
                nbasis = size(basisvals,1);
            else
                nbasis = size(basisvals.bvals.obs,1);
            end
        end
        if ~isempty(data)
            coefs = zeros(nbasis,size(data,2));
        else
            if ~isempty(weights)
                coefs = zeros(nbasis,size(weights,1));
            else
                error(['Cannot determine the dimension of the state ', ...
                       'vector -- please provide one of the ', ...
                       'arguments COEFS or FDOBJ.']);
            end
        end
        nrep = 1;
    end
end

%  ------------------------------------------------------------------------
%  Define struct object LIK containing handles for functions for 
%  evaluating the values of the error sum of squares for the data and 
%  their derivatives.
%  Names of the members of the struct object LIK are:
%  'fn'      'dfdx'    'dfdy'    'dfdp'    'd2fdx2'  'd2fdxdy'
%  'd2fdy2'  'd2fdxdp' 'd2fdydp' 'more'    'bvals'
%  ------------------------------------------------------------------------

%  Level 1 of lik is struct containing functions for evaluating error 
%  sum of squares for fit to data

lik = make_SSElik;

%  Add the member MORE to LIK that defines the transformation of the 
%  process to fit the data.  POS == 0 means no transformation, otherwise an
%  exponential transformation is applied to provide a positive fit.

if ~poslik                              % Map from stats to obs                  
    if isstruct(likfn)                  % All derivatives available analytically
        lik.more = likfn;
        if ~isfield(lik.more, 'more') || isempty(lik.more.more)
            lik.more.more = likmore;
        end
    else                                % Finite-difference for derivatives
        
        

        lik.more = make_findif_ode;
        
        likmoremore.fn   = likfn;
        likmoremore.more = likmore;
        likmoremore.eps  = diffeps;
        
        lik.more.more    = likmoremore;
    end
else                                    % States given on log scale
    if isstruct(likfn)                  % All derivatives available analytically
        lik.more = make_exp;            
        lik.more.more = likfn;
        if ~isfield(lik.more.more, 'more') || isempty(lik.more.more.more)
            lik.more.more.more = likmore;
        end
    else                                % Finite-difference for derivatives

        lik.more = make_findif_ode;
        temp = make_exp;
        lik.more.more.fn = temp.fn;
        
        
        likmoremoremore.fn   = likfn;
        likmoremoremore.more = likmore;
        likmoremoremore.eps  = diffeps;
        
        lik.more.more.more = likmoremoremore;        
    end
        
%    likfn = make_id;   <---- This default value should be set somewhere at
%    beginning? See issue #3
end

%                  <---------------- This section doesn't make sense now?
%  Members MORE and WHICHOBS are referenced as arguments but 
%  usually not used

% lik.more.more     = [];   
% lik.more.whichobs = [];
%                  <--------------- 


%  ------------------------------------------------------------------------
%  Define struct object PROC containing functions for evaluating the
%  penalty term and its derivatives
%  ------------------------------------------------------------------------

%  Level 1 of proc is struct containing functions for evaluating error 
%  sum of squares for fit to derivative of trajectory

proc = make_SSEproc;

%  Define a struct object PROCMORE containing handles to code for
%  evalution of right side of differential equation

if isstruct(fn)
    %  Argument FN is a struct object with members being function handles 
    %  for functions for evaluating the right side of the
    %  differential equation and its derivatives
    procmore      = fn;
    procmore.more = more;  %  contains any additional information req'd
    findif = 0;
else
    %  Argument FN is a functional data object for evaluating only 
    %  the value of the right side of the differential equation.
    %  Required derivatives are computed by differencing.
    %  Set up handles to the differencing code
    if strcmp(class(fn), 'function_handle')
        tempstruct.fn   = fn;
        tempstruct.more = more;
        tempstruct.eps  = diffeps;
        temp.more = tempstruct;
        
 % else if FN is a 'pomp' function. See R code. <------------------- POMP
 
    else
        error('FN is not a function_handle object.'); 
    end
    
    % base level 1 of procmore contains discrete approximation functions
    
    procmore = make_findif_ode; 
    
    if posproc
        %  penalty to be evaluated on log scale, level 2 contains the
        %  logging functions
        
        logtranstemp = make_logtrans;

        tempstruct.fn   = logtranstemp.fn;
        tempstruct.more = temp.more;
%        tempstruct.eps  = diffeps;  <------------ Line 285 Makes redundant
        procmore.more   = tempstruct;
    else
        %  level 2 of procmore contains equation RHS evaluation function
        procmore.more = temp.more;
    end
    %  Add an amount defining differences
    procmore.more.eps  = diffeps;
    findif = 1;
end

%  Add  member MORE to PROC containing:
%    PROCMORE if variables are not constrained to have positive values
%    the struct object returned by function MAKE_LOGTRANS with
%    an additional member containing PROCMORE

if posproc==0 || findif
    %  level 2 of proc is procmore
    proc.more = procmore;
else
    %  level 2 of proc is struct containing logging functions
    proc.more      = make_logtrans;
    %  level 3 of proc is procmore
    proc.more.more = procmore;
end

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
    
    if isempty(quadrature) || isempty(quadrature.qpts)
        rangeval = getbasisrange(basisobj);
        params   = getbasispar(basisobj);
        knots    = [rangeval(1), params, rangeval(2)];
        quadvals = MakeQuadPoints(knots,5);
        qpts     = quadvals(:,1);
    else
        if strcmp(class(quadrature),'double')
            qpts = quadrature;
        else
            if isstruct(quadrature)
                qpts = quadrature.qpts;
            else
                error(['Argument QUADRATURE is neither a matrix nor ', ...
                       'a struct object.']);
            end
        end
    end
    qdims = size(qpts);
    if qdims(1) == 1 || qdims(2) == 1
        qpts = qpts(:);
    else
        error('Quadrature points not a vector.');
    end
    nqpts = length(qpts);
    
    %  basis value array for quadrature points
    
    if ~discrete
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
        proc.more.qpts    = repmat(qpts,nrep,1);
    else
        %  discrete time model
        basismat          = eval_basis(times,basisvals);
        proc.bvals.bvals  = repmat(basismat(1:(n-1),:),nrep,1);
        proc.bvals.dbvals = repmat(basismat(2:n,:),nrep,1);
        proc.more.weights = repmat(1/(n-1),(n-1)*nrep,nbasis*nrep);
        proc.more.qpts    = repmat(times(1:(n-1)),nrep,1);
    end
else
    %  Basis values supplied  as argument
    %  quadrature is ignored
    if discrete && (is.numeric(basisvals) || isempty(basisvals))
        if isempty(basisvals)
            basisvals = Diagonal(size(coefs,1));
        end
        lik.bvals = repmat(basisvals,nrep,1);
        n = size(basisvals,1);
        bvalsstruct.bvals  = repmat(basisvals(1:n-1,:),   nrep,1);
        bvalsstruct.dbvals = repmat(basisvals(2:n,:),     nrep,1);
        proc.bvals = bvalsstruct;
        proc.more.weights = repmat((1/(n-1))*ones(n-1,1),nrep,1);
        proc.more.qpts    = repmat(times(1:n-1),         nrep,1);
    else
        lik.bvals         = repmat(basisvals.bvals.obs,nrep,1);
        tempstruct.bvals  = repmat(basisvals.bvals,    nrep,1);
        tempstruct.dbvals = repmat(basisvals.dbvals,   nrep,1);
        proc.bvals  = tempstruct;
        proc.more.qpts = repmat(basisvals.qpts,1,nrep);
        if ~isempty(basisvals.qwts)
            proc.more.weights = remat(basisvals.qwts,nrep,nrep);
        else
            nqpts = n;
            proc.more.weights = ones(nqpts*nrep,nrep)/nqpts;
        end
    end
end

%  check consistency of data dimensions with dimensions of other arrays

if ~isempty(data)
    if length(size(data)) == 2
        if nrep > 1
            error(['DATA dimensions inconsistent with coefficient ', ...
                   'dimensions.']);
        end
        if size(data,1) ~= length(times)
            error('DATA dimensions inconsistent with length of TIMES.');
        end
    else
        if size(data,2) ~= nrep || size(data,1) ~= length(times)
            error(['DATA dimensions, length of TIMES and COEFS ', ...
                   'dimensions are inconsistent.']);
        end
        data = reshape(data,n*nrep,nvar);
        times = repmat(times,nrep,1);
    end
end

%  ------------------------------------------------------------------------
%  Define weights for observations and add them in member MORE
%    of the LIK object
%  ------------------------------------------------------------------------

if isempty(obsweights)
    if ~isempty(data)
        lik.more.weights = ones(n,nrep);
    else
        lik.more.weights = [];
    end
else
    if length(size(obsweights)) == 3
        weights = ...
            reshape(obsweights,size(obsweights,1),size(obsweights,2));
        if ~isempty(data)
            weights = checkweights(weights, [], data);
        elseif length(size(obsweights)) == 2
            lik.more.weights = ...
            reshape(obsweights,size(obsweights,1)*nrep,size(obsweights,2));
        end
    else
        lik.more.weights = obsweights;
    end  
end

%  ------------------------------------------------------------------------
%  Define vector LAMBDA of smoothing parameters and add as a
%  member of member MORE of the PROC object as weights multiplied by 
%  values in LAMBDA
%  ------------------------------------------------------------------------

if length(lambda) > 1
    proc.more.weights = proc.more.weights*diag(lambda);
else
    proc.more.weights = lambda.*proc.more.weights;
end

%  ------------------------------------------------------------------------

end
