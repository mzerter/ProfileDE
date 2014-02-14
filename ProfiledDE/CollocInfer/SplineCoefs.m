function [f, dfdc, d2fdc2, d2fdcdp, gdata, geqtn] = ...
    SplineCoefs(coefs, times, data, pars, lik, proc)
%  Evaluate the inner optimization objective function and their
%  partial derivatives for a set of coefficient and parameter values
%
%  Arguments:   (All are required)
%  COEFS      ...  A matrix containing coefficients for the expansion.  
%                  It has NVAR columns, and its number of rows is
%                  NBASIS*NVAR.
%  TIMES      ...  The times of observation of those variables that are
%                  measured.  These times are assumed to be common to all
%                  measured variables.
%  DATA       ...  Data matrix.  The number of columns is equal to the
%                  number of variables (including unobserved variables).
%                  The number of columns is equal to the number of 
%                  observation points per variable times the number of
%                  replications.
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

%  Last modified 8 January 2011

if nargin < 6
    error('Not all arguments are supplied.');
end

%  Codefficient array may come in as a vector ... reshape it to be a
%  matrix with NBASIS rows

nbasis = size(lik.bvals,2);
coefs  = reshape(coefs(:),nbasis,numel(coefs)/nbasis);
nvar   = size(coefs,2);

%  values of data fitting functions

x  = lik.bvals*coefs;

% error sum of squares for data
 
fdata  = sum(lik.fn(data, times, x, pars, lik.more));
% disp(['fdata = ', num2str(fdata)])

%  approximate error integral of squares for the penalty term

feqtn = proc.fn(coefs, proc.bvals, pars, proc.more);
% disp(['feqtn = ', num2str(feqtn)])

%  total error sum of squares

f = fdata + feqtn;

%  evaluate the coefficient-gradient if required

if nargout > 1
    gdata = lik.bvals'*lik.dfdx(data,times,x,pars,lik.more);
    geqtn = proc.dfdc(coefs, proc.bvals, pars, proc.more);
    g     = gdata + geqtn;
    dfdc  = g(:);
end

%  evaluate the c-Hessian if required

if nargout > 2
    d2likdx2 = lik.d2fdx2(data, times, x, pars, lik.more);
    [l,m,n] = size(d2likdx2);
    H = cell(m,n);
    for i = 1:m
        for j = 1:n
            H{i,j} = lik.bvals'*diag(d2likdx2(:,i,j))*lik.bvals;
        end
    end
    d2likdc2  = sparse(blocks2mat(H));
    d2procdc2 = sparse(proc.d2fdc2(coefs, proc.bvals, pars, proc.more));
    d2fdc2    = sparse(d2likdc2 + d2procdc2);
end

%  evaluate the c-p cross derivative if required

if nargout > 3
    d2likdxdp = lik.d2fdxdp(data, times, x, pars, lik.more);
    d2likdcdp = zeros(nvar*nbasis,length(pars));
    for i = 1:length(pars)
        H1i = lik.bvals'*squeeze(d2likdxdp(:,:,i));
        d2likdcdp(:,i)  = H1i(:);
    end
    d2procdcdp = proc.d2fdcdp(coefs, proc.bvals, pars, proc.more);
    d2fdcdp    = sparse(d2likdcdp + d2procdcdp);
end

end

