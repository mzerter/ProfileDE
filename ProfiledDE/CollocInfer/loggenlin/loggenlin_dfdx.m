function dfdxval = loggenlin_dfdx(t,x,p,more)
% LOGGENLIN_FN computes values of fits to logged
%  observations at time t with respect to state values in X.
%  Observations are related to variable values x by the equation
%  Fit = xA' + fB' where x is the NT by NX matrix of state values at 
%  time t, and f is a matrix of forcing function values at time t.
%  Nonzero elements of A and B are taken from parameter vector P.

if any(x > 50)
    error(['Probable error: ', ...
           'a value of X exceeds 50 before exponentiation.']);
end

%  put x into data scale

x = exp(x);
p = exp(p);

if nargin < 4,  more = []; end
%  identify dimensions
[nt,nx] = size(x);  
np = length(p);  
%  check argument MORE
more  = checkmore(more,nx,np);
%  Set up mapping matrix A
Amat  = more.mat;
%  check nsub
if ~isempty(more.sub)
    nsub = size(more.sub,1);
else
    nsub = nx;
end
%  set entries in Amat
if ~isempty(more.sub)
    for i=1:nsub
        Amat(more.sub(i,1),more.sub(i,2)) = p(more.sub(i,3));
    end
end
%  compute contribution to fit values from Y
fitval = x * Amat';
if any(fitval <= 0)
    error('A model value is not positive.');
end
[nout, na] = size(Amat);
dfdxval = zeros(nt,nout,nx);
for k=1:nout
    for j=1:nx
        dfdxval(:,k,j) = Amat(k,j)./fitval(:,k);
    end
end
end
