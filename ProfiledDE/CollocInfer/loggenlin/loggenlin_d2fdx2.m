function d2fdx2val = loggenlin_d2fdx2(t,x,p,more)
%  LOGGENLIN_D2FDX2 computes second partial derivatives of values of fits  
%  to logged observations at time t with respect to state values in X.
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
if nargin < 4,  more = []; end
%  identify dimensions
[nt,nx] = size(x);  
np = length(p);  
%  check argument MORE
more  = checkmore(more,nx,np);
%  Set up mapping matrix A
Amat  = more.mat;
%  check argument MORE
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
%  hessian
nout = size(Amat,1);
d2fdx2val = zeros(nt,nout,nx,nx);
for i=1:nt
    fitvali = fitval(i,:);
    fitmati = repmat(fitvali',1,nx);
    temp = Amat./fitmati;
    for j=1:nout
        d2fdx2val(i,j,:,:) = -temp'*temp;
    end
end
%  contributions to fit derivatives from each forcing vector are zero
end
