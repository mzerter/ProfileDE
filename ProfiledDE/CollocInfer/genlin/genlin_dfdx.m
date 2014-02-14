function dfdxval = genlin_dfdx(t,y,p,more)
% GENLIN_DFDX computes partial derivatives of values of fits to 
%  observations at time t with respect to state values in X.
%  Observations are related to variables by the equation
%  Fit = Ay + Bf where y is the column vector of state values at time t,
%  and f is a matrix of forcing function values at time t.
%  Nonzero elements of A and B are taken from parameter vector P.

if nargin < 4,  more = []; end
%  identify dimensions
[nt,ny] = size(y);  np = length(p);  
if ~isempty(more.sub)
    nsub = size(more.sub,1);
else
    nsub = ny;
end
%  check argument MORE
more  = checkmore(more,ny,np);
%  Set up mapping matrix A
Amat  = more.mat;
if ~isempty(more.sub)
    for i=1:nsub
        Amat(more.sub(i,1),more.sub(i,2)) = p(more.sub(i,3));
    end
end
%  compute contribution to fit values from Y
nout = size(more.mat,1);
dfdxval = zeros(nt,nout,ny);
for i=1:nt
    dfdxval(i,:,:) = Amat;
end
%  contributions to fit derivatives from each forcing vector are zero

end
