function d2fdx2val = genlin_d2fdx2(t,y,p,more)
% GENLIN_D2FDX2 computes second partial derivatives of values of fits to 
%  observations at time t with respect to state values in X.
%  Observations are related to variables by the equation
%  Fit = Ay + Bf where y is the column vector of state values at time t,
%  and f is a matrix of forcing function values at time t.
%  Nonzero elements of A and B are taken from parameter vector P.

if nargin < 4,  more = []; end
%  identify dimensions
[nt,ny] = size(y);  np = length(p); 
%  check argument MORE
more  = checkmore(more,ny,np);
%  compute contribution to fit values from Y
nout = size(more.mat,1);
d2fdx2val = zeros(nt,nout,ny,ny);
%  contributions to fit derivatives from each forcing vector are zero

end
