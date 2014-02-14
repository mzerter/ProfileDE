function d2fdxdpval = genlin_d2fdxdp(t,y,p,more)
% GENLIN_D2FDXDP computes second partial derivatives of values of fits to 
%  observations at time t with respect to state values in Y and
%  parameter values in P
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
more = checkmore(more,ny,np);
nout = size(more.mat,1);
%  compute contribution to fit derivatives from Y
d2fdxdpval = zeros(nt,nout,ny,np);
for i = 1:nsub;
    d2fdxdpval(:,more.sub(i,1),more.sub(i,2),more.sub(i,3)) = 1;
end
%  contributions to fit derivatives from each forcing vector are zero

end

