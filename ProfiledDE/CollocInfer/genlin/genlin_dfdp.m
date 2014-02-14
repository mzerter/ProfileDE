function dfdpval = genlin_dfdp(t,y,p,more)
% GENLIN_DFDP computes partial derivatives of values of fits to 
%  observations at time t with respect to parameters.
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
%  compute contributions to fit derivatives with respect to P from Y
nout = size(more.mat,1);
dfdpval = zeros(nt,nout,np);
for i = 1:nsub
    dfdpval(:,more.sub(i,1),more.sub(i,3)) = y(:,more.sub(i,2));
end
%  compute contributions to fit derivatives with respect to P 
%     from each forcing vector
if ~isempty(more.force)
    %  set up matrix of forcing function values at time t
    nf = length(more.force);
    fs = zeros(nt,nf);
    for i = 1:nf
        if is_fd(more.force{i})
            fs(:,i) = eval_fd(t,more.force{i});
        else
            fs(:,i) = more.force{i}(t,more.force.input);
        end
    end
    %  Add contributions to fit from forcing functions
    nsub = size(more.force.sub,1);
    for i = 1:nsub
        dfdpval(:,more.force.sub(i,1),more.force.sub(i,3)) =  ...
        dfdpval(:,more.force.sub(i,1),more.force.sub(i,3)) +  ...
                 fs(:,more.force.sub(i,2));
    end
end

end
