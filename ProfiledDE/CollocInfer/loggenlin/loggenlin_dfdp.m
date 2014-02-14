function dfdpval = loggenlin_dfdp(t,x,p,more)
% LOGGENLIN_DFDP computes partial derivatives of values of fits to 
%  logged observations at time t with respect to parameters.
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
%  compute contributions to fit derivatives with respect to P from Y
fitval = x * Amat';
if any(fitval <= 0)
    error('A model value is not positive.');
end
%  compute di ln Z di theta
[nout, na] = size(Amat);
dfdpval = zeros(nt,np,nout);
for k = 1:nout
    dfdpval(:,more.sub(k,1),more.sub(k,3)) = x(more.sub(k,2))/fitval(:,k);
end
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
    for i = 1:nf
        dfdpval(:,more.force.sub(i,1),more.force.sub(i,3)) =  ...
        dfdpval(:,more.force.sub(i,1),more.force.sub(i,3)) +  ...
             fs(:,more.force.sub(i,2))./fitval;
    end
end

end
