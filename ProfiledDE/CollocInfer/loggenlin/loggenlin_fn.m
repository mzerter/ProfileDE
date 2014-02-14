function fnval = loggenlin_fn(t,x,p,more)  
% LOGGENLIN_FN computes values of fits to logged observations 
%   at time t.
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

[nt,nx] = size(x);
if nargin < 4,  more = [];  end
%  identify dimensions
np      = length(p);  
if ~isempty(more.sub)
    nsub = size(more.sub,1);
else
    nsub = nx;
end
%  check argument MORE
more  = checkmore(more,nx,np);
%  Set up mapping matrix A
Amat  = more.mat;
if ~isempty(more.sub)
    for i=1:nsub
        Amat(more.sub(i,1),more.sub(i,2)) = p(more.sub(i,3));
    end
end

%  compute contribution to fit values from Y

fitval = x * Amat';

%  compute contributions to fit values from each forcing function vector
if  ~isempty(more.force)
    %  set up matrix of forcing function values at time t
    nf = length(more.force);
    fs = zeros(nt,nf);
    for i = 1:nf
        if isa_fd(more.force{i})
            fs(:,i) = eval_fd(t,more.force{i});
        else
            fs(:,i) = more.force{i}(t,more.force_input);
        end
    end
    %  set up mapping matrix B
    Bmat  = more.force.mat;
    nBsub = size(more.force.sub,1);
    for i=1:nBsub
        Bmat(more.force.sub(i,1),more.force.sub(i,2)) = ...
            p(more.force.sub(:,3));
    end
    %  Add contributions to fit from forcing functions
    fitval = fitval + fs * Bmat';
end

%  return fit to data in log scale

if all(fitval > 0)
    fnval = log(fitval);
else
    error('A model value is not positive.');
end

end
