function fnval = genlin_fn(t,y,p,more)  
% GENLIN_FN computes values of fits to observations at time t.
%  Observations are related to variables by the equation
%  Fit = Ay + Bf where y is the column vector of state values at time t,
%  and f is a matrix of forcing function values at time t.
%  Nonzero elements of A and B are taken from parameter vector P.

if nargin < 4,  more = [];  end
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
%  compute row vector of contribution to fit values from Y
fnval = y*Amat';
%  compute contributions to fit values from each forcing function vector
if  ~isempty(more.force)
    %  set up matrix of forcing function values at time t
    nf = length(more.force);
    fs = zeros(nt,nf);
    for i = 1:nf
        if isa_fd(more.force{i})
            fs(:,i) = eval_fd(t,more.force{i});
        else
            fs(:,i) = more.force{i}(t,more.force.input);
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
    fnval = fnval + fs * Bmat';
end

end
