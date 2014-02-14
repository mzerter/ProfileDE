function fnval = genlin_fun_ode(t,y,p,more)
% GENLIN_FUN_ODE computes values of fits to observations at time t.
%  Observations are related to variables by the equation
%  Fit = Ay + Bf where y is the column vector of state values at time t,
%  and f is a matrix of forcing function values at time t.
%  Nonzero elements of A and B are taken from parameter vector P.

if nargin < 4,  more = []; end
%  check argument MORE
n    = size(y,2);
npar = length(p);
more = checkmore(more,n,npar);
%  set up mapping matrix A
[Arow,Acol] = size(more.mat);
Amat  = zeros(Arow,Acol);
nAsub = size(more.sub,1);
for i=1:nAsub
    Amat(more.sub(i,1),more.sub(1,2)) = p(more.sub(i,3));
end
r1 = y * Amat';
%  compute contributions to fit values from each forcing function vector
if  ~isempty(more.force)
    %  set up matrix of forcing function values at time t
    fs = zeros(length(t),length(more.force));
    for i = 1:length(more.force)
        if isa_fd(more.force{i})
            fs(:,i) = eval_fd(t,more.force{i});
        else
            fs(:,i) = more.force{i}(t,more.force_input);
        end
    end
    %  set up mapping matrix B
    [Brow,Bcol] = size(more.force.mat);
    Bmat  = zeros(Brow,Bcol);
    nBsub = size(more.force.sub,1);
    for i=1:nBsub
        Bmat(more.force.sub(i,1),more.force.sub(i,2)) = ...
            p(more.force.sub(:,3));
    end
    %  Add contributions to fit from forcing functions
    r1 = r1 + fs * Bmat';
end
fnval = r1;
end
