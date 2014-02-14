

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Some Utilities
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Vnew = NeweyWest_Var(H, g, maxlag)

V = inv(H);
I = 0*H;
if isempty(maxlag)
    n = size(g,1);
    maxlag = max(5, n^(0.25));
end
if maxlag > 0
    [m,n] = size(g);
    for i = 1:n-1
        for j = (i+1):n
            I(i,j) = Newey_West(g(:,i),g(:,j),maxlag);
            I(j,i) = I(i,j);
        end
    end
end
Vnew = V * ( I+ g'*g) * V;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Covar = ...
    Profile_covariance(pars, times, data, coefs, lik, proc, ...
                       active, in_method, options_in, eps)
                   
if nargin < 10, eps = 1e-6;       end
if nargin <  9, options_in = [];  end
if nargin <  8, in_method = [];   end
if nargin <  7, active = [];      end
if isempty(active)
    active = 1:length(pars);
end
apars = pars(active);
H = zeros(length(apars),length(apars));
sumlik = 0;
g = ProfileDP(apars, pars, times, data, coefs, lik, proc, active, sumlik);
gg = applysum(g);
for i = 1:length(apars)
%     if file.exists('curcoefs.tmp'))file.remove('curcoefs.tmp')end
%     if file.exists('optcoefs.tmp'))file.remove('optcoefs.tmp')end
%     if file.exists('counter.tmp'))file.remove('counter.tmp')end
    tpars = apars;
    tpars(i) = tpars(i) + eps;
    sumlik = 1;
    tf = ProfileErr(tpars, pars, times, data, coefs, lik, proc, ...
                    in_method, options_in, active);
    tg = ProfileDP( tpars, pars, times, data, coefs, lik, proc, ...
                    active, sumlik);
    H(:,i) = (tg-gg)/eps;
%     if file.exists('curcoefs.tmp'))file.remove('curcoefs.tmp')end
%     if file.exists('optcoefs.tmp'))file.remove('optcoefs.tmp')end
%     if file.exists('counter.tmp'))file.remove('counter.tmp')end
end
Covar = NeweyWest_Var( 0.5*(H'+H) , g,5);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sumval = applysum(x)

switch ndims(x)
    case 2
        sumval = sum(x);
    case 3
        sumval = sum(squeeze(sum(x,1),2));
    otherwise
        error('X now two- or three-dimensional.');
end

end

function [newdata, dims] = dataformat(data)

dims = size(data);
if length(dims)>2
    [l,m,n] = size(data);
    newdata = reshape(data,l*m,n);
    dims(1) = length(coefs)/(m*n);
else
    newdata = data;
    dims    = size(coefs);
end

end

function [coefs, nrep] = coefsformat(coefs)

if ndims(coefs) > 2
%     if isempty(colnames)
%         colnames = dimnames(coefs)((3))
%     end
    [nbasis, nrep, nvar] = size(coefs);
    coefs = reshape(coefs,nbasis*nrep,nvar);
else
    nrep = 1;
%     colnames = colnames(coefs)
end

end

function out = blocks2mat(H)  
% Cell array of matrices -> large matrix

out = [];
for i = 1:length(H)
    tout = H{i,1};
    if length(H{i,j}) > 1
        for j = 2:length(H{i,j})
            tout = [tout,H{i,j}];
        end
    end
    if i > 1
        out = [out; tout];
    else
        out = tout;
    end
end

end

% Newey West Calculations, with thanks to Steve Ellner

% GAUSS trimr function: trims n1 rows from the start and n2 rows
%  from the end of a matrix or vector 

function atrim = trimr (a,n1,n2)

da = size(a);
if isempty(da)
    atrim = a((n1+1):(length(a)-n2));
else
    atrim = a((n1+1):(da(1)-n2),:);
end

end

function out = Newey_West(x,y,maxlag)

w=1-(1:maxlag)/(maxlag+1); 
w=w/length(x);
out=mean(x*y);
for i = 1:maxlag
    out = out + w(i)*sum(trimr(x,i,0)*trimr(y,0,i)) + ...
                w(i)*sum(trimr(y,i,0)*trimr(x,0,i));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function weightsmat = checkweights(weights, whichrows, diffs)

[m,n] = size(diffs);
if isempty(whichrows)
    whichrows = 1:m;
end
if isempty(weights)
    weightsmat = ones(m,n);
elseif length(weights) == n
    weightsmat = reshape(weights,n,m)';
elseif length(weights(whichrows)) == m
    weightsmat = reshape(weights,n,m)';
elseif all(size(weights(whichrows,:)) == size(diffs))
    weightsmat = weights(whichrows,:);
else
    error('Dimension of weights does not match that of data');
end

end

%  ------------------------------------------------------------------------

%%%%%%% Squared Process Error %%%%%%%%%%%%%%%%%%%%%%

%  this file contains functions for evaluating the total SSE for ODE
%  and its derivative values:
%  It uses function make_SSElik but substitutes derivative value
%  for data value and quadrature point values for time values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




