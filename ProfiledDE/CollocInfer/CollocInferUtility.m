

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



