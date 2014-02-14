function fnval = ProfileErr_AllPar(pars, times, data, coefs, lik, proc, ...
                                   in_method, options_in, sgn)
% Outer Optimization Objective

if nargin < 9, sgn = 1;          end
if nargin < 8, options_in = [];  end
if nargin < 7, in_method  = [];  end
coefs = reshape(coefs,numel(coefs),1); 
% if file.exists('curcoefs.tmp'))                      
%     altcoefs = read.table('curcoefs.tmp'))               
%     if  ~(length(altcoefs)==length(coefs)) )
%         error(['Variables = curcoefs.tmp do not conform;  ', ...
%                'file exists from previous experiments?']);
%     end
% else
% altcoefs = coefs;
% end
% if file.exists('counter.tmp')                       
%     counter = read.table('counter.tmp')                  
%     niter = counter(nrow(counter),1)
% else
%     counter = reshape([1,0,pars],1,length(pars)+2)
% end
[m,n] = size(lik.bvals);
% altdevals = lik.bvals*reshape(altcoefs, m, length(altcoefs)/n);
% colnames(altdevals) = proc.more.names
f1 = SplineCoefs(coefs, times, data, lik, proc, pars);
f2 = SplineCoefs(as.vector(altcoefs), times, data, lik, proc, pars);
if f2 < f1
    coefs = altcoefs;
end
Ires = inneropt(data, times, pars, coefs, lik, proc, ...
                options_in, in_method);
ncoefs = Ires.coefs;
devals = lik.bvals * ncoefs;
% colnames(devals) = proc.more.names

f  = sum(lik.fn(data, times, devals,    pars, lik.more));
% f2 = sum(lik.fn(data, times, altdevals, pars, lik.more));
% if f <= f2
%     write.table(ncoefs, file='curcoefs.tmp', col.names=FALSE, row.names=FALSE)
%     %      save(ncoefs, niter, pars, file='bestcoefs.Rdata')
% end
% write.table(ncoefs, file='optcoefs.tmp', col.names=FALSE, row.names=FALSE)

% if niter==0)
%     counter(1, 2) = f
%     write.table(counter, file='counter.tmp', col.names=FALSE, row.names=FALSE)
% end

% if niter>=1
%     if f < counter(niter,2)
%         counter = rbind(counter, c(niter+1, f, pars))
%         write.table(counter, file='counter.tmp', col.names=FALSE, row.names=FALSE)
%     end
% end
if ~isempty(lik.report)
    print(f)
end
fnval = sgn*f;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dpval = ProfileDP_AllPar(pars, times, data, coefs, lik, proc, ...
    in_method, options_in, sgn, sumlik)

if nargin < 10, sumlik = 1;          end
if nargin <  9, sgn = 1;             end
if nargin <  8, options_in = [];     end
if nargin <  7, in_method  = [];     end

% if file.exists('optcoefs.tmp')
%     altcoefs = read.table('optcoefs.tmp')
%     if  ~(length(altcoefs)==length(coefs))
%         error('Variables = curcoefs.tmp do not conform; ', ...
%               'file exists from previous experiments?')
%     else
% coefs = altcoefs;
%     end
% end
devals = lik.bvals * coefs;
% colnames(devals) = proc.more.names;
coefs   = reshape(coefs,numel(coefs),1);
[f, dfdc, d2fdc2, d2fdcdp] = ...
    SplineCoefs(coefs, times, data, lik, proc, pars);
dcdp = -ginv(d2fdc2) * d2fdcdp;
if sumlik
    dlikdc = lik.bvals'*lik.dfdx(data, times, devals, pars, lik.more);
    dlikdc = reshape(dlikdc,numel(dlikdc),1);
%     dpterm = lik.dfdp(data, times, devals, pars, lik.more);
    df = dcdp'*dlikdc + ...
         applysum(lik.dfdp(data, times, devals, pars, lik.more));
    df = sgn*df;
    dpval = reshape(df,numel(df),1);
else
    dlikdx = lik.dfdx(data, times, devals, pars, lik.more);
    
    dlikdp = lik.dfdp(data, times, devals, pars, lik.more);
    
    dlikdc = [];
    for i = 1:size(dlikdx,2)
        dlikdc = [dlikdc, diag(dlikdx(:,i))*lik.bvals];
    end
    dpval = dlikdc * dcdp + dlikdp;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = ProfileErr(pars, allpars, times, data, coefs, lik, proc, ...
                        active, in_method, options_in, sgn)

if nargin < 11, sgn = 1;                  end
if nargin < 10, options_in = [];          end
if nargin <  9, in_method  = [];          end
if nargin <  8, active=1:length(allpars); end
allpars(active) = pars;
f = ProfileErr_AllPar(allpars, times, data, coefs, lik, proc, ...
                      in_method, options_in, sgn);
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = ProfileDP(pars, allpars, times, data, coefs, lik, proc, ...
                       active, in_method, options_in, sgn, sumlik)

if nargin < 12, sumlik = 1;               end
if nargin < 11, sgn = 1;                  end
if nargin < 10, options_in = [];          end
if nargin <  9, in_method  = [];          end
if nargin <  8, active=1:length(allpars); end
allpars(active) = pars;
g = ProfileDP_AllPar(allpars, times, data, coefs, lik, proc, ...
                     in_method, options_in, sgn, sumlik);
if sumlik
%     names(g) = names(allpars)
    g = g(active);
else
%     colnames(g) = names(allpars)
    g = g(:,active);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value, gradient] = ...
    ProfileList(pars, allpars, times, data, coefs, lik, proc, ...
                active, in_method, options_in, sgn)

if nargin < 11, sgn = 1;                  end
if nargin < 10, options_in = [];          end
if nargin <  9, in_method  = [];          end
if nargin <  8, active=1:length(allpars); end
value    = ProfileErr(pars, allpars, times, data, coefs, lik, proc, ...
                      in_method, options_in, sgn, active);
gradient = ProfileDP( pars, allpars, times, data, coefs, lik, proc, ...
                      in_method, options_in, sgn, active);

end
