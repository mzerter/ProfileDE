function [res, Dres] = ProfileSSE(pars, times, data, coefs, allpars, ... 
                                   lik, proc, active, ...
                                   in_method, options_in, dcdp, oldpars)                            

if nargin < 12, oldpars    = [];  end
if nargin < 11, dcdp       = [];  end
if nargin < 10, options_in = [];  end
if nargin <  9, in_method  = [];  end
if nargin <  8, active     = [];  end

if isempty(active)
    active=1:length(allpars); 
end

allpars(active) = pars;

[res, Dres] = ProfileSSE_AllPar(allpars, times, data, coefs, lik, proc, ...
                                in_method, options_in, dcdp, oldpars);
Dres = Dres(:,active);

end





