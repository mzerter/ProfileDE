function [value, gradient] = ...
    ProfileErr(pars, times, data, coefs, allpars, lik, proc, ...
               active, in_method, options_in)

if nargin < 10, options_in = [];  end
if nargin <  9, in_method  = [];  end
if nargin <  8, active     = [];  end

if isempty(active)
    active=1:length(allpars); 
end

allpars(active) = pars;

[value, gradient] = ...
    ProfileErr_AllPar(allpars, times, data, coefs, lik, proc, ...
                      in_method, options_in);

gradient = gradient(active);

end
