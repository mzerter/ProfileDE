function [fnval, dfdpval] = ...
            ParsMatch(pars, coefs, proc, allpars, active)  

%  Last modified 8 January 2011

if nargin < 5 
    active=1:length(pars);  
    allpars = pars;
end

allpars(active) = pars;
fnval   = proc.fn(  coefs,proc.bvals,allpars,proc.more);
dfdpval = proc.dfdp(coefs,proc.bvals,allpars,proc.more);
dfdpval = dfdpval(active);

end
