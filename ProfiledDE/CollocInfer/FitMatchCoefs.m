function [fnval, dfdc, d2fdc2] = ...
    FitMatchCoefs(coefs, allcoefs, which, pars, proc)

allcoefs(:,which) = reshape(coefs,size(allcoefs,1),length(which));

fnval = proc.fn(allcoefs,proc.bvals,pars,proc.more);

if nargout > 1
    [n,npars] = size(allcoefs);
    dfdc = proc.dfdc(allcoefs,proc.bvals,pars,proc.more);
    g    = reshape(dfdc,size(allcoefs));
    dfdc = g(:,which);
    dfdc = dfdc(:);
end

if nargout > 2
    [m,n] = size(allcoefs);
    ind = reshape(1:numel(allcoefs),m,n);
    ind = ind(:,which);
    H = proc.d2fdc2(allcoefs,proc.bvals,pars,proc.more);
    d2fdc2 = H(ind,ind);
end

end