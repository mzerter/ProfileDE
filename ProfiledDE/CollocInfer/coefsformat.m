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

