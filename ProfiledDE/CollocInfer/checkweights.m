function weights = checkweights(weights, whichrows, diffs)

[m,n] = size(diffs);

if isempty(whichrows)
    whichrows = 1:m;
end
    
if isempty(weights)
    weights = ones(m,n);
    return;
elseif numel(weights) == n
    weights = repmat(weights(:)',m,1);
    return;
elseif all(size(weights(whichrows,:)) == [m,1])
    weights = repmat(weights(:),1,n);
    return;
elseif all(size(weights(whichrows,:)) == [m,n])
    weights = weights(whichrows,:);
    return;
else
    error('Dimension of weights is not consistent with that of data.');
end
