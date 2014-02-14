function tmpsum = applysum(tmparray,i)
dims = size(tmparray);
ndim = length(dims);
if i <= ndim
    ilim = dims(i);
    tmpsum = zeros(1,ilim);
    for j=1:ilim
        tmpsum(j) = sum(reshape(tmparray,numel(tmparray),1));
    end
else
    error('i exceeds number of dimensions of array');
end