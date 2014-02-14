function [newdata, dims] = dataformat(data, coefs)

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

