%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fnval = SSE(data, times, x, pars, more)

%  last modified 9 November 2013

if ~isfield(more, 'more')
    more.more = [];
end

%  Evaluate the fit to the data:
%    if more.fn is id.fn  in folder id,  x      is returned
%    if more.fn is exp.fn in folder exp, exp(x) is returned
%    other state variable transformations may be added 
%    in the same way

fdevals = more.fn(times, x, pars, more.more);

%  Define the values of the process at observation times

difs = data - fdevals;
difs(isnan(difs)) = 0;

if isfield(more,'which')
    whichobs = more.which;
else
    whichobs = [];
end

weights = checkweights(more.weights, whichobs, difs);

nvar = size(data,2);
fnval = 0;
for j=1:nvar
    indexj = ~isnan(data(:,j));
    if (indexj ~= 0)
        difsj = data(indexj,j) - fdevals(indexj,j);
%         disp(size(weights))
%         disp(size(difsj))
        fnval = fnval + sum(weights(:,j).*difsj.^2);
    end
end

