function fnval = findif_ode_fn(times,y,p,more)

if ~isfield(more, 'more')
    more.more = [];
end
fnval = more.fn(times,y,p,more.more);

end

