function more = checkmore_cvar(more,n)
% checks additional arguments to cvar

if isempty(more)
    more.sub = (1:n)'*ones(1,3);
    % default to diagonal covariance
else
    if isempty(more.sub)
        more.sub = zeros(0,3);
    end
end
if isempty(more.mat)
    if size(more.sub,1) == 0
        more.mat = diag(ones(n,1));
    else
        more.mat = zeros(n,n);
    end
end

end
