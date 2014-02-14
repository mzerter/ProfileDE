function more = checkmore_genlin(more,n)

if isempty(more)
    more.sub = [kronecker(reshape(1:n,n,1),reshape(1,n,1)), ...
        kronecker(reshape(1,n,1),  reshape(1:n,n,1)),1:n^2];
    more.mat = zeros(n,n);
else
    if isfield(more, 'sub')
        if isempty(more.sub)
            more.sub = [kronecker(reshape(1:n,n,1),reshape(1,n,1)), ...
                kronecker(reshape(1,n,1),  reshape(1:n,n,1)),1:n^2];
        else
            if size(more.sub,2)==1 % can just specify which params to use
                more.sub = [kronecker(reshape(1:n,n,1),reshape(1,n,1)), ...
                    kronecker(reshape(1,n,1),  reshape(1:n,n,1)),more.sub];
            end
        end
    end
end

if ~isfield(more, 'mat') || isempty(more.mat)
    more.mat = zeros(n,n);
end

if isfield(more, 'force')
    if ~isempty(more.force)
        m = length(more.force);
        if isempty(more.force.mat)
            more.force.mat = zeros(n,m);
        end
        if isempty(more.force.sub)
            more.force.sub =[kronecker(reshape(1:n,m,1),reshape(1,m,1)), ...
                kronecker(reshape(1,n,1),reshape(1:m,n,1)),1:(n*m)];
        end
        if size(more.force.sub,2)==1 || is_vector(more.force.sub)
            more.force.sub = ...
                [more.force.sub,more.force.sub, ...
                size(more.sub,1)+length(more.force.sub)];
        end
    else  % one force per column, and specify which params
        if size(more.force.sub,2)==2
            more.force.sub = [more.force.sub(:,1), ...
                more.force.sub(:,1), ...
                more.force.sub(:,2)];
        end
        
    end
    
end

end
