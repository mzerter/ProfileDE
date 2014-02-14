function more = checkmore(more,n,npar)
% Checks the struct more for use with the genlin package

if isempty(more)
    error('Argument MORE is empty in a call to GENLIN function.');
else
    if ~isfield(more,'mat') || isempty(more.mat)
        error('MORE.MAT is either not present or empty.');
    end
    if ~isfield(more,'sub') || isempty(more.sub)
        error('MORE.SUB is either not present or empty.');
    end
    [Arow,Acol] = size(more.mat);
    nAsub = size(more.sub,1);
    for i=1:nAsub
        if more.sub(i,1) > Arow || ...
           more.sub(i,2) > Acol || ...
           more.sub(i,3) > npar
            error(['Subscripts in row ', num2str(i), ...
                   ' of more.sub are inconsistent ', ...
                   'with the size of more.mat.']);
        end
    end
    if ~isfield(more,'force'),  more.force = []; end
    if ~isempty(more.force)
        m = length(more.force);
        if ~isfield(more,'force_mat'),  more.force_mat = zeros(n,m);  end
        if ~isfield(more,'force_sub')
            more.force_sub = [kron((1:n)',ones(m,1)) ...
                              kron(ones(n,1),(1:m)')];
        end
        if size(more.force_sub,2) == 1
            more.force_sub = [more.force_sub more.force_sub];
        end
        nBsub = size(more.force.sub,1);
        for i=1:nBsub
            if more.force.sub(i,1) > Brow || ...
                    more.force.sub(i,2) > Bcol || ...
                    more.force.sub(i,3) > npar
                error(['Subscripts in row ', num2str(i), ...
                    ' of more.force.sub are inconsistent ', ...
                    'with the size of more.force.mat.']);
            end
        end
        if ~isfield(more,'force_input'),  more.force_input = [];  end
    end

end