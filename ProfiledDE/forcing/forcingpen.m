function r = forcingpen(coefs,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forcingpen
%
% penalizes a forcing funciton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

quads = getquadvals(more.basis{1});

forces = Make_fdcell(coefs,more.basis);

pens = eval_fdcell(quads(:,1),forces,more.deg);

lambda = more.lambda.*ones(prod(size(forces)),1);

for(i = 1:length(forces))
    pens{i} = diag(sqrt(more.lambda(i)*quads(:,2)))*pens{i};
end

r = cell2mat(pens);

end