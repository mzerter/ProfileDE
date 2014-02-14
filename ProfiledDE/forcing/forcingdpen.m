function r = forcingdpen(coefs,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forcingdpen
%
% derivative of forcing penalty wrt coefs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


quads = getquadvals(more.basis{1});

pens = eval_basis_cell(quads(:,1),more.basis,more.deg);

lambda = more.lambda.*ones(prod(size(pens)),1);

for(i = 1:prod(size(pens)))
    pens{i} = diag(sqrt(more.lambda(i)*quads(:,2)))*pens{i};
end

r = mattdiag_cell(pens,0);

end