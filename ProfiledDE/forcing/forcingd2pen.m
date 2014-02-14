function r = forcingdpen(coefs,more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forcingd2pen
%
% second derivative of forcing penalty wrt coefs
%
% Note that unlike forcingpen and forcingdpen this does not output values
% at quadrature points, but rather provides a full Hessian. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


quads = getquadvals(more.basis{1});

pens = eval_basis_cell(quads(:,1),more.basis,more.deg);

for(i = 1:prod(size(pens)))
    pens{i} = more.lambda(i)*pens{i}'*diag(quads(:,2))*pens{i};
end

r = mattdiag_cell(pens,0);

end