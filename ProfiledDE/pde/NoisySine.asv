%  smooth a noisy sine function using pdetool

addpath('../../fdaM')

xbasis = create_polygonal_basis([0,1],linspace(0,1,21));
ybasis = create_polygonal_basis([-0.1,0.1],linspace(-0.1,0.1,2));

xarg = linspace(0,1,21)';
yarg = linspace(-0.1,0.1,2)';

u0mat = sin(2*pi*xarg)*ones(1,2);
u0mat = u0mat + randn(21,2).*0.2;

noisysinebifd = smooth_bibasis(xarg, yarg, u0mat, xbasis, ybasis);

bifdmat = eval_bifd(xarg, yarg, bifdobj);

surf(bifdmat)

np = size(p,2);
u0 = zeros(np,1);
for i=1:np
    u0(i) = eval_bifd(p(1,i),p(2,i),bifdobj);
end


