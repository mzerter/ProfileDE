function f = SineForcedStrip_f(p,t,u,time)

%  compute values of forcing function at barycenters of triangles

t1ind = t(1,:);
t2ind = t(2,:);
t3ind = t(3,:);

x = (p(1,t1ind) + p(1,t2ind) + p(1,t3ind))./3;
y = (p(2,t1ind) + p(2,t2ind) + p(2,t3ind))./3;

% f = x.*exp(-x.^2/2);
% f = sin(2.*pi.*x);
nt = size(t,2);
f = zeros(1,nt);
ind = 0.4 <= x <= 0.6;
f(ind) = 1;
