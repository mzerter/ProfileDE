%  load boundary polygon

addpath('../fdaM')

load smthmtl

%   set up initial mesh

[p, e, t] = initmesh(dl);
pdemesh(p, e, t)  %  plot the mesh

%  refine the mesh   

p = jigglemesh(p, e, t);
pdemesh(p, e, t)  %  plot the mesh

[p, e, t] = refinemesh(dl, p, e, t);

%  define the parabolic PDE:

%   du       di^2 u
% d --  =  c ------ - a u + f
%   dt       di x^2

d = 1;  c = 1;  a = 0;   f = 0;

%  set up boundary condition  matrix  
%    see  function assemb for details

bcol = [1 0 1 1 '0' '0']';
b = bcol;
for i=2:sum(gd(2,:)), b = [b bcol]; end

%  set up the initial condition matrix

xp = p(1,:)';
yp = p(2,:)';
u0 = interp2(X, Y, Z, xp, yp);
pdeplot(p, e, t, 'xydata', u0, 'colormap', 'hot')
title('\fontsize{16} Income distribution')
  
%  set up the time values

tlist = [0, 10.^(-4:.5:-1)];
tlist = [0, 10.^(-5:0.5:-3)];

tlist = [0,1e-5];

%  solve the equation

u = parabolic(u0, tlist, b, p, e, t, c, a, f, d);

%  image plot

% without flow lines

for j=1:length(tlist)
   pdeplot(p, e, t, 'xydata', u(:,j),...
           'colormap', 'hot')
  title(['\fontsize{16} Wealth at Time ',num2str(tlist(j))])
  pause
end

%  with flow lines


%  The following code sets up a surface plot of the results
%    But I didn't find this at all helpful ... it's too 
%    complex to see what's going on.

%  exchange coordinates to get a better view for surface
%    plotting

[X,Y] = meshgrid(yi,-xi);  %  convert lattice to matrix form

%  interpolate data to grid
%  note:  parabolic is very sensitive to choice of method,
%         only v4 (Version 4) seems to work

method = 'v4';
Z = griddata(x, y, income, X, Y, method);  
Z = griddata(y,-x, inctot, X, Y, method);  
Z = griddata(x, y, incden, X, Y, method);  
Z = griddata(x, y, popden, X, Y, method);  

for j=1:3
  m = gd(2,j);
  temp = gd(3:(m+2),j);
  gd(3:(m+2),j) = gd((m+3):(2*m+2),j);
  gd((m+3):(2*m+2),j) = -temp;
end

dl = decsg(gd, sf, ns);
pdegplot(dl)    %  plot the geometry

%   set up initial mesh

[p, e, t] = initmesh(dl);
pdemesh(p, e, t)  %  plot the mesh

xp = p(1,:)';
yp = p(2,:)';
u0 = interp2(X, Y, Z, xp, yp);
pdeplot(p, e, t, 'xydata', u0, 'colormap', 'hot')
title('\fontsize{16} Initial Distribution Wealth')
u = parabolic(u0, tlist, b, p, e, t, c, a, f, d);

%  image plot

% without flow lines
for j=1:length(tlist)
   pdeplot(p, e, t, 'xydata', u(:,j),...
           'colormap', 'hot')
  title(['\fontsize{16} Wealth at Time ',num2str(tlist(j))])
  pause
end

%  with flow lines
for j=1:length(tlist)
   [ux,uy] = pdegrad(p, t, u(:,j));
   pdeplot(p, e, t, 'xydata', u(:,j),...
           'colormap', 'hot',...
           'flowdata', [-ux;-uy])
  title(['\fontsize{16} Wealth at Time ',num2str(tlist(j))])
  pause
end

%  The following code sets up a surface plot of the results
%    But I didn't find this at all helpful ... it's too 
%    complex to see what's going on.

%  exchange coordinates to get a better view for surface
%    plotting

[X,Y] = meshgrid(yi,-xi);  %  convert lattice to matrix form

%  interpolate data to grid
%  note:  parabolic is very sensitive to choice of method,
%         only v4 (Version 4) seems to work

method = 'v4';
Z = griddata(x, y, income, X, Y, method);  
Z = griddata(y,-x, inctot, X, Y, method);  
Z = griddata(x, y, incden, X, Y, method);  
Z = griddata(x, y, popden, X, Y, method);  

for j=1:3
  m = gd(2,j);
  temp = gd(3:(m+2),j);
  gd(3:(m+2),j) = gd((m+3):(2*m+2),j);
  gd((m+3):(2*m+2),j) = -temp;
end

dl = decsg(gd, sf, ns);
pdegplot(dl)    %  plot the geometry

%   set up initial mesh

[p, e, t] = initmesh(dl);
pdemesh(p, e, t)  %  plot the mesh

xp = p(1,:)';
yp = p(2,:)';
u0 = interp2(X, Y, Z, xp, yp);
pdeplot(p, e, t, 'xydata', u0, 'colormap', 'hot')
title('\fontsize{16} Initial Distribution Wealth')
u = parabolic(u0, tlist, b, p, e, t, c, a, f, d);

%  surface plot

for j=1:length(tlist)
   [ux,uy] = pdegrad(p, t, u(:,j));
   pdeplot(p, e, t, 'xydata', u(:,j), 'xystyle','interp',...
           'zdata',u(:,j),'zstyle','continuous',...
           'colormap', 'hot',...
           'colorbar','off',...
           'flowdata', [-ux;-uy])
  title(['\fontsize{16} Wealth at Time ',num2str(tlist(j))])
  pause
end





