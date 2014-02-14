
load montisland.txt

bdpol = montisland;
npoly = size(bdpol,1);
bdpolmn = mean(montisland);
bdpol = bdpol - ones(npoly,1)*bdpolmn;
indpol = 100:100:npoly;

%  plot the outline of the island

plot(bdpol(:,1),bdpol(:,2))
for i=indpol
   text(bdpol(i,1),bdpol(i,2),num2str(i/100))
end

%  eliminate the Lachine canal

indpol = (1:npoly);
indpol = indpol(indpol< 879 | indpol > 977);
bdpol = bdpol(indpol,:);
npoly = size(bdpol,1);

indpol = 50:50:npoly;
plot(bdpol(:,1),bdpol(:,2))
for i=indpol
   text(bdpol(i,1),bdpol(i,2),num2str(i/50))
end

% eliminate Rene Levesque park in Lachine

indpol = (1:npoly);
indpol = indpol(indpol<= 832 | indpol >= 863);
bdpol = bdpol(indpol,:);
npoly = size(bdpol,1);

%  smooth boundary using a Fourier series decomposition

addpath('c:/Program Files/matlab71/fdaM')

circ = (0:(npoly-1))./(npoly-1);
period = 1;
nbasis = 101;
basisobj = create_fourier_basis([0,1], nbasis, period);

bdpolfd = data2fd(bdpol,circ,basisobj);

plot(bdpolfd)
hold on
plot(circ, bdpol(:,1), 'b-')
plot(circ, bdpol(:,2), 'g-')
hold off

bdpolmat = eval_fd(circ, bdpolfd);
plot(bdpolmat(:,1),bdpolmat(:,2))

%  get a subsample of these points

nbd = 100;  %  number of distinct points on outer boundary
circcrse = linspace(0,1,nbd+1)';
bdpolmat = eval_fd(circcrse, bdpolfd);
plot(bdpolmat(:,1),bdpolmat(:,2), 'o-')

bdlo = min(bdpolmat);
bdhi = max(bdpolmat);

save smthmtl

%  load data

load MtlCens.txt

n = size(MtlCens,1);            %  number of points
x = MtlCens(:,1) - bdpolmn(1);  % centered longitude
y = MtlCens(:,2) - bdpolmn(2);  % centered latitude
income = MtlCens(:,4);   %  average income
inctot = MtlCens(:,4).*MtlCens(:,3)./mean(MtlCens(:,3));
incden = MtlCens(:,4)./MtlCens(:,5);  % income density
popden = MtlCens(:,3)./MtlCens(:,5);  % pop. density

%  get rid of points outside of bounding rectangle

inside = (1:n);
inside = inside(x(inside) >= bdlo(1) & x(inside) <= bdhi(1));
inside = inside(y(inside) >= bdlo(2) & y(inside) <= bdhi(2));
n = length(inside);
x = x(inside);
y = y(inside);

%  eliminate points outside of island by manual identification
%    as well as certain points in Dorval airport and 
%    surrounding industrial area, and in a region in the
%    extreme north-east sector

plot(bdpolmat(:,1),bdpolmat(:,2), '.-', x, y, 'go')
[xout,yout] = ginput;
inside = 1:n;
m = length(xout);
for j=1:m
   xtmp = x(inside);
   ytmp = y(inside);
   indj  = (abs(xtmp-xout(j)) < 5e-3) & ...
           (abs(ytmp-yout(j)) < 5e-3);
   inside = inside(~indj);
end
plot(bdpolmat(:,1),bdpolmat(:,2), 'o-',...
     x(inside), y(inside), 'go')
[xout,yout] = ginput;

%  set up the final set of coordinates for data points

n = length(inside);
x = x(inside);
y = y(inside);

%  define dependent variables for admissible points

income = income(inside);
inctot = inctot(inside);
incden = incden(inside);
popden = popden(inside);

%  define two regions where nobody can go

%  Dorval airport

plot(bdpolmat(:,1),bdpolmat(:,2), '.-',x, y, 'go')
[xout,yout] = ginput;
hold on
plot([xout;xout(1,1)],[yout;yout(1,1)],'ro-')
hold off

m = length(xout);
dorval = [2;m;xout;yout];

%  Water treatment and oil refinery area

plot(bdpolmat(:,1),bdpolmat(:,2), '.-',x, y, 'go')
[xout,yout] = ginput;
hold on
plot([xout;xout(1,1)],[yout;yout(1,1)],'ro-')
hold off

m = length(xout);
eastend = [2;m;xout;yout];

%  plot the boundaries

plot(bdpolmat(:,1),bdpolmat(:,2), '.-')
hold on
m = dorval(2);
plot([dorval(   3: (m+2)  );dorval(  3)], ...
     [dorval((m+3):(2*m+2));dorval(m+3)], 'b.-')
m = eastend(2);
plot([eastend(   3: (m+2)  );eastend(  3)], ...
     [eastend((m+3):(2*m+2));eastend(m+3)], 'b.-')
hold off

%  plot the sampling points with interior boundaries

plot(bdpolmat(:,1),bdpolmat(:,2), '.-', x, y, 'go')
hold on
m = dorval(2);
plot([dorval(   3: (m+2)  );dorval(  3)], ...
     [dorval((m+3):(2*m+2));dorval(m+3)], 'b.-')
m = eastend(2);
plot([eastend(   3: (m+2)  );eastend(  3)], ...
     [eastend((m+3):(2*m+2));eastend(m+3)], 'b.-')
hold off

%  set up the rectangular sampling grid 

ngrid = 201;
xi = linspace(bdlo(1),bdhi(1),ngrid)';
yi = linspace(bdlo(2),bdhi(2),ngrid)';

[X,Y] = meshgrid(xi,yi);  %  convert lattice to matrix form

%  interpolate data to grid
%  note:  parabolic is very sensitive to choice of method,
%         only v4 (Version 4) seems to work

method = 'v4';
Z = griddata(x, y, income, X, Y, method);  
% Z = griddata(x, y, inctot, X, Y, method);  
% Z = griddata(x, y, incden, X, Y, method);  
% Z = griddata(x, y, popden, X, Y, method);  

mesh(X,Y,Z)

%  set up geometry description matrix
%    see function decsg, 5-29 in manual

%  first column is exterior boundary

gdmtl = [2; nbd; bdpolmat(1:nbd,1); bdpolmat(1:nbd,2)];
ngdmtl = size(gdmtl,1);
gd = zeros(ngdmtl,3);
gd(:,1) = gdmtl;

%  add two interior polygons

gd(1:size(dorval, 1),2) = dorval;
gd(1:size(eastend,1),3) = eastend;

%  name space matrix
ns = ['Montreal'; 'Dorval  '; 'Eastend ']';

%  set formula
sf = 'Montreal-(Dorval+Eastend)';

%  set up decomposed geometry

dl = decsg(gd, sf, ns);
pdegplot(dl)    %  plot the geometry

save smthmtl
