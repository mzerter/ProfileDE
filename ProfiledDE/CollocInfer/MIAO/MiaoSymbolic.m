x1 = sym('x1');
x2 = sym('x2');
x3 = sym('x3');

p1 = sym('p1');
p2 = sym('p2');
p3 = sym('p3');
p4 = sym('p4');
p5 = sym('p5');

logten = sym('logten');

r1 = (p1 - p2*x3)/logten;
r2 = (p2*x1*x3 - p3*x2)/(logten*x2);
r3 = (p4*x2 - p5*x3)/(logten*x3);

% p-derivatives

%  r1

dr1dp1 = diff(r1,p1)  % 1/logten
dr1dp2 = diff(r1,p2)  % -x3/logten

%  r2

dr2dp2 = diff(r2,p2)  % (x1*x3)/(logten*x2)
dr2dp3 = diff(r2,p3)  % -1/logten

%  r3

dr3dp1 = diff(r3,p4)  %  x2/(logten*x3)
dr3dp2 = diff(r3,p5)  % -1/logten

%  x-derivatives

%  r1

dr1dx3 = diff(r1,x3)  % -p2/logten

%  r2

dr2dx1 = diff(r2,x1)  % (p2*x3)/(logten*x2)
dr2dx2 = simplify(diff(r2,x2))  % -(p2*x1*x3)/(logten*x2^2)
dr2dx3 = diff(r2,x3)  % (p2*x1)/(logten*x2)

%  r3

dr3dx2 = diff(r3,x2)  % p4/(logten*x3)
dr3dx3 = diff(r3,x3)  % -p5/(logten*x3) - (p4*x2 - p5*x3)/(logten*x3^2)

%  xx-derivatives

%  r2 x1

d2r2dx1x2 = diff(dr2dx1,x2)  % -(p2*x3)/(logten*x2^2)  
d2r2dx1x3 = diff(dr2dx1,x3)  % p2/(logten*x2)

%  r2 x2

d2r2dx2x1 = diff(dr2dx2,x1)  % -(p2*x3)/(logten*x2^2) 
d2r2dx2x2 = simplify(diff(dr2dx2,x2))  % (2*p2*x1*x3)/(logten*x2^3)
d2r2dx2x3 = diff(dr2dx2,x3)  % -(p2*x1)/(logten*x2^2) 

%  r2 x3

d2r2dx3x1 = diff(dr2dx3,x1)  % p2/(logten*x2)  
d2r2dx3x2 = diff(dr2dx3,x2)  % -(p2*x1)/(logten*x2^2)  

%  r3

d2r3dx2x3 = diff(dr3dx2,x3)  % -p4/(logten*x3^2)  

d2r3dx3x2 = diff(dr3dx3,x2)  % -p4/(logten*x3^2)  
d2r3dx3x3 = simplify(diff(dr3dx3,x3))  % (2*p4*x2)/(logten*x3^3)  

%  xp-derivatives

%  r1

d2r1dx3p2 = diff(dr1dx3,p2)  % -1/logten

%  r2

d2r2dx1p2 = diff(dr2dx1,p2)  % x3/(logten*x2)

d2r2dx2p2 = diff(dr2dx2,p2)  % -(x1*x3)/(logten*x2^2)
d2r2dx3p2 = diff(dr2dx3,p2)  % x1/(logten*x2)

%  r3

d2r3dx2p4 = diff(dr3dx2,p4)  % 1/(logten*x3)

d2r3dx3p4 = diff(dr3dx3,p4)  % -x2/(logten*x3^2)







