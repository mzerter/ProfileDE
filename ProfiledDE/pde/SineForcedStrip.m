
%  Set up for a rectange with height 0.2, width 1, centered on 0

% Define the rectangle in pdetool:

pderect([0,1,-0.1,0.1]);

%  exported arrays in Draw step:

% gd =
%     3.0000
%     4.0000
%          0
%     1.0000
%     1.0000
%          0
%    -0.1000
%    -0.1000
%     0.1000
%     0.1000

% sf =
% R1

% ns =
%     82
%     49

%  Exported arrays in Boundary step:
%  Boundary segments numbered starting the bottom and
%  going counter-clockwise.  

% g =
%     2.0000    2.0000    2.0000    2.0000
%          0    1.0000    1.0000         0
%     1.0000    1.0000         0         0
%    -0.1000   -0.1000    0.1000    0.1000
%    -0.1000    0.1000    0.1000   -0.1000
%     1.0000    1.0000    1.0000    1.0000
%          0         0         0         0

% Conditions are:  Neumann, Dirichlet, Neumann, Neumann

% b =
%      1     1     1     1
%      0     1     0     0
%      1     1     1     1
%      1     1     1     1
%     48     1    48    48
%     48     1    48    48
%     48    48    48    48
%     48    48    48    48
%     49    49    49    49
%     48    48    48    48

%  Exported coefficients in PDE step:

% a =
% 1          

% c =
% 0.01       

% d =
% 1.0        

% Mesh after initmesh and one refinemesh:

% 159 points
% 52  edges
% 264 triangles

%  Solve parameters:  time = linspace(0,10,11)

%  size(u) = [159, 11]

