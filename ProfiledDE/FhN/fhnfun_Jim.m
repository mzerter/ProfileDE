function Dy = fhnfun_Jim(t,y,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The FitzHugh-Nagumo equations in vector form
%  re-parameterized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = y(1,:);
R = y(2,:);
a = p(1);  %  phase angle
b = p(2);  %  recovery reaction speed per voltage  reaction speed
c = p(3);  %  voltage  reaction speed
Dy(1,:) = c*(V - V.^3/3 + R);
Dy(2,:) = -b*(sin(a)*R - cos(a) + V/b)/c;

