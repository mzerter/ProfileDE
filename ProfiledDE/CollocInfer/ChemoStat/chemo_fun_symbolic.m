y1 = sym('y1');
y2 = sym('y2');
y3 = sym('y3');
y4 = sym('y4');
y5 = sym('y5');

p1  = sym('p1');
p2  = sym('p2');
p3  = sym('p3');
p4  = sym('p4');
p5  = sym('p5');
p6  = sym('p6');
p7  = sym('p7');
p8  = sym('p8');
p9  = sym('p9');
p10 = sym('p10');
p11 = sym('p11');
p12 = sym('p12');
p13 = sym('p13');
p14 = sym('p14');
p15 = sym('p15');
p16 = sym('p16');

f1 = sym('f1');
f2 = sym('f2');
f3 = sym('f3');
f4 = sym('f4');
f5 = sym('f5');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       define Q functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q  = p3*y2 + p4*y3;
Qs = Q/(p2 + exp(-p1*(Q - p16)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       define five functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = p6*(p5 - y1) - ...
         p12*y2*y1/(p10 + y1) - ...
         p12*y3*y1/(p11 + y1);
f2 = y2*(p9*p12*y1/(p10 + y1) - ...
         p3*p13*(y4 + y5)/(p15 + Qs) - p6);
f3 = y3*(p9*p12*y1/(p11 + y1) - ...
         p4*p13*(y4 + y5)/(p15 + Qs) - p6);
f4 = y4*(p14*p13*Q/(p15 + Qs) - ...
                 (p6 + p7 + p8));
f5 = p8*y4 - (p6 + p7)*y5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       define derivatives of five functions wrt five y's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y1

disp(simplify(diff(f1,'y1')))
disp(simplify(diff(f2,'y1')))
disp(simplify(diff(f3,'y1')))
disp(simplify(diff(f4,'y1')))
disp(simplify(diff(f5,'y1')))

% y2

disp(simplify(diff(f1,'y2')))
disp(simplify(diff(f2,'y2')))
disp(simplify(diff(f3,'y2')))
disp(simplify(diff(f4,'y2')))
disp(simplify(diff(f5,'y2')))

%y3

disp(simplify(diff(f1,'y3')))
disp(simplify(diff(f2,'y3')))
disp(simplify(diff(f3,'y3')))
disp(simplify(diff(f4,'y3')))
disp(simplify(diff(f5,'y3')))

%y4

disp(simplify(diff(f1,'y4')))
disp(simplify(diff(f2,'y4')))
disp(simplify(diff(f3,'y4')))
disp(simplify(diff(f4,'y4')))
disp(simplify(diff(f5,'y4')))

%y5

disp(simplify(diff(f1,'y5')))
disp(simplify(diff(f2,'y5')))
disp(simplify(diff(f3,'y5')))
disp(simplify(diff(f4,'y5')))
disp(simplify(diff(f5,'y5')))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    define second derivatives of five functions wrt five y's 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  y1-y1
disp(simplify(diff(diff(f1,'y1'),'y1')))
disp(simplify(diff(diff(f2,'y1'),'y1')))
disp(simplify(diff(diff(f3,'y1'),'y1')))
disp(simplify(diff(diff(f4,'y1'),'y1')))
disp(simplify(diff(diff(f5,'y1'),'y1')))

%  y2-y1
disp(simplify(diff(diff(f1,'y2'),'y1')))
disp(simplify(diff(diff(f2,'y2'),'y1')))
disp(simplify(diff(diff(f3,'y2'),'y1')))
disp(simplify(diff(diff(f4,'y2'),'y1')))
disp(simplify(diff(diff(f5,'y2'),'y1')))

% y3-y1
disp(simplify(diff(diff(f1,'y3'),'y1')))
disp(simplify(diff(diff(f2,'y3'),'y1')))
disp(simplify(diff(diff(f3,'y3'),'y1')))
disp(simplify(diff(diff(f4,'y3'),'y1')))
disp(simplify(diff(diff(f5,'y3'),'y1')))

% y4-y1
disp(simplify(diff(diff(f1,'y4'),'y1')))
disp(simplify(diff(diff(f2,'y4'),'y1')))
disp(simplify(diff(diff(f3,'y4'),'y1')))
disp(simplify(diff(diff(f4,'y4'),'y1')))
disp(simplify(diff(diff(f5,'y4'),'y1')))

% y5-y1
disp(simplify(diff(diff(f1,'y5'),'y1')))
disp(simplify(diff(diff(f2,'y5'),'y1')))
disp(simplify(diff(diff(f3,'y5'),'y1')))
disp(simplify(diff(diff(f4,'y5'),'y1')))
disp(simplify(diff(diff(f5,'y5'),'y1')))

%  y1-y2
disp(simplify(diff(diff(f1,'y1'),'y2')))
disp(simplify(diff(diff(f2,'y1'),'y2')))
disp(simplify(diff(diff(f3,'y1'),'y2')))
disp(simplify(diff(diff(f4,'y1'),'y2')))
disp(simplify(diff(diff(f5,'y1'),'y2')))

%  y2-y2
disp(simplify(diff(diff(f1,'y2'),'y2')))
disp(simplify(diff(diff(f2,'y2'),'y2')))
disp(simplify(diff(diff(f3,'y2'),'y2')))
disp(simplify(diff(diff(f4,'y2'),'y2')))
disp(simplify(diff(diff(f5,'y2'),'y2')))

%  y3-y2
disp(simplify(diff(diff(f1,'y3'),'y2')))
disp(simplify(diff(diff(f2,'y3'),'y2')))
disp(simplify(diff(diff(f3,'y3'),'y2')))
disp(simplify(diff(diff(f4,'y3'),'y2')))
disp(simplify(diff(diff(f5,'y3'),'y2')))

%  y4-y2
disp(simplify(diff(diff(f1,'y4'),'y2')))
disp(simplify(diff(diff(f2,'y4'),'y2')))
disp(simplify(diff(diff(f3,'y4'),'y2')))
disp(simplify(diff(diff(f4,'y4'),'y2')))
disp(simplify(diff(diff(f5,'y4'),'y2')))

%  y5-y2
disp(simplify(diff(diff(f1,'y5'),'y2')))
disp(simplify(diff(diff(f2,'y5'),'y2')))
disp(simplify(diff(diff(f3,'y5'),'y2')))
disp(simplify(diff(diff(f4,'y5'),'y2')))
disp(simplify(diff(diff(f5,'y5'),'y2')))

%  y1-y3
disp(simplify(diff(diff(f1,'y1'),'y3')))
disp(simplify(diff(diff(f2,'y1'),'y3')))
disp(simplify(diff(diff(f3,'y1'),'y3')))
disp(simplify(diff(diff(f4,'y1'),'y3')))
disp(simplify(diff(diff(f5,'y1'),'y3')))

%  y2-y3
disp(simplify(diff(diff(f1,'y2'),'y3')))
disp(simplify(diff(diff(f2,'y2'),'y3')))
disp(simplify(diff(diff(f3,'y2'),'y3')))
disp(simplify(diff(diff(f4,'y2'),'y3')))
disp(simplify(diff(diff(f5,'y2'),'y3')))

%  y3-y3
disp(simplify(diff(diff(f1,'y3'),'y3')))
disp(simplify(diff(diff(f2,'y3'),'y3')))
disp(simplify(diff(diff(f3,'y3'),'y3')))
disp(simplify(diff(diff(f4,'y3'),'y3')))
disp(simplify(diff(diff(f5,'y3'),'y3')))

%  y4-y3
disp(simplify(diff(diff(f1,'y4'),'y3')))
disp(simplify(diff(diff(f2,'y4'),'y3')))
disp(simplify(diff(diff(f3,'y4'),'y3')))
disp(simplify(diff(diff(f4,'y4'),'y3')))
disp(simplify(diff(diff(f5,'y4'),'y3')))

%  y5-y3
disp(simplify(diff(diff(f1,'y5'),'y3')))
disp(simplify(diff(diff(f2,'y5'),'y3')))
disp(simplify(diff(diff(f3,'y5'),'y3')))
disp(simplify(diff(diff(f4,'y5'),'y3')))
disp(simplify(diff(diff(f5,'y5'),'y3')))

%  y1-y4
disp(simplify(diff(diff(f1,'y1'),'y4')))
disp(simplify(diff(diff(f2,'y1'),'y4')))
disp(simplify(diff(diff(f3,'y1'),'y4')))
disp(simplify(diff(diff(f4,'y1'),'y4')))
disp(simplify(diff(diff(f5,'y1'),'y4')))

%  y2-y4
disp(simplify(diff(diff(f1,'y2'),'y4')))
disp(simplify(diff(diff(f2,'y2'),'y4')))
disp(simplify(diff(diff(f3,'y2'),'y4')))
disp(simplify(diff(diff(f4,'y2'),'y4')))
disp(simplify(diff(diff(f5,'y2'),'y4')))

%  y3-y4
disp(simplify(diff(diff(f1,'y3'),'y4')))
disp(simplify(diff(diff(f2,'y3'),'y4')))
disp(simplify(diff(diff(f3,'y3'),'y4')))
disp(simplify(diff(diff(f4,'y3'),'y4')))
disp(simplify(diff(diff(f5,'y3'),'y4')))

%  y4-y4
disp(simplify(diff(diff(f1,'y4'),'y4')))
disp(simplify(diff(diff(f2,'y4'),'y4')))
disp(simplify(diff(diff(f3,'y4'),'y4')))
disp(simplify(diff(diff(f4,'y4'),'y4')))
disp(simplify(diff(diff(f5,'y4'),'y4')))

%  y5-y4
disp(simplify(diff(diff(f1,'y5'),'y4')))
disp(simplify(diff(diff(f2,'y5'),'y4')))
disp(simplify(diff(diff(f3,'y5'),'y4')))
disp(simplify(diff(diff(f4,'y5'),'y4')))
disp(simplify(diff(diff(f5,'y5'),'y4')))

%  y1-y5
disp(simplify(diff(diff(f1,'y1'),'y5')))
disp(simplify(diff(diff(f2,'y1'),'y5')))
disp(simplify(diff(diff(f3,'y1'),'y5')))
disp(simplify(diff(diff(f4,'y1'),'y5')))
disp(simplify(diff(diff(f5,'y1'),'y5')))

%  y2-y5
disp(simplify(diff(diff(f1,'y2'),'y5')))
disp(simplify(diff(diff(f2,'y2'),'y5')))
disp(simplify(diff(diff(f3,'y2'),'y5')))
disp(simplify(diff(diff(f4,'y2'),'y5')))
disp(simplify(diff(diff(f5,'y2'),'y5')))

%  y3-y5
disp(simplify(diff(diff(f1,'y3'),'y5')))
disp(simplify(diff(diff(f2,'y3'),'y5')))
disp(simplify(diff(diff(f3,'y3'),'y5')))
disp(simplify(diff(diff(f4,'y3'),'y5')))
disp(simplify(diff(diff(f5,'y3'),'y5')))

%  y4-y5
disp(simplify(diff(diff(f1,'y4'),'y5')))
disp(simplify(diff(diff(f2,'y4'),'y5')))
disp(simplify(diff(diff(f3,'y4'),'y5')))
disp(simplify(diff(diff(f4,'y4'),'y5')))
disp(simplify(diff(diff(f5,'y4'),'y5')))

%  y5-y5
disp(simplify(diff(diff(f1,'y5'),'y5')))
disp(simplify(diff(diff(f2,'y5'),'y5')))
disp(simplify(diff(diff(f3,'y5'),'y5')))
disp(simplify(diff(diff(f4,'y5'),'y5')))
disp(simplify(diff(diff(f5,'y5'),'y5')))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   derivatives of 5 f's wrt last 14 p's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p3

disp(simplify(diff(f1,'p3')))
disp(simplify(diff(f2,'p3')))
disp(simplify(diff(f3,'p3')))
disp(simplify(diff(f4,'p3')))
disp(simplify(diff(f5,'p3')))

% p4

disp(simplify(diff(f1,'p4')))
disp(simplify(diff(f2,'p4')))
disp(simplify(diff(f3,'p4')))
disp(simplify(diff(f4,'p4')))
disp(simplify(diff(f5,'p4')))

% p5

disp(simplify(diff(f1,'p5')))
disp(simplify(diff(f2,'p5')))
disp(simplify(diff(f3,'p5')))
disp(simplify(diff(f4,'p5')))
disp(simplify(diff(f5,'p5')))

% p6

disp(simplify(diff(f1,'p6')))
disp(simplify(diff(f2,'p6')))
disp(simplify(diff(f3,'p6')))
disp(simplify(diff(f4,'p6')))
disp(simplify(diff(f5,'p6')))

% p7

disp(simplify(diff(f1,'p7')))
disp(simplify(diff(f2,'p7')))
disp(simplify(diff(f3,'p7')))
disp(simplify(diff(f4,'p7')))
disp(simplify(diff(f5,'p7')))

% p8

disp(simplify(diff(f1,'p8')))
disp(simplify(diff(f2,'p8')))
disp(simplify(diff(f3,'p8')))
disp(simplify(diff(f4,'p8')))
disp(simplify(diff(f5,'p8')))

% p9

disp(simplify(diff(f1,'p9')))
disp(simplify(diff(f2,'p9')))
disp(simplify(diff(f3,'p9')))
disp(simplify(diff(f4,'p9')))
disp(simplify(diff(f5,'p9')))

% p10

disp(simplify(diff(f1,'p10')))
disp(simplify(diff(f2,'p10')))
disp(simplify(diff(f3,'p10')))
disp(simplify(diff(f4,'p10')))
disp(simplify(diff(f5,'p10')))

% p11

disp(simplify(diff(f1,'p11')))
disp(simplify(diff(f2,'p11')))
disp(simplify(diff(f3,'p11')))
disp(simplify(diff(f4,'p11')))
disp(simplify(diff(f5,'p11')))

% p12

disp(simplify(diff(f1,'p12')))
disp(simplify(diff(f2,'p12')))
disp(simplify(diff(f3,'p12')))
disp(simplify(diff(f4,'p12')))
disp(simplify(diff(f5,'p12')))

% p13

disp(simplify(diff(f1,'p13')))
disp(simplify(diff(f2,'p13')))
disp(simplify(diff(f3,'p13')))
disp(simplify(diff(f4,'p13')))
disp(simplify(diff(f5,'p13')))

% p14

disp(simplify(diff(f1,'p14')))
disp(simplify(diff(f2,'p14')))
disp(simplify(diff(f3,'p14')))
disp(simplify(diff(f4,'p14')))
disp(simplify(diff(f5,'p14')))

% p15

disp(simplify(diff(f1,'p15')))
disp(simplify(diff(f2,'p15')))
disp(simplify(diff(f3,'p15')))
disp(simplify(diff(f4,'p15')))
disp(simplify(diff(f5,'p15')))

% p16 

disp(simplify(diff(f1,'p16')))
disp(simplify(diff(f2,'p16')))
disp(simplify(diff(f3,'p16')))
disp(simplify(diff(f4,'p16')))
disp(simplify(diff(f5,'p16')))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     second derivatives of f's wrt y's and wrt last 14 p's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  y1-p3
disp(simplify(diff(diff(f1,'y1'),'p3')))
disp(simplify(diff(diff(f2,'y1'),'p3')))
disp(simplify(diff(diff(f3,'y1'),'p3')))
disp(simplify(diff(diff(f4,'y1'),'p3')))
disp(simplify(diff(diff(f5,'y1'),'p3')))

%  y2-p3
disp(simplify(diff(diff(f1,'y2'),'p3')))
disp(simplify(diff(diff(f2,'y2'),'p3')))
disp(simplify(diff(diff(f3,'y2'),'p3')))
disp(simplify(diff(diff(f4,'y2'),'p3')))
disp(simplify(diff(diff(f5,'y2'),'p3')))

%  y3-p3
disp(simplify(diff(diff(f1,'y3'),'p3')))
disp(simplify(diff(diff(f2,'y3'),'p3')))
disp(simplify(diff(diff(f3,'y3'),'p3')))
disp(simplify(diff(diff(f4,'y3'),'p3')))
disp(simplify(diff(diff(f5,'y3'),'p3')))

%  y4-p3
disp(simplify(diff(diff(f1,'y4'),'p3')))
disp(simplify(diff(diff(f2,'y4'),'p3')))
disp(simplify(diff(diff(f3,'y4'),'p3')))
disp(simplify(diff(diff(f4,'y4'),'p3')))
disp(simplify(diff(diff(f5,'y4'),'p3')))

%  y5-p3
disp(simplify(diff(diff(f1,'y5'),'p3')))
disp(simplify(diff(diff(f2,'y5'),'p3')))
disp(simplify(diff(diff(f3,'y5'),'p3')))
disp(simplify(diff(diff(f4,'y5'),'p3')))
disp(simplify(diff(diff(f5,'y5'),'p3')))



%  y1-p4
disp(simplify(diff(diff(f1,'y1'),'p4')))
disp(simplify(diff(diff(f2,'y1'),'p4')))
disp(simplify(diff(diff(f3,'y1'),'p4')))
disp(simplify(diff(diff(f4,'y1'),'p4')))
disp(simplify(diff(diff(f5,'y1'),'p4')))

%  y2-p4
disp(simplify(diff(diff(f1,'y2'),'p4')))
disp(simplify(diff(diff(f2,'y2'),'p4')))
disp(simplify(diff(diff(f3,'y2'),'p4')))
disp(simplify(diff(diff(f4,'y2'),'p4')))
disp(simplify(diff(diff(f5,'y2'),'p4')))

%  y3-p4
disp(simplify(diff(diff(f1,'y3'),'p4')))
disp(simplify(diff(diff(f2,'y3'),'p4')))
disp(simplify(diff(diff(f3,'y3'),'p4')))
disp(simplify(diff(diff(f4,'y3'),'p4')))
disp(simplify(diff(diff(f5,'y3'),'p4')))

%  y4-p4
disp(simplify(diff(diff(f1,'y4'),'p4')))
disp(simplify(diff(diff(f2,'y4'),'p4')))
disp(simplify(diff(diff(f3,'y4'),'p4')))
disp(simplify(diff(diff(f4,'y4'),'p4')))
disp(simplify(diff(diff(f5,'y4'),'p4')))

%  y5-p4
disp(simplify(diff(diff(f1,'y5'),'p4')))
disp(simplify(diff(diff(f2,'y5'),'p4')))
disp(simplify(diff(diff(f3,'y5'),'p4')))
disp(simplify(diff(diff(f4,'y5'),'p4')))
disp(simplify(diff(diff(f5,'y5'),'p4')))



%  y1-p5
disp(simplify(diff(diff(f1,'y1'),'p5')))
disp(simplify(diff(diff(f2,'y1'),'p5')))
disp(simplify(diff(diff(f3,'y1'),'p5')))
disp(simplify(diff(diff(f4,'y1'),'p5')))
disp(simplify(diff(diff(f5,'y1'),'p5')))

%  y2-p5
disp(simplify(diff(diff(f1,'y2'),'p5')))
disp(simplify(diff(diff(f2,'y2'),'p5')))
disp(simplify(diff(diff(f3,'y2'),'p5')))
disp(simplify(diff(diff(f4,'y2'),'p5')))
disp(simplify(diff(diff(f5,'y2'),'p5')))

%  y3-p5
disp(simplify(diff(diff(f1,'y3'),'p5')))
disp(simplify(diff(diff(f2,'y3'),'p5')))
disp(simplify(diff(diff(f3,'y3'),'p5')))
disp(simplify(diff(diff(f4,'y3'),'p5')))
disp(simplify(diff(diff(f5,'y3'),'p5')))

%  y4-p5
disp(simplify(diff(diff(f1,'y4'),'p5')))
disp(simplify(diff(diff(f2,'y4'),'p5')))
disp(simplify(diff(diff(f3,'y4'),'p5')))
disp(simplify(diff(diff(f4,'y4'),'p5')))
disp(simplify(diff(diff(f5,'y4'),'p5')))

%  y5-p5
disp(simplify(diff(diff(f1,'y5'),'p5')))
disp(simplify(diff(diff(f2,'y5'),'p5')))
disp(simplify(diff(diff(f3,'y5'),'p5')))
disp(simplify(diff(diff(f4,'y5'),'p5')))
disp(simplify(diff(diff(f5,'y5'),'p5')))



%  y1-p6
disp(simplify(diff(diff(f1,'y1'),'p6')))
disp(simplify(diff(diff(f2,'y1'),'p6')))
disp(simplify(diff(diff(f3,'y1'),'p6')))
disp(simplify(diff(diff(f4,'y1'),'p6')))
disp(simplify(diff(diff(f5,'y1'),'p6')))

%  y2-p6
disp(simplify(diff(diff(f1,'y2'),'p6')))
disp(simplify(diff(diff(f2,'y2'),'p6')))
disp(simplify(diff(diff(f3,'y2'),'p6')))
disp(simplify(diff(diff(f4,'y2'),'p6')))
disp(simplify(diff(diff(f5,'y2'),'p6')))

%  y3-p6
disp(simplify(diff(diff(f1,'y3'),'p6')))
disp(simplify(diff(diff(f2,'y3'),'p6')))
disp(simplify(diff(diff(f3,'y3'),'p6')))
disp(simplify(diff(diff(f4,'y3'),'p6')))
disp(simplify(diff(diff(f5,'y3'),'p6')))

%  y4-p6
disp(simplify(diff(diff(f1,'y4'),'p6')))
disp(simplify(diff(diff(f2,'y4'),'p6')))
disp(simplify(diff(diff(f3,'y4'),'p6')))
disp(simplify(diff(diff(f4,'y4'),'p6')))
disp(simplify(diff(diff(f5,'y4'),'p6')))

%  y5-p6
disp(simplify(diff(diff(f1,'y5'),'p6')))
disp(simplify(diff(diff(f2,'y5'),'p6')))
disp(simplify(diff(diff(f3,'y5'),'p6')))
disp(simplify(diff(diff(f4,'y5'),'p6')))
disp(simplify(diff(diff(f5,'y5'),'p6')))



%  y1-p7
disp(simplify(diff(diff(f1,'y1'),'p7')))
disp(simplify(diff(diff(f2,'y1'),'p7')))
disp(simplify(diff(diff(f3,'y1'),'p7')))
disp(simplify(diff(diff(f4,'y1'),'p7')))
disp(simplify(diff(diff(f5,'y1'),'p7')))

%  y2-p7
disp(simplify(diff(diff(f1,'y2'),'p7')))
disp(simplify(diff(diff(f2,'y2'),'p7')))
disp(simplify(diff(diff(f3,'y2'),'p7')))
disp(simplify(diff(diff(f4,'y2'),'p7')))
disp(simplify(diff(diff(f5,'y2'),'p7')))

%  y3-p7
disp(simplify(diff(diff(f1,'y3'),'p7')))
disp(simplify(diff(diff(f2,'y3'),'p7')))
disp(simplify(diff(diff(f3,'y3'),'p7')))
disp(simplify(diff(diff(f4,'y3'),'p7')))
disp(simplify(diff(diff(f5,'y3'),'p7')))

%  y4-p7
disp(simplify(diff(diff(f1,'y4'),'p7')))
disp(simplify(diff(diff(f2,'y4'),'p7')))
disp(simplify(diff(diff(f3,'y4'),'p7')))
disp(simplify(diff(diff(f4,'y4'),'p7')))
disp(simplify(diff(diff(f5,'y4'),'p7')))

%  y5-p7
disp(simplify(diff(diff(f1,'y5'),'p7')))
disp(simplify(diff(diff(f2,'y5'),'p7')))
disp(simplify(diff(diff(f3,'y5'),'p7')))
disp(simplify(diff(diff(f4,'y5'),'p7')))
disp(simplify(diff(diff(f5,'y5'),'p7')))



%  y1-p8
disp(simplify(diff(diff(f1,'y1'),'p8')))
disp(simplify(diff(diff(f2,'y1'),'p8')))
disp(simplify(diff(diff(f3,'y1'),'p8')))
disp(simplify(diff(diff(f4,'y1'),'p8')))
disp(simplify(diff(diff(f5,'y1'),'p8')))

%  y2-p8
disp(simplify(diff(diff(f1,'y2'),'p8')))
disp(simplify(diff(diff(f2,'y2'),'p8')))
disp(simplify(diff(diff(f3,'y2'),'p8')))
disp(simplify(diff(diff(f4,'y2'),'p8')))
disp(simplify(diff(diff(f5,'y2'),'p8')))

%  y3-p8
disp(simplify(diff(diff(f1,'y3'),'p8')))
disp(simplify(diff(diff(f2,'y3'),'p8')))
disp(simplify(diff(diff(f3,'y3'),'p8')))
disp(simplify(diff(diff(f4,'y3'),'p8')))
disp(simplify(diff(diff(f5,'y3'),'p8')))

%  y4-p8
disp(simplify(diff(diff(f1,'y4'),'p8')))
disp(simplify(diff(diff(f2,'y4'),'p8')))
disp(simplify(diff(diff(f3,'y4'),'p8')))
disp(simplify(diff(diff(f4,'y4'),'p8')))
disp(simplify(diff(diff(f5,'y4'),'p8')))

%  y5-p8
disp(simplify(diff(diff(f1,'y5'),'p8')))
disp(simplify(diff(diff(f2,'y5'),'p8')))
disp(simplify(diff(diff(f3,'y5'),'p8')))
disp(simplify(diff(diff(f4,'y5'),'p8')))
disp(simplify(diff(diff(f5,'y5'),'p8')))



%  y1-p9
disp(simplify(diff(diff(f1,'y1'),'p9')))
disp(simplify(diff(diff(f2,'y1'),'p9')))
disp(simplify(diff(diff(f3,'y1'),'p9')))
disp(simplify(diff(diff(f4,'y1'),'p9')))
disp(simplify(diff(diff(f5,'y1'),'p9')))

%  y2-p9
disp(simplify(diff(diff(f1,'y2'),'p9')))
disp(simplify(diff(diff(f2,'y2'),'p9')))
disp(simplify(diff(diff(f3,'y2'),'p9')))
disp(simplify(diff(diff(f4,'y2'),'p9')))
disp(simplify(diff(diff(f5,'y2'),'p9')))

%  y3-p9
disp(simplify(diff(diff(f1,'y3'),'p9')))
disp(simplify(diff(diff(f2,'y3'),'p9')))
disp(simplify(diff(diff(f3,'y3'),'p9')))
disp(simplify(diff(diff(f4,'y3'),'p9')))
disp(simplify(diff(diff(f5,'y3'),'p9')))

%  y4-p9
disp(simplify(diff(diff(f1,'y4'),'p9')))
disp(simplify(diff(diff(f2,'y4'),'p9')))
disp(simplify(diff(diff(f3,'y4'),'p9')))
disp(simplify(diff(diff(f4,'y4'),'p9')))
disp(simplify(diff(diff(f5,'y4'),'p9')))

%  y5-p9
disp(simplify(diff(diff(f1,'y5'),'p9')))
disp(simplify(diff(diff(f2,'y5'),'p9')))
disp(simplify(diff(diff(f3,'y5'),'p9')))
disp(simplify(diff(diff(f4,'y5'),'p9')))
disp(simplify(diff(diff(f5,'y5'),'p9')))



%  y1-p10
disp(simplify(diff(diff(f1,'y1'),'p10')))
disp(simplify(diff(diff(f2,'y1'),'p10')))
disp(simplify(diff(diff(f3,'y1'),'p10')))
disp(simplify(diff(diff(f4,'y1'),'p10')))
disp(simplify(diff(diff(f5,'y1'),'p10')))

%  y2-p10
disp(simplify(diff(diff(f1,'y2'),'p10')))
disp(simplify(diff(diff(f2,'y2'),'p10')))
disp(simplify(diff(diff(f3,'y2'),'p10')))
disp(simplify(diff(diff(f4,'y2'),'p10')))
disp(simplify(diff(diff(f5,'y2'),'p10')))

%  y3-p10
disp(simplify(diff(diff(f1,'y3'),'p10')))
disp(simplify(diff(diff(f2,'y3'),'p10')))
disp(simplify(diff(diff(f3,'y3'),'p10')))
disp(simplify(diff(diff(f4,'y3'),'p10')))
disp(simplify(diff(diff(f5,'y3'),'p10')))

%  y4-p10
disp(simplify(diff(diff(f1,'y4'),'p10')))
disp(simplify(diff(diff(f2,'y4'),'p10')))
disp(simplify(diff(diff(f3,'y4'),'p10')))
disp(simplify(diff(diff(f4,'y4'),'p10')))
disp(simplify(diff(diff(f5,'y4'),'p10')))

%  y5-p10
disp(simplify(diff(diff(f1,'y5'),'p10')))
disp(simplify(diff(diff(f2,'y5'),'p10')))
disp(simplify(diff(diff(f3,'y5'),'p10')))
disp(simplify(diff(diff(f4,'y5'),'p10')))
disp(simplify(diff(diff(f5,'y5'),'p10')))



%  y1-p11
disp(simplify(diff(diff(f1,'y1'),'p11')))
disp(simplify(diff(diff(f2,'y1'),'p11')))
disp(simplify(diff(diff(f3,'y1'),'p11')))
disp(simplify(diff(diff(f4,'y1'),'p11')))
disp(simplify(diff(diff(f5,'y1'),'p11')))

%  y2-p11
disp(simplify(diff(diff(f1,'y2'),'p11')))
disp(simplify(diff(diff(f2,'y2'),'p11')))
disp(simplify(diff(diff(f3,'y2'),'p11')))
disp(simplify(diff(diff(f4,'y2'),'p11')))
disp(simplify(diff(diff(f5,'y2'),'p11')))

%  y3-p11
disp(simplify(diff(diff(f1,'y3'),'p11')))
disp(simplify(diff(diff(f2,'y3'),'p11')))
disp(simplify(diff(diff(f3,'y3'),'p11')))
disp(simplify(diff(diff(f4,'y3'),'p11')))
disp(simplify(diff(diff(f5,'y3'),'p11')))

%  y4-p11
disp(simplify(diff(diff(f1,'y4'),'p11')))
disp(simplify(diff(diff(f2,'y4'),'p11')))
disp(simplify(diff(diff(f3,'y4'),'p11')))
disp(simplify(diff(diff(f4,'y4'),'p11')))
disp(simplify(diff(diff(f5,'y4'),'p11')))

%  y5-p11
disp(simplify(diff(diff(f1,'y5'),'p11')))
disp(simplify(diff(diff(f2,'y5'),'p11')))
disp(simplify(diff(diff(f3,'y5'),'p11')))
disp(simplify(diff(diff(f4,'y5'),'p11')))
disp(simplify(diff(diff(f5,'y5'),'p11')))



%  y1-p12
disp(simplify(diff(diff(f1,'y1'),'p12')))
disp(simplify(diff(diff(f2,'y1'),'p12')))
disp(simplify(diff(diff(f3,'y1'),'p12')))
disp(simplify(diff(diff(f4,'y1'),'p12')))
disp(simplify(diff(diff(f5,'y1'),'p12')))

%  y2-p12
disp(simplify(diff(diff(f1,'y2'),'p12')))
disp(simplify(diff(diff(f2,'y2'),'p12')))
disp(simplify(diff(diff(f3,'y2'),'p12')))
disp(simplify(diff(diff(f4,'y2'),'p12')))
disp(simplify(diff(diff(f5,'y2'),'p12')))

%  y3-p12
disp(simplify(diff(diff(f1,'y3'),'p12')))
disp(simplify(diff(diff(f2,'y3'),'p12')))
disp(simplify(diff(diff(f3,'y3'),'p12')))
disp(simplify(diff(diff(f4,'y3'),'p12')))
disp(simplify(diff(diff(f5,'y3'),'p12')))

%  y4-p12
disp(simplify(diff(diff(f1,'y4'),'p12')))
disp(simplify(diff(diff(f2,'y4'),'p12')))
disp(simplify(diff(diff(f3,'y4'),'p12')))
disp(simplify(diff(diff(f4,'y4'),'p12')))
disp(simplify(diff(diff(f5,'y4'),'p12')))

%  y5-p12
disp(simplify(diff(diff(f1,'y5'),'p12')))
disp(simplify(diff(diff(f2,'y5'),'p12')))
disp(simplify(diff(diff(f3,'y5'),'p12')))
disp(simplify(diff(diff(f4,'y5'),'p12')))
disp(simplify(diff(diff(f5,'y5'),'p12')))



%  y1-p13
disp(simplify(diff(diff(f1,'y1'),'p13')))
disp(simplify(diff(diff(f2,'y1'),'p13')))
disp(simplify(diff(diff(f3,'y1'),'p13')))
disp(simplify(diff(diff(f4,'y1'),'p13')))
disp(simplify(diff(diff(f5,'y1'),'p13')))

%  y2-p13
disp(simplify(diff(diff(f1,'y2'),'p13')))
disp(simplify(diff(diff(f2,'y2'),'p13')))
disp(simplify(diff(diff(f3,'y2'),'p13')))
disp(simplify(diff(diff(f4,'y2'),'p13')))
disp(simplify(diff(diff(f5,'y2'),'p13')))

%  y3-p13
disp(simplify(diff(diff(f1,'y3'),'p13')))
disp(simplify(diff(diff(f2,'y3'),'p13')))
disp(simplify(diff(diff(f3,'y3'),'p13')))
disp(simplify(diff(diff(f4,'y3'),'p13')))
disp(simplify(diff(diff(f5,'y3'),'p13')))

%  y4-p13
disp(simplify(diff(diff(f1,'y4'),'p13')))
disp(simplify(diff(diff(f2,'y4'),'p13')))
disp(simplify(diff(diff(f3,'y4'),'p13')))
disp(simplify(diff(diff(f4,'y4'),'p13')))
disp(simplify(diff(diff(f5,'y4'),'p13')))

%  y5-p13
disp(simplify(diff(diff(f1,'y5'),'p13')))
disp(simplify(diff(diff(f2,'y5'),'p13')))
disp(simplify(diff(diff(f3,'y5'),'p13')))
disp(simplify(diff(diff(f4,'y5'),'p13')))
disp(simplify(diff(diff(f5,'y5'),'p13')))



%  y1-p14
disp(simplify(diff(diff(f1,'y1'),'p14')))
disp(simplify(diff(diff(f2,'y1'),'p14')))
disp(simplify(diff(diff(f3,'y1'),'p14')))
disp(simplify(diff(diff(f4,'y1'),'p14')))
disp(simplify(diff(diff(f5,'y1'),'p14')))

%  y2-p14
disp(simplify(diff(diff(f1,'y2'),'p14')))
disp(simplify(diff(diff(f2,'y2'),'p14')))
disp(simplify(diff(diff(f3,'y2'),'p14')))
disp(simplify(diff(diff(f4,'y2'),'p14')))
disp(simplify(diff(diff(f5,'y2'),'p14')))

%  y3-p14
disp(simplify(diff(diff(f1,'y3'),'p14')))
disp(simplify(diff(diff(f2,'y3'),'p14')))
disp(simplify(diff(diff(f3,'y3'),'p14')))
disp(simplify(diff(diff(f4,'y3'),'p14')))
disp(simplify(diff(diff(f5,'y3'),'p14')))

%  y4-p14
disp(simplify(diff(diff(f1,'y4'),'p14')))
disp(simplify(diff(diff(f2,'y4'),'p14')))
disp(simplify(diff(diff(f3,'y4'),'p14')))
disp(simplify(diff(diff(f4,'y4'),'p14')))
disp(simplify(diff(diff(f5,'y4'),'p14')))

%  y5-p14
disp(simplify(diff(diff(f1,'y5'),'p14')))
disp(simplify(diff(diff(f2,'y5'),'p14')))
disp(simplify(diff(diff(f3,'y5'),'p14')))
disp(simplify(diff(diff(f4,'y5'),'p14')))
disp(simplify(diff(diff(f5,'y5'),'p14')))



%  y1-p15
disp(simplify(diff(diff(f1,'y1'),'p15')))
disp(simplify(diff(diff(f2,'y1'),'p15')))
disp(simplify(diff(diff(f3,'y1'),'p15')))
disp(simplify(diff(diff(f4,'y1'),'p15')))
disp(simplify(diff(diff(f5,'y1'),'p15')))

%  y2-p15
disp(simplify(diff(diff(f1,'y2'),'p15')))
disp(simplify(diff(diff(f2,'y2'),'p15')))
disp(simplify(diff(diff(f3,'y2'),'p15')))
disp(simplify(diff(diff(f4,'y2'),'p15')))
disp(simplify(diff(diff(f5,'y2'),'p15')))

%  y3-p15
disp(simplify(diff(diff(f1,'y3'),'p15')))
disp(simplify(diff(diff(f2,'y3'),'p15')))
disp(simplify(diff(diff(f3,'y3'),'p15')))
disp(simplify(diff(diff(f4,'y3'),'p15')))
disp(simplify(diff(diff(f5,'y3'),'p15')))

%  y4-p15
disp(simplify(diff(diff(f1,'y4'),'p15')))
disp(simplify(diff(diff(f2,'y4'),'p15')))
disp(simplify(diff(diff(f3,'y4'),'p15')))
disp(simplify(diff(diff(f4,'y4'),'p15')))
disp(simplify(diff(diff(f5,'y4'),'p15')))

%  y5-p15
disp(simplify(diff(diff(f1,'y5'),'p15')))
disp(simplify(diff(diff(f2,'y5'),'p15')))
disp(simplify(diff(diff(f3,'y5'),'p15')))
disp(simplify(diff(diff(f4,'y5'),'p15')))
disp(simplify(diff(diff(f5,'y5'),'p15')))



%  y1-p16
disp(simplify(diff(diff(f1,'y1'),'p16')))
disp(simplify(diff(diff(f2,'y1'),'p16')))
disp(simplify(diff(diff(f3,'y1'),'p16')))
disp(simplify(diff(diff(f4,'y1'),'p16')))
disp(simplify(diff(diff(f5,'y1'),'p16')))

%  y2-p16
disp(simplify(diff(diff(f1,'y2'),'p16')))
disp(simplify(diff(diff(f2,'y2'),'p16')))
disp(simplify(diff(diff(f3,'y2'),'p16')))
disp(simplify(diff(diff(f4,'y2'),'p16')))
disp(simplify(diff(diff(f5,'y2'),'p16')))

%  y3-p16
disp(simplify(diff(diff(f1,'y3'),'p16')))
disp(simplify(diff(diff(f2,'y3'),'p16')))
disp(simplify(diff(diff(f3,'y3'),'p16')))
disp(simplify(diff(diff(f4,'y3'),'p16')))
disp(simplify(diff(diff(f5,'y3'),'p16')))

%  y4-p16
disp(simplify(diff(diff(f1,'y4'),'p16')))
disp(simplify(diff(diff(f2,'y4'),'p16')))
disp(simplify(diff(diff(f3,'y4'),'p16')))
disp(simplify(diff(diff(f4,'y4'),'p16')))
disp(simplify(diff(diff(f5,'y4'),'p16')))

%  y5-p16
disp(simplify(diff(diff(f1,'y5'),'p16')))
disp(simplify(diff(diff(f2,'y5'),'p16')))
disp(simplify(diff(diff(f3,'y5'),'p16')))
disp(simplify(diff(diff(f4,'y5'),'p16')))
disp(simplify(diff(diff(f5,'y5'),'p16')))




