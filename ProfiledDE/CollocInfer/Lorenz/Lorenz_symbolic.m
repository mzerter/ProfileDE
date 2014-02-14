y1 = sym('y1');
y2 = sym('y2');
y3 = sym('y3');

p1  = sym('p1');
p2  = sym('p2');
p1  = sym('p3');

f1 = sym('f1');
f2 = sym('f2');
f3 = sym('f3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   define three functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = -p1*y1 + p1*y2;
f2 =  p2*y1 - y2 + y1*y3;
f3 = -p1*y3 + y1*y2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       define derivatives of three functions wrt three y's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y1

disp(simplify(diff(f1,'y1')))
disp(simplify(diff(f2,'y1')))
disp(simplify(diff(f3,'y1')))

% y2

disp(simplify(diff(f1,'y2')))
disp(simplify(diff(f2,'y2')))
disp(simplify(diff(f3,'y2')))

%y3

disp(simplify(diff(f1,'y3')))
disp(simplify(diff(f2,'y3')))
disp(simplify(diff(f3,'y3')))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    define second derivatives of three functions wrt three y's 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  y1-y1
disp(simplify(diff(diff(f1,'y1'),'y1')))
disp(simplify(diff(diff(f2,'y1'),'y1')))
disp(simplify(diff(diff(f3,'y1'),'y1')))

%  y2-y1
disp(simplify(diff(diff(f1,'y2'),'y1')))
disp(simplify(diff(diff(f2,'y2'),'y1')))
disp(simplify(diff(diff(f3,'y2'),'y1')))

% y3-y1
disp(simplify(diff(diff(f1,'y3'),'y1')))
disp(simplify(diff(diff(f2,'y3'),'y1')))
disp(simplify(diff(diff(f3,'y3'),'y1')))


%  y1-y2
disp(simplify(diff(diff(f1,'y1'),'y2')))
disp(simplify(diff(diff(f2,'y1'),'y2')))
disp(simplify(diff(diff(f3,'y1'),'y2')))

%  y2-y2
disp(simplify(diff(diff(f1,'y2'),'y2')))
disp(simplify(diff(diff(f2,'y2'),'y2')))
disp(simplify(diff(diff(f3,'y2'),'y2')))

%  y3-y2
disp(simplify(diff(diff(f1,'y3'),'y2')))
disp(simplify(diff(diff(f2,'y3'),'y2')))
disp(simplify(diff(diff(f3,'y3'),'y2')))


%  y1-y3
disp(simplify(diff(diff(f1,'y1'),'y3')))
disp(simplify(diff(diff(f2,'y1'),'y3')))
disp(simplify(diff(diff(f3,'y1'),'y3')))

%  y2-y3
disp(simplify(diff(diff(f1,'y2'),'y3')))
disp(simplify(diff(diff(f2,'y2'),'y3')))
disp(simplify(diff(diff(f3,'y2'),'y3')))

%  y3-y3
disp(simplify(diff(diff(f1,'y3'),'y3')))
disp(simplify(diff(diff(f2,'y3'),'y3')))
disp(simplify(diff(diff(f3,'y3'),'y3')))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              derivatives of 3 f's wrt 3 p's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p1

disp(simplify(diff(f1,'p1')))
disp(simplify(diff(f2,'p1')))
disp(simplify(diff(f3,'p1')))

% p2

disp(simplify(diff(f1,'p2')))
disp(simplify(diff(f2,'p2')))
disp(simplify(diff(f3,'p2')))

% p3

disp(simplify(diff(f1,'p3')))
disp(simplify(diff(f2,'p3')))
disp(simplify(diff(f3,'p3')))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     second derivatives of 3 f's wrt 3 y's and wrt 3 p's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  y1-p1
disp(simplify(diff(diff(f1,'y1'),'p1')))
disp(simplify(diff(diff(f2,'y1'),'p1')))
disp(simplify(diff(diff(f3,'y1'),'p1')))

%  y2-p1
disp(simplify(diff(diff(f1,'y2'),'p1')))
disp(simplify(diff(diff(f2,'y2'),'p1')))
disp(simplify(diff(diff(f3,'y2'),'p1')))

%  y3-p1
disp(simplify(diff(diff(f1,'y3'),'p1')))
disp(simplify(diff(diff(f2,'y3'),'p1')))
disp(simplify(diff(diff(f3,'y3'),'p1')))



%  y1-p2
disp(simplify(diff(diff(f1,'y1'),'p2')))
disp(simplify(diff(diff(f2,'y1'),'p2')))
disp(simplify(diff(diff(f3,'y1'),'p2')))

%  y2-p2
disp(simplify(diff(diff(f1,'y2'),'p2')))
disp(simplify(diff(diff(f2,'y2'),'p2')))
disp(simplify(diff(diff(f3,'y2'),'p2')))

%  y3-p2
disp(simplify(diff(diff(f1,'y3'),'p2')))
disp(simplify(diff(diff(f2,'y3'),'p2')))
disp(simplify(diff(diff(f3,'y3'),'p2')))



%  y1-p3
disp(simplify(diff(diff(f1,'y1'),'p3')))
disp(simplify(diff(diff(f2,'y1'),'p3')))
disp(simplify(diff(diff(f3,'y1'),'p3')))

%  y2-p3
disp(simplify(diff(diff(f1,'y2'),'p3')))
disp(simplify(diff(diff(f2,'y2'),'p3')))
disp(simplify(diff(diff(f3,'y2'),'p3')))

%  y3-p3
disp(simplify(diff(diff(f1,'y3'),'p3')))
disp(simplify(diff(diff(f2,'y3'),'p3')))
disp(simplify(diff(diff(f3,'y3'),'p3')))


