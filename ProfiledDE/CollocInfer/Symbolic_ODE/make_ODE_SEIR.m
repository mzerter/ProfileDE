function SEIR = make_ODE_SEIR()
%% make_ODE_SEIR Replaces usual make process for SEIR
% Place in CollocInfer/SEIR folder and try out with SEIRdemo.m
% by changing make_SEIR to make_ODE_SEIR

% Create symbolic variables
syms t;
p = sym('p%d_%d',[1,8]);
x = sym('x%d_%d',[1,3]);

% Create a function that describes the ODE using symbolic variables like 
% SEIR below.

% Define B(t,p) and ODE
B =  p(6) + p(7)*sin(2*pi*t) + p(8)*cos(2*pi*t);

f = [-( B*x(1)*(p(2) + x(3)) ) + p(1) - p(3)*x(1);
     B*x(1)*(p(2) + x(3)) - p(4)*x(2) - p(3)*x(2);
     p(4)*x(2) - p(5)*x(3) - p(3)*x(3)];

% Create an instance of ODE  
odeproblem = ODE(f,t,x,p);

% Set the ode to use exp(p) in functions; 
% odeproblem.SetExponentialParameter();

% Compute the derivatives needed for collocinfer
odeproblem.computeDerivatives();

% Get function handles for ODE, what collocinfer is expecting
SEIR = odeproblem.make();

odeproblem.plot([0,10],[1,1,1],[1,1,1,1,1,1,1,1],0);
end
