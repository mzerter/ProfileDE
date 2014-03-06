function FhN = made_ODE_CSTR()
%% make_ODE_CSTR Replaces usual make process for CSTR


% Create symbolic variables, size appropriately
syms t;
p = sym('p%d_%d',[1,3]);
x = sym('x%d_%d',[1,2]);

% Create a function that describes the ODE using symbolic variables like 
% CSTR below.

% Define System Parameters

B_CC = 

% Define ODE


f = [p(3)*(x(1) - x(1).^3/3 + x(2));
     (-1/p(3))*(x(1) - p(1) + p(2) * x(2))];

% Create an instance of ODE  
odeproblem = ODE(f,t,x,p);

% Set the ode to use exp(p) in functions; 
% odeproblem.SetExponentialParameter();

% Compute the derivatives needed for collocinfer
odeproblem.computeDerivatives();

% Get function handles for ODE, what collocinfer is expecting
FhN = odeproblem.make();

odeproblem.plot([0,20],[-1,1],[0.2,0.2,3],[]);
end