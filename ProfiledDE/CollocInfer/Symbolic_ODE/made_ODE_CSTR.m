function FhN = made_ODE_CSTR()
%% make_ODE_CSTR Replaces usual make process for CSTR


% Create symbolic variables, size appropriately
syms t;
p = sym('p%d_%d',[1,3]);
x = sym('x%d_%d',[1,2]);

% Create a function that describes the ODE using symbolic variables like 
% CSTR below.

% Define System Parameters


% INPUTS:
% % Known functions of time.
% 
% % Function of t and F_in
% B_CC = [];
% % Function of F_CO and F_IN
% B_TT = [];
% % Function of t and F_IN
% B_TC = [];
% % Function of F_CO

% These functions have to be defined

% a = [];
syms testfn(t);
f = [(p(1) + 1/(testfn(t)+x(2))).*x(1);0];
aval = testfn(1);
disp(aval);
% Define ODE
% 
% 
% f = [p(3)*(x(1) - x(1).^3/3 + x(2));
%      (-1/p(3))*(x(1) - p(1) + p(2) * x(2))];

% Create an instance of ODE  
odeproblem = ODE(f,t,x,p);

% Set the ode to use exp(p) in functions; 
% odeproblem.SetExponentialParameter();

% Compute the derivatives needed for collocinfer
odeproblem.computeDerivatives();
odeproblem.displayResult();


% Get function handles for ODE, what collocinfer is expecting
 FhN = odeproblem.make();

% odeproblem.plot([0,20],[-1,1],[0.2,0.2,3],[]);
end
%     function x = testfn(t)
%         if t>10
%             x = 1;
%         else x = 2;
%         end;
%     end