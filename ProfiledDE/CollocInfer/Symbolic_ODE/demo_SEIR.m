
clear all; close all; clc;

% Create symbolic variables
syms t;
p = sym('p%d_%d',[1,8]);
x = sym('x%d_%d',[1,3]);

% Create a function that describes the ODE using symbolic variables like 
% SEIR below.

% Define B(t,p) and ODE
Beta =  p(6) + p(7)*sin(2*pi*t) + p(8)*cos(2*pi*t);

f = [-( Beta*x(1)*(p(2) + x(3)) ) + p(1) - p(3)*x(1);
     Beta*x(1)*(p(2) + x(3)) - p(4)*x(2) - p(3)*x(2);
     p(4)*x(2) - p(5)*x(3) - p(3)*x(3)];

% Create an instance of ODE  
odeproblem = ODE(f,t,x,p);

% Set the ode to use exp(p); 
% odeproblem.SetExponentialParameter();

% Compute the derivatives needed for collocinfer
odeproblem.computeDerivatives();

% Display the results
odeproblem.displayResult();

% Get function handles for ODE, what collocinfer in expecting
fn_ODE = odeproblem.make();


% Test the results against existing function, set up a problem and 
% evaluate the functions.

x = [1,2,3;4,5,7];
SEIRpars = [1.0000e+05, 0.0000e+00, 2.0000e-02, 4.5625e+01, ...
            7.3000e+01, 2.4820e-04, 0.0000e+00, 1.9856e-05]; 
t = 0;



% Need to add necessary path here...!!
addpath('../SEIR')

% Here is the usual make for SEIR
fn_SEIR = make_SEIR();

% The beta function is specified in more for the regular make.
beta_fun  = @(t,p,more) p(6) + p(7)*sin(2*pi*t) + p(8)*cos(2*pi*t);
beta_dfdp = @(t,p,more) [ones(length(t),1), sin(2*pi*t), cos(2*pi*t)];
beta_ind  = [6,7,8];

betamore.beta_fun  = beta_fun;
betamore.beta_dfdp = beta_dfdp;
betamore.beta_ind  = beta_ind;


%Compare output, can do this for .dfdx ect..
disp('Value for regular make')
fn_SEIR.dfdx(t,x,SEIRpars,betamore)
disp('Value for make using ODE')
fn_ODE.dfdx(t,x,SEIRpars,betamore)



