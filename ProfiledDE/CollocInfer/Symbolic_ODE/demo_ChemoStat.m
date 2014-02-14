clear all; close all; clc;

% Create symbolic variables
syms t;
p = sym('p%d_%d',[1,16]);
x = sym('x%d_%d',[1,5]);

syms Q Qs;

Q = p(3)*x(2) + p(4)*x(3);

Qs = Q*exp(10*(Q - p(16)))/(1 + exp(10*(Q - p(16)))) + ...
     p(16)/(1 + exp(10*(Q - p(16))));

% Create a function that describes the ODE using symbolic variables like 
% Chemo example below.

f = [ p(6)*(p(5) - x(1)) - p(12)*x(2)*x(1)/(p(10) + x(1)) - p(12)*x(3)*x(1)/(p(11) + x(1));
 x(2)*(p(9)*p(12)*x(1)/(p(10) + x(1)) - p(3)*p(13)*(x(4) + x(5))/(p(15) + Qs) - p(6));
 x(3)*(p(9)*p(12)*x(1)/(p(11) + x(1)) - p(4)*p(13)*(x(4) + x(5))/(p(15) + Qs) - p(6));
 x(4)*(p(14)*p(13)*Q/(p(15) + Qs) - (p(6) + p(7) + p(8)));
 p(8)*x(4) - (p(6) + p(7))*x(5)];


odeproblem = ODE(f,t,x,p);

% Set the ode to use exp(p); 
odeproblem.SetExponentialParameter();


% Compute the derivatives needed for collocinfer
odeproblem.computeDerivatives();

% Display the results
odeproblem.displayResult();

% Test the results against existing function, set up a problem and 
% evaluate the functions.

x = [1,2,3,4,5;4,5,7,8,10];
p = [1:16];
t = 0;
more = 0;

fn = odeproblem.make();

% Need to add necessary path here...!!
addpath('../ChemoStat')

fn_chemo = make_chemo();

fn_chemo.fn(t,x,p,more) 

fn.fn(t,x,p,more)

