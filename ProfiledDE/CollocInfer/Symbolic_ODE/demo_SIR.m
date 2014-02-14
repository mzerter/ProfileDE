% Declare symbolic variables

syms t;
p = sym('p%d_%d',[1,3]);
x = sym('x%d_%d',[1,3]);

% Create a function that describes the ODE using symbolic variables like 
% SIR below.

f = [- p(1)*x(1)*x(2) - p(3)*x(1) + p(3);
     p(1)*x(1)*x(2) - (p(2) + p(3)) * x(2);
     p(2)*x(2) - p(3)*x(3)];
 
% Create an instance of ODE  
odeproblem = ODE(f,t,x,p);

% Set the ode to use exp(p); 
odeproblem.SetExponentialParameter();

odeproblem.computeDerivatives();

% Display the results
odeproblem.displayResult();

x = [1,2,3;4,5,7];
p = [1,2,7];
t = 0;
more = 0;

fn = odeproblem.make();

% Need to add necessary path here...!!
addpath('../SIR')

fn_SIR = make_SIR();

fn_SIR.dfdx(t,x,p,more)
fn.dfdx(t,x,p,more)

