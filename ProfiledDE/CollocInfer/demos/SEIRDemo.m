%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This demonstration illustrates the use of log transformations with
% CollocInfer. 
%  We use the SEIR equations, various modifications of which are used to model
%  the spread of disease in epidemeiology.  The letters stand for:
%
%  Susceptibles:  the large subpopulation of those who have not been infected
%  Exposed:       the smaller subpopulation of those who have been infected,
%                 but are not yet infectious
%  Infectious:    those who have already been infected, and are now also 
%                 able to infect others
%  Recovered:     these were infected, but are now recovered and are, usually,
%                 immune to further infections.
%
%  The equations are used to model measles infections in the Canadian province
%  of Ontario over a two-year period.  The data analyzed here are derived from
%  a larger dataset analyzed in greater detail by Hooker, Ellner, De Vargas
%  Roditi and Earn (2011).  
%
%  The basic equations are here modified to all for a seasonally varying
%  infection rate, a realistic feature for modeling childhood diseases where
%  the infection rate increases during the school year.  The function that
%  provides the seasonal variability is beta(t).
%
% The equations are as follows.  Each equation has been organized so that
% the first term on the right defines the exponential decline of the 
% subpopulation if it were left to itself, and the remaining terms, called
% "forcing" terms, define the inputs from the other variables and external
% sources that cause the population to increase (or decline further).
%
% Sdot = -nu*S         + mu - beta(t)*S*(I+i) (Susceptibles)
% Edot = -(sigma+nu)*E + beta(t)*S*(I+i)      (Exposed)
% Idot = -(gamma+nu)*I + sigma*E              (Infectious)
% Rdot = -nu*R         + gamma*I               (Recovered)
% 
% Here beta(t) is parameterized by a sinusoidal function plus a constant,
%      beta(t) = p_1 + p_2 sin(2 pi t) + p_3 cos(2 pi t)
% where the period is one year.
%
% Other parameters in the equations are
% mu    - the arrival rate of new individuals into the population at large
% nu    - the death rate in the population at large
% i     - an external source of infected individuals
% sigma - the rate of movement from Exposed to Infectious.
% gamma - the rate of recovery from infection.
%
% Because the number of susceptibles is usually many times the number of
% infected individuals, S is on a quite different scale than E and I. 
% In these data, S varies around about 450,000, but E and I vary between
% about 400 and 9000. 
% Moreover, all of these variables cannot in preinciple be negative.  Both 
% computational efficiency and positivity argue for fitting the logarithm of the 
% data, and also modelling the logarithm of each variable. 
% This requires a simple modification of the SEIR equations above, consisting
% of dividing the right side of the equations by the exponential of the 
% respective variables.
% The CollocInfer package contains a number of useful tools that let you make 
% this transition without needing to re-code all of your differential equations.

%  add paths to required functions

addpath('../../fdaM')
addpath('SEIR')
addpath('SSE')
addpath('id')
addpath('exp')
addpath('findif')
addpath('logtrans')
      
% Get some data and parameters

%  data

fid  = fopen('SEIRdata.txt','rt');
SEIRdata = fscanf(fid,'%f');
fclose(fid);

%  observation times

SEIRtimes = linspace(0,5,261)';

%  parameter values

SEIRpars = [1.0000e+05, 0.0000e+00, 2.0000e-02, 4.5625e+01, ...
            7.3000e+01, 2.4820e-04, 0.0000e+00, 1.9856e-05]; 

%  variable names
        
SEIRvarnames = ['S'; 'E'; 'I'];

%  parameter names

SEIRparnames = ['mu   '; 'i    '; 'nu   '; 'sigma'; ...
                'gamma'; 'b0   '; 'b1   '; 'b2   '];

% Now augment the data for I  so that S and E measurements are listed as NA

data = [nan.*ones(length(SEIRdata),2),SEIRdata];

% Calculate the natural logs of the observations  
% Note: using common logs would have made the results easier for clients to
% interpret, and would have required dividing the SEIR right sides above by
% 10 to the power of the variable values rather than vy their exponentials.

logdata = log(data);

% Function make.SEIR is distributed with the package, and defines 
% the right sides of the differential equations as well as various of their
% partial derivatives required during the computation.   The code can be 
% viewed by typing "make.SEIR" into the command window in R.

SEIRfn = make_SEIR;

% Define the time-varying infection rate beta  as a sinusoid with period
% one year plus a constant

beta_fun  = @(t,p,more) p(6) + p(7)*sin(2*pi*t) + p(8)*cos(2*pi*t);

% Define as well its partial derivatives with respect to the coefficients

beta_dfdp = @(t,p,more) [ones(length(t),1), sin(2*pi*t), cos(2*pi*t)];

%  indices of the coefficients in the parameter vector SEIRpars

beta_ind  = [6,7,8];

%  Bundle these two functions as well as the names of their coefficients
%  into a struct object called betamore

betamore.beta_fun  = beta_fun;
betamore.beta_dfdp = beta_dfdp;
betamore.beta_ind  = beta_ind;

% Set up a B-spline functional basis object to represent the state vector

rr     = [0,5];                  %  the range of observations times
knots  = rr(1):2/52:rr(2);       %  knots at 52 equally spaced values
norder = 3;                      %  order of the B-spline basis functions,
                                 %  in this case piece-wise quadratic
nbasis = length(knots)+norder-2; %  the number of basis functions

%  set up the basis object

bbasis = create_bspline_basis(rr,nbasis,norder,knots);

% To get an initial estimate of the states we smooth the observed I 
% component and set the other coefficients to zero.  

% smooth the log I values

DEfd = smooth_basis(SEIRtimes,logdata(:,3),fdPar(bbasis,1,0.1));

% plot the smooth plus data

plotfit_fd(logdata(:,3),SEIRtimes,DEfd)

% Augment the coefficient vector by two columns.  The first column
% contains all 13's, which is about the logarithm of the susceptible 
% population size, and the second contains all 7's, about the estimated
% size of the exposed subpopulation.  The ratio of the two population sizes
% is here approximated by exp(13-7), or about 400.

coefs0 = [ones(nbasis,2), getcoef(DEfd)];
coefs0(:,1) = coefs0(:,1).*13;
coefs0(:,2) = coefs0(:,2).*7;

% set up the functional data object for the three variables

DEfd = fd(coefs0,bbasis);

% We use the usual least squares measures of the fidelity of the estimated
% functions S(t), E(t) and I(t) to the data and to the differential equation.
% This means that we can save the user the effort of setting up the lik and
% proc objects by calling function LS.setup, which output these objects.

% We will want to represent the state variables on the log scale so that they 
% remain always positive. To do this, we set the 'posproc' argument of 
% LS.setup to 1. We will also compare the log I(t) to the log of the data 
% directly, and therefore set argument 'poslik' to 0, the default.

%  Here, too, we specify the level of emphasis on fitting each differential
%  equation by setting the corresponding smoothing or bandwidth parameter
%  lambda.  The following two statements specify the lambda to be 100 for
%  S, and 1 for E and I.

SEIRlambda = ones(1,3);
SEIRlambda(1) = 100;

%  struct object betamore is now passed into LS.setup in order to make
%  it available as a functional parameter defined by its three coefficients

%  run LS.setup

[SEIRlik, SEIRproc] = LS_setup(SEIRfn, SEIRtimes, [], [], SEIRlambda, ...
                               DEfd, betamore, [], [], [], 0, 1);
                  
%  these control options modifiy the default settings for Matlab
%  function fminunc called within the inner optimization loop.

control_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'Display', 'iter');

%  The next task is to improve the initial coefficient estimates 13 and 7
%  defined above for S and E, respectively by using function FitMatchOpt
%  that optimize coefficients for unobserved variables.  The 'which' argument
%  of this function indicates which variables's coefficients are to be
%  optimized with respect to the fit to the observed data for I.
%  Function proc.time is used here to show the elapsed time, the time
%  titled 'user' is the current elapsed session time in seconds.
%  The optimization function to be used is function nlminb.

tic;
coefs1 = FitMatchOpt(coefs0, [1,2], SEIRpars, SEIRproc, control_in);
toc

% Comment:  Optimization went much better with the LargeScale
% setting turned off rather then when turned on.  When turned off,
% Matlab function uses the quasi-Newton method and does not use
% the Hessian matrix.

% Let's have a look at the three functions and the fit to the I data

DEfd1 = fd(coefs1(:,3),bbasis);

plotfit_fd(log(SEIRdata),SEIRtimes,DEfd1)

% We can now run an initial smooth using the estimated coefficients as 
%  starting points. 

%  The starting coefficient matrix for the inner optimizations must be
%  declared global and then defined at this point.  This is required
%  because inner optimizations should be started with the coefficients
%  as optimized on the last inner optimization, and this one of a few
%  techniques for achieving this.  The R package has even more trouble with
%  this in the sense of needing to write the optimized coefficients to a 
%  file after each inner optimization, and then retrieving prior to a new
%  optimization.  The Matlab strategy used here does all this in memory,
%  which is considerably faster.

global INNEROPT_COEFS0

%  Initialize the inner optimization with the coefficients produced by
%  the FitMatchOpt optimization

INNEROPT_COEFS0 = coefs1; 

%  Modify the default optimization setting to display results for
%  each iteration

control_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'Display', 'off');

%  Now carry out an inner optimization to refine the FitMatchOpt results

tic;
coefs2 = inneropt(SEIRtimes, logdata, SEIRpars, SEIRlik, SEIRproc, [], ...
                  control_in);
toc;

% Comment:  Again optimization went better with the LargeScale
% setting turned off rather then when turned on.  Nevertheless, it
% still failed to optimize within the default limit of 400 iterations.

%  define the new function estimates

DEfd2 = fd(coefs2(:,3),bbasis);

% Has this changed much?  Plot the new function estimates along with the old

plot(DEfd2)
hold on
plot(DEfd1)
hold off

% All these analyses have been using the initial parameter estimates in
% order to get starting values for the coefficients defining the three
% functions.  Now we are ready to optimize the fit to the data with
% respect to the parameters in SEIRpars.  Here, however, we only optimize
% the values of the 'i' parameter defining the external input of infectives
% and the three coefficients defining the time-varying infection rate beta. 
% These are selected using the argument 'active' of outer optimization
% function outeropt.

% The computation time involved in this step is substantial (around 10 minutes
% on the computer used while preparing these notes).   You may
% want to disable output buffering by clicking the popup menu item 
% 'Buffered output'that you see when you click on the 'Misc' item at the
% top of the command window.

active = [2, 6, 7, 8];

%  Define the inner and outer optimization settings

control_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'Display', 'off');
                  
control_out = optimset('LargeScale', 'on', 'GradObj', 'on', ...
                       'Hessian', 'off', ...
                       'Display', 'iter');

tic;
[SEIRpars3, coefs3, value, gradient] = ...
    outeropt(SEIRtimes, logdata, coefs2, SEIRpars, SEIRlik, SEIRproc, active, ...
             [], control_in, [], control_out);
toc

%  Comments:  
%             Outer LargeScale off terminated in 5 its., and at 3.55656
%             Outer LargeScale on  terminated in 9 its., and at 3.54875
%     Whether Inner LargeScale was on or off made no difference to the
%     final function values.

%  initialize the inner optimization 

control_in = optimset('LargeScale', 'on', 'GradObj', 'on', ...
                      'Hessian', 'on', 'Diagnostics', 'off', ...
                      'Display', 'off');
                  
control_out = optimset('Jacobian', 'on', 'MaxIter', 500, ...
                       'Display', 'iter');

INNEROPT_COEFS0 = coefs2; 

tic;
Profile_LS_struct = ...
    Profile_LS(SEIRfn, SEIRtimes, logdata, coefs2, SEIRpars, [], SEIRlambda, ...
               DEfd, betamore, [], [], active, ...
               [], control_in, [], control_out, [], [], 1);
toc

%  The Matlab function lsqnonlin, which uses a Gauss-Newton method, does
%  not do as well as function fminunc called by outeropt, since it fails
%  to converge in the allowed number of iterations, only reaches a
%  final criterion value of 3.56699, and takes longer to get there.

% Some plots

%  display the initial and estimated parameter values

[SEIRpars(active); SEIRpars3(active)]

% Plot the estimated trajectories

DEfd3 = fd(coefs3,bbasis);
plot(DEfd3)

% Let's compare this to the data

plotfit_fd(logdata(:,3),SEIRtimes,DEfd3(3))
ylabel('Fit to Data')

% We can also look at the discrepancy between the estimated trajectory and the
% differential equation

traj = eval_fd(SEIRtimes,DEfd3);     % estimated trajectory

dtraj = eval_fd(SEIRtimes,DEfd3,1);  % derivative of the estimated trajectory

ftraj = SEIRproc.more.fn(SEIRtimes,traj,SEIRpars3,SEIRproc.more.more);   % Trajectory predicted by ODE

plot(SEIRtimes,dtraj)
ylabel('SEIR derivatives' )

plot(SEIRtimes,ftraj)

plot(SEIRtimes,dtraj-ftraj)
ylabel('Fit to Model')

%% The alternative is to exponentiate the state before we compare to original data.
% This can take a very long time and is only recommended if you really need to do
% it, or have a couple of hours to wait. 

[SEIRlik2, SEIRproc2] = LS_setup(SEIRfn, SEIRtimes, [], [], [], DEfd, betamore, ...
                         [], [], [], 1, 1);
                     
INNEROPT_COEFS0 = coefs3; 

%  Now carry out an inner optimization to refine the FitMatchOpt results

%  Modify the default optimization setting to display results for
%  each iteration

control_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'on', 'Diagnostics', 'off', ...
                      'Display', 'iter', 'TolFun', 1e-3);

coefs4 = inneropt(SEIRtimes, data, SEIRpars3, SEIRlik2, SEIRproc2, [], control_in);

control_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'Display', 'off');
                  
control_out = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                       'Hessian', 'off', ...
                       'Display', 'iter');
                   
[SEIRpars5, coefs5, value, gradient] = ...
    outeropt(SEIRtimes, data, coefs4, SEIRpars3, SEIRlik2, SEIRproc2, active, ...
             [], control_in, [], control_out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some more basic setup operations
%
% Here we go through the steps necessary to set up the SEIR equations manually. 
% This allows us several options for dealing with the positivity of the 
% state vector. 
%
% 1. We can ignore it and hope for the best. 
%
% 2. We can take a log transform of the ODE and then exponentiate the solutions
% to compare to the data
%
% 3. We can take a log transform of the ODE and compare this to the log data. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% First of all, we need the values of the basis at the observation times and the 
% quadrature points. 

qpts = 0.5*(knots(1:(length(knots)-1)) + knots(2:length(knots)));

bvals.obs = eval_basis(SEIRtimes,bbasis);

bvals.bvals  = eval_basis(qpts,bbasis);
bvals.dbvals = eval_basis(qpts,bbasis,1);

%%% This proc object is just the standard squared error deviation
% from the right hand side of a differential equation

sproc               = make_SSEproc;
sproc.bvals         = bvals;
sproc.more          = make_SEIR;
sproc.more.more     = betamore;
sproc.more.qpts     = qpts;
sproc.more.weights  = ones(length(qpts),3)*diag([1e2,1e0,1e0]);
sproc.more.names    = SEIRvarnames;
sproc.more.parnames = SEIRparnames;

%%% However, ODEs are often much more numerically stable if represented on a
% log scale. The make_logtrans function will take the right hand side functions
% and derivatives defined for any differential equation and provide the equivalent
% system for the log state. Note that this does affect the way  you represent
% your state when considering the observations. 

lsproc                = make_SSEproc;
lsproc.bvals          = bvals;
lsproc.more           = make_logtrans;
lsproc.more.more      = make_SEIR;
lsproc.more.more.more = betamore;
lsproc.more.qpts      = qpts;
lsproc.more.weights   = ones(length(qpts),3)*diag([1e2,1e0,1e0]);

%%% Lik objects, this is the standard squared error. 

slik               = make_SSElik;
slik.bvals         = eval_basis(SEIRtimes,bbasis);
slik.more          = make_id;
slik.more.weights  = ones(size(data));
slik.more.names    = SEIRvarnames;
slik.more.parnames = SEIRparnames;

% Log transform transformation for the lik object. For this, we note that we 
% have represented the trajectory on the log scale and will need to transform 
% back, this makes the numerics much harder and it can take a very long time to 
% converge.  

lslik               = make_logstate_lik;
lslik.bvals         = slik.bvals;
lslik.more.weights  = slik.more.weights;
lslik.more          = slik;
lslik.more.parnames = SEIRparnames;

% Numerically things work much better on the log scale

dfd = smooth_basis(SEIRtimes,logdata(:,3),bbasis);

coefs0 = zeros(nbasis,3);
coefs0(:,3) = getcoef(dfd);
 
coefs1 = FitMatchOpt(coefs, which, SEIRpars3, lsproc, control_in);

INNEROPT_COEFS0 = coefs1; 

control_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'on', 'Diagnostics', 'off', ...
                      'Display', 'iter');

coefs2 = inneropt(SEIRtimes, logdata, SEIRpars, lslik, lsproc, [], control_in);

res3 = outeropt(SEIRtimes,data,res2.coefs,SEIRpars,lslik,lsproc);

%% Or just log the observations, this is faster. 

INNEROPT_COEFS0 = coefs3; 

control_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'on', 'Diagnostics', 'off', ...
                      'Display', 'iter');

coefs2 = inneropt(SEIRtimes, logdata, SEIRpars3, lslike, lsproc, [], control_in);

control_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'Display', 'off');
                  
control_out = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                       'Hessian', 'off', ...
                       'Display', 'iter');
                   
[SEIRpars3, coefs3, value, gradient] = ...
    outeropt(SEIRtimes, logdata, coefs2, SEIRpars, slik, lsproc, active, ...
             [], control_in, [], control_out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

