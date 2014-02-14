% SIR model for Bombay 1906 dataset
%
% Author: Mathieu Zerter  
% Date: Oct 1, 2013     
% Partner:  
% 
% 
% SIRs of S (say for a sort function): 
% 
%    >> grades = [ 95 52 95 80 71 85 78 12 93 99 ]; 
%    >> sorted_array = sort( array ) 
%    ans = 
%           [ 12 52 71 78 80 85 93 95 95 99 ] 
% 
% 
%% Add Paths

clear;

addpath('../../fdaM')
addpath('SIR')
addpath('SSE')
addpath('id')
addpath('logtrans')


%% Load Data

% Load Data
load ('BOMBAY_1906.MAT', '-mat')
% load ('SIRdata.mat', '-mat')

% Split into data vector and time vector

n  = size(BOMBAY_1906,1);

%  define the data matrix, putting NaN's where thereis no data

SIR_Data = zeros(n,3);
SIR_Data(:,1) = NaN;
SIR_Data(:,2) = BOMBAY_1906;
SIR_Data(:,3) = NaN;

%% Setup 

% Time points

t_max  = 30;
t_span = [0,t_max];

SIR_Time = linspace(0,t_max,n);

%  Define the single data value for S.  The population of Bombay was about
%  one million at that point.  Problem: it seems pretty unreasonable to 
%  suppose that the suceptible population was as large as that.  The plague
%  only infects people living with rats.

% Basis Functions

knots  = (0:3:t_max)';
norder = 4;
nbasis = length(knots) + norder - 2;
lambda = 1e-4;

% What does this do apart from name variables?

SIR_fdnames    = cell(1,3);                                  
SIR_fdnames{3} = cell(1,2);
SIR_fdnames{3}{1} = 'Variables';
SIR_fdnames{3}{2} = cell(1,3);
SIR_fdnames{3}{2}{1} = 'S';
SIR_fdnames{3}{2}{2} = 'I';
SIR_fdnames{3}{2}{3} = 'R';

SIR_basis = create_bspline_basis(t_span, nbasis, norder, knots);

% Set up a functional parameter object.  When the first argument is a 
% basis, the functional data object defaults to all coefficients zer0.

SIR_fdPar = fdPar(SIR_basis, 2, lambda);

% Smooth infected cases

SIR_Ifd = smooth_basis(SIR_Time, SIR_Data(:,2), ...
                             SIR_fdPar, ones(n,1), ...
                             SIR_fdnames);
  
plotfit_fd(SIR_Data(:,2), SIR_Time, SIR_Ifd)
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Infected (I)')
title('\fontsize{16} Smoothing of number of infected data')
                
%  set initial coefficient values

SRcoefs = zeros(nbasis,1);

Scoefs = SRcoefs;
Rcoefs = SRcoefs;
Icoefs = getcoef(SIR_Ifd);

SIR_coefs0 = [Scoefs, Icoefs, Rcoefs];

% Make Functions

SIR_fn = make_SIR;

%  Run the model with these parameters and initial conditions  
%  N, SIR_data(1,2), 0

%  set initial parameter values and N

N = 2.4e4;
beta0 = 5.3e-5;
nu0   = 0.94;
mu0   = 1e-10;
R0    = N*beta0/nu0;

SIR_pars0 = [beta0, nu0, mu0]';

%   convert to logarithmic scale

logbeta0 = log(beta0);
lognu0   = log(nu0);
logmu0   = log(mu0);

SIR_pars0 = [logbeta0, lognu0, logmu0]';

SIR_Data0 = [2.4e4, SIR_Data(1,2), 0]';

[tODE, yODE] = ode45(@SIR_fn.fn, t_span, SIR_Data0, [], SIR_pars0, []);

figure(1)
subplot(3,1,1)
plot(tODE,yODE(:,1),'-')
subplot(3,1,2)
plot(tODE,yODE(:,2),'-',SIR_Time,SIR_Data(:,2),'o')
subplot(3,1,3)
plot(tODE,yODE(:,3),'-')

%  smooth the data to get coefficients for S and R

n = length(tODE);
SIR_Sfd = smooth_basis(tODE, yODE(:,1), SIR_fdPar);
SIR_Rfd = smooth_basis(tODE, yODE(:,3), SIR_fdPar);

SIR_logSfd = smooth_basis(tODE,     log(yODE(:,1)),     SIR_fdPar);
SIR_logIfd = smooth_basis(SIR_Time, log(SIR_Data(:,2)), SIR_fdPar);
SIR_logRfd = smooth_basis(tODE,     log(yODE(:,3)+1),   SIR_fdPar);

SIR_coefs0(:,1) = getcoef(SIR_logSfd);
SIR_coefs0(:,2) = getcoef(SIR_logIfd);
SIR_coefs0(:,3) = getcoef(SIR_logRfd);

%  Set individual lambda values

SIR_lambda = 1e-2;


%% RUN EXAMPLE

global INNEROPT_COEFS0

INNEROPT_COEFS0 = SIR_coefs0; 
      
SIR_struct = Profile_LS(fn, SIR_Time, SIR_Data, ...
                            SIR_coefs0, SIR_pars0, ...
                            SIR_basis, SIR_lambda);

SIR_coefs = SIR_struct.coefs; 
SIR_pars  = SIR_struct.pars;

% SIR_coefs0 = SIR_coefs;

%% DISPLAY RESULTS

beta_est = SIR_pars(1);
nu_est   = SIR_pars(2);
R0_est   = beta_est*N/nu_est;
disp('Parameter estimates:')
disp(['Estimate beta = ',num2str(beta_est)])
disp(['Estimate nu   = ',num2str(nu_est)])
disp(['Estimate R0   = ',num2str(R0_est)])

SIR_fd = fd(SIR_coefs,SIR_basis);

figure(2)
subplot(3,1,1)
plotfit_fd(SIR_Data(:,1),SIR_Time,SIR_fd(1));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Susceptible (S)')
subplot(3,1,2)
plotfit_fd(SIR_Data(:,2),SIR_Time,SIR_fd(2));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Infected (I)')
subplot(3,1,3)
plotfit_fd(SIR_Data(:,3),SIR_Time,SIR_fd(3));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Recovered (r)')

%%  use FitMatchOpt to set up initial coefficients

%  This code requires a previous run of Profile_LS in order to define
%  the proc argument of FitMatchOpt
%  Otherwise, an explicit execution of LS_setup is required to define
%  the proc argument

coefs0 = zeros(nbasis,3);
coefs0(:,2) = getcoef(SIR_Ifd);
coefs0 = FitMatchOpt(coefs0, [1,3], SIR_pars, SIR_struct.proc);

SIR_fd0 = fd(coefs0,SIR_basis);

figure(2)
subplot(3,1,1)
plotfit_fd(SIR_Data(:,1),SIR_Time,SIR_fd0(1));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Susceptible (S)')
subplot(3,1,2)
plotfit_fd(SIR_Data(:,2),SIR_Time,SIR_fd0(2));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Infected (I)')
subplot(3,1,3)
plotfit_fd(SIR_Data(:,3),SIR_Time,SIR_fd0(3));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Recovered (r)')

SIR_lambda = 1e0;
SIR_coefs0 = coefs0;

INNEROPT_COEFS0 = SIR_coefs0; 

SIR_struct = Profile_LS(fn,SIR_Time,SIR_Data,SIR_coefs0,...
                            SIR_pars0,SIR_basis,SIR_lambda);

SIR_coefs = SIR_struct.coefs; 
SIR_pars  = SIR_struct.pars;


%% DISPLAY RESULTS

beta_est = SIR_pars(1);
nu_est   = SIR_pars(2);
R0_est   = beta_est*N/nu_est;
disp('Parameter estimates:')
disp(['Estimate beta = ',num2str(beta_est)])
disp(['Estimate nu   = ',num2str(nu_est)])
disp(['Estimate R0   = ',num2str(R0_est)])

SIR_fd = fd(SIR_coefs,SIR_basis);

figure(2)
subplot(3,1,1)
plotfit_fd(SIR_Data(:,1),SIR_Time,SIR_fd(1));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Susceptible (S)')
subplot(3,1,2)
plotfit_fd(SIR_Data(:,2),SIR_Time,SIR_fd(2));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Infected (I)')
subplot(3,1,3)
plotfit_fd(SIR_Data(:,3),SIR_Time,SIR_fd(3));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Recovered (r)')

%%  Fit log data 

logdata = log(SIR_Data);

% Smooth log infected cases

SIR_logIfd = smooth_basis(SIR_Time, logdata(:,2), SIR_fdPar);

plotfit_fd(logdata(:,2), SIR_Time, SIR_logIfd)
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} log Infected (I)')
title('\fontsize{16} Smoothing of log ofnumber of infected data')

%  set poslik and posproc

poslik  = 1;  %  because the fit is to log data
posproc = 1;  %  because penalty term is on log of path

%  run LS.setup

[SIR_lik, SIR_proc] = ...
    LS_setup(SIR_fn, SIR_Time, [], [], SIR_lambda, ...
             SIR_logIfd, [], [], [], [], poslik, posproc);
                  
%  The next task is to improve the initial coefficient estimates 13 and 7
%  defined above for S and E, respectively by using function FitMatchOpt
%  that optimize coefficients for unobserved variables.  The 'which' 
%  argument of this function indicates which variables's coefficients are 
%  to be optimized with respect to the fit to the observed data for I.
%  Function proc.time is used here to show the elapsed time, the time
%  titled 'user' is the current elapsed session time in seconds.
%  The optimization function to be used is function nlminb.

%  set up the initial coefficients

SIR_coefs0(:,1) = getcoef(SIR_logSfd);
SIR_coefs0(:,2) = getcoef(SIR_logIfd);
SIR_coefs0(:,3) = getcoef(SIR_logRfd);

%  optimize the coefficients given the initial parameters

options_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'Display', 'iter');

which = 1:3;

tic;
SIR_coefs1 = FitMatchOpt(SIR_coefs0, which, SIR_pars0, ...
                            SIR_proc, options_in);
toc

SIR_fd1 = fd(SIR_coefs1, SIR_basis);

plot(SIR_fd1)

%  0.21 seconds, 93 iterations
% We can now run an initial smooth using the estimated coefficients as 
%  starting points. 

%  The starting coefficient matrix for the inner optimizations must be
%  declared global and then defined at this point.  This is required
%  because inner optimizations should be started with the coefficients
%  as optimized on the last inner optimization, and this one of a few
%  techniques for achieving this.  

%  Initialize the inner optimization with the coefficients produced by
%  the FitMatchOpt optimization

global INNEROPT_COEFS0

INNEROPT_COEFS0 = SIR_coefs1; 

%  Modify the default optimization setting to display results for
%  each iteration

options_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'Display', 'iter');

%  Now carry out an inner optimization to refine the FitMatchOpt results

%  set up variables for inside inneropt

tic;
SIR_coefs2 = inneropt(SIR_Time, logdata, SIR_pars0, SIR_lik, ...
                         SIR_proc, [], options_in);
toc;

%  0.5 seconds,  320 iterations

DEfd2 = fd(SIR_coefs2(:,2),SIR_basis);

% Has this changed much?  Plot the new function estimates along with the 
% old

plot(DEfd2)
hold on
plot(DEfd1)
hold off

% All these analyses have been using the initial parameter estimates in
% order to get starting values for the coefficients defining the three
% functions.  Now we are ready to optimize the fit to the data with
% respect to the parameters in SIR_pars.  Here, however, we only optimize
% the values of the 'i' parameter defining the external input of infectives
% and the three coefficients defining the time-varying infection rate beta. 
% These are selected using the argument 'active' of outer optimization
% function outeropt.

%  Define the inner and outer optimization settings

options_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'Display', 'off');
                  
options_out = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                       'Hessian', 'off', ...
                       'Display', 'iter');

tic;
[SIR_pars3, SIR_coefs3, value, gradient] = ...
    outeropt(SIR_Time, logdata, SIR_coefs2, SIR_pars0, ...
             SIR_lik, SIR_proc, ...
             [], [], options_in, [], options_out);
toc


%% RUN EXAMPLE

logdata(:,1) = log(1e6);

active = 1:2;

SIR_lambda = 1e2;

% Make Functions

SIR_fn = make_SIR;

fn       = SIR_fn;
times    = SIR_Time;
data     = logdata;
coefs    = SIR_coefs0;
allpars  = SIR_pars0;
basisobj = SIR_basis;
lambda   = SIR_lambda;

fdobj       = [];
more        = [];
obsweights  = [];
quadrature  = [];
options_out = [];
out_method  = [];
options_in  = []; 
in_method   = [];
poslik      = 0;
posproc     = 0;
discrete    = 0;
diffeps     = 0;

global INNEROPT_COEFS0

INNEROPT_COEFS0 = SIR_coefs0; 

tic;
SIR_struct = Profile_LS(SIR_fn, SIR_Time, logdata, SIR_coefs0, ...
                        SIR_pars0, SIR_basis,  SIR_lambda, ...
                        [], [], [], [], active);
toc
                            
% 0.3 secs, 6 iterations
                            
SIR_coefs4 = SIR_struct.coefs; 
SIR_pars   = SIR_struct.pars;
SIR_fd4    = fd(SIR_coefs4,SIR_basis);

%% DISPLAY RESULTS

beta_est = exp(SIR_pars(1));
nu_est   = exp(SIR_pars(2));
mu_est   = exp(SIR_pars(3));
R0_est   = N*beta_est/nu0;

disp('Parameter estimates:')
disp(['Estimate beta = ',num2str(beta_est)])
disp(['Estimate nu   = ',num2str(nu_est)])
disp(['Estimate mu   = ',num2str(mu_est)])
disp(['Estimate R0   = ',num2str(R0_est)])

timefine = linspace(0,30,201);
SIR_mat4 = eval_fd(timefine, SIR_fd4);

figure(2)
subplot(3,1,1)
plot(timefine, SIR_mat4(:,1));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} log10 Susceptible (S)')
subplot(3,1,2)
plot(timefine, SIR_mat4(:,2), '-', SIR_Time, logdata(:,2), 'o');
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} log10 Infected (I)')
subplot(3,1,3)
plot(timefine, SIR_mat4(:,3));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} log10 Recovered (R)')

figure(3)
subplot(3,1,1)
plot(timefine, exp(SIR_mat4(:,1)));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Susceptible (S)')
subplot(3,1,2)
plot(timefine, exp(SIR_mat4(:,2)), '-', ...
     SIR_Time, SIR_Data(:,2), 'o');
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Infected (I)')
subplot(3,1,3)
plot(timefine, exp(SIR_mat4(:,3)));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Recovered (R)')

%  Conclusion:  Manipulating N and setting all of S to N have no impact on
%  the fit to the I data, and the effect of changing lambda over 1e-2 to
%  1e2 is miniscule on the fit.  Basically, there is no information in I
%  about the level of S and its derivative can be change to accommodate
%  all of these changes.

%  However, forcing the parameters to be positive by exponentiation was
%  very effective, and should always be used.

