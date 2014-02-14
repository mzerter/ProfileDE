% SIV model for flu dataset of Miao, et al (2010)
%
% 
% 
%% Add Paths

clear;

addpath('../../fdaM')
addpath('Miao')
addpath('SSE')
addpath('id')

%% Load Data

% Load Data
MiaoData = xlsread('vl.xls');

% Split into data vector and time vector

Earlyind = find(MiaoData(:,1) <= 5);
Earlyind = Earlyind([1:35,37:60]);  % drop outlier value of 3.5 at time 2.5

SIVlog_Time = MiaoData(Earlyind,1);

n  = length(Earlyind);  %  59

SIVlog_Data = zeros(n,3);
SIVlog_Data(:,1) = NaN;
SIVlog_Data(:,2) = NaN;
% SIVlog_Data(1,1) = 5.8e5;
SIVlog_Data(1,1) = log10(5.8) + 5;  %  log S version
% SIVlog_Data(1,2) = 0;
% SIVlog_Data(1,2) = -2;  %  log I version
SIVlog_Data(:,3) = MiaoData(Earlyind,2);


%% Setup 

% Time points

t_max  = 5;
t_span = [0,t_max];

% Basis Functions

knots  = (0:0.1:t_max)';
norder = 4;
nbasis = length(knots) + norder - 2;  % 53

% What does this do apart from name variables?

SIVlog_fdnames    = cell(1,3);                                  
SIVlog_fdnames{3} = cell(1,2);
SIVlog_fdnames{3}{1} = 'Variables';
SIVlog_fdnames{3}{2} = cell(1,3);
SIVlog_fdnames{3}{2}{1} = 'S';
SIVlog_fdnames{3}{2}{2} = 'I';
SIVlog_fdnames{3}{2}{3} = 'log10 V';

SIVlog_basis = create_bspline_basis(t_span, nbasis, norder, knots);

% Set up a functional parameter object.  When the first argument is a 
% basis, the functional data object defaults to all coefficients zer0.

lambda = 1e-0;
SIVlog_fdPar = fdPar(SIVlog_basis, 2, lambda);

% Smooth Basis

SIVlog_Vfd = smooth_basis(SIVlog_Time, SIVlog_Data(:,3), ...
                             SIVlog_fdPar, ones(n,1), ...
                             SIVlog_fdnames);

figure(1)
plotfit_fd(SIVlog_Data(:,3), SIVlog_Time, SIVlog_Vfd)
xlabel('\fontsize{13} Days')
ylabel('\fontsize{13} Viral load (V)')
title('\fontsize{16} Smoothing of EID50 viral load')
                
%  set initial parameter values

rho   = 0;
beta  = 2.4e-6;
delta = 6.0e-1;
pival = 1e2;
c     = 4.2;

%  parameters to be estimated

active = [2,3,5]';

%  set up initial parameter values

SIVlog_pars0 = [rho,beta,delta,pival,c];

%  set up log of starting values

SIVlog_pars0 = log(SIVlog_pars0);

%  set up initial coefficients

SIcoefs = zeros(nbasis,1);

Scoefs = SIcoefs;
% Scoefs(1,1) = 5.8e5;
Scoefs(1,1) = log10(5.8) + 5;  %  log S version 

Icoefs = SIcoefs;

Vcoefs = getcoef(SIVlog_Vfd);

SIVlog_coefs0 = [Scoefs, Icoefs, Vcoefs];

% Make Functions

SIVlog_fn = make_SIVlog;

%  Set individual lambda values

SIVlog_lambda = 1e-2;

%%  Use FitMatchOpt to optimize coefficients given initial parameters

%  set poslik and posproc

poslik  = 0;  %  because data are already logged
posproc = 0;  %  because ode equations are for logged data

%  run LS.setup

[SIVlog_lik, SIVlog_proc] = ...
    LS_setup(SIVlog_fn, SIVlog_Time, [], [], SIVlog_lambda, ...
             SIVlog_Vfd, [], [], [], [], poslik, posproc);
                  
options_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'Display', 'iter');

%  The next task is to improve the initial coefficient estimates 13 and 7
%  defined above for S and E, respectively by using function FitMatchOpt
%  that optimize coefficients for unobserved variables.  The 'which' 
%  argument of this function indicates which variables's coefficients are 
%  to be optimized with respect to the fit to the observed data for I.
%  Function proc.time is used here to show the elapsed time, the time
%  titled 'user' is the current elapsed session time in seconds.
%  The optimization function to be used is function nlminb.

which = [1,2];

coefs0 = SIVlog_coefs0;
pars   = SIVlog_pars0;
proc   = SIVlog_proc;

tic;
SIVlog_coefs1 = FitMatchOpt(SIVlog_coefs0, which, SIVlog_pars0, ...
                            SIVlog_proc, options_in);
toc

%  0.2 seconds, 18 iterations

% Let's have a look at the three functions and the fit to the I data

DEfd1 = fd(SIVlog_coefs1(:,3),SIVlog_basis);

figure(1)
plotfit_fd(SIVlog_Data(:,3),SIVlog_Time,DEfd1)
xlabel('\fontsize{13} Day')
ylabel('\fontsize{13} log10 viral load V')

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


%  Initialize the inner optimization with the coefficients produced by
%  the FitMatchOpt optimization

global INNEROPT_COEFS0

INNEROPT_COEFS0 = SIVlog_coefs1; 

%  Modify the default optimization setting to display results for
%  each iteration

options_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'Display', 'iter');

%  Now carry out an inner optimization to refine the FitMatchOpt results

%  set up variables for inside inneropt

tic;
SIVlog_coefs2 = inneropt(SIVlog_Time, SIVlog_Data, SIVlog_pars0, SIVlog_lik, ...
                         SIVlog_proc, [], options_in);
toc;

%  0.1 seconds,  6 iterations

DEfd2 = fd(SIVlog_coefs2(:,3),SIVlog_basis);

% Has this changed much?  Plot the new function estimates along with the 
% old

plot(DEfd2)
hold on
plot(DEfd1)
hold off

% All these analyses have been using the initial parameter estimates in
% order to get starting values for the coefficients defining the three
% functions.  Now we are ready to optimize the fit to the data with
% respect to the parameters in SIVlog_pars.  Here, however, we only optimize
% the values of the 'i' parameter defining the external input of infectives
% and the three coefficients defining the time-varying infection rate beta. 
% These are selected using the argument 'active' of outer optimization
% function outeropt.

% The computation time involved in this step is substantial (around 10 
% minutes on the computer used while preparing these notes).   You may
% want to disable output buffering by clicking the popup menu item 
% 'Buffered output'that you see when you click on the 'Misc' item at the
% top of the command window.

%  Define the inner and outer optimization settings

options_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'Display', 'off');
                  
options_out = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                       'Hessian', 'off', ...
                       'Display', 'iter');

tic;
[SIVlog_pars3, SIVlog_coefs3, value, gradient] = ...
    outeropt(SIVlog_Time, SIVlog_Data, SIVlog_coefs2, SIVlog_pars0, ...
             SIVlog_lik, SIVlog_proc, ...
             active, [], options_in, [], options_out);
toc


%% RUN EXAMPLE

SIVlog_lambda = 1e0;

global INNEROPT_COEFS0

INNEROPT_COEFS0 = SIVlog_coefs3; 

tic;
SIVlog_struct = Profile_LS(fn, SIVlog_Time,   SIVlog_Data, ...
                                SIVlog_coefs3, SIVlog_pars0, ...
                                SIVlog_basis,  SIVlog_lambda, ...
                                [], [], [], [], active);
toc
                            
% 5.4 secs, 19 iterations
                            
SIVlog_coefs4 = SIVlog_struct.coefs; 
SIVlog_pars   = SIVlog_struct.pars;
SIVlog_fd4    = fd(SIVlog_coefs4,SIVlog_basis);

% SIVlog_coefs0 = SIVlog_coefs;

%% DISPLAY RESULTS

rho_est   = SIVlog_pars(1);
beta_est  = SIVlog_pars(2);
delta_est = SIVlog_pars(3);
pi_est    = SIVlog_pars(4);
c_est     = SIVlog_pars(5);

disp('Parameter estimates:')
disp(['Estimate rho   = ',num2str(rho_est)])
disp(['Estimate beta  = ',num2str(beta_est)])
disp(['Estimate delta = ',num2str(delta_est)])
disp(['Estimate pi    = ',num2str(pi_est)])
disp(['Estimate c     = ',num2str(c_est)])

SIVlog_fd4 = fd(SIVlog_coefs4,SIVlog_basis);

figure(2)
subplot(3,1,1)
plotfit_fd(SIVlog_Data(:,1),SIVlog_Time,SIVlog_fd4(1));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} log10 Susceptible (S)')
subplot(3,1,2)
plotfit_fd(SIVlog_Data(:,2),SIVlog_Time,SIVlog_fd4(2));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} log10 Infected (I)')
subplot(3,1,3)
plotfit_fd(SIVlog_Data(:,3),SIVlog_Time,SIVlog_fd4(3));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} log10 Viral load (V)')

timefine = linspace(0,5,201);
SIV_vec4 = eval_fd(timefine, SIVlog_fd4);
SIV_fd4  = smooth_basis(timefine, 10.^SIV_vec4, ...
                        fdPar(SIVlog_basis,2,1e-8));

figure(3)
subplot(3,1,1)
plot(timefine,SIV_vec4(:,1),'-');
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Susceptible (S)')
subplot(3,1,2)
plot(timefine,SIV_vec4(:,1),'-');
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Infected (I)')
subplot(3,1,3)
plot(SIVlog_Time,10.^SIVlog_Data(:,3),'o', ...
     timefine,SIV_vec4(:,3),'b-');
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Viral load (V)')

%%  use FitMatchOpt to set up initial coefficients

%  This code requires a previous run of Profile_LS in order to define
%  the proc argument of FitMatchOpt
%  Otherwise, an explicit execution of LS_setup is required to define
%  the proc argument

coefs0 = zeros(nbasis,3);
coefs0(:,3) = getcoef(SIVlog_Vfd);
coefs0 = FitMatchOpt(coefs0, [1,2], SIVlog_pars, SIVlog_struct.proc);

SIVlog_fd0 = fd(coefs0,SIVlog_basis);

figure(2)
subplot(3,1,1)
plotfit_fd(SIVlog_Data(:,1),SIVlog_Time,SIVlog_fd0(1));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Susceptible (S)')
subplot(3,1,2)
plotfit_fd(SIVlog_Data(:,2),SIVlog_Time,SIVlog_fd0(2));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Infected (I)')
subplot(3,1,3)
plotfit_fd(SIVlog_Data(:,3),SIVlog_Time,SIVlog_fd0(3));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Recovered (r)')

SIVlog_lambda = 1e0;
SIVlog_coefs0 = coefs0;

INNEROPT_COEFS0 = SIVlog_coefs0; 

SIVlog_struct = Profile_LS(fn,SIVlog_Time,SIVlog_Data,SIVlog_coefs0,...
                            SIVlog_pars0,SIVlog_basis,SIVlog_lambda);

SIVlog_coefs = SIVlog_struct.coefs; 
SIVlog_pars  = SIVlog_struct.pars;


%% DISPLAY RESULTS

beta_est = SIVlog_pars(1);
nu_est   = SIVlog_pars(2);
R0_est   = beta_est*N/nu_est;
disp('Parameter estimates:')
disp(['Estimate beta = ',num2str(beta_est)])
disp(['Estimate nu   = ',num2str(nu_est)])
disp(['Estimate R0   = ',num2str(R0_est)])

SIVlog_fd = fd(SIVlog_coefs,SIVlog_basis);

figure(2)
subplot(3,1,1)
plotfit_fd(SIVlog_Data(:,1),SIVlog_Time,SIVlog_fd(1));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Susceptible (S)')
subplot(3,1,2)
plotfit_fd(SIVlog_Data(:,2),SIVlog_Time,SIVlog_fd(2));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Infected (I)')
subplot(3,1,3)
plotfit_fd(SIVlog_Data(:,3),SIVlog_Time,SIVlog_fd(3));
xlabel('\fontsize{13} Weeks')
ylabel('\fontsize{13} Recovered (r)')



                   


