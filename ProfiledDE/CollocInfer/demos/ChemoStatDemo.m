%  add paths to required functions

addpath('../../fdaM')
addpath('fhn')
addpath('SSE')
addpath('id')
addpath('logtrans')
addpath('logstate_lik')
addpath('genlin')
addpath('findif')
addpath('loggenlin')
addpath('ChemoStat')

% The Chemostat equations represent a four-species Chemostat plus the  
% resource of Nitrogen. There are two species of Algae with varying
% defenses against Rotifers. The Rotifers themselves are divided into two 
% class -- breeding and senescent, although these two are very tightly 
% coupled.
%
% A full description of these equations can be found in the user manual. 
% The five state variables for the equations are
% 
% N  - nitrogen content in the Chemostat
% C1 - Algal type 1
% C2 - Algal type 2 
% B  - Breeding Rotifers
% S  - Senescent Rotifers
%
% The system has 16 parameters, also described in the user manual. Notable
% features include that only the sums C1+C2 and B+S can be observed.  
% Further, an unknown fraction of each is counted at each time. This  
% requires us to set up a model for the observation process along with the  
% ODE.

%  if this .mat file is set up, go to line 172

load ChemoStatDemoSetup

% First we load up some data

load ChemoData.txt

% The first two of these parameters give the fractions of Algae and 
% Rotifers that are counted. The remaining parameters are all positive and  
% using their logged values is helpful. 

ChemoTime = 0:100;

ChemoPars = [1.000e+02 1.000e+00 5.000e-02 1.000e+00 1.600e+02 ...
             3.000e-01 5.500e-02 4.000e-01 2.025e-03 4.400e+00 ...
             2.200e+00 2.700e+02 1.100e-02 1.700e+02 1.500e-01 ...
             1.500e-01];
         
ChemoParNames = ['a1    '; ...    
                 'a2    '; ...    
                 'p1    '; ...    
                 'p2    '; ...    
                 'NI    '; ... 
                 'delta '; ...         
                 'm     '; ...
                 'lambda'; ...
                 'XC    '; ...   
                 'KC1   '; ...    
                 'KC2   '; ...    
                 'rho   '; ...      
                 'G     '; ...   
                 'XB    '; ...    
                 'KB    '; ... 
                 'Qstar '];

ChemoVarNames = ['N '; 'C1'; 'C2'; 'B '; 'S '];

% Parameters 'p1' and 'p2' represent relative palatability of the two algal
% clones, as such only one can be estimated and we fix p2 = 0. 

ChemoLogPars = log(ChemoPars);

ChemoLogData = log(ChemoData);

active = [1:2,5,7:16];     

% We'll choose a fairly large value of lambda. 

lambda = ones(5,1);  

% We need some basis functions

rng        = [0,100];
knots      = rng(1):0.5:rng(2); 
% knots  = [0:3:6, 6.5:0.5:13.5, 14:4:50, 50.5:0.5:59.5, 60:4:100];
nbasis     = length(knots)+2;
ChemoBasis = create_bspline_basis(rng,nbasis,4,knots);

% We will also have to set up the basis matrices manually. 

%  basis values at observation times

bvals_obs = eval_basis(ChemoTime,ChemoBasis);

%  basis values at quadrature points

mids = [0,(knots(1:(length(knots)-1)) + 0.25),100];
bvals_proc.bvals  = eval_basis(mids,ChemoBasis);
bvals_proc.dbvals = eval_basis(mids,ChemoBasis,1);

% We can now set up the proc object. We will want to take a log 
% transformation of the state here for numerical stability. 
% In general it is better to do finite differencing AFTER the 
% log transformation rather than before it. 

proc                = make_SSEproc;          % Sum of squared errors
proc.bvals          = bvals_proc;            % Basis values

proc.more           = make_findif_ode;       % Finite differencing
Chemotrans          = make_id;
logTransStruct.fn   = Chemotrans.fn;
logTransStruct.eps  = 1e-8;
proc.more.qpts      = mids;                  % Quadrature points
proc.more.weights   = ones(5,1).*lambda;     % Quadrature weights 
proc.more.names     = ChemoVarNames;         % Variable names
proc.more.parnames  = ChemoParNames;         % Parameter names

proc.more.more      = logTransStruct;        % Log transform
 
proc.more.more.more.fn = @chemo_fun;         % ODE function

%  display the four-level cascade of functions in proc

disp('proc: SSE functions for evaluating fit to DE')
disp(proc)                %  SSE functions for evaluating fit to DE
disp('proc.more: findif   functions for difference-derivatives')
disp(proc.more)           %  findif   functions for difference-derivatives
disp('proc.more.more: logtrans functions')
disp(proc.more.more)      %  logtrans functions
disp('proc.more.more.more: chemo functions for right side of DE')
disp(proc.more.more.more) %  chemo functions for right side of DE

% For the lik object we need to both represent the linear combination 
%  transform and we need to model the observation process. 

% First to represent the observation process, we can use the genlin
% functions. These produce a linear combination of the the states 
% (they can be used in proc objects for linear systems, too).  

lik = make_logstate_lik;

temp_lik      = make_SSElik;
temp_lik.more = make_loggenlin;

% Genlin requires a more object with two elements. The 'mat' element
% gives a template for the matrix defining the linear combination. This is
% all zeros 2x5 in our case for the two observations from five states. 
% The 'sub' element specifies which elements of the parameters should be
% substituted into the mat element. 'sub' should be a kx3 matrix, each
% row defines the row (1) and column (2) of 'mat' to use and the element
% of the parameter vector (3) to add to it. 

temp_lik.more.more.mat = zeros(2,5);
temp_lik.more.more.sub = [1,2,1; ...
                          1,3,1; ...
                          2,4,2; ...
                          2,5,2];
temp_lik.more.weights  = [100,1];

% Finally, we tell CollocInfer that the trajectories are represented on
% the log scale, and these must be exponentiated before comparing them 
% to the data. 

lik.more  = temp_lik;
lik.bvals = bvals_obs;

%  Display the three-level cascade of functions in lik, with the 
%  fourth level being the mat and sub matrices required by genlin fns.

disp('lik: logstate functions')
disp(lik)              %  logstate functions 
disp('lik.more: SSE functions for evaluating fit')
disp(lik.more)          %  SSE functions for evaluating fit
disp('lik.more.more: loggenlin functions fitting composite observations')
disp(lik.more.more)     %  genlin functions fitting composite observations
disp('lik.more.more.more: mat and sub matrices for genlin functions')
disp(lik.more.more.more)  %  mat and sub matrices for genlin functions

% Now let's try running this

% Because we don't have direct observations of any state, we'll use a 
% starting smooth obtained from generating some ODE solutions

%  define observation vector at time 0

y0 = 0.05.*exp(randn(2,1))*log([2, 0.1, 0.4, 0.2, 0.1]);

%  test loggenlin functions

lik.more.more.fn(     0, y0, logpars, lik.more.more.more)
lik.more.more.dfdx(   0, y0, logpars, lik.more.more.more)
lik.more.more.dfdp(   0, y0, logpars, lik.more.more.more)
lik.more.more.d2fdx2( 0, y0, logpars, lik.more.more.more)
lik.more.more.d2fdxdp(0, y0, logpars, lik.more.more.more)

%  test chemo_ode function

chemo_ode(0, y0', logpars)

y0 = log([2, 0.1, 0.4, 0.2, 0.1]);

%  approximate the solution to the equation

[odetrajt, odetrajy] = ode45(@chemo_ode,ChemoTime,y0,[],ChemoLogPars);

plot(odetrajt, odetrajy)

DEfd = smooth_basis(ChemoTime, odetrajy, ...
                    fdPar(ChemoBasis,int2Lfd(2),1e-6));
                
coefs0 = getcoef(DEfd);

subplot(1,1,1)
for i=1:5
%     subplot(5,1,i)
    plot(coefs0(:,i))
    pause
end

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

%  initialize the inner optimization with the coefficients produced by
%  the FitMatchOpt optimization

INNEROPT_COEFS0 = coefs0; 

% Now, with parameters fixed, we'll estimate coefficients. 

[f, grad] = SplineCoefs(coefs0, ChemoTime, ChemoLogData, ChemoLogPars, ...
                        lik, proc);

grad = reshape(grad,nbasis,5); 
subplot(1,1,1)
for i=1:5
%     subplot(5,1,i)
    plot(grad(:,i))
    pause
end

%  call to inneropt

control_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'Display', 'iter');

INNEROPT_COEFS0 = coefs0; 


tic;
coefs1 = inneropt(ChemoTime, ChemoData, ChemoLogPars, lik, proc, ...
                  [], control_in);
toc

% We'll for the trajectory and also the appropriate sum of exponentiated
% states to compare to the data. 

traj    = lik.bvals*coefs1;
obstraj = lik.more.more.fn(ChemoTime,exp(traj),logpars,lik.more.more.more);

% Plot these against the data

subplot(2,1,1)
plot(ChemoTime, obstraj(:,1), '-', ChemoTime, ChemoData(:,1), 'bo' )
ylabel('Chlamy')
xlabel('')
subplot(2,1,2)
plot(ChemoTime, obstraj(:,2), '-', ChemoTime, ChemoData(:,2), 'bo' )
ylabel('Brachionus')
xlabel('days')

% Now we can continue with the outer optimization


control_in = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'Display', 'off', ...
                      'MaxIter', 100);


control_out = optimset('LargeScale', 'off', 'GradObj', 'on', ...
                      'Hessian', 'off', 'Diagnostics', 'off', ...
                      'MaxIter', 100, ...
                      'Display', 'iter');

INNEROPT_COEFS0 = coefs1; 

tic;
[pars_opt2, coefs2, value, gradient] = ...
         outeropt(ChemoTime, ChemoData, coefs1, logpars, ...
                  lik, proc, active, ...
                  [], control_in, [], control_out);
toc

% We'll extract the resulting parameters and coefficients. 

npars = pars_opt2;
coefs2 = reshape(coefs2,size(coefs1));

% And obtain an estimated trajectory and the exponentiated sum to comprare
% to the data. 

traj  = lik.bvals*coefs2;
ptraj = lik.more.more.fn(ChemoTime,exp(traj),npars,lik.more.more.more);

% Lets have a look at how much we changed our parameters on the original 
% scale. 

new.pars = npars;
new.pars(3:16) = exp(new.pars(3:16));

disp([log10(ChemoPars)', log10(new.pars)'])
disp(log10(new.pars))
disp(log10(new.pars./ChemoPars))

% Now we can produce a set of diagnostic plots. 

% Firstly, a representation of the trajectory compared to the data. 

subplot(2,1,1)
plot(ChemoTime,ptraj(:,1),'b-', ChemoTime,ChemoData(:,1), 'bo')
ylabel('Chlamy')
xlabel('')
subplot(2,1,2)
plot(ChemoTime,ptraj(:,2),'b-', ChemoTime,ChemoData(:,2), 'bo')
ylabel('Brachionus')
xlabel('days')


% Now we'll plot both the derivative of the trajectory and the value of the
% differential equation right hand side at each point. This represents the 
% fit to the model. 

traj2  = proc.bvals.bvals*coefs2;
dtraj2 = proc.bvals.dbvals*coefs2;

ftraj2 = proc.more.fn(proc2.more.qpts,traj2,npars,proc.more.more);


for i = 1:5
  subplot(5,1,i)  
  plot(mids,dtraj2(:,i))
  xlabel('')
  ylabel(ChemoVarnames(i))
  lines(mids,ftraj2(:,i))
end
legend('Smooth','Model','location','NorthWest')

% Solving the differential equation from the estiamted initial condition
% of the trajectory allows us to compare the qualitative behavior of 
% our estimate to that of the differential equation. 

y0 = traj(1,:);
odetraj = ode45(@chemo.ode,y0,ChemoTime,[],npars);


subplot(2,1,1)
plot(ChemoTime,traj)
ylabel('')
title('Reconstructed Trajectories')
legend(ChemoVarnames,'location','NorthEast')
subplot(2,1,2)
plot(ChemoTime,odetraj(:,2:6))
ylabel('')
title('ODE Solution')

% We can also compare the pattern of observations predicted by the 
% differential equation and that estimated by our methods. 

otraj = lik.more.more.fn(ChemoTime,exp(odetraj(:,2:6)),npars, ...
                         lik.more.more.more);

subplot(2,1,1)
plot(ChemoTime,ptraj)
xlabel('days')
ylabel('')
title('Predicted Observations -- Smooth')
plot(ChemoTime,ChemoData,'.')
legend(['Algae','Rotifers'],'location','NorthEast')
subplot(2,1,2)
plot(ChemoTime,otraj)
xlabel('days')
ylabel('')
title('Predicted Observations -- ODE')
legend(['Algae','Rotifers'],'location','NorthEast')

