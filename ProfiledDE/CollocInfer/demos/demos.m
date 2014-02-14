          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %        FitzHugh-Nagumo Example         %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  add paths to required functions

addpath('../../fdaM')
addpath('fhn')
addpath('SSE')
addpath('id')
addpath('findif')
addpath('logtrans')

% Obtain some pre-generated data

load FhNdata.txt
FhNtimes = (0:0.5:20)';
FhNpars  = [0.2, 0.2, 3.0];

% Now we want to set up a basis expansion; this makes use of
% routines in the fda library

knots  = (0:0.5:20)';
norder = 3;
nbasis = length(knots) + norder - 2;
range  = [0,20];

bbasis = create_bspline_basis(range, nbasis, norder, knots);

% Collocation points for quadrature at midpoints between knots

qpts = knots(1:(length(knots)-1))+diff(knots)/2;  

% We'll start off by creating a smooth of each state variable to get 
% initial values for the parameters.

lambda = 1;
Lfdobj = int2Lfd(1);  %  penalize the size of the first derivative
bfdPar = fdPar(bbasis,Lfdobj,lambda);

%  Set up variable names, using the case names to identify variables
%  See P. 48 of Ramsay, Hooker and Graves (2009)

%  first set up cell array of length 2 for variable names

FhNvarnames    = cell(1,2);
FhNvarnames{1} = 'Variables';  %  the generic name
FhNvarnames{2} = ['V';'R'];    %  character array of specific names

%  Now set up the cell array of length 3 for the fdnames object

FhNfdnames = cell(1,3);
FhNfdnames{1} = 'Time';
FhNfdnames{2} = FhNvarnames;
FhNfdnames{3} = 'Variables';  %  Also the generic name

%  now smooth the data to get a preliminary estimate the paths

DEfd   = smooth_basis(FhNtimes,FhNdata,bfdPar,[],FhNfdnames);

%  extract the coefficient matrix to provide a preliminary estimate

FhNcoefs0 = getcoef(DEfd);

%  The first three rows should be:

%    -0.1969    0.8180
%    -0.0443    0.8838
%     0.6719    0.8102
% 
%  ------- Option 1: pre-defined set-up for squared error  ------------

%  the penalty parameter for the Profiling process:

lambda = 1000;

% First just optimize the coefficients to improve on their initial
%   estimates 

%  This requires an initial coefficient matrix estimate, which is
%  communicated to the inner optimization by Matlab as the values of
%  the global variable INNEROPT_COEFS0

global INNEROPT_COEFS0
INNEROPT_COEFS0 = FhNcoefs0;

%  define the struct object fhn to define functions for evaluating
%  the fits to the data and their derivatives

fhn = make_fhn;

%  Proceed with the new data smoothing

IDEfd1 = Smooth_LS(fhn, FhNtimes, FhNdata, ...
                   FhNcoefs0, FhNpars, bbasis, lambda, DEfd);

% Let's have a look at this

plotfit_fd(FhNdata, FhNtimes, IDEfd1);

% Now we can do the profiling in the easiest way:  using function
% ProfileLS that automatically defines the least squares data-fitting
% and differential equation fitting methods.

LS_struct1 = Profile_LS(fhn, FhNtimes, FhNdata, FhNcoefs0, ...
                        FhNpars, bbasis, lambda);

in_method = [];            % Inner Optimization
control_in.reltol=1e-8;
control_in.maxit=100;
control_in.trace=0;  % Optimization control
control_in.LargeScale = 'on';
control_in.GradObj='on';
control_in.Hessian='on';
control_in.Display='off';

out_method = [];            % Outer Optimization
control_out = optimset('Jacobian', 'on', 'MaxIter', 100, ...
                       'Display', 'iter', ...
                       'TolFun', 1e-6);

LS_struct1 = Profile_LS(fhn, FhNtimes, FhNdata, FhNcoefs0, ...
                        FhNpars, bbasis, lambda, [], [], [], [], [], ...
                        in_method, control_in, out_method, control_out);


% Look at the parameter estimates:

disp(['Parameters: ', num2str(LS_struct1.pars)])

%  display the fits to the data

ODEfd1 = fd(LS_struct1.coefs, bbasis);

plotfit_fd(FhNdata,FhNtimes,ODEfd1)

%  display one-half the error sum of squares of the data fit 

SSEhalf = sum(LS_struct1.res.^2)/2;

disp(['SSE/2 = ',num2str(SSEhalf)])

% SSE/2 = 2.4432

% If we only coded the FitzHugh Nagumo funciton with no derivatives we 
% would just have make_fhn.fn, in which case we default to estimating the 
% derivatives by finite differencing.  This achieved by replacing the first
% argument by the value of member "fn" that only evaluates the fit.

IDEfd1a	= Smooth_LS(fhn.fn, FhNtimes, FhNdata, FhNcoefs0, FhNpars, ...
                    bbasis, lambda);

% Now we can do the profiling with differencing

LS_struct1a = Profile_LS(fhn.fn, FhNtimes, FhNdata, FhNcoefs0, ...
                         FhNpars, bbasis, lambda);

disp(['Parameters: ', num2str(LS_struct1a.pars)])

% Parameters: 0.22951      0.3056      2.7411

%  -------  Option 2:  set-up functions for proc and lik objects -------
%                and then call optimizers

%  set up the lik and proc objects associated with least squares fit
%  measures

[lik, proc] = LS_setup(fhn,FhNtimes,FhNcoefs0,bbasis,lambda);

% Now we can get initial parameter estimates from 'gradient matching'

npars = ParsMatchOpt(FhNpars,FhNcoefs0,proc);

disp(['Parameters: ', num2str(npars)])

% Parameters: 0.27885    0.085152      1.3974

% Smoothing can be done more specifically with 'inneropt'

INNEROPT_COEFS0 = FhNcoefs0;

Ires2 = inneropt(FhNtimes,FhNdata,npars,lik,proc);

% And we can also do profiling by calling the outer optimization
% function 'outeropt'

npars = outeropt(FhNtimes, FhNdata, FhNcoefs0, FhNpars, lik, proc);

method_out = [];
control_out = optimset('GradObj', 'on', ...
                       'Hessian', 'on', ...
                       'MaxIter', 100, ...
                       'display', 'iter', ...
                       'TolFun',  1e-6);

npars = outeropt(FhNtimes, FhNdata, FhNcoefs0, FhNpars, lik, proc, ...
                 [], [], [], method_out, control_out);

disp(['Parameters: ', num2str(npars)])

%  ----------- Option 3:  set everything up manually  ----------------

%% lik object

% we are using squared error
lik2              = make_SSElik;    
% values of the basis at observation times
lik2.bvals        = eval_basis(FhNtimes,bbasis);  
% identity transformation
lik2.more         = make_id;           
% only use data for V
lik2.more.weights = [1,0]; 
%  Needed by Matlab but not by R
lik2.more.more    = [];

%% proc object

FhNvarnames = ['V', 'R'];
FhNparnames = ['a', 'b', 'c'];
% we are using squared error
proc2 = make_SSEproc;                    
% values and derivative of basis at collocation points
proc2.bvals.bvals   = eval_basis(qpts,bbasis);    
proc2.bvals.dbvals  = eval_basis(qpts,bbasis,1);  
% FitzHugh-Nagumo right hand side
proc2.more          = make_fhn;                   
% State variable names
proc2.more.names    = FhNvarnames;            
% Parameter names
proc2.more.parnames = FhNparnames;         
% Collocation or quadrature points
proc2.more.qpts     = qpts;                    
% Weight relative to observations of each collocation point.
proc2.more.weights  = [1000,1000];        
proc2.more.more     = [];                                          
proc2.more.eps      = 1e-6;

% We can also specify optimization methods

in_method = [];            % Inner Optimization
control_in.reltol=1e-8;
control_in.maxit=100;
control_in.trace=0;  % Optimization control
control_in.LargeScale = 'on';
control_in.GradObj='on';
control_in.Hessian='on';

out_method = [];            % Outer Optimization
control_out.trace=6;
control_out.maxit=100;
control_out.reltol=1e-8; % Optimization control
control_out.LargeScale = 'on';
control_out.GradObj='on';
control_out.Hessian='on';
                                                                                   
% Now we will also remove the 'R' part of the data (ie, assume we did not 
% measure it) and set the corresponding coefficients to zero

new.FhNdata = FhNdata;
new.FhNdata(:,2) = nan;

new.coefs = FhNcoefs0;
new.coefs(:,2) = 0;

% We can now fix the smooth or 'V' and try and choose 'R' to make the 
% differential equation match as well as possible.

new.coefs2 = FitMatchOpt(new.coefs, 2, FhNpars, proc2, control_in);

% And we can call the same inner optimization as above, this time with our
% own control parameters and method

INNEROPT_COEFS0 = FhNcoefs0;

in_method = [];

control_in.Display = 'iter';
control_in.trace   = 2;

new.coefs3 = inneropt(FhNtimes,new.FhNdata,FhNpars,lik2,proc2, ...
                      in_method,control_in);

% And we can also do profiling with specified controls

method_out = [];
control_out = optimset('GradObj', 'on', ...
                       'Hessian', 'on', ...
                       'MaxIter', 100, ...
                       'display', 'iter', ...
                       'TolFun',  1e-6);

npars = outeropt(FhNtimes, new.FhNdata, new.coefs3, FhNpars, lik2, ...
                 proc2, [], [], [], method_out, control_out);

disp(['Parameters: ', num2str(npars)])

% We can also use finite differencing to calculate derivatives of the 
% right hand side of the FitzHugh-Nagumo equations, this is defined by 
% modifying the proc object.

proc2a = make_SSEproc;                    
proc2a.bvals.bvals  = eval_basis(qpts,bbasis);    
proc2a.bvals.dbvals = eval_basis(qpts,bbasis,1);  
proc2a.more = make_findif_ode;                   
proc2a.more.names = FhNvarnames;            
proc2a.more.parnames = FhNparnames;         
proc2a.more.qpts = qpts;                    
proc2a.more.weights  = [1000,1000];        
Fhnfn = make_fhn;                           
proc2a.more.more.fn = Fhnfn.fn; 
proc2a.more.more.reltol=1e-6; 
proc2a.more.more.more = [];                        
proc2a.more.more.eps = 1e-3;                                                      

in_method = [];
new.coefs3 = inneropt(FhNtimes,new.FhNdata,FhNpars,lik2,proc2a,...
                      in_method,control_in);

% And we can also do profiling with specified controls

method_out = [];
control_out = optimset('GradObj', 'on', ...
                       'Hessian', 'on', ...
                       'MaxIter', 100, ...
                       'display', 'iter', ...
                       'TolFun',  1e-6);
                   
nparsa = outeropt(FhNtimes,new.FhNdata,new.coefs3,FhNpars, ...
                  lik2,proc2a, [], [], [], out_method,control_out);
              
disp(['Parameters: ', num2str(nparsa)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This demonstration illustrates the use of log transformations with
% CollocInfer. We use the SEIR equations with a seasonally varying
% infection rate for this purpose.
%
% The equations are given as
%
% Sdot = mu - beta(t)*S*(I+i) - nu*S        (Susceptibles)
% Edot = beta(t)*S*(I+i) - (sigma+nu)*E     (Exposed)
% Idot = sigma*E - (gamma+nu)*I              (Infectious)
% 
% Here beta(t) - the infection rate - is parameterized by a sinusoidal 
% function plus a constant.
%
% Traditionally, there is an additional state
%
% Rdot = gamma*I - nu*R
%
% However, since we only observe I, R contributes nothing to the data fit, 
% and we have removed it from the system. 
%
% Other parameters are
% i     - a visiting process
% nu    - death rate
% sigma - the rate of movement from Exposed to Infectious.
% gamma - the rate of recovery from infection.
%
% It is generally more stable to solve the equations for the log states 
% rather than the states themselves. CollocInfer contains a number of 
% useful tools that let you make this transition without needing to re-code 
% all of your differential equations.

%  add paths to required functions

addpath('../../fdaM')
addpath('SEIR')
addpath('SSE')
addpath('id')
addpath('exp')
addpath('findif')
addpath('logtrans')
      
% %%%% Get some data and parameters
% 
% SEIRpars = SEIRpars
% 

fid  = fopen('SEIRdata.txt','rt');
SEIRdata = fscanf(fid,'%f');
fclose(fid);

SEIRtimes = linspace(0,5,261)';

SEIRpars = [1.0000e+05, 0.0000e+00, 2.0000e-02, 4.5625e+01, ...
            7.3000e+01, 2.4820e-04, 0.0000e+00, 1.9856e-05]; 

SEIRvarnames = ['S'; 'E'; 'I'];

SEIRparnames = ['mu   '; 'i    '; 'nu   '; 'sigma'; 'gamma'; ...
                'b0   '; 'b1   '; 'b2   '];

%%%% Now format the data so that S and E measurements are listed as nan

data = [nan.*ones(length(SEIRdata),2),SEIRdata];

% We'll also look at the log observations

logdata = log(data);

%%%% define the right side evaluation function

SEIRfn = make_SEIR;

%%%% A couple of functions to define the infection rate

beta_fun  = @(t,p,more) p(6) + p(7)*sin(2*pi*t) + p(8)*cos(2*pi*t);
beta_dfdp = @(t,p,more) [ones(length(t),1), sin(2*pi*t), cos(2*pi*t)]; 
beta_ind  = [6,7,8];

betamore.beta_fun  = beta_fun;
betamore.beta_dfdp = beta_dfdp;
betamore.beta_ind  = beta_ind;

% Create a collocation basis to represent the state vector

rr     = [0,5];
knots  = rr(1):2/52:rr(2);
% knots  = [0,2.5,5];
norder = 3;
nbasis = length(knots)+norder-2;

bbasis = create_bspline_basis(rr,nbasis,norder,knots);

% To get an initial estimate of the states we smooth the observed I component
% and set the other coefficients to zero.  

DEfd = smooth_basis(SEIRtimes,logdata(:,3),fdPar(bbasis,1,0.1));

plotfit_fd(log(SEIRdata),SEIRtimes,DEfd)

coefs = [zeros(nbasis,2),getcoef(DEfd)];
DEfd = fd(coefs,bbasis);

% We will want to represent the state variables on the log scale so that they remain
% always positive. To do this, we set the 'posproc' component of LS_setup to 1. We
% will also compare the logstate to the log of the data directly, and therefore set
% 'poslik' to 0.

% We call LS_setup first and use the outputted lik and proc objects to pull the 
% coefficients for the unobserved state variables into line with the 
% differential equation. 

lambda = [100,1,1];
[lik,proc] = LS_setup(SEIRfn, SEIRtimes, [], [], lambda, DEfd,...
                      betamore, [], [], [], [], 1, 0);

allcoefs = coefs;
coefs = allcoefs(:,1:2);
coefs = coefs(:);

[fnval, dfdc, d2fdc2] = ...
    FitMatchCoefs(coefs, allcoefs, [1,2], SEIRpars, proc);

Hdiag = diag(d2fdc2);
Hdiag(1:10)

coefs1 = FitMatchOpt(coefs,1:2,SEIRpars,proc);

% Let's have a look at the result

DEfd1 = fd(coefs1(:,3),bbasis);

plot(DEfd1)

plotfit_fd(log(SEIRdata),SEIRtimes,DEfd1)

% We can now run an initial smooth using the estimated coefficients as starting points. 

global INNEROPT_COEFS0

coefs2 = inneropt(SEIRtimes, logdata, SEIRpars, lik, proc);

%  bring coefs2 from R for debugging purposes

fid  = fopen('SEIRcoefs2.txt','rt');
SEIRcoefs2 = fscanf(fid,'%f');
fclose(fid);
SEIRcoefs2 = reshape(SEIRcoefs2,length(SEIRcoefs2)/3,3);

% Has this changed much?

DEfd2 = fd(coefs2(:,3),bbasis);

plot(DEfd2)
hold on
plot(DEfd1)
hold off

%  check out function ProfileErr_Allpar

active = [2, 6, 7, 8];

[value, gradient] = ...
    ProfileErr(SEIRpars(active), SEIRtimes, logdata, coefs, ...
               SEIRpars, lik, proc, active);
                  
% And we can call the optimizing functions, optimizing only with 
% respect to parameter 'i', and the three coefficients for 'beta'. 

[SEIRpars3, coefs3, value, gradient] = ...
    outeropt(SEIRtimes, logdata, coefs2, SEIRpars, lik, proc, active);

control_in = optimset('LargeScale', 'on', 'GradObj', 'on', ...
                      'Hessian', 'on', 'Diagnostics', 'off', ...
                      'Display', 'iter');
control_out = optimset('LargeScale', 'off', 'GradObj', 'off', ...
                       'Hessian', 'off', ...
                       'Display', 'iter');
                   
[SEIRpars3, coefs3, value, gradient] = ...
    outeropt(SEIRtimes, logdata, coefs2, SEIRpars, lik, proc, active, ...
             [], control_in, [], control_out);

% Some plots

%% First, the estimated trajectories

DEfd3 = fd(res3.coefs,bbasis);
plot(DEfd3)

% Let's compare this to the data

plotfit_fd(logdata(:,3),SEIRtimes,DEfd3(3))
ylabel('Fit to Data')

% We can also look at the discrepancy between the estimated trajectory and the
% differential equation

traj = eval_fd(SEIRtimes,DEfd3);     % estimated trajectory

dtraj = eval_fd(SEIRtimes,DEfd3,1);  % derivative of the estimated trajectory

ftraj = proc.more.fn(SEIRtimes,traj,res3.pars,proc.more.more);   % Trajectory predicted by ODE


plot(SEIRtimes,dtraj)
ylabel('SEIR derivatives' )
plot(SEIRtimes,ftraj)

plot(SEIRtimes,dtraj-ftraj)
ylabel('Fit to Model')

%% The alternative is to exponentiate the state before we compare to original data.
% This can take a very long time and is only recommended if you really need to do
% it, or have a couple of hours to wait. 

objs2 = LS_setup(SEIRfn,SEIRtimes,[],[],[],DEfd,betamore,data, ...
                 1,1,[100,1,1]);

lik2  = objs2.lik;
proc2 = objs2.proc;

res2 = inneropt(SEIRtimes,data,res3.pars,proc2);

res3 = outeropt(SEIRtimes,data,res3.coefs,res3.pars,lik2,proc2);

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

bvals.obs = eval_basis(times,bbasis);

bvals.bvals  = eval_basis(qpts,bbasis);
bvals.dbvals = eval_basis(qpts,bbasis,1);

%%% This proc object is just the standard squared error deviation
% from the right hand side of a differential equation

sproc = make_SSEproc;
sproc.bvals = bvals;
sproc.more = make_SEIR;
sproc.more.more = betamore;
sproc.more.qpts = qpts;
sproc.more.weights = ones(length(qpts),3)*diag(c(1e2,1e0,1e0));
sproc.more.names = SEIRvarnames;
sproc.more.parnames = SEIRparnames;

%%% However, ODEs are often much more numerically stable if represented on a
% log scale. The make_logtrans function will take the right hand side functions
% and derivatives defined for any differential equation and provide the equivalent
% system for the log state. Note that this does affect the way  you represent
% your state when considering the observations. 

lsproc = make_SSEproc;
lsproc.bvals = bvals;
lsproc.more = make_logtrans;
lsproc.more.more = make_SEIR;
lsproc.more.more.more = betamore;
lsproc.more.qpts = qpts;
lsproc.more.weights = ones(length(qpts),3)*diag([1e2,1e0,1e0]);

%%% Lik objects, this is the standard squared error. 

slik = make_SSElik;
slik.bvals = eval_basis(times,bbasis);
slik.more = make_id;
slik.more.weights = ones(size(data));
slik.more.names = SEIRvarnames;
slik.more.parnames = SEIRparnames;

% Log transform transformation for the lik object. For this, we note that we 
% have represented the trajectory on the log scale and will need to transform 
% back, this makes the numerics much harder and it can take a very long time to 
% converge.  

lslik = make_logstate.lik;
lslik.bvals = slik.bvals;
lslik.more.weights = slik.more.weights;
lslik.more = slik;
lslik.more.parnames = SEIRparnames;

% Numerically things work much better on the log scale

dfd = data2fd(logdata(:,3),times,bbasis);

coefs = zeros(nbasis,3);
coefs(:,3) = dfd.coefs;
 
res = FitMatchOpt(coefs,1:2,lsproc,SEIRpars,[]);

res2 = inneropt(SEIRtimes,logdata,SEIRpars,lsproc);

res3 = outeropt(SEIRtimes,data,res2.coefs,SEIRpars,lslik,lsproc);

%% Or just log the observations, this is faster. 

res2 = inneropt(SEIRtimes,logdata,SEIRpars,lsproc);

res3 = outeropt(SEIRtimes,logdata,res2.coefs,SEIRpars,slik,lsproc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Groundwater level measurements as a function of rainfall.  
%
% The data were collected in North Vancouver near a site of a 
% mud slide that caused two fatalities in 2004.
%
% The groundwater levels are a relatively smooth process with a
% rather high signal to noise ratio.  Rainfall events, however,
% are sporadic and can be considered as events with times and 
% intensities, called a marked point process in statistical theory.
%
% Groundwater level is modeled by a first order linear differential
% equation forced by rainfall.
%
% A constant coefficient version of the equation is
% 
% ..\dot{Gend(t) = - \beta G(t) + \alpha R(t - \delta) + \mu..
%
% and a variable coefficient version is
%
% ..\dot{Gend(t) = - \beta(t) G(t) + \alpha(t) R(t - \delta) + \mu(t)..
%
% Coefficient \beta controls the speed of the response
% of groundwater to a rainfall event.  
%
% Coefficient \alpha controls the size of the response
% of groundwater to a rainfall event.   The gain in the system
% is defined as K = \alpha / \beta.
%
% Constant \mu is an adjustment to the overall level of
% groundwater required because the origin of the measurements
% can be viewed as arbitrary or irrelevant.
%
% Lag constant \delta allows a time lapse before groundwater
% begins to respond to a rainfall event.  It was estimated by
% inspecting plots to be three hours for the purposes of these
% analyses, and is therefore not actually estimated by these
% analyses themselves.
%
% The following code allows for each of the remaining three
% parameters to be either constant, or themselves functions 
% of time.  That is,
%
% Although the constant model does a fair job of 
% describing the data, allowing for a mild time-variation in
% these parameters substantially improved the fit.  We used
% five order 4 B-splines to model this time dependency.  
%
% Last modified 26 April 2010
%
%% Setting up the data 
%

%  add paths to required functions

addpath('../../fdaM')
addpath('fhn')
addpath('SSE')
addpath('id')

% Input the data to be analyzed.

load groundwater
load rainfall

%  rainfall should be lagged by about three hours

%Ninput = length(rainfall)
%lag    = 3
%rainfall = c(rep(0,2*lag), rainfall(1:(Ninput-2*lag)))

%  a subset of the data are selected, running from day 79 into day 92
%  if all the data are used, R has to allocate additional storage

% index = 1597:3024  
%index = 1907:2221  %  here we select a smaller set of data 
%yobs = as.matrix(groundwater(index))
%zobs = as.matrix(rainfall(index))
%N = length(yobs)
%tobs = (index - index(1))


yobs = NSgroundwater;
zobs = NSrainfall;
tobs = NStimes;

N    = length(tobs);
rangeval = [0,N-1];

%  plot the groundwater data

subplot(2,1,1)
plot(tobs, yobs, '.')  
xlabel('Time (hours)'); 
ylabel('Groundwater level (metres)')
subplot(2,1,2)
plot(tobs, zobs, '.')  
xlabel('Time (hours)') 
ylabel('Hourly rainfall (millimetres)')

%  set up cumulative rainfall

%  Set up a basis for rainfall itself: 
%     An order 1 or piece-wise constant spline

norder = 1;
nbasis = N + norder - 2;
rainbasis  = create_bspline_basis(rangeval, nbasis, norder);
rainfd = smooth_basis(tobs, zobs, rainbasis);

%  plot the data

plotfit_fd(zobs, tobs, rainfd, N)

%% The basis for representing groundwater: G(t)

% Define the knot sequence.  A knot is placed at each mid-hour
% point, giving maximum capacity to fit the data.  Variable
% rangeval is a vector of length 2 containing the lower and
% upper time limits of the observations, and is in file 
% NorthShoreData.mat

knots  = [rangeval(1), ...
          linspace(rangeval(1)+0.5, rangeval(2)-0.5, N-1), ...
           rangeval(2)];
norder   = 3;
nbasis   = length(knots) + norder - 2;
basisobj = create_bspline_basis(rangeval, nbasis, norder, knots);
 
% The constant coefficient basis is required for constant beta and alpha.  

conbasis = create_constant_basis(rangeval);

%% Definition of functional basis objects for the coefficient functions

% Set up the basis for the coefficient function \beta(t) or
% \beta(G) determining the speed of response and recovery
%
% select the command that specifies the argument of \beta
%
% This basis is for \beta(t) as a function of time.  You may
% change the number of basis functions that are used.
%
% Here we use a constant basis for each coefficient, but at end
% we can do that analysis again with five order 5 B-spline basis functions.

nbetabasis = 1;
if nbetabasis > 1
    betabasis = create_bspline_basis(rangeval, nbetabasis);
else
    betabasis = conbasis;
end

% Set up the basis for the coefficient function \alpha(t)
% multiplying rainfall

nalphabasis = 1;
if nalphabasis > 1
    alphabasis = create_bspline_basis(rangeval, nalphabasis);
else
    alphabasis = conbasis;
end

%  store the bases and rainfall object in a list 'more'

more.betabasis  = betabasis;
more.alphabasis = alphabasis;
more.rainfd     = rainfd;

% Definition of the required named list containing the right side
%  evaluation functions.  

NSfn = make_NS;

%% Prepare the constant coefficient analysis

% Set up the functional parameter object that defines the
% roughness penalty.  Here an important choice must be made:
% how large should the smoothing parameter \lambda be?
% One may be advised to begin with a lowish value, such as
% \lambda = 1, complete the analysis with this value, then
% increase the value to something like \lambda = 100, using
% as starting values for the parameter array pars at this stage
% the array newpars computed at the profiled estimation analysis
% for the previous smaller value.  

lambdaDE = 1e0;  % A good value for the initial analysis.     
penorder = 1;    % The penalty is first order
GfdPar   = fdPar(basisobj,penorder,lambdaDE);     

% Smooth the data with this roughness penalty to get initial
% values for coefficients defining G(t).

DEfd0 = smooth_basis(tobs, yobs, GfdPar);

% plot the data and the fit

plotfit_fd(yobs, tobs, DEfd0) 
xlabel('Time (hours)') 
ylabel('Groundwater level (metres)') 

% extract the coefficient array

coefs0 = getcoef(DEfd0);

% Initial parameter values for constant bases

pars0 = zeros( nbetabasis + nalphabasis, 1); 

%pars0 = matrix(c(rep(res1.pars(1),nbetabasis), 
%                 rep(res1.pars(1),nalphabasis)), 
%                 nbetabasis + nalphabasis, 1) 

%% Profiled estimation step

control_out.trace  = 6;
control_out.reltol = 1e-6;

% This is the main optimization for profiled estimation.  
% It uses initial values for parameters and coefficients that
% are set up above.  This requires a lot of computation, and
% will take a number of minutes per iteration.  Find something
% else to do in the meantime.

lambda = 1e0;

%  estimate the parameters using differencing

res0 = Profile_LS(make_NS.fn, tobs, yobs, coefs0, pars0, basisobj, ... 
                  lambda, more, 'ProfileGN', control_out);

res0.pars;  %  display the parameter estimates

res0.outer.result.value;  % display the final function value


%  estimate the parameters using the derivative evaluation functions

lambda = 1e0;

res1 = Profile_LS(NSfn, tobs, yobs, coefs0, pars0, basisobj, lambda, ...
                  more, 'ProfileGN', control_out);

res1.pars;  %  display the parameter estimates

res1.outer.result.value;  % display the final function value

%  set up the functional data object for the fit

DEfd1 = fd(res1.coefs, basisobj);

%  increase lambda a few times, starting each time with 
%  parameters estimated from the last value

lambda = 1e2;

res2 = Profile_LS(NSfn, tobs, yobs, res1.coefs, res1.pars, basisobj, ...
                  lambda, elsemore, 'ProfileGN', control_out);

res2.pars  %  display the parameter estimates

DEfd2 = fd(res2.coefs, basisobj);

lambda = 1e4;

res3 = Profile_LS(NSfn, tobs, yobs, res2.coefs, res2.pars, basisobj, ...
                  lambda, elsemore, 'ProfileGN', control_out);

res3.pars;  %  display the parameter estimates

DEfd3 = fd(res3.coefs, basisobj);

lambda = 1e6;

res4 = Profile_LS(NSfn, tobs, yobs, res3.coefs, res3.pars, basisobj, ...
                   lambda, elsemore, 'ProfileGN', control_out);

res4.pars  %  display the parameter estimates

DEfd4 = fd(res4.coefs, basisobj);

lambda = 1e8;

res5 = Profile_LS(NSfn, tobs, yobs, res4.coefs, res4.pars, basisobj, ...
                  lambda, elsemore, 'ProfileGN', control_out);

res5.pars  %  display the parameter estimates

DEfd5 = fd(res5.coefs, basisobj);

%  compare fits to the data

subplot(1,1,1)
plotfit_fd(yobs, tobs, DEfd1) 
xlabel('Time (hours)') 
ylabel('Groundwater level (metres)')
lines(DEfd2)
lines(DEfd3)
lines(DEfd4)
lines(DEfd5)

%  By the time we get to res4, the fit to the data is essentially
%  a soluton to the differential equation.  We see that it has
%  rather less capacity to fit the data, but still captures some
%  of the shapes in the data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Henon map examples
%
% This example demonstrates the use of CollocInfer on a discrete-time system.
% The Henon Map is a classical dynamical system that exhibits chaos.
% It's equations are given as
%
% x(t+1) = 1 - a*x(t)^2 + y(t)
% y(t+1) = b*x(t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Data Generation  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hpars = [1.4,0.3];

ntimes = 200;

x = [-1,1];
X = zeros(ntimes+20,2);
X(1,:) = x;

for i = 2:(ntimes+20)
    X(i,:) = make_Henon.ode(i,X(i-1,:),hpars,[]); 
end

X = X(20+1:ntimes,:);

Y = X + 0.05*randn(ntimes,2);

t = 1:ntimes;

coefs = Y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Optimization Control  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

control.trace = 0;
control.maxit = 1000;
control.maxtry = 10;
control.reltol = 1e-6;
control.meth = [];

control_in = control;
control_in.reltol = 1e-12;

control_out = control;
control_out.trace = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Optimization       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hpars2 = [1.3,0.4];          % Perturbed parameters
lambda = 10000;

%%% SSE for discrete process%%%%

Ires1 = Smooth_LS(make_Henon,t,Y,coefs,hpars2,[], ...
                    lambda,[],control_in,0,0,1);
  
Ores1 = Profile_LS(make_Henon,t,Y,coefs,hpars2,[], ...
                   lambda,[],[],control_in,control_out,0,0,0,1);

%%% ProfileErr with SSEproc %%%%

profile_obj = LS_setup(make_Henon,[],lambda,coefs,[],0,0,1);
lik  = profile_obj.lik;
proc = profile_obj.proc;
   
Ires2 = inneropt(t,Y,hpars2,lik,proc,[],control_in);

Ores2 = outeropt(t,Y,coefs,hpars2,lik,proc,[],[],control_in,control_out);
    
%%% Dproc %%%%%%%%%%%%

var = [1,0.01];

Ires3 = Smooth_multinorm(make_Henon,Y,t,hpars2,coefs,[], ...
                         var,[],control_in,0,0,1);

Ores3 = Profile_multinorm(make_Henon,Y,t,hpars2,coefs, ...
            [],var,[],[],[],[], ...
            [],control_in,control_out,1e-6,[],0,0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Example Diagnostics -- Learning the FitzHugh-Nagumo Equations
%
% Demonstration estimation of function functions. This also provides
% a template for random-effects treatment within profiling.

%  add paths to required functions

addpath('../../fdaM')
addpath('fhn')
addpath('SSE')
addpath('id')

% First Create Some Data

t = 0:0.05:20;

pars = [0.2,0.2,3];

x0 = [-1,1];
y = ode45(@make_fhn.fn.ode,x0,t,[],pars);
y = y(:,2:3);

data = y + 0.2*randn(401,2);

% Now a basis object

knots = 0:0.2:20;
norder = 3;
nbasis = length(knots) + norder - 2;
range = [0,20];

bbasis = create_bspline_basis(range,nbasis,norder,knots);


% Initial values for coefficients will be obtained by smoothing

DEfd = data2fd(data,t,bbasis);
coefs = DEfd.coefs;

% Usual meta-parameters; quadrature points, weights and knots

lambda = [100,100];
qpts = knots;
qwts = rep(1/length(knots),length(knots));

qwts = qwts*t(lambda);
weights = ones(size(data));


% Now I define a measurement process log likelihood along with some
% additional features: in this case it's squared error. 

likmore = make_id;
likmore.weights = weights;

lik = make_SSElik;
lik.more = likmore;
lik.bvals = eval_basis(t,bbasis);

% Proc is a process log likelihood -- in this case treated as squared
% discrepancy from the ODE definition. 

procmore = make_genlin;
procmore.names = varnames;
procmore.parnames = parnames;
procmore.more.mat=zeros(2,2);
procmore.more.sub=[1,1,1;1,2,2;2,1,3;2,2,4];

procmore.weights = qwts;
procmore.qpts = qpts;

proc = make_SSEproc;
proc.more = procmore;
proc.bvals.bvals  = eval_basis(procmore.qpts,bbasis,0);
proc.bvals.dbvals = eval_basis(procmore.qpts,bbasis,1);

spars = [0,1,-1,0];

Ires = inneropt(t,data,spars,lik,proc,[]);

Ores = outeropt(t,data,Ires.coefs,spars,lik,proc);

traj  = proc.bvals.bvals  * Ores.coefs;
dtraj = proc.bvals.dbvals * Ores.coefs;
ftraj = dtraj - ...
    proc.more.fn(proc.more.qpts,dtraj,Ores.pars,proc.more.more);

m = 0;
for i = 1:2
    for j = 1:2
        m = m + 1;
        subplot(2,2,m)
        plot(traj(:,i),ftraj(:,j))
    end
end

%% Now we estimate some forcing functions

fbasis = create_bspline_basis(range,23,4);

dproc = make_SSEproc;
dproc.more = make_diagnostics;
dproc.more.qpts = procmore.qpts;
dproc.more.weights = procmore.weights;
dproc.more.more = procmore;
dproc.more.more.p = Ores.pars;
dproc.more.more.which = 1:2;
dproc.more.more.psi = eval_basis(procmore.qpts,fbasis);
dproc.bvals.bvals   = eval_basis(procmore.qpts,bbasis,0);
dproc.bvals.dbvals  = eval_basis(procmore.qpts,bbasis,1);

dpars = zeros(2*fbasis.nbasis,1);  
  
dOres = outeropt(t,data,Ires.coefs,dpars,lik,dproc); 

% Trajectories

force = dproc.more.more.psi * reshape(dOres.par,fbasis.nbasis,2);
traj  = dproc.bvals.bvals   * dOres.coefs;

plot(traj(:,1),force(:,1))



