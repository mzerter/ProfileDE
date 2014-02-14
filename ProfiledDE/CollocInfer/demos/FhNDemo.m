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

control_out = optimset('GradObj', 'on', ...
                       'Hessian', 'on', ...
                       'MaxIter', 100, ...
                       'display', 'iter', ...
                       'TolFun',  1e-6);

npars = outeropt(FhNtimes, FhNdata, FhNcoefs0, FhNpars, lik, proc, ...
                 [], [], [], [], control_out);

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
                       'Hessian', 'off', ...
                       'MaxIter', 100, ...
                       'display', 'iter', ...
                       'TolFun',  1e-6);

npars = outeropt(FhNtimes, new.FhNdata, new.coefs3, FhNpars, lik2, ...
                 proc2, [], [], [], method_out, control_out);

disp(['Parameters: ', num2str(npars)])

% We can also use finite differencing to calculate derivatives of the 
% right hand side of the FitzHugh-Nagumo equations, this is defined by 
% modifying the proc object.

proc2a               = make_SSEproc;                    
proc2a.bvals.bvals   = eval_basis(qpts,bbasis);    
proc2a.bvals.dbvals  = eval_basis(qpts,bbasis,1);  
proc2a.more          = make_findif_ode;                   
proc2a.more.names    = FhNvarnames;            
proc2a.more.parnames = FhNparnames;         
proc2a.more.qpts     = qpts;                    
proc2a.more.weights  = [1000,1000]; 

Fhnfn                   = make_fhn;                           
proc2a.more.more.fn     = Fhnfn.fn; 
proc2a.more.more.reltol = 1e-6; 
proc2a.more.more.more   = [];                        
proc2a.more.more.eps    = 1e-3;                                                      

in_method = [];
coefs3 = inneropt(FhNtimes,new.FhNdata,FhNpars,lik2,proc2a,...
                      in_method,control_in);

% And we can also do profiling with specified controls

method_out = [];
control_out = optimset('GradObj', 'on', ...
                       'Hessian', 'off', ...
                       'MaxIter', 100, ...
                       'display', 'iter', ...
                       'TolFun',  1e-6);
                   
nparsa = outeropt(FhNtimes,new.FhNdata,new.coefs3,FhNpars, ...
                  lik2,proc2a, [], [], [], out_method,control_out);
              
disp(['Parameters: ', num2str(nparsa)])
