%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       FitzHugh-Nagumo Example       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  add paths to required functions

addpath('../../fdaM')
addpath('fhn')
addpath('SSE')
addpath('id')

% Obtain some pre-generated data

load FHN_data.txt

% Now we want to set up a basis expansion; this makes use of
% routines in the fda library

knots  = (0:0.5:20)';
norder = 3;
nbasis = length(knots) + norder - 2;
range  = [0,20];

bbasis = create_bspline_basis(range, nbasis, norder, knots);
	
% We'll start off by creating a smooth of each state variable to get initial
% values for the parameters.

bfdPar = fdPar(bbasis,int2Lfd(1),1);
fdnames = cell(1,3);                                  
fdnames{3} = ['V','R'];
DEfd = smooth_basis(FhNtimes,FhNdata,bfdPar,fdnames);

coefs = getcoefs(DEfd);

%%%% Option 1: pre-defined set-up for squared error

lambda = 1000;

% First just optimize the coefficients

IDEfd1 = Smooth_LS(@make_fhn, FhNdata, FhNtimes, FhNpars,...
                   coefs, bbasis, lambda, ...
                   [], [], [], [], []);

% Let's have a look at this

plotfit_fd(FhNtimes, FhNdata, IDEfd1);

% Now we can do the profiling

fhn = make_fhn;
LS_struct1 = Profile_LS(fd,FhNtimes,FhNdata,coefs,FhNpars, ...
                         bbasis,lambda)

% And look at the result

ODEfd1 = fd(LS_struct1.coefs,bbasis);
plotfit_fd(FhNtimes,FhNdata,ODEfd1)
  
% If we only coded the FitzHugh Nagumo funciton with no derivatives we would just
% have make_fhn.fn, in which case we default to estimating the derivatives by
% finite differencing

Ires1a	= Smooth_LS(fhn.fn,FhNdata,FhNtimes,FhNpars,coefs,bbasis,lambda,
                    []);

% Now we can do the profiling

Ores1a = Profile_LS(fhn.fn,FhNdata,FhNtimes,FhNpars,coefs, ...
                     bbasis,lambda);
  
%%%% Option 2:  set-up functions for proc and lik objects and then call
% optimizers

profile_obj = LS.setup(FhNpars,fhn,lambda,FhNtimes, ...
                       coefs,bbasis);
lik  = profile_obj.lik;
proc = profile_obj.proc;

% Now we can get initial parameter estimates from "gradient matching"

pres = ParsMatchOpt(FhNpars,coefs,proc);
npars = pres.pars;

% Smoothing can be done more specifically with 'inneropt'

Ires2 = inneropt(FhNdata,FhNtimes,npars,coefs,lik,proc);

% And we can also do profiling

Ores2 = outeropt(FhNdata,FhNtimes,npars,coefs,lik,proc);

%%%% Option 3:  set everything up manually

%% lik object

lik2 = make_SSElik;                      % we are using squared error
lik2.bvals = eval_basis(FhNtimes,bbasis);  % values of the basis at observation times
lik2.more = make_id;                     % identity transformation
lik2.more.weights = [1,0];                % only use data for V

%% proc object

qpts = knots(1:(length(knots)-1))+diff(knots)/2;  % Collocation points at midpoints
                                                 % between knots

proc2 = make_SSEproc                    % we are using squared error
proc2.bvals.bvals  = eval_basis(qpts,bbasis);    % values and derivative of
proc2.bvals.dbvals = eval_basis(qpts,bbasis,1));  % basis at collocation points
proc2.more = make_fhn;                   % FitzHugh-Nagumo right hand side
proc2.more.names = FhNvarnames;            % State variable names
proc2.more.parnames = FhNparnames;         % Parameter names
proc2.more.qpts = qpts;                    % Collocation points
proc2.more.weights  = [1000,1000];        % Weight relative to observations of
                                          % each collocation point.
                                          
% We can also specify optimization methods

in.meth = [];            % Inner Optimization
control.in.reltol=1e-8;
control.in.maxit=100;
control.in.trace=0;  % Optimization control

out.meth = [];            % Outer Optimization
control.out.trace=6;
control.out.maxit=100;
control.out.reltol=1e-8; % Optimization control
                                          
                                          
% Now we will also remove the 'R' part of the data (ie, assume we did not measure
% it) and set the corresponding coefficients to zero

new.FhNdata = FhNdata;
new.FhNdata(:,2) = nan;

new.coefs = coefs;
new.coefs(:,2) = 0;

% We can now fix the smooth or 'V' and try and choose 'R' to make the differential
% equation match as well as possible.

fres = FitMatchOpt(new.coefs,2,FhNpars,proc2);
new.coefs2 = fres.coefs;

% And we can call the same inner optimization as above, this time with our
% own control parameters and method

Ires3 = inneropt(new.FhNdata,FhNtimes,FhNpars,new.coefs2,lik2,proc2, ...
                 in.meth,control.in);

% And we can also do profiling with specified controls

Ores3 = outeropt(new.FhNdata,FhNtimes,FhNpars,new.coefs2,lik2,proc2,
  in.meth,out.meth,control.in,control.out);

% We can also use finite differencing to calculate derivatives of the right hand side of the
% FitzHugh-Nagumo equations, this is defined by modifying the proc object.

proc2a = make_SSEproc                    % we are using squared error
proc2a.bvals.bvals  = eval_basis(qpts,bbasis);    % values and derivative of
proc2a.bvals.dbvals = eval_basis(qpts,bbasis,1));  % basis at collocation points
proc2a.more = make_findif.ode;                   % FitzHugh-Nagumo right hand side
proc2a.more.names = FhNvarnames;            % State variable names
proc2a.more.parnames = FhNparnames;         % Parameter names
proc2a.more.qpts = qpts;                    % Collocation points
proc2a.more.weights  = [1000,1000];        % Weight relative to observations of
                                           % each collocation point.
proc2a.more.more.fn = make_fhn.fn; 
proc2a.more.more.eps=1e-6; % Tell findif that the function
                                                      % to difference is fhn and the
                                                      % difference stepsize


Ires3a = inneropt(new.FhNdata,FhNtimes,FhNpars,new.coefs2,lik2,proc2a,...
                  in.meth,control.in);

% And we can also do profiling with specified controls

Ores3a = outeropt(new.FhNdata,FhNtimes,FhNpars,new.coefs2, ...
  lik2,proc2a,in.meth,out.meth,control.in,control.out);

