%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Spread of disease SIR Example       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Last modified 12 August 2013

%  add paths to required functions, assuming this code is executed
%  in the CollocInfer folder in ProfiledDE

addpath('../../fdaM')
addpath('SIR')
addpath('SSE')
addpath('id')

%  set time span in weeks

tmax  = 5;
tspan = [0,tmax];

% Obtain some pre-generated data

load SIRdata
nobs = size(SIRdata,1);
SIRtimes = linspace(0,tmax,nobs)';

% Now we want to set up a basis expansion; this makes use of
% routines in the fda library

knots  = (0:0.25:tmax)';
norder = 4;
nbasis = length(knots) + norder - 2;
SIRrng = [0,tmax];

SIRbasis = create_bspline_basis(SIRrng, nbasis, norder, knots);
	
lambda = 1e-4;
SIRfdPar = fdPar(SIRbasis, 2, lambda);

%  smooth the data using FDA function smooth_basis

SIRfdnames = cell(1,3);                                  
SIRfdnames{3} = cell(1,2);
SIRfdnames{3}{1} = 'SIR Variables';
SIRfdnames{3}{2} = cell(1,3);
SIRfdnames{3}{2}{1} = 'S';
SIRfdnames{3}{2}{2} = 'I';
SIRfdnames{3}{2}{3} = 'R';

SIR_xfd = smooth_basis(SIRtimes, SIRdata, SIRfdPar, ones(nobs,1), ...
                       SIRfdnames);
                   
%  set the initial coefficints from this smooth as well as the 
%  values of lambda for ProfiledDE
                   
SIRcoefs0 = getcoef(SIR_xfd);

%  plot the fit arising from the initial smooth

SIRfd0 = fd(SIRcoefs0, SIRbasis);
plotfit_fd(SIRdata, SIRtimes, SIRfd0)

%  define initial parameters

N = 1e3;
R0 = 5;
gammaval = 1e0;  %  recovery rate
muval    = 0e0;  %  change in susceptible population
betaval  = (gammaval + muval)*R0/N;  %  infection rate
SIRpars0 = [betaval, gammaval, muval];

% Make the SIR functions

fn = make_SIR;

%  set lambda values

SIRlambda = 1e0*[1,1,1]';

%  set global parameter containing coefficienjts

global INNEROPT_COEFS0
INNEROPT_COEFS0 = SIRcoefs0; 

%  Now do the profiling estimation

SIRstruct = Profile_LS(fn,SIRtimes,SIRdata,SIRcoefs0,SIRpars0, ...
                          SIRbasis,SIRlambda);

SIRcoefs = SIRstruct.coefs; 
SIRpars  = SIRstruct.pars;

disp(['Parameter estimates: ',num2str(SIRpars)])

% And look at the result

SIRfd = fd(SIRcoefs,SIRbasis);
plotfit_fd(SIRdata,SIRtimes,SIRfd)
  
%  try estimation without the first variable

SIRdata1 = SIRdata;
SIRdata1(:,1) = NaN;
SIRdata1( 1: 3,1) = SIRdata( 1: 3,1);
SIRdata1(19:21,1) = SIRdata(19:21,1);
SIRcoefs1 = SIRcoefs0;
SIRpars1  = SIRpars0;

SIRstruct1 = Profile_LS(fn,SIRtimes,SIRdata1,SIRcoefs1,SIRpars1, ...
                        SIRbasis,SIRlambda);
 
SIRcoefs1 = SIRstruct1.coefs;
SIRpars1  = SIRstruct1.pars;
 
disp(['Parameter estimates: ',num2str(SIRpars1)])

% And look at the result

SIRfd1 = fd(SIRcoefs1,SIRbasis,SIRfdnames);
plotfit_fd(SIRdata,SIRtimes,SIRfd1)

%  try estimation without the first and last variable

SIRdata2 = SIRdata1;
SIRdata2(:,3) = NaN;

SIRstruct2 = Profile_LS(fn,SIRtimes,SIRdata2,SIRcoefs0,SIRpars0, ...
                          SIRbasis,SIRlambda);
 
SIRcoefs2 = SIRstruct2.coefs;
SIRpars2  = SIRstruct2.pars;

disp(['Parameter estimates: ',num2str(SIRpars2)])

% And look at the result

SIRfd2 = fd(SIRcoefs2,SIRbasis,SIRfdnames);
plotfit_fd(SIRdata,SIRtimes,SIRfd2)

