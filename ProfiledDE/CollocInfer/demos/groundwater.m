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

fid  = fopen('NSdata.txt','rt');
NSdata = reshape(fscanf(fid,'%f'),3,315)';
fclose(fid);

yobs = fscanf(fid,'%f');
zobs = fscanf(fid,'%f');
tobs = 0:314;

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

