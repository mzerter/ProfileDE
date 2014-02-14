% Store Current Folder so that it can be reset at the end of this script.
clear;
currentFolder = pwd;




% CollocInfer Path 

% This should be changed!! :
cd('/Users/Jim/Documents/MATLAB/ProfiledDE/CollocInfer/')





%  add paths to required functions

addpath('../../fdaM')


addpath('SEIR')
addpath('SSE')
addpath('id')
addpath('exp')
addpath('findif')
addpath('logtrans')
addpath('logstate_lik')
addpath('../')

% Get some data and parameters

%  data

fid  = fopen('SEIRdata.txt','rt');
SEIRdata = fscanf(fid,'%f');
fclose(fid);

nobs = size(SEIRdata,1);




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

SEIRfdPar = fdPar(bbasis, 1, 0.1);

SEIR_xfd = smooth_basis(SEIRtimes, logdata(:,3), SEIRfdPar, ones(nobs,1));

C = getcoef(SEIR_xfd);


coefs = zeros(length(C),3);

coefs(:,1) = 13;
coefs(:,2) = 7;
coefs(:,3) = C;

DEfd = fd(coefs,bbasis);


SEIRlambda = ones(1,3);
SEIRlambda(1) = 100;


active = [2, 6, 7, 8];

%  struct object betamore is now passed into LS.setup in order to make
%  it available as a functional parameter defined by its three coefficients

%  run LS.setup

[SEIRlik, SEIRproc] = LS_setup(SEIRfn, SEIRtimes, [], [], SEIRlambda, ...
                               DEfd, betamore, [], [], [], 0, 1);
                  

%  initialize the inner optimization 

control_in = optimset('LargeScale', 'on', 'GradObj', 'on', ...
                      'Hessian', 'on', 'Diagnostics', 'off', ...
                      'Display', 'off');
                  
control_out = optimset('Jacobian', 'on', 'MaxIter', 100, ...
                       'Display', 'iter');

                   
global INNEROPT_COEFS0;                   
INNEROPT_COEFS0 = coefs; 

tic;
Profile_LS_struct = ...
    Profile_LS(SEIRfn, SEIRtimes, logdata, coefs, SEIRpars,bbasis, SEIRlambda, ...
               DEfd, betamore, [], [],active , ...
               [], control_in, [], control_out, [], 0, 1);
toc;


% RESULTS:

%  Iteration  Func-count     f(x)          step          optimality   CG-iterations
%      0          1         4.10265                      1.63e+03
%      1          2         4.09148             10       1.51e+03            0
%      2          3         4.08907              1       1.49e+03            0
%      3          4         4.08709              1       1.47e+03            0
%      4          5         4.08534              1       1.46e+03            0
%      5          6         4.08354              1       1.44e+03            0
%      6          7         4.08185              1       1.42e+03            0
%      7          8         4.08025              1       1.41e+03            0
%      8          9         4.07869              1       1.39e+03            0
%      9         10         4.07745              1       1.38e+03            0
%     10         11         4.07601              1       1.36e+03            0
%     11         12         4.07459              1       1.35e+03            0
%     12         13         4.07321              1       1.33e+03            0
%     13         14         4.07187              1       1.32e+03            0
%     14         15         4.07058              1        1.3e+03            0
%     15         16         4.06933              1       1.29e+03            0
%     16         17         4.06813              1       1.27e+03            0
%     17         18         4.06697              1       1.26e+03            0
%     18         19         4.06586              1       1.25e+03            0
%     19         20         4.06479              1       1.23e+03            0
%     20         21         4.06377              1       1.22e+03            0
%     21         22         4.06259              1       1.21e+03            0
%     22         23         4.06185              1       1.19e+03            0
%     23         24          4.0608              1       1.18e+03            0
%     24         25         4.05996              1       1.17e+03            0
%     25         26         4.05918              1       1.16e+03            0
%     26         27         4.05844              1       1.14e+03            0
%     27         28         4.05775              1       1.13e+03            0
%     28         29          4.0571              1       1.12e+03            0
%     29         30         4.05648       0.954582       1.11e+03            0
%     30         31         4.05585       0.626954        1.1e+03            0
%     31         32         4.05525        0.39014       1.09e+03            0
%     32         33         4.05457       0.216111       1.09e+03            0
%     33         34         4.05388       0.098805       1.08e+03            0
%     34         35         4.05318      0.0189606       1.08e+03            0
%     35         36         4.05247      0.0364387       1.07e+03            0
%     36         37         4.05176       0.074983       1.07e+03            0
%     37         38         4.05105       0.102011       1.07e+03            0
%     38         39         4.05035       0.120081       1.07e+03            0
%     39         40         4.04964       0.133605       1.06e+03            0
%     40         41         4.04895       0.141227       1.06e+03            0
%     41         42         4.04825       0.146187       1.06e+03            0
%     42         43         4.04756       0.151713       1.05e+03            0
%     43         44         4.04688       0.153924       1.05e+03            0
%     44         45          4.0462       0.153834       1.05e+03            0
%     45         46         4.04552       0.156384       1.05e+03            0
%     46         47         4.04485       0.156353       1.04e+03            0
%     47         48         4.04419       0.156033       1.04e+03            0
%     48         49         4.04353       0.155824       1.04e+03            0
%     49         50         4.04288       0.155532       1.04e+03            0
%     50         51         4.04223       0.155343       1.03e+03            0
%     51         52         4.04158       0.154593       1.03e+03            0
%     52         53         4.04094       0.154141       1.03e+03            0
%     53         54         4.04031       0.151148       1.03e+03            0
%     54         55         4.03968       0.151719       1.03e+03            0
%     55         56         4.03906       0.151417       1.02e+03            0
%     56         57         4.03844       0.150935       1.02e+03            0
%     57         58         4.03783       0.150391       1.02e+03            0
%     58         59         4.03722       0.149814       1.02e+03            0
%     59         60         4.03661        0.14921       1.01e+03            0
%     60         61         4.03601       0.148582       1.01e+03            0
%     61         62         4.03542       0.147168       1.01e+03            0
%     62         63         4.03483       0.145343       1.01e+03            0
%     63         64         4.03424       0.147829       1.01e+03            0
%     64         65          4.0337       0.143972          1e+03            0
%     65         66         4.03309       0.145307          1e+03            0
%     66         67         4.03255        0.14421            998            0
%     67         68         4.03195       0.144506            996            0
%     68         69         4.03142       0.141905            994            0
%     69         70         4.03083       0.141376            992            0
%     70         71          4.0303       0.141728            989            0
%     71         72         4.02972       0.140075            988            0
%     72         73          4.0292       0.138438            985            0
%     73         74         4.02866       0.137587            983            0
%     74         75         4.02813       0.140023            981            0
%     75         76         4.02759       0.138825            979            0
%     76         77         4.02706       0.135922            977            0
%     77         78         4.02653       0.135391            975            0
%     78         79         4.02601       0.134917            973            0
%     79         80         4.02549       0.133139            971            0
%     80         81         4.02498       0.133386            969            0
%     81         82         4.02447       0.133026            967            0
%     82         83         4.02396       0.134589            965            0
%     83         84         4.02346       0.131277            963            0
%     84         85         4.02296       0.129353            961            0
%     85         86         4.02246       0.130787            959            0
%     86         87         4.02197        0.12997            957            0
%     87         88         4.02148       0.128834            955            0
%     88         89           4.021       0.128693            953            0
%     89         90         4.02052       0.128264            951            0
%     90         91         4.02004       0.127787            950            0
%     91         92         4.01956       0.127285            948            0
%     92         93         4.01909       0.124931            946            0
%     93         94         4.01862       0.125328            944            0
%     94         95         4.01816       0.124992            942            0
%     95         96          4.0177       0.122785            940            0
%     96         97         4.01724       0.123216            938            0
%     97         98         4.01679       0.122931            936            0
%     98         99         4.01633       0.124329            935            0
%     99        100         4.01589       0.121239            933            0
%    100        101         4.01544       0.121508            931            0
%    101        102         4.01496       0.120562            930            0
% 
% Solver stopped prematurely.
% 
% lsqnonlin stopped because it exceeded the iteration limit,
% options.MaxIter = 100 (the selected value).
% 
% Elapsed time is 42.721555 seconds.


% Profile_LS_struct.pars(active)
% 
% ans =
% 
%           30.1713988506761      0.000211496315654376      8.10022457089302e-07      1.66037379288473e-05
% 
% 



format long g;

cd(currentFolder);
