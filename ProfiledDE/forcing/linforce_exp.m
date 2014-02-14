odefn       = @genlinode;        % Function for ODE solver (exact)

fn.fn       = @forcingfun;        % RHS function

fn.dfdx     = @forcingdfdx;       % Derivative wrt inputs (Jacobian)
fn.dfdp     = @forcingdfdp;       % Derviative wrt parameters

fn.d2fdx2   = @forcingd2fdx2;     % Hessian wrt inputs
fn.d2fdxdp  = @forcingd2fdxdp;    % Cross derivatives wrt inputs and parameters
fn.d2fdp2   = @forcingd2fdp2;     % Hessian wrt parameters

fn.d3fdx2dp = @forcingd3fdx2dp;   % Third derivative wrt inputs, inputs, pars 
fn.d3fdx3   = @forcingd3fdx3;     % Third derivative wrt inputs
fn.d3fdxdp2 = @forcingd3fdxdp2;   % Third derivative wrt inputs, pars and pars

fn_extras.fn      = @genlinfun;    % Original function
fn_extras.dfdx    = @genlindfdx;   % First derivative
fn_extras.d2fdx2  = @genlind2fdx2; % Second derivative
fn_extras.d3fdx3  = @genlind3fdx3; % Third derivative

fn_extras.extras  = [];        % Original information to fn_extras.fn. 

more = [];

y0 = [-1,1];                      
pars = [-1; 2; -1; 1]; 
A = reshape(pars,2,2);


sigma = 0.5;                    

jitter = 0.2;                   
nrep = 50;


tspan = 0:0.05:20;  
obs_pts{1} = 1:length(tspan);       
obs_pts{2} = 1:length(tspan);      

lambdas = 10^(-2:4);     
lambda0 = 1;        

nknots = 401;
norder = 3;
nquad = 5;     

range = [min(tspan),max(tspan)];

knots_cell = cell(1,length(y0));
knots_cell(:) = {linspace(range(1),range(2),401)};

basis_cell = cell(1,length(y0)); 
Lfd_cell = cell(1,length(y0));

nbasis = zeros(length(y0),1);

bigknots = knots_cell{1};               
nbasis(1) = length(knots_cell{1}) + norder - 2;          

for(i = 2:length(y0))
    bigknots = [bigknots knots_cell{i}];
    nbasis(i) = length(knots_cell{i}) + norder -2;
end

quadvals = MakeQuadPoints(bigknots,nquad);   

for(i = 1:length(y0))
    basis_cell{i} = MakeBasis(range,nbasis(i),norder,...  
        knots_cell{i},quadvals,1);                        
    Lfd_cell{i} = fdPar(basis_cell{i},1,lambda0);         
end

odeopts = odeset('RelTol',1e-13);
[full_time,full_path] = ode45(odefn,tspan,y0,odeopts,[pars; 1; 0],more);


for(g = 1:nrep)

    Ycell = path;
    for(i = 1:length(path))
        Ycell{i} = path{i} + sigma*randn(size(path{i}));
    end

    [DEfd,forces] = linforceest(basis_cell,basis_cell,A,1:2,10000,...
        0.0001,2,Tcell,Ycell,wts);
    
    