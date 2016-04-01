%{
This program sets up and calls the functions that
solve the Golosov-Lucas model.

Written by:     Victoria Gregory
Date:           3/28/2016
%}

clear all;
clc;
dbstop if error;
cd '/Users/victoriagregory/Dropbox/MenuCostsModel/CollocationModel/VictoriaCode/Golosov-Lucas'

%% Settings

% What to solve for
options.solvexL     = 'Y';      % Solve for p and L given a Y 
options.solveeq     = 'Y';      % Solve equilibrium
options.solveKS     = 'Y';      % Solve Krussel-Smith

% Tolerances, iterations
options.Nbell       = 2;        % Number of Bellman (Contraction) iterations
options.Nnewt       = 15;       % Maximum number of Newton steps
options.tolc        = 1e-8;     % Tolerance on value functions
options.tolgolden   = 1e-6;     % Tolerance for golden search
options.itermaxL    = 5000;     % Maximum iterations to find stationary dist L
options.tolL        = 1e-11;    % Tolerance on L

% Set-up for state space
glob.n          = [10,5];       % Number of nodes in each dimension
glob.nf         = [300,5];      % Number of points for x and nu in histogram L
glob.curv       = 1;            % Grid curvature for x on (0,1] (1 is no curvature)
glob.spliorder  = [3,1];        % Order of splines (always use linear if shocks are discrete (not AR1))
glob.xmin       = exp(-0.4);    % Lower bound on x
glob.xmax       = exp(0.5);     % Upper bound on x

% NOTE (VG): resulting k grid will be n(1)+spliorder(1)-1
% Creating the cubic spline space adds 3-1=2 points.

% Model parameters
param.rho       = 0.04;
param.gamma     = 2;
param.epsilon   = 7;
param.alpha     = 6;
param.eta       = 0.55;
param.sigmanu   = sqrt(0.011);
param.k         = 0.0025;
param.mu        = 0.0064;
param.sigmam    = sqrt(0);
param.R         = param.rho + param.mu;
param.beta      = exp(-param.R);

% Print / plot 
options.print       = 'Y';      % Print out c-solution convergence

%% Setup problem
fprintf('Setup\n');
[param,glob]    = setup_ss(param,glob,options);      
fprintf('Setup complete\n');

%% Solve only x and L for a given cbar
switch options.solvexL
    case 'Y'
        cbar                = 0.35;      % Conjectured value of cbar    
        options.cresult     = [];     % Holds previous solution for c. Empty in this case.
        eq                  = solve_xL(cbar,param,glob,options);  
        fprintf('cbarin = %1.2f,\tcbarout = %1.2f\n',cbar,eq.cbar);
end

%% Solve equilibrium
switch options.solveeq
    case 'Y'
        options.tolcbar     = 0.0001;           % Tolerance on cbar
        options.cbarlb      = 0.1;              % cbar lower bound
        options.cbarub      = 1;                % cbar upper boud
        options.itermaxcbar = 30;               % Max iterations of bisection
        options.eqplot      = 'Y'; 
        options.eqprint     = 'Y'; 
        options.print       = 'N';
        options.Loadc       = 'Y';              % For new guess of cbar use old c as starting guess
        options.plotSD      = 'N';              % If Y plot steady state distribution
        eq                  = solve_eq(param,glob,options); 
end
