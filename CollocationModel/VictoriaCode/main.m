%{
This program sets up and calls the functions that
solve a simple menu cost model (taken from Terry 
(2008)).

Written by:     Victoria Gregory
Date:           2/26/2016
%}

clear all;
clc;
dbstop if error;
cd '/Users/victoriagregory/Dropbox/MenuCostsModel/CollocationModel/VictoriaCode'

%% Settings

% What to solve for
options.solvepL     = 'N';      % Solve for p and L given a Y 
options.solveeq     = 'N';      % Solve equilibrium
options.solveKS     = 'Y';      % Solve Krussel-Smith

% Tolerances, iterations
options.Nbell       = 2;        % Number of Bellman (Contraction) iterations
options.Nnewt       = 15;       % Maximum number of Newton steps
options.tolc        = 1e-8;     % Tolerance on value functions
options.tolgolden   = 1e-6;     % Tolerance for golden search
options.itermaxL    = 5000;     % Maximum iterations to find stationary dist L
options.tolL        = 1e-11;    % Tolerance on L

% Set-up for state space
glob.n          = [30,5];        % Number of nodes in each dimension
glob.nf         = [300,5];    % Number of points for p and a in histogram L
glob.curv       = 1;            % Grid curvature for p/P on (0,1] (1 is no curvature)
glob.spliorder  = [3,1];        % Order of splines (always use linear if shocks are discrete (not AR1))
glob.pmin       = 0.75;         % Lower bound on p
glob.pmax       = 1.50;         % Upper bound on p

% NOTE (VG): resulting k grid will be n(1)+spliorder(1)-1
% Creating the cubic spline space adds 3-1=2 points.

% Model parameters
param.beta      = 0.99;     % discount factor
param.delta     = 0.3;      % relative weight of labor-to-consumption in utility
param.sigma     = 1;        % risk aversion coefficient
param.phi       = 0.5;      % inveser labour supply elasticity
param.theta     = 5;        % elasticity of substitution
param.alpha     = 2/3;      % returns to labour
param.rhoa      = 0.35;     % persistence of productivity
param.sigmazeta = 0.225;    % stddev of productivity shocks
param.Phi       = 0.156;    % menu cost in labour units
param.mu        = 0.006;    % s.s money growth
param.rhom       = 0.37;     % persistence of money growth
param.sigmaeps  = 0.0048;   % stddev of money growth shocks
param.tauc      = 0.005;    % tolerance for forecasting rule
param.n         = 5000;     % number of firms
param.T         = 96;       % simulation length
param.S         = 25;       % simulations for computing forecasting coeffs
param.s         = 100;      % simulations for moment computations

% Print / plot 
options.print       = 'Y';      % Print out c-solution convergence

%% Setup problem
fprintf('Setup\n');
[param,glob]    = setup_ss(param,glob,options);      
fprintf('Setup complete\n');


%% Solve only p and L for a given output Y
switch options.solvepL
    case 'Y'
        Y                   = 1.01;    % Conjectured value of Y    
        options.cresult     = [];   % Holds previous solution for c. Empty in this case.
        eq                  = solve_pL(Y,param,glob,options);  
        fprintf('Yin = %1.2f,\tYout = %1.2f\n',Y,eq.Y);
end

%plot(glob.sf(1:50,1)./(eq.Pa),eq.v.vf(1:50))
% out=funbas(glob.fspace,glob.sf)*eq.c;
% plot(glob.sf(350:400,1)./(eq.Pa),eq.v.vf(350:400),glob.sf(350:400,1)./(eq.Pa),out(350:400))
% legend('RHS','LHS')

%% Solve equilibrium
switch options.solveeq
    case 'Y'
        options.tolY        = 0.0001;           % Tolerance on output
        options.Ylb         = 0.1;              % Output lower bound
        options.Yub         = 5;               % Output upper boud
        options.itermaxY    = 30;               % Max iterations of bisection
        options.eqplot      = 'Y'; 
        options.eqprint     = 'Y'; 
        options.print       = 'N';
        options.Loadc       = 'Y';              % For new guess of p use old c as starting guess
        options.plotSD      = 'N';              % If Y plot steady state distribution
        eq                  = solve_eq_menucost(param,glob,options); 
end

%% Set up for Krussel-Smith

% State space
glob.n          = [20,5,3,3];       % Number of nodes in each dimension
glob.nf         = [300,8,6,6];      % Number of points for p and a in histogram L
glob.curv       = 1;                % Grid curvature for p/P on (0,1] (1 is no curvature)
glob.spliorder  = [3,1,1,1];        % Order of splines (always use linear if shocks are discrete (not AR1))

% Law of motion - initial guesses
cKS.b0     = 0.001;
cKS.b1     = 0.5;
cKS.b2     = 0.1;

%% Setup problem
fprintf('Setup\n');
[param,glob]    = setup_ks(cKS,param,glob,options);      
fprintf('Setup complete\n');

%% For now, just going to solve the problem on the expanded state space

switch options.solveKS
    case 'Y'
    options.cresult        = [];   % Holds previous solution for c. Empty in this case.
    [c,v]                  = solve_KS(cKS,param,glob,options);
end

outc=funbas(glob.fspace,glob.s)*c(end/3+1:2*end/3);
outk=funbas(glob.fspace,glob.s)*c(1:end/3);
vf = max(v.vk,v.vc);
vfs = max(outc,outk);
plot(glob.pgrid,vf(1:glob.n(1)+2),glob.pgrid,vfs(1:glob.n(1)+2));

