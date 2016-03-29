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
options.solvepL     = 'Y';      % Solve for p and L given a Y 
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
glob.nf         = [300,5];      % Number of points for p and nu in histogram L
glob.curv       = 1;            % Grid curvature for p on (0,1] (1 is no curvature)
glob.spliorder  = [3,1];        % Order of splines (always use linear if shocks are discrete (not AR1))
glob.pmin       = 0.75;         % Lower bound on p
glob.pmax       = 1.50;         % Upper bound on p

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
