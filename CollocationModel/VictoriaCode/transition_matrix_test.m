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

% Law of motion
b0     = 0.015;
b1     = 0.3;
b2     = 0.25;

% Sizes of grids
Na     = 5;
Ndm    = 5;
Ny     = 5;

% Folding a, y, Dm into a single var to get entire state space in one go
Nvar = [Na; ...
    Ndm; ...
    Ny];
muvar = [0; ...
    param.mu*(1-param.rhom); ...
    b0 + b2*param.mu*(1-param.rhom)];
Avar = [param.rhoa, 0, 0; ...
    0, param.rhom, 0; ...
    0, b2*param.rhom, b1];
Svar = [param.sigmazeta; ...
    param.sigmaeps; ...
    b2*param.sigmaeps];

% Create grid and transition matrix for A, delta M, Y
[grid,Transition_aydm,~]=tauchenvar(Nvar,muvar,Avar,Svar);
grid = exp(grid);
