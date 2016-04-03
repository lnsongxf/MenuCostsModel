function [c,v,cKS,R2,sim] = solve_KS(cKS,eq,param,glob,options)
%SOLVE_KS Implements the Krussel Smith Algorithm for the model
%-------------------------------------------------
%   Simulates the aggregate states, solves the model at each step, 
%   
%   INPUTS
%   - cKS       = parameters for the law of motion (lnY = b0 + b1lnY-1 + b2 Dm)
%   - eq        = ???
%   - param     = model parameters
%   - glob      = global variables, including those used in the approximating functions
%   - options   = options, including continuous vs. discrete Markov, Tauchen vs. Rouwenhurst, 
%
%   OUTPUTS
%   - At        = Simulated path for aggregate productivity, A
%   - pt        = Simulated path for price level, p
%   - Kt        = Simulated path for aggregate capital, K
%-------------------------------------------------
%
%   NEED TO FIX THE TIMING OF THE MONEY SHOCKS!!!!
%
%
%% Run setup file again with new cKS params
cKS0 = cKS;
glob = setup_agg(param,glob,cKS0,options);


%% Unpack
ns          = size(glob.s,1); 
T           = options.T;

%% Initialise guesses (if val.cresult has an old guess in it, use that)
cKold       = zeros(ns,1);
cCold       = zeros(ns,1);
cEold       = zeros(ns,1);
cold        = [cKold;cCold;cEold];


%% Solve value function problem
% Define equilibrium output and price
glob.P      = eq.P;     %%%% DEBUG
glob.Y      = eq.Y;     %%%% DEBUG
[c,v]       = solve_cKS(cold,cKS0,param,glob,options);   % Get new val func approx coefficients

%% Set up simulations

cKS      = zeros(3,options.KSsim);
% rng(219);
for s=1:options.KSsim
    [cKS(:,s),R2(s,1),sim] = simulate_KS(c,v,eq,param,glob,options);
    s
end
 
    
end














