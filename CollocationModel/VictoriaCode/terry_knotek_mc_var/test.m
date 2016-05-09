clear all;
clc
% Model parameters
param.beta      = 0.99;     % discount factor
param.delta     = 0.5;      % relative weight of labor-to-consumption in utility
param.sigma     = 1;        % risk aversion coefficient
param.phi       = 0.5;      % inveser labour supply elasticity
param.theta     = 5;        % elasticity of substitution
param.alpha     = 2/3;      % returns to labour
param.rhoa      = 0.35;     % persistence of productivity
param.sigmazeta = 0.225;    % stddev of productivity shocks
param.Phi       = 0.156;    % menu cost in labour units
param.mu        = 0.006;    % s.s money growth
param.rhom       = 0.37;    % persistence of money growth
param.sigmaeps  = 0.0048;   % stddev of money growth shocks
param.tauc      = 0.005;    % tolerance for forecasting rule
param.n         = 5000;     % number of firms
param.T         = 96;       % simulation length
param.S         = 25;       % simulations for computing forecasting coeffs
param.s         = 100;      % simulations for moment computations

% Law of motion - initial guesses
cKS.b0     = -0.03;
cKS.b1     = 0.75;
cKS.b2     = 0.25;

% VAR matrices
A1         = [cKS.b0 + cKS.b2*param.mu*(1-param.rhom); param.mu*(1-param.rhom)];
A2         = [cKS.b1, cKS.b2*param.rhom; 0, param.rhom];
Sigma      = [param.sigmaeps^2*cKS.b2^2 param.sigmaeps^2*cKS.b2; param.sigmaeps^2*cKS.b2 param.sigmaeps^2];

% construct Markov chain
[Pr_mat,Pr_mat_key,zbar] = fn_var_to_markov(eye(2),A1,A2,Sigma,[3; 3],1000,2);    
mgrid = unique(Pr_mat_key(2,:))';
%Ygrid = exp(unique(Pr_mat_key(1,:)))';  
Ygrid = unique(Pr_mat_key(1,:))';  

% run some tests: simulate straight from the VAR
T = 10000;
eps = normrnd(0,param.sigmaeps,[T,1]);
state = zeros(2,T);
state(:,1) = [Ygrid(2) mgrid(2)]';
for t=2:T
    state(:,t) = A1 + A2*state(:,t-1) + [cKS.b2 1]'.*eps(t-1);
end

% simulate from the Markov chain
[path, path_ind] = simul_markov(1:length(Pr_mat_key),Pr_mat,T);
state_markov = zeros(2,T);
for t=1:T
    state_markov(:,t) = Pr_mat_key(:,path(t));
end

% run regressions on Markov chain data
X = [ones(T-1,1) state_markov(1,1:end-1)' state_markov(2,1:end-1)'];
Y1 = state_markov(1,2:end)';
beta1 = regress(Y1,X)';
resid1 = Y1-X*beta1';
Y2 = state_markov(2,2:end)';
beta2 = regress(Y2,X)';
resid2 = Y2-X*beta2';
A1hat = [beta1(1); beta2(1)];
A2hat = [beta1(2:end); beta2(2:end)];
Sigmahat = cov(resid1,resid2);
