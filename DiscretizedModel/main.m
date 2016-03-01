%--------------------------------------------------------
%
% This file computes the firm pricing model with menu costs from Terry &
% Knotek (2008)
%
% James Graham
% 2/18/2016
%
%---------------------------------------------------------
% Add CompEcon package
% p = genpath('E:\Dropbox\Economics\Matlab_codes\CompEcon');
p = genpath('C:\Users\James\Dropbox\Economics\Matlab_codes\CompEcon');
addpath(p);

% cd('E:\Dropbox\Economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\MenuCostsModel\DiscretizedModel')
cd('C:\Users\James\Dropbox\economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\MenuCostsModel\DiscretizedModel')
close all
clear
clc


%% parameters

parms.beta      = 0.99;     % discount factor
parms.delta     = 0.3;      % relative weight of labor-to-consumption in utility
parms.sigma     = 1;        % risk aversion coefficient
parms.phi       = 0.5;      % inveser labour supply elasticity
parms.theta     = 5;        % elasticity of substitution
parms.alpha     = 2/3;      % returns to labour
parms.rhoa      = 0.35;     % persistence of productivity
parms.sigmazeta = 0.225;    % stddev of productivity shocks
parms.Phi       = 0.156;    % menu cost in labour units
parms.mu        = 0.006;    % s.s money growth
parms.rhom       = 0.37;     % persistence of money growth
parms.sigmaeps  = 0.0048;   % stddev of money growth shocks
parms.tauc      = 0.005;    % tolerance for forecasting rule
parms.n         = 5000;     % number of firms
parms.T         = 96;       % simulation length
parms.S         = 25;       % simulations for computing forecasting coeffs
parms.s         = 100;      % simulations for moment computations

%% Create gridspace
% Number of grids for each variable
parms.Npp = 101;  %349;
parms.Na  = 3;
parms.Ny  = 3;
parms.Ndm = 3;
parms.Npi = 3;

% Real price grid
parms.pPmin = 0.75;
parms.pPmax = 1.25;
parms.pPgrid = linspace(parms.pPmin,parms.pPmax,parms.Npp);

% Output, money growth, and inflation grids
% NOTE: construct a VAR, then used VAR-Tauchen to construct grids. Note
% that for many parameterizations inflation is trending down...
b0 = 0.001;
b1 = 0.5;
b2 = 0.1;

% %Productivity grid
[parms.trans_a, agrid] = tauchen(parms.Na,parms.rhoa,parms.sigmazeta,0);
parms.agrid = exp(agrid);   % Take [ln(a_t)] back to non-logged vals 

Nvar = [parms.Ndm;...
    parms.Ny; ...
    parms.Npi];
muvar = [parms.mu*(1-parms.rhom); ...
    b0 + b2*parms.mu*(1-parms.rhom); ...
    -b0 + (1-b2)*parms.mu*(1-parms.rhom)];
Avar = [parms.rhom, 0, 0; ...
    b2*parms.rhom, b1, 0; ...
    (1-b2)*parms.rhom, (1-b1), 0];
Svar = [parms.sigmaeps; ...
    b2*parms.sigmaeps; ...
    (1-b2)*parms.sigmaeps];
[agggrid,parms.trans_agg,~]=tauchenvar(Nvar,muvar,Avar,Svar);
parms.agggrid = exp(agggrid);   % Take [D(ln(M_t)), ln(Y_t), ln(P_t/P_t-1)] back to non-logged vals  

parms.grid = [repmat(parms.agrid,1,parms.Ndm*parms.Ny*parms.Npi); kron(parms.agggrid,ones(1,parms.Na)) ];
Trans_admypi = kron(parms.trans_agg,parms.trans_a);

% Normalize so rows sum to one
s=size(Trans_admypi, 2);
Trans_admypi = spdiags(sum(Trans_admypi,2), 0, s, s)\Trans_admypi;
parms.trans = Trans_admypi;   % Make sure transiton probs are Prob(t->t+1)

%% Value function iteration
%--------------------------
% Need to speed up an insane amount if gonna get up to 349 grid points for
% real prices...
%
% Get rid of loops, vectorize as much as possible...
%
%--------------------------
tic 
[V_final, Vc_final, Vk_final] = val_iter(parms, 0.1);
toc

% Plot value function
idx = ((parms.Na*parms.Ny*parms.Ndm*parms.Npi) - 1)/2 + 1;

figure
plot( parms.pPgrid,squeeze(V_final(idx,:)) )
xlabel('Real price')
ylabel('Value')
title('Value function at mean value of state variables')


%% Compute stationary distribution
dmss = median(parms.grid(2,:));
Yss = median(parms.grid(3,:));
piss = median(parms.grid(4,:));

pPfun = @(pP,a) pricefunc(parms,pP,a,dmss,Yss,piss,Vk_final,Vc_final,V_final);
% ^^^ SOMETHING NOT WORKING HERE....
% FINISH UP TRANSLATING EIGENVECTOR METHOD TO THIS PROBLEM
Qpp = zeros(parms.Npp*parms.Na,parms.Npp);
row = 1;
for a = 1:parms.Na
    for pP = 1:parms.Npp
        lower = find(pPfun(pP,a) >= parms.pPgrid);
        lower = lower(end);
        upper = find(pPfun(pP,a) < parms.pPgrid);
        upper = upper(1);
        
        if isempty(upper)
            Qpp(row,upper) = 1;
        else
            Qpp(row,lower) = (parms.pPgrid(upper) - pPfun(pP,a))/...
                (parms.pPgrid(upper) - parms.pPgrid(lower));
            Qpp(row,upper) = (pPfun(pP,a) - parms.pPgrid(lower))/...
                (parms.pPgrid(upper) - parms.pPgrid(lower));
        end
        row = row + 1;
    end
end

Qa = kron(parms.trans_a,ones(parms.Npp,1));
col = 1;
Q = zeros(parms.Npp*parms.Na,parms.Npp*parms.Na);
for a = 1:parms.Na
    for pP = 1:parms.Npp
        Q(:,col) = Qpp(:,pP).*Qa(:,a);        
        col = col + 1;
    end
end
    
% Pertub the Q matrix so eigenvalues are unique
eta = min(nonzeros(Q))/(2*parms.Npp);
index = find(Q == 0);
Q(index) = eta;

for i = 1:size(Q,1)
    Q(i,:) = Q(i,:)/(sum(Q(i,:)));
end

% Find eigenvector ==> stationary distribution
[V,D,W] = eig(Q);
V = real(V);
D = real(D);
W = real(W);
index = find(diag(D) > 0.999999);








%% Krussel Smith Step...
% TO DO: simulate aggregate states, OLS step, find new coefficients, put
% into larger loop around the value function iteration step. Create a
% markov chain simulation stand-alone code.

% simulate 5000 firms for 96 periods
T = 50;
F = 500;

sim_a = nan(T,F);
sim_a(1,:) = median(parms.agrid);
pPsim = nan(T,F);
sim_agg = nan(3,T);
sim_agg(:,1) = median(parms.grid(2:end,:),2);

for t = 1: T
    agg_idx = find(sim_agg(:,t) == parms.agggrid,);
    sim_a(t+1,f) = randsample(parms.agrid,1,true,parms.trans_a(agg_idx,:));
%     pPout = newprice(parms,pP,statenum,Vk_final,Vc_final,V_final);

    for f = 1:F
        a_idx = find(sim_a(t,f) == parms.agrid);
        sim_a(t+1,f) = randsample(parms.agrid,1,true,parms.trans_a(a_idx,:));
        
    end
            
end





