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
p = genpath('E:\Dropbox\Economics\Matlab_codes\CompEcon');
% p = genpath('C:\Users\James\Dropbox\Economics\Matlab_codes\CompEcon');
addpath(p);

cd('E:\Dropbox\Economics\2015_2016_material\AdvMacro_Midrigan\TermPaper')
% cd('C:\Users\James\Dropbox\economics\2015_2016_material\AdvMacro_Midrigan\TermPaper')
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


% Folding a, Dm, y, pi into a single var to get entire state space in one go
% Nvar = [parms.Na; ...
%     parms.Ndm; ...
%     parms.Ny; ...
%     parms.Npi];
% muvar = [0; ...
%     parms.mu*(1-parms.rhom); ...
%     b0 + b2*parms.mu*(1-parms.rhom); ...
%     -b0 + (1-b2)*parms.mu*(1-parms.rhom)];
% Avar = [parms.rhoa, 0, 0, 0; ...
%     0, parms.rhom, 0, 0; ...
%     0, b2*parms.rhom, b1, 0; ...
%     0, (1-b2)*parms.rhom, (1-b1), 0];
% Svar = [parms.sigmazeta; ...
%     parms.sigmaeps; ...
%     b2*parms.sigmaeps; ...
%     (1-b2)*parms.sigmaeps];

% [grid,Transition_admypi,~]=tauchenvar(Nvar,muvar,Avar,Svar);
% grid(1,:) = exp(grid(1,:));   % take ln(a_t) back to a_t
% grid(3,:) = exp(grid(3,:));   % take ln(Y_t) back to Y_t
% grid(4,:) = exp(grid(4,:));   % take ln(P_t/P_t-1) back to P_t/P_t-1
% parms.grid = grid;
% 
% agrid = grid(1,1:parms.Na);
% dmgrid = grid(2,1:parms.Na:parms.Na*parms.Ndm);
% Ygrid = grid(3,1:parms.Na*parms.Ndm:parms.Na*parms.Ndm*parms.Ny);
% pigrid = grid(4,1:parms.Na*parms.Ndm*parms.Ny:end);   
% 
% parms.trans = Transition_admypi';   % Make sure transiton probs are Prob(t->t+1)


%% Value function iteration
%--------------------------
% Need to speed up an insane amount if gonna get up to 349 grid points for
% real prices...
%
% Get rid of loops, vectorize as much as possible...
%
% GET RID OF LOOPS
%--------------------------
tic 
% Initialize value function
T = 500;
Vk(:,:,T) = nan(parms.Na*parms.Ndm*parms.Ny*parms.Npi, parms.Npp);
Vc(:,:,T) = nan(parms.Na*parms.Ndm*parms.Ny*parms.Npi, parms.Npp);
V(:,:,T+1) = nan(parms.Na*parms.Ndm*parms.Ny*parms.Npi, parms.Npp);


% Load previous solutionss
V(:,:,1) = eye(parms.Na*parms.Ny*parms.Ndm*parms.Npi, parms.Npp);
% tmp = load(['FinalV_Na' num2str(parms.Na) ...
%     '_Ndm' num2str(parms.Ndm) ...
%     '_Ny' num2str(parms.Ny) ... 
%     '_Npi' num2str(parms.Npi) ...
%     '_Npp' num2str(parms.Npp)]);
% V(:,:,1) = tmp.V_final;

% pre-allocate real profit functions for speed

for j = 1:parms.Npp
    % realprofit(flag,parms,pP,a,Y)
    profK(:,j)  = realprofit('K',parms,parms.pPgrid(j),parms.grid(1,:),parms.grid(3,:));
end
profC = realprofit('C',parms,parms.pPgrid,parms.grid(1,:),parms.grid(3,:));


% % Pre-find matrix indicies when changing
% for j = 1:parms.Npp
%     tmp = abs( repmat(parms.pPgrid,parms.Na*parms.Ndm*parms.Ny*parms.Npi,1) - ...
%         repmat(parms.pPgrid(1,j)./parms.grid(4, :)',1,parms.Npp) ); % how far from grid?
%     [~, idx_pP] = min(tmp,[],2); % index of closest values
% end


h = waitbar(0,'Value function iteration in progress...');
for t = 1:T
    waitbar(t / T)
    
    % Precompute expectation for speed
%     E_V = nan(parms.Na*parms.Ndm*parms.Ny*parms.Npi,parms.Npp);
%     for i = 1:parms.Na*parms.Ndm*parms.Ny*parms.Npi
%         E_V(i,j) = V(i,idx_pP(i),t);  % Move to new pP grid point in V0
%     end
    
    Vc(:,:,t)   = repmat(valfun('C', parms, parms.pPgrid, V(:,:,t),profC), 1, parms.Npp);
    
    for j = 1:parms.Npp 
        Vk(:,j,t)   = valfun('K', parms, parms.pPgrid(j), V(:,:,t), profK(:,j));        
    end
    
    V(:,:,t+1) = bsxfun(@max,Vk(:,:,t),Vc(:,:,t));  % bsxfun = binary operations on matrices

    disp(['Norm = ' num2str(norm(V(:,:,t+1) - V(:,:,t))) ]) 
end
close(h) 

V_final = V(:,:,end);
Vc_final = Vc(:,:,end);
Vk_final = Vk(:,:,end);

toc 


% Save the final value fuction 
save(['FinalV_Na' num2str(parms.Na) ...
    '_Ny' num2str(parms.Ny) ... 
    '_Ndm' num2str(parms.Ndm) ...
    '_Npi' num2str(parms.Npi) ...
    '_Npp' num2str(parms.Npp) '.mat'], 'V_final');

% NOTE: MUCH faster now. 
% 30,000 grid points ~ 1 minute
% 60,000 grid points ~ 2 minutes.
% 160,000 grid points ~ 9 minutes

% Can make even faster by eliminating the 'j' loop.

% Precompute expectations for additional speed? Check accuracy, seems to be
% something wrong when we do that...

% Look for other speed improvements...
%
%% Plot value function
idx = ((parms.Na*parms.Ny*parms.Ndm*parms.Npi) - 1)/2 + 1;

figure
plot( parms.pPgrid,squeeze(V(idx,:,end)) )
xlabel('Real price')
ylabel('Value')
title('Value function at mean value of state variables')



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





