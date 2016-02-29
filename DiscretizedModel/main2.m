%--------------------------------------------------------
%
% This file computes the firm pricing model with menu costs from Terry &
% Knotek (2008)
% 
% Updated to get rid of inflation as a state variable. Only have 4 state
% variables now.
% 
% James Graham
% 2/18/2016
%
%---------------------------------------------------------
% Add CompEcon package
% p = genpath('E:\Dropbox\Economics\Matlab_codes\CompEcon');
p = genpath('C:\Users\James\Dropbox\Economics\Matlab_codes\CompEcon');
addpath(p);

% cd('E:\Dropbox\Economics\2015_2016_material\AdvMacro_Midrigan\TermPaper')
cd('C:\Users\James\Dropbox\economics\2015_2016_material\AdvMacro_Midrigan\TermPaper')
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
parms.Npp = 51;  %349;
parms.Na  = 3;
parms.Ny  = 3;
parms.Ndm = 3;

% Real price grid
parms.pPmin = 0.75;
parms.pPmax = 1.25;
parms.pPgrid = linspace(parms.pPmin,parms.pPmax,parms.Npp);

% Productivity grid 
% % [Transition_a, agrid] = tauchen(Na,parms.rhoa,parms.sigmazeta,0);
% % agrid2 = exp(agrid);

% Output, money growth, and inflation grids
% NOTE: construct a VAR, then used VAR-Tauchen to construct grids. Note
% that for many parameterizations inflation is trending down...

b0 = 0.001;
b1 = 0.5;
b2 = 0.1;

% Folding a, y, Dm into a single var to get entire state space in one go
Nvar = [parms.Na; ...
    parms.Ndm; ...
    parms.Ny];
muvar = [0; ...
    parms.mu*(1-parms.rhom); ...
    b0 + b2*parms.mu*(1-parms.rhom)];
Avar = [parms.rhoa, 0, 0; ...
    0, parms.rhom, 0; ...
    0, b2*parms.rhom, b1];
Svar = [parms.sigmazeta; ...
    parms.sigmaeps; ...
    b2*parms.sigmaeps];


[grid,Transition_aydm,~]=tauchenvar(Nvar,muvar,Avar,Svar);
parms.grid = exp(grid);   % take grid points from ln(a), ln(Y), Dln(M) --> a, Y, DM  

agrid = grid(1,1:parms.Na);
dmgrid = grid(2,1:parms.Na:parms.Na*parms.Ndm);
Ygrid = grid(3,1:parms.Na*parms.Ndm:parms.Na*parms.Ndm*parms.Ny);

parms.trans = Transition_aydm';   % Make sure transiton probs are Prob(t->t+1)


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
T = 100;
Vk(:,:,T) = nan(parms.Na*parms.Ny*parms.Ndm, parms.Npp);
Vc(:,:,T) = nan(parms.Na*parms.Ny*parms.Ndm, parms.Npp);
V(:,:,T+1) = nan(parms.Na*parms.Ny*parms.Ndm, parms.Npp);

% Load previous solutionss
% V(:,:,1) = eye(parms.Na*parms.Ny*parms.Ndm, parms.Npp);
tmp = load(['FinalV2_Na' num2str(parms.Na) ...
    '_Ny' num2str(parms.Ny) ... 
    '_Ndm' num2str(parms.Ndm) ...
    '_Npp' num2str(parms.Npp)]);
V(:,:,1) = tmp.V2_final;

% pre-allocate real profit functions for speed
for j = 1:parms.Npp
    profK(:,j)  = realprofit2('K',parms,parms.pPgrid(j));
end
profC = realprofit2('C',parms,parms.pPgrid);

% % Pre-find matrix indicies when changing
% for j = 1:parms.Npp
%     tmp = abs( repmat(parms.pPgrid,parms.Na*parms.Ndm*parms.Ny*parms.Npi,1) - ...
%         repmat(parms.pPgrid(1,j)./parms.grid(4, :)',1,parms.Npp) ); % how far from grid?
%     [~, idx_pP] = min(tmp,[],2); % index of closest values
% end


h = waitbar(0,'Value function iteration in progress...');
for t = 1:T
    waitbar(t / T)
    
    for j = 1:parms.Npp
        Vk(:,j,t)   = valfun2('K', parms, parms.pPgrid(j), V(:,:,t), profK(:,j));
    end
    
    Vc(:,:,t)   = repmat(valfun2('C', parms, parms.pPgrid, V(:,:,t),profC), 1, parms.Npp);
    

    
    V(:,:,t+1) = bsxfun(@max,Vk(:,:,t),Vc(:,:,t));  % bsxfun = binary operations on matrices
    
    disp(['Norm = ' num2str(norm(V(:,:,t+1) - V(:,:,t))) ])
end

close(h) 



toc 

% Save the final value fuction 
V2_final = V(:,:,end);
save(['FinalV2_Na' num2str(parms.Na) ...
    '_Ny' num2str(parms.Ny) ... 
    '_Ndm' num2str(parms.Ndm) ...
    '_Npp' num2str(parms.Npp) '.mat'], 'V2_final');

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
idx = ((parms.Na*parms.Ny*parms.Ndm) - 1)/2 + 1;

figure
plot( parms.pPgrid, squeeze(V(idx,:,end)) )
xlabel('Real price')
ylabel('Value')
title('Value function at mean value of state variables')



%% Krussel Smith Step...

% TO DO: This step! 


