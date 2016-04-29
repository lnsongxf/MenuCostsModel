%------------------------------
%
%   Main file for Calvo model solved via collocation. Based on similar code
%   for the menu costs model in Terry and Knotek II. 
%
%   James Graham
%   4/22/2016
%
%   Based on original code by Simon Mongey (NYU, 2015)

%% 

cd('E:\Dropbox\Economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\MenuCostsModel\CollocationModel_Calvo_James')
p = genpath('E:\Dropbox\Economics\Matlab_codes\CompEcon');
addpath(p);

% cd('C:\Users\James\Dropbox\economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\MenuCostsModel\CollocationModel_Calvo_James')
% p = genpath('C:\Users\James\Dropbox\Economics\Matlab_codes\CompEcon');
% addpath(p);

close all
clear
clc

%% Set all options

% Things to do
options.solvecL       = 'Y';      % Solve only coefficients c and stationary dist L 
options.solveeq       = 'Y';      % Solve equilibrium (not if agg uncertainty)

options.polfun      = 'Y';      % 
options.solveKS     = 'Y';      % Solve Krussel-Smith
options.sim         = 'Y';      % Solve simulation
options.irf         = 'N';      % Solve IRFs

% Model options 
options.discmethod     = 'R';      % If 'T' use Tauchen, if 'R' use Rouwenhurst

% Compute stationary distribution?
options.stationarydist  ='N';   % Don't compute stationary distribution for aggregate uncertainty case
options.solveKS = 'Y';          % Solve Krussel Smith step


% Tolerances, iterations
options.Nbell       = 2;        % Number of Bellman (Contraction) iterations
options.Nnewt       = 15;       % Maximum number of Newton steps
options.tolc        = 1e-8;     % Tolerance on value functions
options.tolgolden   = 1e-6;     % Tolerance for golden search
options.itermaxL    = 5000;     % Maximum iterations to find stationary dist L
options.tolL        = 1e-11;    % Tolerance on L
options.tolYeq      = 1e-6;    % Tolerance for eqm Y in no agg uncertainty 
options.tolYks      = 1e-2;    % Tolerance for eqm Y in KS step
options.KSsim       = 5;        % Number of KS simulations

% For computation of equilibrium
options.Ylb         = 0.5;              % Output lower bound
options.Yub         = 10;               % Output upper boud
% options.Nfirms      = 5000;             % Number of firms for simulation


% Print / plot 
options.print       = 'Y';      % Print out c-solution convergence
options.eqprint     = 'N';      % Print out equilibrium convergence steps
options.plotSD      = 'N';      % Plot stationary distribution while solving equilibrium
options.fontsize    = 14;       % Plot fontsize
options.fignum      = 888;


%% Statespace parameters
glob.n          = [10,5,3,3];   % Number of nodes in each dimension: [Np,Na,Nm,Ny]
glob.nf         = [500,5];   % Number of points for pP and a in histogram L
glob.curv       = 1;           % Curvature for pP (1 is no curvature, <1 is curvature)
glob.spliorder  = [3,1,1,1];   % Order of splines (use linear if exogenous vars are discrete (not AR1))
glob.pPmin       = 0.75;       % Lower bound on real price
glob.pPmax       = 1.50;         % 25;        % Upper bound on real price

%% Model parameters
% NOTE: delta=0.3 seems to work fine, but delta=0.352 gets closer to the
% paper's plotted policy function. However, delta=0.352 doesn't seem to be
% stable when solving the KS algorithm step...


param.beta      = 0.99;     % discount factor
param.delta     = 0.5;      % 0.5333;   % relative weight of labor-to-consumption in utility
param.sigma     = 1;        % risk aversion coefficient
param.phielas   = 0.5;      % inveser labour supply elasticity
param.theta     = 5;        % elasticity of substitution
param.alpha     = 2/3;      % returns to labour
param.rhoa      = 0.35;     % persistence of productivity
param.sigzeta   = 0.225;      % stddev of productivity shocks
param.Phicost   = 0.156;    % menu cost in labour units
param.mu        = 0.006;    % s.s money growth
param.rhom      = 0.37;     % persistence of money growth
param.sigmaeps  = 0.0048;   % stddev of money growth shocks
param.lambda    = 1-0.3070; % Calvo probability of keeping price
options.MC      = 'N';      % No menu cost in the Calvo model

 %% NO AGGREGATE UNCERTAINTY

% Setup no aggregate uncertainty problem
glob = setup_noagg(param,glob,options);

%% Solve value function approx coefficients and stationary distribution for a given output Y

if strcmp(options.solvecL,'Y');
    options.plotpolicyfun = 'N';      % If Y, plot policy functions
    Y                   = 0.8416;  % Conjectured value of output, Y
    options.cresult     = [];   % Holds previous solution for c. Empty in this case.
    eq                  = solve_cL(Y,param,glob,options);
    glob.c              = eq.c;
%     fprintf('Yin = %1.2f,\tYout = %1.2f\n',Y,eq.Y);
%     fprintf('Pin = %1.2f,\tPout = %1.2f\n',1,eq.P);
    fprintf('--------------------------------------\n');
end

%% Solve equilibrium
if strcmp(options.solveeq,'Y');
    options.Ylb             = 0.1;              % Output lower bound
    options.Yub             = 5;               % Output upper boud
    options.itermaxp        = 50;               % Max iterations of bisection
    options.eqplot          = 'Y';
    options.eqprint         = 'Y';
    options.print           = 'N';
    options.Loadc           = 'Y';              % For new guess of p use old c as starting guess
    options.plotSD          = 'N';              % If Y plot steady state distribution
    options.plotpolicyfun   = 'N';            % If Y, plot policy functions
    eq                      = solve_eq(param,glob,options);
end

figure;
plot(glob.pPgridf,eq.LpP(1:end/2),'-o','color','k')
% hold all
% plot(glob.pPgridf,eq.LpP(end/2+1:end),'-o','color','r')
grid on
% legend('Keep price','Change price')
xlabel('Real price','fontsize',options.fontsize)
ylabel('Density','fontsize',options.fontsize)
set(gca,'fontsize',options.fontsize)


%% Stationary model moments

% Output: mean and standard deviation
momCalvo.y_mean = eq.Y;
momCalvo.y_std  = 0;

% inflation: mean and standard deviation
momCalvo.pi_mean = exp(param.mu);
momCalvo.pi_std  = 0;

% standard deviation of individual (real) prices
% (not controlling for weights)
momCalvo.price_mean = eq.L'*eq.v.pPdist;
momCalvo.price_std = sqrt(eq.L'*(eq.v.pPdist - momCalvo.price_mean).^2);

% price_changes
price_change = eq.v.pPdist - [glob.sf(:,1); glob.sf(:,1)];
nonzero_price_change = abs(price_change(price_change~=0));
nonzero_orig_price = [glob.sf(:,1); glob.sf(:,1)];
nonzero_orig_price = nonzero_orig_price(price_change~=0);
nonzero_L = eq.L(price_change~=0)/sum(eq.L(price_change~=0));   % Drop 0 change, then re-normalize dist

% mean absolute size of price change
momCalvo.avg_abs_price_change = nonzero_L'*(nonzero_price_change./nonzero_orig_price);
momCalvo.std_abs_price_change = sqrt(nonzero_L'*...
    (nonzero_price_change./nonzero_orig_price - momCalvo.avg_abs_price_change).^2);

% mean price increase, sd new prices
momCalvo.price_increase_ind      = (momCalvo.price_change>0);
momCalvo.dist_pr_inc             = eq.L.*momCalvo.price_increase_ind;
momCalvo.dist_pr_inc             = momCalvo.dist_pr_inc./sum(momCalvo.dist_pr_inc,1); % Normalize columns to sum to 1

momCalvo.log_price_increase      = log(momCalvo.price_change.*momCalvo.price_increase_ind);
momCalvo.pr_log_price_increase   = momCalvo.log_price_increase.*momCalvo.dist_pr_inc;
momCalvo.avg_log_price_increase  = nansum(momCalvo.pr_log_price_increase,1);

% frequency of price change
momCalvo.dist_price_change       = eq.L.*(1-eq.v.ind);
momCalvo.freq_price_change       = sum(momCalvo.dist_price_change,1);
momCalvo.avg_duration_prices     = 1/momCalvo.freq_price_change;




%% AGGREGATE UNCERTAINTY

close all
glob.damp           = 0;      % Dampening parameter for updating
options.burn        = 10;       % Burn in period for simulations
options.simplot     = 'N';      % Plot simulations in real time
options.eqprint     = 'N';      % 
options.seed        = 'Y';      % Ensures same simulation path each time
options.irf         = 'N';      % 
options.T           = 96;       % Length of simulation period
options.T_KSiter    = 25;       % simulations for computing forecasting coeffs
options.tolKS       = 1e-4;     % Convergence criterion for norm of CKS params
options.KSsim       = 1;        % Number of simulations to average over for KS 
cKS0                = [-0.0164; 0.8863; 0.9784]; % Initialize KS coeffs

for KSiter = 1:options.T_KSiter
    
    [c,v,cKS,R2,sim]    = solve_KS(cKS0,eq,param,glob,options);
    cKSnew              = glob.damp*cKS0 + (1-glob.damp)*mean(cKS,2);
    conv                = norm(cKSnew-cKS0)/norm(cKS0);
    
    fprintf('----------------\n') 
    fprintf('%2i. D(cKS) = %2.4f \n',KSiter,conv);
    fprintf('R^2 = %2.4f \n',mean(R2));
    fprintf('b0 = %2.4f \n',cKSnew(1));
    fprintf('b1 = %2.4f \n',cKSnew(2));
    fprintf('b2 = %2.4f \n',cKSnew(3));
    fprintf('----------------\n') 
    
    if conv < options.tolKS     % && mean(R2) > 0.7
        cKS0 = cKSnew;
        fprintf('----------------\n')
        fprintf('Solved KS step\n')
        fprintf('----------------\n')
        break
    end
    
    cKS0 = cKSnew;
end

cKS = cKS0;

%% Simulate and compute model moments

options.simplot     = 'Y';
options.irf         = 'N';
glob = setup_agg(param,glob,cKS,options);
% Set starting point for Y = mean of Y from KS simulation 
eq.Y = mean(sim.Yt); 
[~,~,sim] = simulate_KS(c,v,eq,param,glob,options);


% Output: mean and standard deviation
momA.y_mean = mean(sim.Yt(options.burn:end));
momA.y_std  = std(sim.Yt(options.burn:end));

% inflation: mean and standard deviation
pi      = sim.Pt(2:options.T)./sim.Pt(1:options.T-1);
momA.pi_mean = mean(pi(options.burn:end));
momA.pi_std  = std(pi(options.burn:end));

% standard deviation of individual (real) prices
% (not controlling for weights)
momA.price_std = std(sim.pPdist(:,2:end),1);

% price_changes
momA.price_change = sim.pPdist - sim.p_state;

% mean absolute size of price change
tmp1 = sim.Lt'*(abs(momA.price_change)./sim.p_state);
momA.avg_abs_price_change = nanmean(diag(tmp1))*100;
tmp2 = sim.Lt'*...
    (abs(momA.price_change)./sim.p_state - repmat(diag(tmp1)',2500,1)).^2;
momA.std_abs_price_change = nanmean(sqrt(diag(tmp2)))*100;

% mean price increase, sd new prices
momA.price_increase_ind      = (momA.price_change>0);
momA.dist_pr_inc             = sim.Lt.*momA.price_increase_ind;
momA.dist_pr_inc             = bsxfun(@rdivide, momA.dist_pr_inc, sum(momA.dist_pr_inc,1)); % Normalize columns to sum to 1

momA.log_price_increase      = log(momA.price_change.*momA.price_increase_ind);
momA.pr_log_price_increase   = momA.log_price_increase.*momA.dist_pr_inc;
momA.avg_log_price_increase  = nansum(momA.pr_log_price_increase,1);
mean(momA.avg_log_price_increase(20:end))

% frequency of price change
momA.dist_price_change       = bsxfun(@times,sim.Lt,(1-sim.ind));
momA.freq_price_change       = sum(momA.dist_price_change,1);
momA.mean_freq_price_change  = mean(momA.freq_price_change(20:options.T));
momA.avg_duration_prices     = 1/momA.mean_freq_price_change;


% standard deviation of new prices
momA.avg_price_t             = sum(bsxfun(@times,sim.Lt,log(sim.p_state)),1);
momA.prices_increased        = sim.p_state.*momA.price_increase_ind;
momA.prices_increased(prices_increased==0)=NaN; 
momA.prices_increased        = log(momA.prices_increased);
momA.log_dev_prices_inc      = bsxfun(@minus,momA.prices_increased,momA.avg_price_t);
momA.sd_new_prices           = nanstd(momA.log_dev_prices_inc,1);
mean(momA.sd_new_prices(20:end))


% Frequency of price change = 70%, seems high. Implied annual inflation is
% ~2%, which is probably about right. 


%% Plot IRFs

options.irf         = 'Y';
options.simplot     = 'N';  
options.T = 50; 
glob = setup_agg(param,glob,cKS,options);

% Get stationary distribution in KS
eq.Y = mean(sim.Yt); 
eq.L = mean(sim.Lt,2);
eq.pi = mean(sim.Pt(2:end)./sim.Pt(1:end-1));
[~,~,sim] = simulate_KS(c,v,eq,param,glob,options);

% Re-run using stationary eqm values from KS as starting points
options.simplot     = 'Y';  
eq.Y = sim.Yt(end); 
eq.L = sim.Lt(:,end);
tmp = mean(sim.Pt(2:end)./sim.Pt(1:end-1));
eq.pi = tmp(end);
[~,~,sim] = simulate_KS(c,v,eq,param,glob,options);


figure;
subplot(1,3,1)
plot((sim.DMt(1:20) - exp(param.mu))/exp(param.mu)*100, 'linewidth',2);
grid on
title('\Delta M_t IRF')
ylabel('% deviation from trend')
subplot(1,3,2)
plot((sim.Yt(1:20)-eq.Y)/eq.Y*100, 'linewidth',2);
grid on
title('Y_t IRF')
ylabel('% deviation from trend')
subplot(1,3,3)
pisim = sim.Pt(2:20)./sim.Pt(1:20-1);
plot( 400*(pisim - eq.pi)/eq.pi, 'linewidth',2);
grid on
title('\pi_t IRF')
ylabel('% deviation from trend (Annualized)')

save('calvo.mat','sim')

calvo = load('calvo')

load('TEMP')