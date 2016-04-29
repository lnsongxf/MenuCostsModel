%------------------------------
%   Includes aggregate uncertainty in this file 
%
%   Main file for collocation solution to the menu costs model in Terry and
%   Knotek II (2008). 
%
%   James Graham
%   3/6/2016
%
%   Based on original code by Simon Mongey (NYU, 2015)
%--------------------------------
%
%   THINGS TO DO
%   - NOTE: I think one reason things are going wrong is that we're not
%   accounting for steady state money growth in the no-aggregate
%   uncertainty case. Since money is growth at rate mu, need to adjust the
%   price tomorrow by that rate...
%   - TOTALLY REWRITE SIMULATIONS. Need to actually track people over time,
%   using 5000 firms. This way, get an actual distribution of firms rather
%   than the distribution implied by the grid (which is computed inside the value function
%   solver). Simon's method for this doesn't seem to work. He doesn't
%   actually compute the equilibrium at each period, t. And he doesn't
%   allow for aggregate uncertainty. 
%   - Once simulations are rewritten, then can code the IRFs. Follows
%   exactly the same code as the simulation except with deterministic
%   money growth shocks  that follow the money IRF path. 
%   - Add option for Tauchen productivity grid
%   - Possibly figure out Rouwenhurst VAR?
%   - Compute IRFs using Terry's version (doesn't require continuous shocks?)
%   - Compute model moments
%   - Tidy up code significantly. Get rid of rednundant files. Add more
%   comments
%   - Figure out what the correct value of delta should be

%% 

cd('E:\Dropbox\Economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\MenuCostsModel\CollocationModel_AggUnc_James')
p = genpath('E:\Dropbox\Economics\Matlab_codes\CompEcon');
addpath(p);

% cd('C:\Users\James\Dropbox\economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\MenuCostsModel\CollocationModel_AggUnc_James')
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

% [NOT NEEDED HERE - PERHAPS CHANGE OPTIONS]
options.MC          = 'Y';        % If 'N' then menu costs are zero

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
param.lambda    = 0.75;     % Calvo probability of keeping price

 %% NO AGGREGATE UNCERTAINTY

% Setup no aggregate uncertainty problem
glob = setup_noagg(param,glob,options);

%% Solve value function approx coefficients and stationary distribution for a given output Y

if strcmp(options.solvecL,'Y');
    options.plotpolicyfun = 'Y';      % If Y, plot policy functions
    Y                   = 0.8416;  % Conjectured value of output, Y
    options.cresult     = [];   % Holds previous solution for c. Empty in this case.
    eq                  = solve_cL(Y,param,glob,options);
    glob.c              = eq.c;
%     fprintf('Yin = %1.2f,\tYout = %1.2f\n',Y,eq.Y);
%     fprintf('Pin = %1.2f,\tPout = %1.2f\n',1,eq.P);
    fprintf('--------------------------------------');
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
    options.plotSD          = 'Y';              % If Y plot steady state distribution
    options.plotpolicyfun   = 'N';            % If Y, plot policy functions
    eq                      = solve_eq(param,glob,options);
end

figure;
plot(glob.pPgridf,eq.LpP,'-o','color','k')
grid on
xlabel('Real price','fontsize',options.fontsize)
ylabel('Density','fontsize',options.fontsize)
set(gca,'fontsize',options.fontsize)


%% Reproduce Figure 1 of GS(2007)
tmpgridnum      = 200;
a_plot          = nodeunif(tmpgridnum,exp(-0.5),exp(0.5));
s_plot          = gridmake(1,a_plot);

% set up state space
glob.Phi_A      = splibas(glob.agrid0,0,glob.spliorder(2),s_plot(:,2));        % Used in Bellman / Newton computing expected values

% begin by plotting the middle line: price firm would pick if it can costlessly adjust: v.Pc
options.MC      = 'N';      % Switch menu cost off
tmp             = size(s_plot,1)/size(glob.P,2);
glob.Emat       = kron(glob.P,speye(tmp));
a_mid           = solve_valfunc_noagg(eq.c,s_plot,eq.Y,param,glob,options);

% for each point in a, find price (lower and upper bound) at which firm is 
% indifferent between changing and keeping
options.MC      = 'Y';  % turn the menu cost back on
pP_low          = zeros(1,length(a_plot));
pP_upp          = zeros(1,length(a_plot));


pPl             = 500;
for n=1:length(a_plot)
    % for the given level of v, set up the state space: lower bound
    pP_plot_low = nodeunif(pPl,min(glob.pPgrid),a_mid.pPstar(n));    % v_mid.Xc(n)
    pP_plot_upp = nodeunif(pPl,a_mid.pPstar(n),max(glob.pPgrid));      % v_mid.Xc(n)
    pP_plot     = [pP_plot_low; pP_plot_upp(2:end)];
    pP_plot_upp = pP_plot_upp(2:end);
    s_plot      = gridmake(pP_plot,a_plot(n));
    glob.Phi_A  = splibas(glob.agrid0,0,glob.spliorder(2),s_plot(:,2));

    % compute lower/upper bounds on price
    v           = solve_valfunc_noagg(eq.c,s_plot,eq.Y,param,glob,options,1);
    dist        = abs(v.vC - v.vK);
    [~,I_low]   = min(dist(1:pPl));
    pP_low(n)   = pP_plot_low(I_low);
    [~,I_upp]   = min(dist(pPl+1:end));
    pP_upp(n)   = pP_plot_upp(I_upp);
end

% Make figure
x = log(a_plot);
index = 1:1:length(x);
y1 = log(pP_low);
y2 = log(pP_upp);
figure
for i=1:length(y1)-1
    h1 = fill(x(index([i i:i+1 i+1])),...        % Plot each polygon in turn
        [y1(i) y2(i:i+1) y1(i+1)],...
        [0.75,0.75,0.75],'EdgeColor','none');
hold on
end
p1 = plot(log(a_plot),log(a_mid.pPstar),'LineWidth',2,'linestyle','--','Color','k');
hold on;
p2 = plot(log(a_plot),log(pP_low),'LineWidth',2,'Color','k');
hold on;
p3 = plot(log(a_plot),log(pP_upp),'LineWidth',2,'Color','k');
ylim([-0.3,0.3])
set(gca,'fontsize',options.fontsize)
legend([p1,p2],'Optimal real price','Adjustment threshholds')
xlabel('Log Productivity','fontsize',options.fontsize)
ylabel('Log Real Price','fontsize',options.fontsize)

%% Stationary model moments

% Output: mean and standard deviation
momNA.y_mean = eq.Y;
momNA.y_std  = 0;

% inflation: mean and standard deviation
momNA.pi_mean = exp(param.mu);
momNA.pi_std  = 0;

% standard deviation of individual (real) prices
% (not controlling for weights)
momNA.price_mean = eq.L'*eq.v.pPdist;
momNA.price_std = sqrt(eq.L'*(eq.v.pPdist - momNA.price_mean).^2);

% price_changes
momNA.price_change = eq.v.pPdist - glob.sf(:,1);
nonzero_price_change = momNA.price_change(momNA.price_change~=0);
nonzero_price_change = abs(nonzero_price_change);
nonzero_orig_price = glob.sf(momNA.price_change~=0,1);
nonzero_L = eq.L(momNA.price_change~=0)/sum(eq.L(momNA.price_change~=0));   % Drop 0 change, then re-normalize dist

% mean absolute size of price change
momNA.avg_abs_price_change = nonzero_L'*(nonzero_price_change./nonzero_orig_price);
momNA.std_abs_price_change = sqrt(nonzero_L'*...
    (nonzero_price_change./nonzero_orig_price - momNA.avg_abs_price_change).^2);

% mean price increase, sd new prices
momNA.price_increase_ind      = (momNA.price_change>0);
momNA.dist_pr_inc             = eq.L.*momNA.price_increase_ind;
momNA.dist_pr_inc             = momNA.dist_pr_inc./sum(momNA.dist_pr_inc,1); % Normalize columns to sum to 1

momNA.log_price_increase      = log(momNA.price_change.*momNA.price_increase_ind);
momNA.pr_log_price_increase   = momNA.log_price_increase.*momNA.dist_pr_inc;
momNA.avg_log_price_increase  = nansum(momNA.pr_log_price_increase,1);

% frequency of price change
momNA.dist_price_change       = eq.L.*(1-eq.v.ind);
momNA.freq_price_change       = sum(momNA.dist_price_change,1);
momNA.avg_duration_prices     = 1/momNA.freq_price_change;

% standard deviation of new prices
% avg_price_t             = sum(bsxfun(@times,sim.Lt,log(sim.p_state)),1);
% prices_increased        = sim.p_state.*price_increase_ind;
% prices_increased(prices_increased==0)=NaN; 
% prices_increased        = log(prices_increased);
% log_dev_prices_inc      = bsxfun(@minus,prices_increased,avg_price_t);
% sd_new_prices           = nanstd(log_dev_prices_inc,1);
% mean(sd_new_prices(20:end))





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
cKS0                = [-0.05; 0.7; 0.25]; % Initialize KS coeffs
% cKS0                = [0.5; 0.1; 0.1]; % Initialize KS coeffs

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
for t = 2:options.T
    nonzero_price_change = momA.price_change(momA.price_change(:,t)~=0,t);
    nonzero_price_change = abs(nonzero_price_change);
    nonzero_orig_price = glob.sf(momA.price_change(:,t)~=0,1);
    nonzero_L = sim.Lt(momA.price_change(:,t)~=0,t)/...
        sum(sim.Lt(momA.price_change(:,t)~=0,t));   % Drop 0 change, then re-normalize dist
    momA.avg_abs_price_change(t) = nonzero_L'*(nonzero_price_change./nonzero_orig_price);
    momA.std_abs_price_change(t) = sqrt(nonzero_L'*...
        (nonzero_price_change./nonzero_orig_price - momA.avg_abs_price_change(t)).^2);
end
momA.avg_abs_price_change = mean(momA.avg_abs_price_change(2:end));    
momA.std_abs_price_change = mean(momA.std_abs_price_change(2:end));

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
options.T           = 50;
% Get stationary distribution in KS
eq.Y = mean(sim.Yt); 
eq.L = mean(sim.Lt,2);
eq.pi = mean(sim.Pt(2:end)./sim.Pt(1:end-1));
glob = setup_agg(param,glob,cKS,options);
[~,~,sim] = simulate_KS(c,v,eq,param,glob,options);

% Re-run using stationary eqm values from KS as starting points
options.simplot     = 'Y';  
eq.Y = sim.Yt(end); 
eq.L = sim.Lt(:,end);
tmp = mean(sim.Pt(2:end)./sim.Pt(1:end-1));
eq.pi = tmp(1);
[~,~,sim] = simulate_KS(c,v,eq,param,glob,options);


figure(102);
subplot(1,3,1)
plot((sim.DMt(1:20) - exp(param.mu))/exp(param.mu)*100, 'linewidth',2,'color','k');
grid on
title('Money growth')
ylabel('% deviation from trend')
set(gca,'fontsize',options.fontsize)
subplot(1,3,2)
plot((sim.Yt(1:20)-eq.Y)/eq.Y*100, 'linewidth',2,'color','k');
grid on
title('Output gap')
ylabel('% deviation from trend')
set(gca,'fontsize',options.fontsize)
subplot(1,3,3)
pisim = log(sim.Pt(2:20)./sim.Pt(1:20-1));
plot( 400*(pisim - pisim(1)), 'linewidth',2,'color','k');
grid on
title('Inflation')
ylabel('% deviation from trend (Annualized)')
set(gca,'fontsize',options.fontsize)

%% Plot against Calvo model

calvo = load('calvo');

figure(103);
subplot(1,3,1)
plot((sim.DMt(1:20) - exp(param.mu))/exp(param.mu)*100, 'linewidth',2,'color','k');
grid on
title('Money growth')
ylabel('% deviation from trend')
set(gca,'fontsize',options.fontsize)
subplot(1,3,2)
plot((calvo.sim.Yt(1:20)-calvo.sim.Yt(1))/calvo.sim.Yt(1)*100, 'linewidth',2,'color','r');
hold on
plot((sim.Yt(1:20)-eq.Y)/eq.Y*100, 'linewidth',2,'color','k','linestyle','--');
grid on
title('Output gap')
ylabel('% deviation from trend')
set(gca,'fontsize',options.fontsize)
subplot(1,3,3)
pisim = log(sim.Pt(2:20)./sim.Pt(1:20-1));
calvopisim = log(calvo.sim.Pt(2:20)./calvo.sim.Pt(1:20-1));
plot( 400*(calvopisim - calvopisim(1)), 'linewidth',2,'color','r');
hold on
plot( 400*(pisim - pisim(1)), 'linewidth',2,'color','k','linestyle','--');
grid on
legend('Calvo','Menu costs')
title('Inflation')
ylabel('% deviation from trend (Annualized)')
set(gca,'fontsize',options.fontsize)

