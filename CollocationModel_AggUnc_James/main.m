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
% Add CompEcon package
% p = genpath('E:\Dropbox\Economics\Matlab_codes\CompEcon');
p = genpath('C:\Users\James\Dropbox\Economics\Matlab_codes\CompEcon');
addpath(p);

% cd('E:\Dropbox\Economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\MenuCostsModel\CollocationModel_AggUnc_James')
cd('C:\Users\James\Dropbox\economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\MenuCostsModel\CollocationModel_AggUnc_James')

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
options.fontsize    = 12;       % Plot fontsize
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


 %% NO AGGREGATE UNCERTAINTY

% Setup no aggregate uncertainty problem
fprintf('Setup\n');
glob = setup_noagg(param,glob,options);
fprintf('Setup complete\n');

%% Solve value function approx coefficients and stationary distribution for a given output Y

if strcmp(options.solvecL,'Y');
    options.plotpolicyfun = 'Y';      % If Y, plot policy functions
    Y                   = 0.9;  % Conjectured value of output, Y
    options.cresult     = [];   % Holds previous solution for c. Empty in this case.
    eq                  = solve_cL(Y,param,glob,options);
    glob.c              = eq.c;
    fprintf('Yin = %1.2f,\tYout = %1.2f\n',Y,eq.Y);
    fprintf('Pin = %1.2f,\tPout = %1.2f\n',1,eq.P);
    fprintf('--------------------------------------');
end

%% Solve equilibrium
if strcmp(options.solveeq,'Y');
    options.tolp            = 0.0001;           % Tolerance on price
    options.Ylb             = 0.1;              % Output lower bound
    options.Yub             = 5;               % Output upper boud
    options.itermaxp        = 30;               % Max iterations of bisection
    options.eqplot          = 'Y';
    options.eqprint         = 'Y';
    options.print           = 'N';
    options.Loadc           = 'Y';              % For new guess of p use old c as starting guess
    options.plotSD          = 'Y';              % If Y plot steady state distribution
    options.plotpolicyfun   = 'N';            % If Y, plot policy functions
    options.fontsize        = 12;
    eq                      = solve_eq(param,glob,options);
end

save TEMP

%% Reproduce Figure 1 of GS(2007)
%{ 
Na           = 100;
NpP          = 100;
agridlongtmp = nodeunif(Na,min(glob.agrid),max(glob.agrid));  % Adds curvature
pPgridlongtmp = nodeunif(NpP,min(glob.pPgrid),max(glob.pPgrid));  % Adds curvature
s_eval       = gridmake(pPgridlongtmp,agridlongtmp);

% Can interpolate for ANY state vector if function is given
interp_vK    = funfitxy(glob.fspace,glob.sf,eq.v.vK); 
interp_vC    = funfitxy(glob.fspace,glob.sf,eq.v.vC); 
interp_funcs = funeval([interp_vK, interp_vC],glob.fspace,s_eval);
vK           = reshape(interp_funcs(:,1), NpP, Na);
vC           = reshape(interp_funcs(:,2), NpP, Na);
vmax = bsxfun(@max,vK,vC);
figure
for i = 1:length(agridlongtmp)
plot(pPgridlongtmp,vmax(:,i))
hold on
end


ind          = (interp_funcs(:,1) < interp_funcs(:,2));
ind          = double(reshape(ind, NpP, Na));

% krn = [1 -1];
for aa = 1:Na
%    changes = conv(krn,ind(:,aa));
%    idx = find(changes==-1,1,'first');          % These are 1 --> 0 transitions (active to inactive)
   idx = find(ind(:,aa)==0,1,'first');          % These are 1 --> 0 transitions (active to inactive)
%     if idx > 100
%        idx = 100;
%    end
   upperbound(aa) = pPgridlongtmp(idx); 
   idx = find(ind(:,aa)==0,1,'last');          % These are 0 --> 1 transitions (inactive to active)
   lowerbound(aa) = pPgridlongtmp(idx);   
end


% Plot the optimal price with no menu cost
cE = eq.c(2*end/3+1:end);   
B               = menufun_noagg('bounds',glob.sf,[],[],eq.Y,param,glob,options);
options.MC      = 'N';      % Turn off the menu cost
glob.Phi_A      = glob.Phi_Af;
obj             = @(pPstar)valfunc_noagg('C',cE,glob.sf,pPstar,eq.Y,param,glob,options);
pPstar          = goldenx(obj,B(:,1),B(:,2));
interp          = funfitxy(glob.fspace,glob.sf,pPstar); 
pPstar          = funeval(interp,glob.fspace,s_eval);
pPstar          = reshape(pPstar,NpP,Na);

figure(777)
plot(log(agridlongtmp), log(pPstar(1,:)),'linestyle','--')
hold all
plot(log(agridlongtmp), log(upperbound),'color','k','linewidth',2)
hold all
plot(log(agridlongtmp), log(lowerbound),'color','k','linewidth',2)
xlabel('Log productivity, a')
ylabel('Log relative price, p/P')
legend('Optimal','Upper bound','Lower bound')

%}

%% AGGREGATE UNCERTAINTY

close all
glob.damp           = 0.5;
options.burn        = 10;
options.simplot     = 'N';
options.eqprint     = 'N';
options.seed        = 'N';      % Ensures same simulation path each time
options.T           = 96; 
options.T_KSiter    = 25;       % simulations for computing forecasting coeffs
options.tolKS       = 0.01;    %1e-2;
cKS0                 = [0.001; 0.5; 0.1]; % Initialize KS coeffs
options.KSsim       = 1;  
for KSiter = 1:options.T_KSiter
    
    [c,v,cKS,R2,sim] = solve_KS(cKS0,eq,param,glob,options);
    
    cKSnew = glob.damp*cKS0 + (1-glob.damp)*mean(cKS,2);
    
    fprintf('----------------\n') 
    fprintf('%2i. D(cKS) = %2.4f \n',KSiter,norm(cKSnew-cKS0)/norm(cKS0));
    fprintf('R^2 = %2.4f \n',mean(R2));
    fprintf('b0 = %2.4f \n',cKSnew(1));
    fprintf('b1 = %2.4f \n',cKSnew(2));
    fprintf('b2 = %2.4f \n',cKSnew(3));
    fprintf('----------------\n') 
    
    if norm(cKSnew-cKS)<options.tolKS
        cKS0 = cKSnew;
        fprintf('----------------\n')
        fprintf('Solved KS step\n')
        fprintf('----------------\n')
        break
    end
    
    cKS0 = cKSnew;
end

cKS = cKS0;

%% Simulate using the KS method

options.simplot     = 'Y';
options.irf         =  'N';
options.irf         = 'N';
glob = setup_agg(param,glob,cKS,options);
% Set starting point for Y = mean of Y from KS simulation 
eq.Y = mean(sim.Yt); 
[~,~,sim] = simulate_KS(c,v,eq,param,glob,options);

sim.Yt           = Yt;
sim.Pt           = Pt;   
sim.Lt           = Lt;
sim.ind          = ind;   
sim.p_state      = p_state;
sim.pPdist       = pPdist;
sim.DMt          = DMt;
sim.Mt           = Mt; 

% inflation: mean and standard deviation
pi      = sim.Pt(2:options.T)./sim.Pt(1:options.T-1)-1;
pi_mean = mean(pi(options.burn:end));
pi_std  = std(pi(options.burn:end));

% standard deviation of individual (real) prices
% (not controlling for weights)
price_std = std(sim.pPdist(:,2:end),1);

% price_changes
price_change = sim.pPdist - sim.p_state;

% mean price increase, sd new prices
price_increase_ind      = (price_change>0);

dist_pr_inc             = sim.Lt.*price_increase_ind;
dist_pr_inc             = bsxfun(@rdivide, dist_pr_inc, sum(dist_pr_inc,1)); % Normalize columns to sum to 1

log_price_increase      = log(price_change.*price_increase_ind);
pr_log_price_increase   = log_price_increase.*dist_pr_inc;
avg_log_price_increase  = nansum(pr_log_price_increase,1);
mean(avg_log_price_increase(20:end))

% frequency of price change
dist_price_change       = bsxfun(@times,sim.Lt,sim.ind);
freq_price_change       = sum(dist_price_change,1);
mean_freq_price_change  = mean(freq_price_change(20:options.T));

% standard deviation of new prices
avg_price_t             = sum(bsxfun(@times,sim.Lt,log(sim.p_state)),1);
prices_increased        = sim.p_state.*price_increase_ind;
prices_increased(prices_increased==0)=NaN; 
prices_increased        = log(prices_increased);
log_dev_prices_inc      = bsxfun(@minus,prices_increased,avg_price_t);
sd_new_prices           = nanstd(log_dev_prices_inc,1);
mean(sd_new_prices(20:end))


% Frequency of price change = 70%, seems high. Implied annual inflation is
% ~2%, which is probably about right. 


%% Plot IRFs

options.irf         = 'Y';
options.simplot     = 'Y';
% Set starting point for Y = mean of Y from KS simulation 
eq.Y = mean(sim.Yt); 
eq.L = mean(sim.Lt,2);
eq.pi = mean(sim.Pt(2:end)./sim.Pt(1:end-1));

options.T = 300; 
[~,~,sim] = simulate_KS(c,v,eq,param,glob,options);

figure;
subplot(3,1,1)
plot((sim.DMt(1:20) - exp(param.mu))/exp(param.mu)*100);
title('\Delta M_t IRF')
ylabel('% deviation from trend')
subplot(3,1,2)
plot((sim.Yt(1:20)-eq.Y)/eq.Y*100);
title('Y_t IRF')
ylabel('% deviation from trend')
subplot(3,1,3)
% plot( Pt(1:end)./[eq.pi; Pt(1:end-1)]-eq.pi);
  plot( (sim.Pt(2:20)./sim.Pt(1:20-1)-eq.pi)/eq.pi*100);
title('\pi_t IRF')
ylabel('% deviation from trend')
