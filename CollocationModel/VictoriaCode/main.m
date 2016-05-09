%{
This program sets up and calls the functions that
solve a simple menu cost model (taken from Terry 
(2008)).

Written by:     Victoria Gregory
Date:           2/26/2016
%}

clear all;
clc;
dbstop if error;
cd '/Users/victoriagregory/Dropbox/MenuCostsModel/CollocationModel/VictoriaCode'

%% Settings

% What to solve for/do
options.solvepL     = 'Y';      % Solve for p and L given a Y 
options.solveeq     = 'Y';      % Solve equilibrium
options.solveKS     = 'Y';      % Solve Krussel-Smith
options.GLfig       = 'Y';      % (s,S) bands from Golosov-Lucas
options.moments     = 'Y';      % Compute model moments
options.IRF         = 'Y';      % Impulse responses

% Tolerances, iterations
options.Nbell       = 2;        % Number of Bellman (Contraction) iterations
options.Nnewt       = 15;       % Maximum number of Newton steps
options.tolc        = 1e-8;     % Tolerance on value functions
options.tolgolden   = 1e-6;     % Tolerance for golden search
options.itermaxL    = 5000;     % Maximum iterations to find stationary dist L
options.tolL        = 1e-11;    % Tolerance on L

% Set-up for state space
glob.n          = [10,5];        % Number of nodes in each dimension
glob.nf         = [500,5];    % Number of points for p and a in histogram L
glob.curv       = 1;            % Grid curvature for p/P on (0,1] (1 is no curvature)
glob.spliorder  = [3,1];        % Order of splines (always use linear if shocks are discrete (not AR1))
glob.pmin       = 0.75;         % Lower bound on p
glob.pmax       = 1.50;         % Upper bound on p

% NOTE (VG): resulting k grid will be n(1)+spliorder(1)-1
% Creating the cubic spline space adds 3-1=2 points.

% Model parameters
param.beta      = 0.99;     % discount factor
param.delta     = 0.5;      % relative weight of labor-to-consumption in utility
param.sigma     = 1;        % risk aversion coefficient
param.phi       = 0.5;      % inveser labour supply elasticity
param.theta     = 5;        % elasticity of substitution
param.alpha     = 2/3;      % returns to labour
param.rhoa      = 0.35;     % persistence of productivity
param.sigmazeta = 0.225;    % stddev of productivity shocks
%param.sigmazeta = 0.350;    % stddev of productivity shocks
%param.Phi       = 0.3;    % menu cost in labour units
param.Phi       = 0.156;      % menu cost in labour units
param.mu        = 0.006;    % s.s money growth
param.rhom       = 0.37;     % persistence of money growth
param.sigmaeps  = 0.0048;   % stddev of money growth shocks
%param.sigmaeps  = 0.003;   % stddev of money growth shocks
param.tauc      = 0.005;    % tolerance for forecasting rule
param.n         = 5000;     % number of firms
param.T         = 96;       % simulation length
param.S         = 25;       % simulations for computing forecasting coeffs
param.s         = 100;      % simulations for moment computations

% Print / plot 
options.print       = 'Y';      % Print out c-solution convergence

%% Setup problem
fprintf('Setup\n');
[param,glob]    = setup_ss(param,glob,options);      
fprintf('Setup complete\n');


%% Solve only p and L for a given output Y
switch options.solvepL
    case 'Y'
        Y                   = 0.8;    % Conjectured value of Y    
        options.cresult     = [];   % Holds previous solution for c. Empty in this case.
        eq                  = solve_pL(Y,param,glob,options);  
        fprintf('Yin = %1.2f,\tYout = %1.2f\n',Y,eq.G_Y);
eq.L'*eq.v.Pp
end
% plot(glob.sf(1:50,1)./(eq.Pa),eq.v.vf(1:50))
% out=funbas(glob.fspace,glob.sf)*eq.c;
% plot(glob.sf(350:400,1)./(eq.Pa),eq.v.vf(350:400),glob.sf(350:400,1)./(eq.Pa),out(350:400))
% legend('RHS','LHS')

%% Solve equilibrium
switch options.solveeq
    case 'Y'
        options.tolY        = 0.0001;           % Tolerance on output
        options.Ylb         = 0.1;              % Output lower bound
        options.Yub         = 5;               % Output upper boud
        options.itermaxY    = 30;               % Max iterations of bisection
        options.eqplot      = 'Y'; 
        options.eqprint     = 'Y'; 
        options.print       = 'N';
        options.Loadc       = 'Y';              % For new guess of p use old c as starting guess
        options.plotSD      = 'N';              % If Y plot steady state distribution
        eq                  = solve_eq_menucost(param,glob,options); 
        
        figure
    L_reshape = reshape(eq.L,glob.nf(1),glob.nf(2));
    density = sum(L_reshape,2);
    plot(glob.pgridf,density)
    plot(glob.pgridf,eq.v.vf(1:500),glob.pgridf,eq.v.vf(501:1000),glob.pgridf,eq.v.vf(1001:1500),glob.pgridf,eq.v.vf(1501:2000))

end

%% Set up for Krussel-Smith

% State space
glob.n          = [glob.n(1),glob.n(2),3,3];                    % Number of nodes in each dimension
glob.nf         = [glob.nf(1),glob.nf(2),6,6];                  % Number of points for p and a in histogram L
glob.curv       = 1;                                            % Grid curvature for p/P on (0,1] (1 is no curvature)
glob.spliorder  = [glob.spliorder(1),glob.spliorder(2),1,1];    % Order of splines (always use linear if shocks are discrete (not AR1))

% Law of motion - initial guesses
cKS.b0     = 0.001;
cKS.b1     = 0.5;
cKS.b2     = 0.1;


% Options
options.print       = 'Y';
options.eqprint     = 'N';  % for market clearing step in simulation
options.toly        = 0.01; % tolerance on output Y
options.T           = 96;   % simulation length
options.S           = 5;   % number of simulations


%% Solve Krussel-Smith problem

switch options.solveKS
    case 'Y'
    
    for itercKS = 1:10     
        % Solve problem, simulate, etc.
        options.cresult             = [];   % Holds previous solution for c. Empty in this case.
        [c,v,KS_coeffs,R2,paths]    = solve_KS(cKS,eq,param,glob,options);

        % update Krusell-Smith coefficients guess
        mean_coeffs            = mean(KS_coeffs,1)
        d_b0                   = norm(mean_coeffs(1)-cKS.b0)/norm(cKS.b0);  
        d_b1                   = norm(mean_coeffs(2)-cKS.b1)/norm(cKS.b1); 
        d_b2                   = norm(mean_coeffs(3)-cKS.b2)/norm(cKS.b2);
        cKS.b0                 = mean_coeffs(1);
        cKS.b1                 = mean_coeffs(2);
        cKS.b2                 = mean_coeffs(3);

        % print and check for convergence
        fprintf('d_b0 = %1.4f\n',d_b0);
        fprintf('d_b1 = %1.4f\n',d_b1);
        fprintf('d_b2 = %1.4f\n',d_b2);
        if d_b0+d_b1+d_b2<0.01,break,end
    end
end

% re-run the setup
fprintf('Setup\n');
[param,glob]    = setup_ks(cKS,param,glob,options);      
fprintf('Setup complete\n');

% plot
outc=funbas(glob.fspace,glob.sf)*c(end/3+1:2*end/3);
outk=funbas(glob.fspace,glob.sf)*c(1:end/3);
vf = max(v.vk,v.vc);
vfs = max(outc,outk);
plot(glob.pgridf,vfs(1:glob.nf(1)));

save temp;

%% 

load temp;

m_mean = mean(paths.mt(20:options.T));
Y_mean = mean(paths.Y(20:options.T));

%% The graph: James' way

% load temp
% 
% Na           = 100;
% NpP          = 100;
% agridlongtmp = nodeunif(Na,min(glob.agrid),max(glob.agrid));  % Adds curvature
% pPgridlongtmp = nodeunif(NpP,min(glob.pgrid),max(glob.pgrid));  % Adds curvature
% m_mean = mean(paths.mt(20:options.T));
% Y_mean = mean(paths.Y(20:options.T));
% s_eval       = gridmake(pPgridlongtmp,agridlongtmp,Y_mean,m_mean);
% 
% % Can interpolate for ANY state vector if function is given
% interp_vK    = funfitxy(glob.fspace,glob.s,v.vk); 
% interp_vC    = funfitxy(glob.fspace,glob.s,v.vc); 
% interp_funcs = funeval([interp_vK, interp_vC],glob.fspace,s_eval);
% vK           = reshape(interp_funcs(:,1), NpP, Na);
% vC           = reshape(interp_funcs(:,2), NpP, Na);
% vmax = bsxfun(@max,vK,vC);
% 
% 
% ind          = (interp_funcs(:,1) < interp_funcs(:,2));
% ind          = double(reshape(ind, NpP, Na));
% 
% % krn = [1 -1];
% for aa = 1:Na
% %    changes = conv(krn,ind(:,aa));
% %    idx = find(changes==-1,1,'first');          % These are 1 --> 0 transitions (active to inactive)
%    idx = find(ind(:,aa)==0,1,'first');          % These are 1 --> 0 transitions (active to inactive)
% %     if idx > 100
% %        idx = 100;
% %    end
%    upperbound(aa) = pPgridlongtmp(idx); 
%    idx = find(ind(:,aa)==0,1,'last');          % These are 0 --> 1 transitions (inactive to active)
%    lowerbound(aa) = pPgridlongtmp(idx);   
% end
% 
% 
% % Plot the optimal price with no menu cost
% cE = c(2*end/3+1:end);   
%     glob.Phi_A  = splibas(glob.agrid0,0,glob.spliorder(2),s_eval(:,2));
%     glob.Phi_Y  = splibas(glob.Ygrid0,0,glob.spliorder(3),s_eval(:,3));
%     glob.Phi_m  = splibas(glob.mgrid0,0,glob.spliorder(4),s_eval(:,4));
%     
%     % compute lower bound on price
%     v2       = solve_valfuncKS(c,s_eval,param,glob,options);
% interp          = funfitxy(glob.fspace,s_eval,v2.Pc); 
% pPstar          = funeval(interp,glob.fspace,s_eval);
% pPstar          = reshape(pPstar,NpP,Na);
% 
% figure(777)
% plot(log(agridlongtmp), log(pPstar(1,:)),'linestyle','--')
% hold all
% plot(log(agridlongtmp), log(upperbound),'color','k','linewidth',2)
% hold all
% plot(log(agridlongtmp), log(lowerbound),'color','k','linewidth',2)
% xlabel('Log productivity, a')
% ylabel('Log relative price, p/P')
% legend('Optimal','Upper bound','Lower bound')


%% model-generated moments

switch options.moments
    case 'Y'

    % inflation: mean and standard deviation
    pi      = paths.P(2:options.T)./paths.P(1:options.T-1)-1;
    pi_mean = mean(pi(20:end));
    pi_std  = std(pi(20:end));
    fprintf('mean inflation: %1.4f\n',pi_mean);
    fprintf('std. dev. inflation: %1.4f\n',pi_std);

    % price_changes
    st = gridmake(glob.pgridf,glob.agridf);
    price_change = bsxfun(@minus, paths.pol, st(:,1));

    % mean price increase, sd new prices
    price_increase_ind      = (price_change>0);
    dist_pr_inc             = paths.L.*price_increase_ind;
    dist_pr_inc_norm        = bsxfun(@rdivide, dist_pr_inc, sum(dist_pr_inc,1));
    log_price_increase      = log(price_change.*price_increase_ind);
    pr_log_price_increase   = log_price_increase.*dist_pr_inc_norm;
    avg_log_price_increase  = nansum(pr_log_price_increase,1);
    mean_log_price_stat     = mean(avg_log_price_increase(20:end));

    % frequency of price change
    dist_price_change       = paths.L.*paths.I;
    freq_price_change       = sum(dist_price_change,1);
    mean_freq_price_change  = mean(freq_price_change(20:options.T));
    fprintf('frequency of price change (monthly): %1.4f\n',mean_freq_price_change/3);
    fprintf('monthly duration: %1.4f\n',-1/log(1-mean_freq_price_change/3));
    
    % proportion of price changes that are increases
    dist_price_increase     = paths.L.*price_increase_ind;
    freq_price_increase     = sum(dist_price_increase,1);
    prop_increase           = freq_price_increase./freq_price_change;
    prop_increase_stat      = mean(prop_increase(20:options.T));
    fprintf('proportion increase: %1.4f\n',prop_increase_stat);

    % standard deviation of new prices
    avg_price_t             = sum(bsxfun(@times,paths.L,log(st(:,1))),1);
    prices_increased        = bsxfun(@times, st(:,1), price_increase_ind);
    prices_increased(prices_increased==0)=NaN; 
    prices_increased        = log(prices_increased);
    log_dev_prices_inc      = bsxfun(@minus,prices_increased,avg_price_t);
    sd_new_prices           = nanstd(log_dev_prices_inc,1);
    sd_new_prices_stat      = mean(sd_new_prices(20:end));
    fprintf('std. dev. new prices: %1.4f\n',sd_new_prices_stat);

    % average absolute size of price changes
    price_change_pct        = abs(100.*bsxfun(@rdivide,price_change,st(:,1)));
    price_change_ind        = (price_change~=0);
    dist_pr_chg             = paths.L.*price_change_ind;
    dist_pr_chg_norm        = bsxfun(@rdivide, dist_pr_chg, sum(dist_pr_chg,1));
    price_change_pct_i      = price_change_pct.*price_change_ind;
    price_change_pct_i(price_change_pct_i==0)=NaN; 
    pr_price_change_pct_i   = price_change_pct_i.*dist_pr_chg_norm;
    avg_abs_price_change    = nansum(pr_price_change_pct_i,1);
    avg_abs_size_stat       = mean(avg_abs_price_change(20:end));
    fprintf('average absolute size of change: %1.4f\n',avg_abs_size_stat);
    
    % std. dev. of absolute sizes
    price_change_dev        = (price_change_pct - avg_abs_size_stat).^2;
    price_change_pct_i_dev  = price_change_dev.*price_change_ind;
    price_change_pct_i_dev(price_change_pct_i_dev==0)=NaN; 
    pr_price_change_pct_i_dev   = price_change_pct_i_dev.*dist_pr_chg_norm;
    var_abs_price_change    = nansum(pr_price_change_pct_i_dev,1);
    std_abs_size_stat       = sqrt(mean(var_abs_price_change(20:end)));
    fprintf('std.dev absolute size of change: %1.4f\n',std_abs_size_stat);
    
    % distribution of price changes
    price_change_pct_all    = 100.*bsxfun(@rdivide,price_change,st(:,1));
    changes = randsample(price_change_pct_all(:,end),10000,true,paths.L(:,end));
    ksdensity(changes(changes~=0),'bandwidth',1.2)
    fig_data_KT.changes = changes;
end

%% impluse responses

switch options.IRF
    case 'Y'

    % draw money growth shocks: note lagged timing here to conform with
    % timing of function simulate_KS.m
    mt          = zeros(1,options.T+1);
    mt(1)       = param.mu;
    mt(2)       = param.mu;
    mt(3)       = mt(1) + param.sigmaeps;
    for t = 4:options.T+1;
        mt(t) = param.mu*(1-param.rhom) + param.rhom*mt(t-1);
        mt(t) = max(min(mt(t),max(glob.mgrid)),min(glob.mgrid));
    end 
    Minit       = 0;
    M_sim       = zeros(1,options.T+1);
    M_sim(1)    = Minit + mt(1);
    M_sim(2)    = Minit + mt(1);
    for t = 3:options.T+1
        M_sim(t) = mt(t) + M_sim(t-1);
    end
    M_sim       = exp(M_sim);

    % compute steady state L, Y
    [~,~,paths_IRF] = simulate_KS(mt,M_sim,c,v,cKS,eq,param,glob,options);
    eq.L            = paths_IRF.L(:,end);
    eq.Y            = paths_IRF.Y(end);
    eq.v.Is         = paths_IRF.I(:,end);
    eq.v.Pp         = paths_IRF.pol(:,end);
    [~,~,paths_IRF] = simulate_KS(mt,M_sim,c,v,cKS,eq,param,glob,options);
    Y_sim           = paths_IRF.Y;
    P_sim           = paths_IRF.P;

    % compute repricing rate along paths
    reprice = sum(paths_IRF.L.*paths_IRF.I,1);

    %% create figure

    figure;

    % output gap
    subplot(3,1,1)
    dev_Y = 100*(paths_IRF.Y - Y_sim(1))./Y_sim(1);
    plot(-1:16,dev_Y(1:18),'LineWidth',2.5,'Color','k')
    grid on;
    xlabel('Quarters after shock','Interpreter','latex')
    ylabel('Deviation from Trend (\%)','Interpreter','latex')
    title('Output','Interpreter','latex','FontSize',20)
    %set(gca,'Ylim',[0 .6])
    set(gca,'Xlim',[-1 16])

    % inflation
    subplot(3,1,2);
    pi = paths_IRF.P(2:options.T-1)./paths_IRF.P(1:options.T-2)-1;
    dev_pi = 100*(pi-param.mu)/param.mu;
    plot(-1:16,dev_pi(1:18),'LineWidth',2.5,'Color','k')
    grid on;
    xlabel('Quarters after shock','Interpreter','latex')
    ylabel('Deviation from Trend (\%)','Interpreter','latex')
    title('Inflation','Interpreter','latex','FontSize',20)
    set(gca,'Xlim',[-1 16])
    
    % aggregate labor
    subplot(3,1,3)
    dev_L = 100*(paths_IRF.Labor - paths_IRF.Labor(end))./paths_IRF.Labor(end);
    plot(-1:16,dev_L(1:18),'LineWidth',2.5,'Color','k')
    grid on;
    xlabel('Quarters after shock','Interpreter','latex')
    ylabel('Deviation from Trend (\%)','Interpreter','latex')
    title('Employment','Interpreter','latex','FontSize',20)
    set(gca,'Xlim',[-1 16])
    %set(gca,'Ylim',[0 0.8])
    hold off;
    %saveas(gcf, '/Users/victoriagregory/Dropbox/Midrigan/TermPaper/irfs2_KT.eps', 'psc2');

end

% save data for figure
fig_data_KT.dev_Y = dev_Y;
fig_data_KT.dev_pi = dev_pi;
fig_data_KT.dev_L = dev_L;

%% Monte Carlo simulations for hazard rates

% first generate equilibrium paths for money, output
% start from "stationary" distibution from above
% take same shocks as original simulation
mt              = paths.mt;
M_sim           = paths.M;
[~,~,paths_equ] = simulate_KS(mt,M_sim,c,v,cKS,eq,param,glob,options);

% explicitly simulate firms
N        = 10000;
Y        = paths_equ.Y;
paths_MC = simulate_MC(N,Y,mt,c,v,cKS,eq,param,glob,options);

% extract the durations
durations = [];
for i = 1:N
    policy_inv = abs(paths_MC.policy(i,:) - 1);
    [pol n]    = RunLength(policy_inv);
    durs       = pol.*n;
    durs       = durs(1:end-1);     % cut off last spell to avoid censoring
    durs       = durs(durs~=0);
    durations  = [durations; durs'];
end

% empirical hazard function
[xi hazard survival] = hazard_discrete(durations);
figure;
plot(xi(1:end-1),hazard(1:end-1)/3)

% save hazard data for comparison
haz_data.paths      = paths_MC;
haz_data.durations  = durations;
haz_data.xi         = xi;
haz_data.hazard     = hazard;
haz_data.surival    = survival;
fname = sprintf('hazards_sp%1.3f_Phi%1.3f_sm%1.2e_th%1.0f.mat', param.sigmazeta, param.Phi, param.sigmaeps, param.theta);
cd Hazards
save(fname, 'haz_data')
cd ..

%% price trajectory of one guy

figure;
plot((2:49)./4,log(paths_MC.price_opt(71,1:48).*paths.P(1:48)),'LineWidth',2,'color',[0,0,0]+0.5)
hold on;
plot((1:48)./4,log(paths_MC.price(71,1:48).*paths.P(1:48)),'LineWidth',2.5,'Color','k')
hold on;
plot((1:48)./4,log(paths.P(1:48)),'LineWidth',2.5,'color',[0,0,0]+0.5,'LineStyle','--')
grid on;
xlabel('Years','FontSize',16,'Interpreter','latex')
ylabel('Log Price','FontSize',16,'Interpreter','latex')
h = legend('Optimal Price','Price','Price level');
set(h,'Interpreter','latex','FontSize',14,'Location','Northwest')
set(gca,'Xlim',[0 12])
hold off;
%saveas(gcf, '/Users/victoriagregory/Dropbox/Midrigan/TermPaper/price_sim_KT.eps', 'psc2');
%plot(1:48,log(paths_MC.price(3,1:48).*paths.P(1:48)),2:49,log(paths_MC.price_opt(3,1:48).*paths.P(1:48)),1:48,log(paths.P(1:48)))

fig_data_KT.price_opt = log(paths_MC.price_opt(71,1:48).*paths.P(1:48));
fig_data_KT.price     = log(paths_MC.price(71,1:48).*paths.P(1:48));
fig_data_KT.price_lev = log(paths.P(1:48));

%% Plot value function and "stationary distribution" 

m_mean = mean(paths.mt(20:options.T));
Y_mean = mean(paths.Y(20:options.T));
asize  = 50;
psize  = 100;

% create transition matrix over finer a grid
[Pa_plot,a_plot,Pssa]  = setup_MarkovZ(asize,param.sigmazeta,param.rhoa,1);
QA_plot                = kron(Pa_plot,ones(psize,1));

% states at which to plot
p_plot = nodeunif(psize,min(glob.pgrid),max(glob.pgrid));

%a_plot = nodeunif(50,min(glob.agrid),max(glob.agrid));
s_plot = gridmake(p_plot,a_plot,Y_mean,m_mean);

% set up state space
glob.Phi_A  = splibas(glob.agrid0,0,glob.spliorder(2),s_plot(:,2));
glob.Phi_Y  = splibas(glob.Ygrid0,0,glob.spliorder(3),s_plot(:,3));
glob.Phi_m  = splibas(glob.mgrid0,0,glob.spliorder(4),s_plot(:,4));

% compute policy and value functions
v_plot = solve_valfuncKS(c,s_plot,param,glob,options);
v_max  = max(v_plot.vk, v_plot.vc);
v_max  = reshape(v_max,[psize asize]);

% find "stationary distribution"
Pp              = min(v_plot.Pp,max(p_plot)).*1/(exp(param.mu));     
fspaceergp      = fundef({'spli',p_plot,0,1});
QP              = funbas(fspaceergp,Pp); 
Q               = dprod(QA_plot,QP); 
L               = ones(size(Q,1),1);
L               = L/sum(L);
for itL = (1:options.itermaxL);
    Lnew    = Q'*L;  
    dL      = norm(Lnew-L)/norm(L);  
    if (dL<options.tolL),break,end;
    if mod(itL,100)==0 
        if strcmp(options.print,'Y')
            fprintf('dL:\t%1.3e\n',dL);
        end
    end
    L       = Lnew;
end
L  = reshape(L,[psize asize]);
%contourf(a_plot,p_plot,L,15)

% smooth out stationary dist.
% a_plot_fine = linspace(min(a_plot),max(a_plot),700);
% p_plot_fine = linspace(min(p_plot),max(p_plot),1000);
% L_smooth = interp2(a_plot,p_plot,L,a_plot_fine,p_plot_fine','spline');

% compute some conditional densities (conditional on productivity)
L_sums = sum(L,1);
L_cond = bsxfun(@rdivide, L, L_sums);

% density of just real price
density = sum(L,2);
%plot(p_plot,density)

% save data
fig_data_KT.p_plot = p_plot;
fig_data_KT.a_plot = a_plot;
fig_data_KT.L      = L;
fig_data_KT.L_cond = L_cond;
fig_data_KT.v_max  = v_max;


%% Replicate Figure 1 in Golosov-Lucas (2007)

m_mean = mean(paths.mt(20:options.T));
Y_mean = mean(paths.Y(20:options.T));

% states at which to plot:
% simulation means of m and Y
a_plot = nodeunif(100,min(glob.agrid),max(glob.agrid));
s_plot = gridmake(1,a_plot,Y_mean,m_mean);

% set up state space
glob.Phi_A  = splibas(glob.agrid0,0,glob.spliorder(2),s_plot(:,2));
glob.Phi_Y  = splibas(glob.Ygrid0,0,glob.spliorder(3),s_plot(:,3));
glob.Phi_m  = splibas(glob.mgrid0,0,glob.spliorder(4),s_plot(:,4));

% begin by plotting the middle line:
% price firm would pick if it can
% costlessly adjust: v.Pc
param.Phi       = 0;
v_mid           = solve_valfuncKS(c,s_plot,param,glob,options);

% for each point in a, find price (lower and upper bound)
% at which firm is indifferent between changing and keeping
param.Phi       = 0.156;     % turn the menu cost back on
p_low = zeros(1,length(a_plot));
p_upp = zeros(1,length(a_plot));

for a=1:length(a_plot)

    % for the given level of a, set up the state space: lower bound
    pl               = 1000;
    p_plot_low       = nodeunif(pl,min(glob.pgrid),v_mid.Pc(a));
    p_plot_upp       = nodeunif(pl,v_mid.Pc(a),max(glob.pgrid));
    p_plot           = [p_plot_low; p_plot_upp(2:end)];
    p_plot_upp       = p_plot_upp(2:end);
    s_plot      = gridmake(p_plot,a_plot(a),Y_mean,m_mean);
    glob.Phi_A  = splibas(glob.agrid0,0,glob.spliorder(2),s_plot(:,2));
    glob.Phi_Y  = splibas(glob.Ygrid0,0,glob.spliorder(3),s_plot(:,3));
    glob.Phi_m  = splibas(glob.mgrid0,0,glob.spliorder(4),s_plot(:,4));

    % compute lower/upper bounds on price
    v           = solve_valfuncKS(c,s_plot,param,glob,options);
    dist        = abs(v.vc - v.vk);
    [~,I_low]   = min(dist(1:pl));
    p_low(a)    = p_plot_low(I_low);
    [~,I_upp]   = min(dist(pl+1:end));
    p_upp(a)    = p_plot_upp(I_upp);

end

% save data
fig_data_KT.a_plot2 = a_plot;
fig_data_KT.p_mid   = v_mid.Pc;
fig_data_KT.p_low   = p_low;
fig_data_KT.p_upp   = p_upp;

save('fig_data_KT.mat','fig_data_KT');

