%{
This program sets up and calls the functions that
solve the Golosov-Lucas model.

Written by:     Victoria Gregory
Date:           3/28/2016
%}

clear all;
clc;
dbstop if error;
cd '/Users/victoriagregory/Dropbox/MenuCostsModel/CollocationModel/VictoriaCode/Golosov-Lucas'

%% Settings

% What to solve for
options.solvexL     = 'Y';      % Solve for p and L given a Y 
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
glob.n          = [10,5];       % Number of nodes in each dimension
glob.nf         = [300,5];      % Number of points for x and nu in histogram L
glob.curv       = 1;            % Grid curvature for x on (0,1] (1 is no curvature)
glob.spliorder  = [3,1];        % Order of splines (always use linear if shocks are discrete (not AR1))
glob.xmin       = exp(-0.4);    % Lower bound on x
glob.xmax       = exp(0.5);     % Upper bound on x

% NOTE (VG): resulting k grid will be n(1)+spliorder(1)-1
% Creating the cubic spline space adds 3-1=2 points.

% Model parameters
param.rho       = 0.04;
param.gamma     = 2;
param.epsilon   = 7;
param.alpha     = 6;
param.eta       = 0.55;
param.sigmanu   = sqrt(0.011);
%param.sigmanu   = 0.18;
param.k         = 0.0025;
param.mu        = 0.0064;
param.sigmam    = 0.0062;
%param.sigmam    = 0.003;
param.R         = param.rho + param.mu;
param.beta      = exp(-param.R);

% Print / plot 
options.print       = 'Y';      % Print out c-solution convergence

%% Setup problem
fprintf('Setup\n');
[param,glob]    = setup_ss(param,glob,options);      
fprintf('Setup complete\n');

%% Solve only x and L for a given cbar
switch options.solvexL
    case 'Y'
        cbar                = 0.3849;      % Conjectured value of cbar    
        options.cresult     = [];     % Holds previous solution for c. Empty in this case.
        eq                  = solve_xL(cbar,param,glob,options);  
        fprintf('cbarin = %1.2f,\tcbarout = %1.2f\n',cbar,eq.cbar);
end

%% Solve equilibrium
switch options.solveeq
    case 'Y'
        options.tolcbar     = 0.0001;           % Tolerance on cbar
        options.cbarlb      = 0.1;              % cbar lower bound
        options.cbarub      = 0.6;                % cbar upper boud
        options.itermaxcbar = 30;               % Max iterations of bisection
        options.eqplot      = 'N'; 
        options.eqprint     = 'Y'; 
        options.print       = 'N';
        options.Loadc       = 'Y';              % For new guess of cbar use old c as starting guess
        options.plotSD      = 'N';              % If Y plot steady state distribution
        eq                  = solve_eq(param,glob,options); 
end

% cbar = eq.cbar;
% eq2 = solve_valfunc_GL(eq.c,glob.sf,cbar,param,glob,options);


%% Plot stationary distribution and value functions

figure;
subplot(1,2,1)
L_reshape = reshape(eq.L,glob.nf(1),glob.nf(2));
density = sum(L_reshape,2);
plot(glob.xgridf,density)

subplot(1,2,2)
valfun = max(eq.v.vk,eq.v.vc);
valfun_reshape = reshape(valfun,glob.nf(1),glob.nf(2));
plot(glob.xgridf,valfun_reshape);

%% Compute moments from stationary dist.

switch options.moments
    case 'Y'

    % price changes
    price_change = bsxfun(@minus, eq.v.Xp, glob.sf(1:length(eq.L),1));

    % mean price increase
    price_increase_ind      = (price_change>0);
    dist_pr_inc             = eq.L.*price_increase_ind;
    dist_pr_inc_norm        = bsxfun(@rdivide, dist_pr_inc, sum(dist_pr_inc,1));
    log_price_increase      = log(price_change.*price_increase_ind);
    pr_log_price_increase   = log_price_increase.*dist_pr_inc_norm;
    avg_log_price_increase  = nansum(pr_log_price_increase,1);

    % frequency of price change
    dist_price_change       = bsxfun(@times,eq.L,eq.v.Is);
    freq_price_change       = sum(dist_price_change,1);

    % standard deviation of new prices
    avg_price_t             = sum(bsxfun(@times,eq.L,glob.sf(1:length(eq.L),1)),1);
    prices_increased        = glob.sf(1:length(eq.L),1).*price_increase_ind;
    prices_increased(prices_increased==0)=NaN; 
    prices_increased        = log(prices_increased);
    log_dev_prices_inc      = bsxfun(@minus,prices_increased,avg_price_t);
    sd_new_prices           = nanstd(log_dev_prices_inc,1);

end
%% Replicate Figure 1

% switch options.GLfig
%     case 'Y'
% 
%     nu_plot = nodeunif(100,exp(-0.5),exp(0.5));
%     s_plot = gridmake(1,nu_plot);
% 
%     % set up state space
%     glob.Phi_nu  = splibas(glob.nugrid0,0,glob.spliorder(2),s_plot(:,2));
% 
%     % begin by plotting the middle line:
%     % price firm would pick if it can
%     % costlessly adjust: v.Pc
%     param.k         = 0;
%     v_mid           = solve_valfunc_GL(eq.c,s_plot,eq.cbar,param,glob,options,1);
% 
%     % for each point in a, find price (lower and upper bound)
%     % at which firm is indifferent between changing and keeping
%     param.k       = 0.0025;     % turn the menu cost back on
%     x_low = zeros(1,length(nu_plot));
%     x_upp = zeros(1,length(nu_plot));
% 
%     for n=1:length(nu_plot)
% 
%         % for the given level of a, set up the state space: lower bound
%         xl               = 1000;
%         x_plot_low       = nodeunif(xl,min(glob.xgrid),v_mid.Xc(n));
%         x_plot_upp       = nodeunif(xl,v_mid.Xc(n),max(glob.xgrid));
%         x_plot           = [x_plot_low; x_plot_upp(2:end)];
%         x_plot_upp       = x_plot_upp(2:end);
%         s_plot          = gridmake(x_plot,nu_plot(n));
%         glob.Phi_nu     = splibas(glob.nugrid0,0,glob.spliorder(2),s_plot(:,2));
% 
%         % compute lower/upper bounds on price
%         v           = solve_valfunc_GL(eq.c,s_plot,eq.cbar,param,glob,options,1);
%         dist        = abs(v.vc - v.vk);
%         [~,I_low]   = min(dist(1:xl));
%         x_low(n)    = x_plot_low(I_low);
%         [~,I_upp]   = min(dist(xl+1:end));
%         x_upp(n)    = x_plot_upp(I_upp);
% 
%     end
% 
%     % make figure
%     figure;
%     plot(log(nu_plot),log(v_mid.Xc),'--','Color','k')
%     hold on;
%     plot(log(nu_plot),log(x_low),'LineWidth',2,'Color','b')
%     hold on;
%     plot(log(nu_plot),log(x_upp),'LineWidth',2,'Color','b')
%     %set(gca,'Xlim',[-0.5 0.5])
%     %set(gca,'Ylim',[-0.3 0.4])
%     grid on;
%     xlabel('Log Productivity')
%     ylabel('Log Real Price')
% 
% end

%% set up for Krussel-Smith

% State space
glob.n          = [glob.n(1),glob.n(2),6];                   % Number of nodes in each dimension
glob.nf         = [glob.nf(1),glob.nf(2),12];                % Number of points for x, nu, cbar in histogram L
glob.curv       = 1;                                         % Grid curvature for x on (0,1] (1 is no curvature)
glob.spliorder  = [glob.spliorder(1),glob.spliorder(2),1];   % Order of splines (always use linear if shocks are discrete (not AR1))
glob.Ne         = 50;                                        % Number of points for money shocks

% Law of motion - initial guesses
cKS.b0     = -0.105;     % constant
cKS.b1     = 0.9;        % coeff on cbar_{t-1}
cKS.b2     = 1;          % coeff on pi_t
% cKS.b0     = -0.6524;     % constant
% cKS.b1     = 0.3179;        % coeff on cbar_{t-1}
% cKS.b2     = .1816;          % coeff on pi_t

% Print / plot 
options.print       = 'Y';
options.tolcagg     = 0.0001;
options.T           = 300;
options.Tburn       = 20;
options.KSit        = 10;
options.KStol       = 0.001;
options.eqprint     = 'N';

%% Solve Krussel-Smith problem

for itercKS = 1:options.KSit 

    fprintf('----------- Simulation Number:%2.0f -----------\n',itercKS);
    % Solve, simulate, etc
    options.cresult = [];   % Holds previous solution for c. Empty in this case.
    [c,v,KS_coeffs,R2,paths]  = solve_KS(cKS,eq,param,glob,options);
    
    % update Krusell-Smith coefficients guess
    d_b0                   = norm(KS_coeffs(1)-cKS.b0)/norm(cKS.b0);  
    d_b1                   = norm(KS_coeffs(2)-cKS.b1)/norm(cKS.b1); 
    d_b2                   = norm(KS_coeffs(3)-cKS.b2)/norm(cKS.b2);
    cKS.b0                 = KS_coeffs(1);
    cKS.b1                 = KS_coeffs(2);
    cKS.b2                 = KS_coeffs(3);

    % print and check for convergence
    fprintf('b0:\t%2.6f\tb1:\t%2.6f\tb2:%2.6f\n',KS_coeffs(1),KS_coeffs(2),KS_coeffs(3));
    fprintf('norm = %1.4f\n',d_b0+d_b1+d_b2);
    if d_b0+d_b1+d_b2<options.KStol,break,end
end

% re-run setup after solving Krussel-Smith
fprintf('Setup\n');
[param,glob]    = setup_ks(cKS,param,glob,options);      
fprintf('Setup complete\n');

save temp;


%% IRF

load temp;

    switch options.IRF
        case 'Y'

    T           = 30;
    options.T   = T;
    %pi_sim     = ones(T,1).*param.mu;
    %pi_sim(2)  = param.mu + 2*param.sigmam;
    epi_sim     = ones(T,1).*exp(param.mu);         % shock series
    %epi_sim(2)  = epi_sim(2)*1.0125;    
    epi_sim(2)  = exp(param.mu+param.sigmam);
    pi_sim      = log(epi_sim);
    eq.L        = mean(paths.L(:,100:end),2);       % starting points for simulation
    eq.cbar     = mean(paths.C(100:end));
    [~,~,paths_IRF] = simulate_KS(pi_sim,c,v,cKS,eq,param,glob,options);
    eq.L            = paths_IRF.L(:,end);
    [~,~,paths_IRF] = simulate_KS(pi_sim,c,v,cKS,eq,param,glob,options);

    % compute stuff
    pct_dev_C = 100*(paths_IRF.C - paths_IRF.C(end))./paths_IRF.C(end);
    pct_dev_Y = 100*(paths_IRF.Y - paths_IRF.Y(end))./paths_IRF.Y(end);
    pct_dev_l = 100*(paths_IRF.l - paths_IRF.l(end))./paths_IRF.l(end);
    pct_dev_U = 100*(paths_IRF.U - paths_IRF.U(end))./paths_IRF.U(end);
    inflation = paths_IRF.P(2:end)./paths_IRF.P(1:end-1) - 1;
    pct_dev_i = 100*(inflation - inflation(end))./inflation(end);
    pct_dev_i = [0 pct_dev_i];
    
    % make plot
    figure;

    % consumption
    subplot(3,1,1)
    plot(-1:16,pct_dev_C(1:18),'LineWidth',2.5,'Color','k')
    hold on;
    plot(-1:16,zeros(1,18),'LineStyle','--','Color','k')
    grid on;
    xlabel('Quarters after shock','Interpreter','latex')
    ylabel('Deviation from Trend (\%)','Interpreter','latex')
    title('Consumption','Interpreter','latex','FontSize',20)
    set(gca,'Xlim',[-1 16])

    % labor
    subplot(3,1,2);
    plot(-1:16,pct_dev_l(1:18),'LineWidth',2.5,'Color','k')
    hold on;
    plot(-1:16,zeros(1,18),'LineStyle','--','Color','k')
    grid on;
    xlabel('Quarters after shock','Interpreter','latex')
    ylabel('Deviation from Trend (\%)','Interpreter','latex')
    title('Employment','Interpreter','latex','FontSize',20)
    set(gca,'Xlim',[-1 16])
    
    % re-pricing rate
%     subplot(3,1,3)
%     plot(-1:16,pct_dev_U(1:18),'LineWidth',2.5,'Color','k')
%     hold on;
%     plot(-1:16,zeros(1,18),'LineStyle','--','Color','k')
%     grid on;
%     xlabel('Quarters after shock','Interpreter','latex')
%     ylabel('Deviation from Trend (\%)','Interpreter','latex')
%     title('Re-Pricing Rate','Interpreter','latex','FontSize',20)
%     set(gca,'Xlim',[-1 16])
%     hold off;
    
    % inflation
    subplot(3,1,3)
    plot(-1:16,pct_dev_i(1:18),'LineWidth',2.5,'Color','k')
    hold on;
    plot(-1:16,zeros(1,18),'LineStyle','--','Color','k')
    grid on;
    xlabel('Quarters after shock','Interpreter','latex')
    ylabel('Deviation from Trend (\%)','Interpreter','latex')
    title('Inflation','Interpreter','latex','FontSize',20)
    set(gca,'Xlim',[-1 16])
    hold off;

    % save data for figure
    fig_data_GL.dev_C = pct_dev_C;
    fig_data_GL.dev_l = pct_dev_l;
    fig_data_GL.dev_i = pct_dev_i;
    
end

%% moments from model with aggregate shocks

switch options.moments
    case 'Y'

    % inflation: mean and standard deviation
%     pi      = exp(paths.logpi);
%     pi_mean = mean(pi(20:end));
%     pi_std  = std(pi(20:end));
%     fprintf('mean inflation: %1.4f\n',pi_mean);
%     fprintf('std. dev. inflation: %1.4f\n',pi_std);
    
    % price_changes
    st = gridmake(glob.xgridf,glob.nugridf);
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
%     price_change_pct_all    = 100.*bsxfun(@rdivide,price_change,st(:,1));
%     changes = randsample(price_change_pct_all(:,end),10000,true,paths.L(:,end));
%     ksdensity(changes(changes~=0), 'bandwidth',1.2)
%     fig_data_GL.changes = changes;
%     hold off;
end
% 

%% Monte Carlo simulations for hazard rates

% first generate equilibrium paths for consumption
options.T       = 300;
pi_sim          = paths.logpi;
C               = paths.C;
[~,~,paths_equ] = simulate_KS(pi_sim,c,v,cKS,eq,param,glob,options)

% explicitly simulate firms
options.T= 100;
N        = 10000;
C        = paths_equ.C;
paths_MC = simulate_MC(N,C,pi_sim,c,v,cKS,eq,param,glob,options);

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
figure
plot(xi(1:8),hazard(1:8)./3)
hold off;

% save hazard data for comparison
haz_data.paths      = paths_MC;
haz_data.durations  = durations;
haz_data.xi         = xi;
haz_data.hazard     = hazard;
haz_data.surival    = survival;
fname = sprintf('hazards_sp%1.3f_k%1.4f_sm%1.2e_eps%1.0f.mat', param.sigmanu, param.k, param.sigmam, param.epsilon);
cd Hazards
save(fname, 'haz_data')
cd ..

%% plot time series of one guy

wage = param.alpha.*param.R.*paths_MC.price_level;
figure;
plot((2:49)./4,log(paths_MC.price_opt(94,1:48).*wage(1:48)),'LineWidth',2,'color',[0,0,0]+0.5)
hold on;
plot((1:48)./4,log(paths_MC.price(94,1:48).*wage(1:48)),'LineWidth',2.5,'Color','k')
hold on;
plot((1:48)./4,log(paths_MC.P(1:48)),'LineWidth',2.5,'color',[0,0,0]+0.5,'LineStyle','--')
grid on;
xlabel('Years','FontSize',16,'Interpreter','latex')
ylabel('Log Price','FontSize',16,'Interpreter','latex')
h = legend('Optimal Price','Price','Price level');
set(h,'Interpreter','latex','FontSize',14,'Location','Northwest')
set(gca,'Xlim',[0 12])
hold off;
%saveas(gcf, '/Users/victoriagregory/Dropbox/Midrigan/TermPaper/price_sim_GL.eps', 'psc2');

fig_data_GL.price_opt = log(paths_MC.price_opt(94,1:48).*wage(1:48));
fig_data_GL.price     = log(paths_MC.price(94,1:48).*wage(1:48));
fig_data_GL.wage      = log(paths_MC.P(1:48));

%% Plot value function and "stationary distribution"

C_mean  = mean(paths.C(20:options.T));
nusize  = 50;
xsize  = 100;

% create transition matrix over finer a grid
[Pnu_plot,nu_plot,Pssnu]  = setup_MarkovZ(nusize,param.sigmanu,(1-param.eta),1);
QNu_plot                  = kron(Pnu_plot,ones(xsize,1));

% states at which to plot
x_plot = nodeunif(xsize,min(glob.xgrid),max(glob.xgrid));

%a_plot = nodeunif(50,min(glob.agrid),max(glob.agrid));
s_plot = gridmake(x_plot,nu_plot,C_mean);

% set up state space
glob.Phi_nu  = splibas(glob.nugrid0,0,glob.spliorder(2),s_plot(:,2));
glob.Phi_c  = splibas(glob.cgrid0,0,glob.spliorder(3),s_plot(:,3));

% compute policy and value functions
v_plot = solve_valfuncKS(c,s_plot,param,glob,options);
v_max  = max(v_plot.vk, v_plot.vc);
v_max  = reshape(v_max,[xsize nusize]);

% find "stationary distribution"
Xp              = min(v_plot.Xp,max(x_plot)).*1/(exp(param.mu));     
fspaceergx      = fundef({'spli',x_plot,0,1});
QX              = funbas(fspaceergx,Xp); 
Q               = dprod(QNu_plot,QX); 
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
L  = reshape(L,[xsize nusize]);

% compute some conditional densities (conditional on productivity)
L_sums = sum(L,1);
L_cond = bsxfun(@rdivide, L, L_sums);

% save data
fig_data_GL.x_plot  = x_plot;
fig_data_GL.nu_plot = nu_plot;
fig_data_GL.L       = L;
fig_data_GL.L_cond  = L_cond;
fig_data_GL.v_max   = v_max;

%% Replicate Figure 1 in Golosov-Lucas (2007)

C_mean  = mean(paths.C(20:options.T));
nu_plot = nodeunif(100,exp(-0.5),exp(0.5));
s_plot  = gridmake(1,nu_plot,C_mean);

% set up state space
glob.Phi_nu  = splibas(glob.nugrid0,0,glob.spliorder(2),s_plot(:,2));
glob.Phi_c   = splibas(glob.cgrid0,0,glob.spliorder(3),s_plot(:,3));

% begin by plotting the middle line:
% price firm would pick if it can
% costlessly adjust: v.Pc
param.k         = 0;
v_mid           = solve_valfuncKS(c,s_plot,param,glob,options);

% for each point in a, find price (lower and upper bound)
% at which firm is indifferent between changing and keeping
param.k       = 0.0025;     % turn the menu cost back on
x_low = zeros(1,length(nu_plot));
x_upp = zeros(1,length(nu_plot));

for n=1:length(nu_plot)

    % for the given level of a, set up the state space: lower bound
    xl               = 1000;
    x_plot_low       = nodeunif(xl,min(glob.xgrid),v_mid.Xc(n));
    x_plot_upp       = nodeunif(xl,v_mid.Xc(n),max(glob.xgrid));
    x_plot           = [x_plot_low; x_plot_upp(2:end)];
    x_plot_upp       = x_plot_upp(2:end);
    s_plot          = gridmake(x_plot,nu_plot(n),C_mean);
    glob.Phi_nu     = splibas(glob.nugrid0,0,glob.spliorder(2),s_plot(:,2));
    glob.Phi_c      = splibas(glob.cgrid0,0,glob.spliorder(3),s_plot(:,3));

    % compute lower/upper bounds on price
    v           = solve_valfuncKS(c,s_plot,param,glob,options);
    dist        = abs(v.vc - v.vk);
    [~,I_low]   = min(dist(1:xl));
    x_low(n)    = x_plot_low(I_low);
    [~,I_upp]   = min(dist(xl+1:end));
    x_upp(n)    = x_plot_upp(I_upp);

end

% save data
fig_data_GL.nu_plot2 = nu_plot;
fig_data_GL.x_mid   = v_mid.Xc;
fig_data_GL.x_low   = x_low;
fig_data_GL.x_upp   = x_upp;

save('fig_data_GL.mat','fig_data_GL');
