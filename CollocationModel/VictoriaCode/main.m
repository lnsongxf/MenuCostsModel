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
% cd '/Users/victoriagregory/Dropbox/MenuCostsModel/CollocationModel/VictoriaCode'
cd('E:\Dropbox\Economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\MenuCostsModel\CollocationModel\VictoriaCode')
%% Settings

% What to solve for
options.solvepL     = 'Y';      % Solve for p and L given a Y 
options.solveeq     = 'Y';      % Solve equilibrium
options.solveKS     = 'Y';      % Solve Krussel-Smith

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
param.Phi       = 0.156;    % menu cost in labour units
param.mu        = 0.006;    % s.s money growth
param.rhom       = 0.37;     % persistence of money growth
param.sigmaeps  = 0.0048;   % stddev of money growth shocks
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
        Y                   = 0.9;    % Conjectured value of Y    
        options.cresult     = [];   % Holds previous solution for c. Empty in this case.
        eq                  = solve_pL(Y,param,glob,options);  
        fprintf('Yin = %1.2f,\tYout = %1.2f\n',Y,eq.Y);
end
eq.L'*eq.v.Pp
% plot(glob.sf(1:50,1)./(eq.Pa),eq.v.vf(1:50))
% out=funbas(glob.fspace,glob.sf)*eq.c;
% plot(glob.sf(350:400,1)./(eq.Pa),eq.v.vf(350:400),glob.sf(350:400,1)./(eq.Pa),out(350:400))
% legend('RHS','LHS')
plot(glob.pgridf,eq.v.vf(1:500),glob.pgridf,eq.v.vf(501:1000),glob.pgridf,eq.v.vf(1001:1500),glob.pgridf,eq.v.vf(1501:2000))


% valF = funbas(glob.fspace,glob.sf)*eq.c;
% valF = reshape(valF, length(glob.pgridf), length(glob.agridf));
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% plot(glob.pgridf, valF)
% xlabel('Real price','fontsize',12)
% ylabel('Value','fontsize',12)
% set(gca, 'fontsize', 12)
% legend('a_1','a_2','a_3','a_4','a_5')



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
end

L_reshape = reshape(eq.L,glob.nf(1),glob.nf(2));
density = sum(L_reshape,2);
plot(density)

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

%% Replicate Figure 1 in Golosov-Lucas (2007)

load temp;

% aggegate states at which to plot:
% simulation means of m and Y
m_mean = mean(paths.mt(20:options.T));
Y_mean = mean(paths.Y(20:options.T));
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

% make figure
figure;
plot(log(a_plot),log(v_mid.Pc),'--','Color','k')
hold on;
plot(log(a_plot),log(p_low),'LineWidth',2,'Color','b')
hold on;
plot(log(a_plot),log(p_upp),'LineWidth',2,'Color','b')
set(gca,'Xlim',[-0.5 0.5])
set(gca,'Ylim',[-0.3 0.4])
grid on;
xlabel('Log Productivity')
ylabel('Log Real Price')

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


%% to do: impluse responses, moments, GL parameterization

%% model-generated moments

% inflation: mean and standard deviation
pi      = paths.P(2:options.T)./paths.P(1:options.T-1)-1;
pi_mean = mean(pi(20:end));
pi_std  = std(pi(20:end));

% standard deviation of optimal prices
% (not controlling for weights)
sd_prices = std(paths.pol(:,2:end),1);

% price_changes
price_change = paths.pol - paths.ps;

% mean price increase, sd new prices
price_increase_ind      = (price_change>0);
dist_pr_inc             = paths.L.*price_increase_ind;
dist_pr_inc_norm        = bsxfun(@rdivide, dist_pr_inc, sum(dist_pr_inc,1));
log_price_increase      = log(price_change.*price_increase_ind);
pr_log_price_increase   = log_price_increase.*dist_pr_inc_norm;
avg_log_price_increase  = nansum(pr_log_price_increase,1);
mean(avg_log_price_increase(20:end))

% frequency of price change
dist_price_change       = bsxfun(@times,paths.L,paths.I);
freq_price_change       = sum(dist_price_change,1);
mean_freq_price_change  = mean(freq_price_change(20:options.T))

% standard deviation of new prices
avg_price_t             = sum(bsxfun(@times,paths.L,log(paths.ps)),1);
prices_increased        = paths.ps.*price_increase_ind;
prices_increased(prices_increased==0)=NaN; 
prices_increased        = log(prices_increased);
log_dev_prices_inc      = bsxfun(@minus,prices_increased,avg_price_t);
sd_new_prices           = nanstd(log_dev_prices_inc,1);
mean(sd_new_prices(20:end))


%% impluse responses

% mean distribution over original simulation
mean_L      = mean(paths.L(:,20:end),2);
Y_mean      = mean(paths.Y(20:end));

% distribution over idiosyncratic states
L_sim       = zeros(length(mean_L),options.T);
L_sim(:,1)  = mean_L;

% policy functions
pol_sim     = zeros(length(mean_L),options.T-1);
I_sim       = zeros(length(mean_L),options.T-1);

% initial price states
p_state     = zeros(length(mean_L),options.T-1);

% draw money growth shocks
mt          = zeros(1,options.T);
mt(1)       = param.mu;
mt(2)       = mt(1);
mt(3)       = mt(2) + param.sigmaeps;
%rng(222);
for t = 4:options.T;
    mt(t) = param.mu*(1-param.rhom) + param.rhom*mt(t-1);
    % keep from going outside the grid points:
    mt(t) = max(min(mt(t),max(glob.mgrid)),min(glob.mgrid));
end 

% vectors for price level, money, output
Minit       = 0;
M_sim       = zeros(1,options.T);
M_sim(1)    = Minit + mt(1);
for t = 2:options.T
    M_sim(t) = mt(t) + M_sim(t-1);
end
M_sim       = exp(M_sim);
P_sim       = zeros(1,options.T+1);
Y_sim       = zeros(1,options.T+1);
P_sim(1)    = (1/Y_mean)*M_sim(1);
Y_sim(1)    = M_sim(1)/P_sim(1);

% simulate

%t=2;
tictic = tic;
for t = 2:options.T 

    ylb    = 0.5*Y_sim(t-1);
    yub    = 1.5*Y_sim(t-1);
%   Yin    = (1/2)*(ylb+yub);
    Pin    = P_sim(t-1);
    Pout   = Pin;

    for itery = 1:50;

        Yin    = (1/2)*(ylb+yub);

        % 1. define state vector
        st     = gridmake(glob.pgridf,glob.agridf,Yin,mt(t));

        % 2. for Y guess, deflate price distribution
        pi     = Pout/P_sim(t-1);

        % 3. create basis matrix for continuation values
        glob.Phi_A  = splibas(glob.agrid0,0,glob.spliorder(2),st(:,2));
        glob.Phi_Y  = splibas(glob.Ygrid0,0,glob.spliorder(3),st(:,3));
        glob.Phi_m  = splibas(glob.mgrid0,0,glob.spliorder(4),st(:,4));

        % 4. compute real price distribution from policy functions
        v  = solve_valfuncKS(c,[st(:,1).*(1/pi) st(:,2:end)],param,glob,options);

        % 5. Multiply by P_t guess to get nominal price distribution
        nom_prices  = v.Pp.*Pin;

        % 6. Use CES to aggregate price level
        Pout        = (L_sim(:,t-1)'*(nom_prices.^(1-param.theta)))^(1/(1-param.theta));

        % 7. Market clearing
        Yout        = M_sim(t)./Pout;

        % 8. Check convergence of Y
        down        = (Yin>Yout); 
        up          = (Yin<Yout);
        ylb         = up*Yin + down*ylb;
        yub         = up*yub + down*Yin;
        if strcmp(options.eqprint,'Y') 
             fprintf('%2i. yin:\t%2.6f\tyout:\t%2.6f\tt:%2.1f\n',itery,Yin,Yout,toc(tictic));
        end
%            if abs(Yout-Yin)<options.toly;fprintf('Solved\n');break;end;
        if abs(Yout-Yin)<options.toly;break;end;
%             Yin = 0.9*Yin + 0.1*Yout;

    end

    Y_sim(t)        = Yout;
    P_sim(t)        = Pout;
    pol_sim(:,t)    = v.Pp;
    I_sim(:,t)      = v.Is;
    p_state(:,t)    = st(:,1).*(1/pi);
%   Pin - Pout

    % 9. Update distributions (is this right?)
    fspaceerg     = fundef({'spli',glob.pgridf,0,1});
    Pp            = max(min(v.Pp,max(glob.pgridf)),min(glob.pgridf));
    Qp            = funbas(fspaceerg,Pp);
    L_sim(:,t)    = dprod(glob.QA,Qp)'*L_sim(:,t-1);

end

% output gap
dev_Y = (Y_sim - Y_mean)./Y_mean;
plot(dev_Y(1:end-1))

% inflation
pi = P_sim(2:options.T-1)./P_sim(1:options.T-2);
pi2 = paths.P(2:options.T-1)./paths.P(1:options.T-2);
pi2_mean = mean(pi2(20:end));
plot((pi(1:12)-pi2_mean)/pi2_mean)
