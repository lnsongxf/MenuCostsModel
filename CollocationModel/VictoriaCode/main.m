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
glob.nf         = [300,5];    % Number of points for p and a in histogram L
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
        Y                   = 1.01;    % Conjectured value of Y    
        options.cresult     = [];   % Holds previous solution for c. Empty in this case.
        eq                  = solve_pL(Y,param,glob,options);  
        fprintf('Yin = %1.2f,\tYout = %1.2f\n',Y,eq.Y);
end
eq.L'*eq.v.Pp
%plot(glob.sf(1:50,1)./(eq.Pa),eq.v.vf(1:50))
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
        options.eqplot      = 'N'; 
        options.eqprint     = 'Y'; 
        options.print       = 'N';
        options.Loadc       = 'Y';              % For new guess of p use old c as starting guess
        options.plotSD      = 'N';              % If Y plot steady state distribution
        eq                  = solve_eq_menucost(param,glob,options); 
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
options.T           = 75;   % simulation length
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
param.Phi       = 0.156; 
p_low = zeros(1,length(a_plot));
p_upp = zeros(1,length(a_plot));

for a=1:length(a_plot)
    
    % for the given level of a, set up the state space: lower bound
    p_plot      = nodeunif(1000,min(glob.pgrid),v_mid.Pc(a));
    s_plot      = gridmake(p_plot,a_plot(a),Y_mean,m_mean);
    glob.Phi_A  = splibas(glob.agrid0,0,glob.spliorder(2),s_plot(:,2));
    glob.Phi_Y  = splibas(glob.Ygrid0,0,glob.spliorder(3),s_plot(:,3));
    glob.Phi_m  = splibas(glob.mgrid0,0,glob.spliorder(4),s_plot(:,4));
    
    % compute lower bound on price
    v_low       = solve_valfuncKS(c,s_plot,param,glob,options);
    dist_low    = abs(v_low.vc - v_low.vk);
    [~,I] = min(dist_low);
    p_low(a) = p_plot(I);
    
    % for the given level of a, set up the state space: upper bound
    p_plot      = nodeunif(1000,v_mid.Pc(a),max(glob.pgrid));
    s_plot      = gridmake(p_plot,a_plot(a),Y_mean,m_mean);
    
    % compute upper bound on price
    v_upp       = solve_valfuncKS(c,s_plot,param,glob,options);
    dist_upp    = abs(v_upp.vc - v_upp.vk);
    [~,I] = min(dist_upp);
    p_upp(a) = p_plot(I);
    
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

