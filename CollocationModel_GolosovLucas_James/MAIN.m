%------------------------------
%   Includes aggregate uncertainty in this file 
%
%   Main file for collocation solution to the menu costs model in Golosov and Lucas (2007). 
%
%   James Graham
%
%   (This is modified code based on original code by Simon Mongey (NYU,
%   2015))

%% 
% Add CompEcon package
% 
% cd('C:\Users\James\Dropbox\economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\MenuCostsModel\CollocationModel_GolosovLucas_James')
% p = genpath('C:\Users\James\Dropbox\Economics\Matlab_codes\CompEcon');
% addpath(p);

cd('E:\Dropbox\Economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\MenuCostsModel\CollocationModel_GolosovLucas_James')
p = genpath('E:\Dropbox\Economics\Matlab_codes\CompEcon');
addpath(p);

clear
clc

%% Set all options

% Things to do
options.solvecL       = 'Y';      % Solve only c and L 
options.solveeq       = 'Y';      % Solve equilibrium (not if agg uncertainty)

options.polfun      = 'Y';      % 
options.solveKS     = 'Y';      % Solve Krussel-Smith
options.sim         = 'Y';      % Solve simulation
options.IRF         = 'N';

% Model options 
options.discmethod  = 'R';      % If 'T' use Tauchen, if 'R' use Rouwenhurst
options.MC          = 'Y';        % If 'N' then menu costs are zero

% Compute stationary distribution?
options.stationarydist  ='N';   % Don't compute stationary distribution for aggregate uncertainty case
options.solveKS = 'Y';          % Solve Krussel Smith step


% Tolerances, iterations
options.Nbell       = 5;        % Number of Bellman (Contraction) iterations
options.Nnewt       = 25;       % Maximum number of Newton steps
options.tolc        = 1e-8;     % Tolerance on value functions
options.tolgolden   = 1e-8;     % Tolerance for golden search
options.itermaxL    = 5000;     % Maximum iterations to find stationary dist L
options.tolL        = 1e-8;    % Tolerance on L
options.tolYeq      = 1e-8;    % Tolerance for eqm Y in no agg uncertainty 
options.tolYks      = 1e-2;    % Tolerance for eqm Y in KS step

% For computation of equilibrium
options.damp        = 0.5;      %  Dampening parameter computation of eqm values
options.T           = 100;       % simulation length
options.Nfirms      = 5000;             % Number of firms for simulation


% Print / plot 
options.print       = 'Y';      % Print out c-solution convergence
options.eqprint     = 'Y';      % Print out equilibrium convergence steps
options.plotSD      = 'N';      % Plot stationary distribution while solving equilibrium
options.plotpolicyfun = 'N';      % If Y, plot policy functions
options.fontsize    = 12;       % Plot fontsize
options.fignum      = 999;

%% Statespace parameters
glob.n          = [10,5,3,3];   % Number of nodes in each dimension: [Np,Nv,Nm,Ny]
glob.nf         = [300,5];   % Number of points for pP and a in histogram L
glob.curv       = 1;           % Curvature for pP (1 is no curvature, <1 is curvature)
glob.spliorder  = [3,1,1,1];   % Order of splines (use linear if exogenous vars are discrete (not AR1))
glob.pPmin       = exp(-0.4);       % Lower bound on real price
glob.pPmax       = exp(0.5);        % Upper bound on real price

%% Model parameters
param.gamma     = 2;                % risk aversion
param.epsilon   = 7;                % elasticity of substition
param.alpha     = 6;                % disutility of labour
param.eta       = 0.55;             % 
param.rhov      = (1-param.eta);    % AR(1) coefficient for productivity
param.sigv      = sqrt(0.011);      % Std dev of productivity process
param.k         = 0.0025;           % menu cost 
param.mu        = 0.0064;           % Quarterly inflation rate: mu = 0.0064
param.sigm      = 0.0062;           % Std dev of money growth shocks
param.rho       = 0.04;             % Annual discount rate: rho = 0.04
param.R         = param.rho + param.mu;  % Stationary interest rate
param.beta      = exp(-param.R);    

glob.piw                = param.mu;   % set the steady state inflation rate

%% NO AGGREGATE UNCERTAINTY

% Setup no aggregate uncertainty problem
glob = setup_noagg(param,glob,options);

%% Solve value function approx coefficients and stationary distribution for a given output Y
if strcmp(options.solvecL,'Y');
    options.plotpolicyfun   = 'N';              % If Y, plot policy functions
    Y                       = 0.5;  % Conjectured value of output, Y
    options.cresult         = [];   % Holds previous solution for c. Empty in this case.
    eq                      = solve_cL(Y,param,glob,options);
    glob.c                  = eq.c;
    fprintf('Yin = %1.2f,\tYout = %1.2f\n',Y,eq.Y);
    fprintf('--------------------------------------');
end

%% Solve equilibrium
if strcmp(options.solveeq,'Y');
    options.Yinit           = 1;
    options.itermaxp        = 30;               % Max iterations of bisection
    options.eqplot          = 'Y';
    options.eqprint         = 'Y';
    options.print           = 'N';
    options.Loadc           = 'Y';              % For new guess of p use old c as starting guess
    options.plotSD          = 'N';              % If Y plot steady state distribution
    options.fontsize        = 12;
    options.plotpolicyfun   = 'N';              % If Y, plot policy functions
    eq                      = solve_eq(param,glob,options);
end

% Compute momnets
ind_c = 1 - eq.v.ind;
freqpricechange = eq.L'*ind_c;


% vC = reshape(eq.v.vC,glob.nf(1),glob.nf(2));
% vK = reshape(eq.v.vK,glob.nf(1),glob.nf(2));
% figure
% subplot(1,2,1)
% plot(max(vC, vK))
% subplot(1,2,1)
% plot(max(vC, vK))

save TEMP


%% Reproduce Figure 1 of GS(2007)
tmpgridnum      = 100;
v_plot          = nodeunif(tmpgridnum,exp(-0.5),exp(0.5));
s_plot          = gridmake(1,v_plot);

% set up state space
glob.Phi_V      = splibas(glob.vgrid0,0,glob.spliorder(2),s_plot(:,2));        % Used in Bellman / Newton computing expected values
glob.Phi_V      = splibas(glob.vgrid0,0,glob.spliorder(2),s_plot(:,2));        % Used in Bellman / Newton computing expected values

% begin by plotting the middle line: price firm would pick if it can costlessly adjust: v.Pc
options.MC      = 'N';      % Switch menu cost off
tmp             = size(s_plot,1)/size(glob.P,2);
glob.Emat       = kron(glob.P,speye(tmp));
v_mid           = solve_valfunc_noagg(eq.c,s_plot,eq.Y,param,glob,options);

% for each point in a, find price (lower and upper bound) at which firm is 
% indifferent between changing and keeping
options.MC      = 'Y';  % turn the menu cost back on
pP_low          = zeros(1,length(v_plot));
pP_upp          = zeros(1,length(v_plot));


pPl             = 500;
for n=1:length(v_plot)
    % for the given level of v, set up the state space: lower bound
    pP_plot_low = nodeunif(pPl,min(glob.pPgrid),v_mid.pPstar(n));    % v_mid.Xc(n)
    pP_plot_upp = nodeunif(pPl,v_mid.pPstar(n),max(glob.pPgrid));      % v_mid.Xc(n)
    pP_plot     = [pP_plot_low; pP_plot_upp(2:end)];
    pP_plot_upp = pP_plot_upp(2:end);
    s_plot      = gridmake(pP_plot,v_plot(n));
    glob.Phi_V  = splibas(glob.vgrid0,0,glob.spliorder(2),s_plot(:,2));

    % compute lower/upper bounds on price
    v           = solve_valfunc_noagg(eq.c,s_plot,eq.Y,param,glob,options,1);
    dist        = abs(v.vC - v.vK);
    [~,I_low]   = min(dist(1:pPl));
    pP_low(n)   = pP_plot_low(I_low);
    [~,I_upp]   = min(dist(pPl+1:end));
    pP_upp(n)   = pP_plot_upp(I_upp);
end

% make figure
figure;
plot(log(v_plot),log(v_mid.pPstar),'--','Color','k')
hold on;
plot(log(v_plot),log(pP_low),'LineWidth',2,'Color','b')
hold on;
plot(log(v_plot),log(pP_upp),'LineWidth',2,'Color','b')

grid on;
xlabel('Log Productivity','fontsize',12)
ylabel('Log Real Price','fontsize',12)


%% Compute impulse responses the GL way
% Note: these impulses are created by hitting the money supply with an
% unexpected shock in period 1, and then agents have perfect foresight over
% the transition path that follows. 

% Reset the model states, parameters, etc
glob = setup_noagg(param,glob,options);

options.numtrans    = 40;              % Number of transition periods
options.damp        = 0.5;
options.IRF         = 'Y';
options.plotSD      = 'N';              % If Y plot steady state distribution
options.eqprint     = 'Y';
options.tolY        = 1e-4;

solve_IRF(eq,param,glob,options)





%% AGGREGATE UNCERTAINTY
param.cKS0 = [0.001; 0.5; 0.1];  % Initial guess for state evolution equation
cKS        = param.cKS0;

% Setup aggregate uncertainty problem (for a particular set of cKS params).
fprintf('Setup Aggregate Uncertainty\n');
glob = setup_agg(param,glob,cKS,options);
fprintf('Setup complete\n');

%% UP TO HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% Solve Krussel-Smith step
close all
glob.damp           = 0.5;
options.burn        = 20;
options.simplot     = 'Y';
options.eqprint     = 'N';
options.seed        = 'Y';      % Ensures same simulation path each time
options.T           = 96; 
options.T_KSiter    = 25;       % simulations for computing forecasting coeffs
options.tolKS       = 1e-2;

for KSiter = 1:options.T_KSiter
    
    % Run setup file again with new cKS params
    glob = setup_agg(param,glob,cKS,options);

    [DMt,Yt,~] = solve_KS(cKS,eq,param,glob,options);
    DMt         = DMt(options.burn+1:end);
    Yt          = Yt(options.burn+1:end);

    % Regression step
    Xt = [ones(length(DMt)-1,1), log(Yt(1:end-1)), log(DMt(2:end))];
    beta = (Xt'*Xt)^(-1)*(Xt'*log(Yt(2:end)));
    resid = log(Yt(2:end)) - Xt*beta;
    Rsquared = 1 - sum(resid.^2)/sum( (log(Yt(2:end)) - mean(log(Yt(2:end)))).^2);
    % Updating coefficients
    cKSnew = glob.damp*cKS + (1-glob.damp)*beta;
    
    fprintf('----------------\n') 
    fprintf('%2i. D(cKS) = %2.4f \n',KSiter,norm(cKSnew-cKS)/norm(cKS));
    fprintf('%2i. R^2 = %2.4f \n',KSiter,Rsquared);
    fprintf('----------------\n') 
    
    if norm(cKSnew-cKS)<options.tolKS
        cKS = cKSnew;
        fprintf('----------------\n')        
        fprintf('Solved KS step\n')  
        fprintf('----------------\n')
        break
    end
    
    cKS = cKSnew;

end

glob.cKS   = cKS;

%% Simulate using the same method as KS

switch options.sim
    case 'Y'
        % Initial vector of firms (can try different things)
        % 1. All pP's, average productivity
        s01                 = repmat(gridmake(glob.pPgrid,glob.agrid(round(glob.Na/2))),500,1);
        % 2. Low pP, all productivities - Growth
        s02                 = repmat(gridmake(glob.pPgridf(10),glob.agrid),1,1);
        % 3. Low pP, all productivities - Growth
%         s04                 = repmat(gridmake(Kss,1*exp(-glob.sige)),5000,1);
        s04                 = repmat(gridmake(min(glob.pPgrid),min(glob.agridf)),10000,1);
        % Choose s0
        s0                  = s01; 
     
        % Options
        options.simsolve    = 'N';          % Solve problem each period (EXACT)
                                            % If N, then linearly
                                            % interpolate policies (FAST) -
        options.constprod   = 'N';          % Keep productivity constant
        options.simplot     = 'Y';
        options.IRF         = 'N';          % Fix as 'N'

        % Solve
        options.AC          = 'Y';
        sim                 = solve_simulation(s0,eq,1,param,glob,options); 
end






%% Simulate model (no uncertainty first?)  
% Simon's version of this doesn't seem to solve for equilibrium... How is
% it really simulating? 
switch options.sim
    case 'Y'
        % Initial vector of firms (can try different things)
        % 1. All pP's, average productivity
        s01                 = repmat(gridmake(glob.pPgrid,glob.agrid(round(glob.Na/2))),500,1);
        % 2. Low pP, all productivities - Growth
        s02                 = repmat(gridmake(glob.pPgridf(10),glob.agrid),1,1);
        % 3. Low pP, all productivities - Growth
%         s04                 = repmat(gridmake(Kss,1*exp(-glob.sige)),5000,1);
        s04                 = repmat(gridmake(min(glob.pPgrid),min(glob.agridf)),10000,1);
        % Choose s0
        s0                  = s01; 
     
        % Options
        options.simsolve    = 'N';          % Solve problem each period (EXACT)
                                            % If N, then linearly
                                            % interpolate policies (FAST) -
        options.constprod   = 'N';          % Keep productivity constant
        options.simplot     = 'Y';
        options.IRF         = 'N';          % Fix as 'N'

        % Solve
        options.AC          = 'Y';
        sim                 = solve_simulation(s0,eq,1,param,glob,options); 
end



%% IRFs for one time monetary shock.

options.Mshockson = 'N'

options.Mshockson = 'Y'




%% Simulate the model 
options.T_simnum         = 10;    %100;      % simulations for moment computations
glob = setup_agg(param,glob,cKS,options);
options.simplot     = 'Y';
options.seed        = 'N';
options.Mshockson   = 'Y';


DMt = nan(options.T - options.burn, options.T_simnum);
Yt  = nan(options.T - options.burn, options.T_simnum);
Pt  = nan(options.T - options.burn, options.T_simnum);
pit = nan(options.T - options.burn - 1, options.T_simnum);  
    
for s = 1:options.T_simnum
        
    [simDM,simY,simP] = solve_GIRF(eq,param,glob,options);
    DMt(:,s)         = simDM(glob.burn+1:end);
    Yt(:,s)          = simY(glob.burn+1:end);
    Pt(:,s)          = simP(glob.burn+1:end);
    pit(:,s)         = Pt(2:end,s)./Pt(1:end-1,s);  
            
    % Average size of price changes   
    % Average duration between price changes    
    % inflation volatility
    % Output volatility
    
end


% Compute moments
moments.mean.DM = mean(DMt,1);
moments.mean.Y = mean(Yt,1);
moments.mean.pi = mean(pit,1);

moments.var.DM = var(DMt,1);
moments.var.Y = var(Yt,1);
moments.var.pi = var(pit,1);


%% TO DO 
% - Simulate model.
% - Compute model moments
% - Compute model IRFs
% - Find correct value of delta



%% PREVIOUS KS CODE, SIMULATION CODE, IRF CODE
%{


%% Solve Krussell-Smith
% load TEMP3

switch options.solveKS
    case 'Y'
        NA              = 3;
        NM              = 3;
        NY              = 3;
        NpP              = 3;
        glob.spliorder  = [3,1,3,1];  % For old model: [k,K,z,A]
            
        
        
        %____________________________________________________
        % Nodes for k and z
        kgrid0          = glob.kgrid0;
        zgrid0          = glob.zgrid0;
        %____________________________________________________
        % Nodes for A
%         rhoA            = glob.rhoz;                     glob.rhoA   = rhoA;  
%         sigA            = glob.sige;                     glob.sigA   = sigA;  
        rhoA            = 0.999;                     
        glob.rhoA   = rhoA;     %%%% DEBUG
        sigA            = 0.005;                     
        glob.sigA   = sigA;     %%%% DEBUG
        pAlb            = 0.05;
        pAub            = 1-pAlb;
        logAlb          = norminv(pAlb,0,sigA/sqrt(1-rhoA^2));
        logAub          = norminv(pAub,0,sigA/sqrt(1-rhoA^2));
        Alb             = exp(logAlb);
        Aub             = exp(logAub);
        Agrid0          = nodeunif(NA,Alb,Aub);     glob.Agrid0 = Agrid0; 
%         Agrid0          = [0.999,1,1.001]';           glob.Agrid0 = Agrid0;    %%%% DEBUG
        %____________________________________________________
        % Nodes for K
        Kmin            = 0.70*eq.Ka;
        Kmax            = 1.30*eq.Ka;
        Kgrid0          = nodeunif(NK,Kmin,Kmax);               glob.Kgrid0 = Kgrid0;
        Kgrid0          = [0.50*eq.Ka,eq.Ka,1.50*eq.Ka]';       glob.Kgrid0 = Kgrid0;   %%%% DEBUG
%         Kgrid0          = [0.999*eq.Ka;1.001*eq.Ka]';  glob.Kgrid0 = Kgrid0; %%%% DEBUG  
%         %____________________________________________________
        % Function space
        fspace      = fundef({'spli',kgrid0,0,glob.spliorder(1)},...
                            {'spli',Kgrid0,0,glob.spliorder(2)},...
                            {'spli',zgrid0,0,glob.spliorder(3)},...
                            {'spli',Agrid0,0,glob.spliorder(4)});
        sgrid       = funnode(fspace); 
        s           = gridmake(sgrid);      glob.s      = s;
        Ns          = size(s,1);            glob.Ns     = Ns;
        kgrid       = unique(s(:,1));       glob.kgrid  = kgrid;
        Kgrid       = unique(s(:,2));       glob.Kgrid  = Kgrid;
        zgrid       = unique(s(:,3));       glob.zgrid  = zgrid;
        Agrid       = unique(s(:,4));       glob.Agrid  = Agrid;
        %__________________________________________________________________
        % First guess of coefficients - capital
        cKS.gamma_p     = [log(eq.p);0;0;0];
        cKS.gamma_K     = [log(eq.Ka);0;0;0];
        %__________________________________________________________________
        % Compute expectations matrix
        Ne_z        = glob.Ne1;
        pvec_z      = nodeunif(Ne_z,glob.plb,1-glob.plb);
        e_z         = norminv(pvec_z,0,glob.sige);
        w_z         = normpdf(e_z,0,glob.sige);
        w_z         = w_z/sum(w_z);
        %__________________________________________________________________
        Ne_A        = 5;
        pvec_A      = nodeunif(Ne_A,0.01,0.99);
        e_A         = norminv(pvec_A,0,sigA);
        w_A         = normpdf(e_A,0,sigA);
        w_A         = w_A/sum(w_A); 
%         Ne_A        = 1;                    %%%% DEBUG
%         e_A         = zeros(Ne_A,1);        %%%% DEBUG
%         w_A         = ones(Ne_A,1)/Ne_A;    %%%% DEBUG  
        %__________________________________________________________________
        iNs         = ones(Ns,1);
        gfun        = @(z,e,rho,zgrid) max(min(exp(rho*log(z)+e),max(zgrid)),min(zgrid));
        %__________________________________________________________________
        w           = kron(w_z,w_A);
        e_A         = kron(ones(Ne_z,1),e_A); 
        e_z         = kron(e_z,ones(Ne_A,1));
        Ne          = Ne_A*Ne_z;
        iNe         = ones(Ne,1);   
        %__________________________________________________________________
        % Evolution of k
        k           = s(:,1);
        %__________________________________________________________________
        % Evolution of K
        % Note: Uses the conjectured function for capital
        K           = s(:,2);
        A           = s(:,4);
        X           = [ones(size(K)),log(K),log(A),log(A).*log(K)];
        Kp          = exp(X*cKS.gamma_K);
        Kp          = max(min(Kp,max(Kgrid)),min(Kgrid));
        %__________________________________________________________________
        % Evolution of A
        A           = s(:,4);
        gA          = gfun(kron(A,iNe),kron(ones(Ns,1),e_A),rhoA,Agrid); 
        %__________________________________________________________________
        % Evolution of z
        z           = s(:,3);
        gz          = gfun(kron(z,iNe),kron(ones(Ns,1),e_z),glob.rhoz,zgrid);
        %__________________________________________________________________
        % Compute expectations matrix 
        Phi_k       = splibas(kgrid0,0,glob.spliorder(1),kron(k,iNe));
        Phi_K       = splibas(Kgrid0,0,glob.spliorder(2),kron(Kp,iNe)); 
        Phi_z       = splibas(zgrid0,0,glob.spliorder(3),gz); 
        Phi_A       = splibas(Agrid0,0,glob.spliorder(4),gA); 
        Phi_zA      = dprod(Phi_A,Phi_z);
        Phi_KzA     = dprod(Phi_zA,Phi_K);
        Phi_kKzA    = dprod(Phi_KzA,Phi_k);
        Ikronw      = kron(speye(Ns),w');
        glob.Emat   = Ikronw*Phi_kKzA;      % Checks comformability 
        % Save components of Emat, then put together given a guess of
        % coefficients
        glob.Emat_Ikronw    = Ikronw;
        glob.Emat_Phi_zA    = Phi_zA;
        glob.Emat_Phi_k     = Phi_k; 
        glob.iNe            = iNe;
        %__________________________________________________________________
        % One time only basis matrices used in computation of continuation
        % value
        Phi_k           = splibas(kgrid0,0,glob.spliorder(1),k);
        Phi_K           = splibas(Kgrid0,0,glob.spliorder(2),K);
        Phi_z           = splibas(zgrid0,0,glob.spliorder(3),z);
        Phi_A           = splibas(Agrid0,0,glob.spliorder(4),A); 
        Phi_zA          = dprod(Phi_A,Phi_z);
        glob.Phi_KzA    = dprod(Phi_zA,Phi_K);          % For contuinuation
        glob.Phi        = dprod(glob.Phi_KzA,Phi_k);    % For jacobian
        %__________________________________________________________________
        % KS algorithm
        options.Nbell   = 3;
        options.print   = 'Y';
        options.tolp    = 0.00001;
        options.eqprint = 'Y';
        %__________________________________________________________________
        % Find best initial capital under constant shocks
%         options.Krun    = 'Y';
%         options.T       = 50;
%         [At,pt,Kt]      = solve_KS(cKS,eq,param,glob,options); 
%         pstar           = pt(find(pt>0,1,'last'));
%         Kstar           = Kt(find(Kt>0,1,'last'));
%         cKS.gamma_p     = [log(pstar);0;0;0];
%         cKS.gamma_K     = [log(Kstar);0;0;0];
%         keyboard
        %__________________________________________________________________
        % Algorithm
        options.Krun = 'N';
        options.T    = 40;
        for itercKS = 1:10; 
            %______________________________________________________________
            % Solve
            [At,pt,Kt]      = solve_KS(cKS,eq,param,glob,options); 
            %______________________________________________________________
            % Tidy up time-series
            Tburn           = 10;
            At              = At(Tburn:end-1);
            pt              = pt(Tburn:end-1);
            Ktp1            = Kt(Tburn+1:end);
            Kt              = Kt(Tburn:end-1);
            %______________________________________________________________
            % Regressions
            % 1. Forecasting rule for Ktp1
            Xt              = [ones(size(Kt)),log(Kt),log(At),log(Kt).*log(At)];
            [gamma_K,~,~,~,stats_K] = regress(log(Ktp1),Xt);
            dgamma_p = norm(gamma_p-cKS.gamma_p)/norm(cKS.gamma_p);  
            % 2. Forecasting rule for pt
            Xt              = [ones(size(Kt)),log(Kt),log(At),log(Kt).*log(At)];
            [gamma_p,~,~,~,stats_p] = regress(log(pt),Xt);
            dgamma_K = norm(gamma_K-cKS.gamma_K)/norm(cKS.gamma_K);
            %______________________________________________________________
            % Print
            fprintf('r_K = %1.3f\t dgamma_K = %1.4f\n',stats_K(1),dgamma_K);
            fprintf('r_p = %1.3f\t dgamma_p = %1.4f\n',stats_p(1),dgamma_p);
            if dgamma_K+dgamma_p<0.001,
                return 
            end
        end
    
end


















%% Policy functions
% % % load TEMP
% % % switch options.polfun
% % %     case 'Y'
% % %         %__________________________________________________________________
% % %         % 1. Varying a for 3 levels of pP
% % %         pPbar    = eq.L'*glob.sf(:,1);
% % %         pPplot   = [.5*pPbar, pPbar, 1.25*pPbar]';
% % %         aplot   = glob.agridf;
% % %         pPpolA   = [];
% % %         for i = 1:numel(pPplot);
% % %             s           = gridmake(pPplot(i),aplot);
% % %             glob.Phi_A  = splibas(glob.agrid0,0,glob.spliorder(2),s(:,2));
% % %             
% % %             v           = solve_valfunc(eq.c,s,eq.Y,param,glob,options,1);
% % %             pPpolA(:,i)  = v.pPstar;
% % %         end
% % %         %__________________________________________________________________
% % %         % 2. Varying pP for 3 levels of a
% % %         abar    = eq.L'*glob.sf(:,2);
% % %         aplot   = [.5*abar, abar, 1.25*abar]';
% % %         pPplot   = glob.pPgridf;
% % %         pPpolP   = [];
% % %         for i = 1:numel(aplot);
% % %             s           = gridmake(pPplot,aplot(i));
% % %             glob.Phi_A  = splibas(glob.zgrid0,0,glob.spliorder(2),s(:,2));
% % %             v           = solve_valfunc(eq.c,s,eq.Y,param,glob,options,1);
% % %             pPpolP(:,i)  = v.pPstar;
% % %         end
% % %         %__________________________________________________________________
% % %         figure(round(1000*rand));
% % %         subplot(2,2,1);plot(glob.pPgridf,pPpolP);title('1A. Real price');xlabel('Real price');legend('Low a','Med a','High a');grid on;
% % % %         subplot(2,2,2);plot(glob.kgridf,ipolK);title('1B. Investment');xlabel('Capital');legend('Low Z','Med Z','High Z');grid on;
% % %         subplot(2,2,3);plot(glob.agridf,pPpolA);title('2A. Real price');xlabel('Productivity');legend('Low pP','Med pP','High pP');grid on;
% % % %         subplot(2,2,4);plot(glob.zgridf,ipolZ);title('2B. Investment');xlabel('Productivity');legend('Low K','Med K','High K');grid on;
% % % end
% % % 
% % % % Interesting: Constprod + s0

%% IRF
% % % % Note: An IRF is only possible from the following (pP,a)
% % % % a must be the unconditional mean of the productivity process
% % % % pP must be the long-run price associated with that a under no-shocks
% % % load TEMP
% % % % options.solveIRF = 'Y';
% % % switch options.solveIRF
% % %     case 'Y'
% % %         if ~strcmp(options.AR1,'Y')
% % %             return
% % %         end;
% % %         %__________________________________________________________________
% % %         Ez                  = glob.Psszf'*glob.zgridf; 
% % %         %__________________________________________________________________
% % %         % 1. Solve for long-run capital for z=E[z]=1 under no shocks
% % %         options.simsolve    = 'N';
% % %         options.constprod   = 'Y';
% % %         options.IRF         = 'N';
% % %         options.simplot     = 'N';
% % %         T                   = 200;
% % %         s0                  = gridmake(median(glob.kgridf),1);
% % %         sim                 = solve_simulation(s0,T,eq,1,param,glob,options);  
% % %         Kss                 = sim.Kt(end,:)'; 
% % %         %__________________________________________________________________
% % %         % 2. IRF - Start at [Kss,1] and solve with AR(1) shock in period 1
% % %         options.simsolve    = 'Y';
% % %         options.constprod   = 'N';
% % %         options.IRF         = 'Y';
% % %         options.simplot     = 'Y';
% % %         T                   = 100;
% % %         s0                  = [Kss,1];
% % %         sim                 = solve_simulation(s0,T,eq,1,param,glob,options); 
% % % end

%% Simulate model
switch options.sim
    case 'Y'
        %__________________________________________________________________
        % Initial vector of firms (can try different things)
        % 1. All a's average productivity
        s01  = repmat(gridmake(glob.pPgrid,glob.agrid(round(glob.Na/2))),500,1);
        % 2. Small a, all productivities - Growth
        s02  = repmat(gridmake(glob.pPgridf(10),glob.agrid),1,1);
        % 3. Small a, all productivities - Growth
%         s04                 = repmat(gridmake(Kss,1*exp(-glob.sige)),5000,1);
        s04  = repmat(gridmake(min(glob.pPgrid),min(glob.agridf)),10000,1);
        % 3. One individual point
        pPbar = eq.L'*glob.sf(:,1);
        abar = eq.L'*glob.sf(:,2);
        s03  = [pPbar,abar];
        % Choose s0
        s0   = s04; 
        %__________________________________________________________________
        % Number of periods
        T    = 400;
        % Options
        options.simsolve    = 'N';          % Solve problem each period (EXACT)
                                            % If N, then linearly
                                            % interpolate policies (FAST) -
        options.constprod   = 'N';          % Keep productivity constant
        options.simplot     = 'Y';
        options.IRF         = 'N';          % Fix as 'N'
        %__________________________________________________________________
        % Solve
        options.MC          = 'Y';
        sim                 = solve_simulation(s0,T,eq,1,param,glob,options); 
end




%}
