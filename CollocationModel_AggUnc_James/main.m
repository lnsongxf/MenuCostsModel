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

%%
% Add CompEcon package
% p = genpath('E:\Dropbox\Economics\Matlab_codes\CompEcon');
p = genpath('C:\Users\James\Dropbox\Economics\Matlab_codes\CompEcon');
addpath(p);

% cd('E:\Dropbox\Economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\Collocation_AggUnc_James')
cd('C:\Users\James\Dropbox\Economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\Collocation_AggUnc_James')

close all
clear
clc

%% Set all options

% Things to do
options.solvecL       = 'Y';      % Solve only c and L 
options.solveeq       = 'Y';      % Solve equilibrium
options.plotpolicyfun = 'Y';            % If Y, plot policy functions

options.polfun      = 'Y';
options.sim         = 'Y';      % Solve simulation
options.solveIRF    = 'N';
% options.solveKS     = 'N';      % Solve Krussel-Smith

% Model options 
options.AR1            = 'N';      % Approx continuous AR1 process, if 'N', choose discretization
options.discmethod     = 'R';      % If 'T' use Tauchen, if 'R' use Rouwenhurst

% [NOT NEEDED HERE - PERHAPS CHANGE OPTIONS]
options.MC          = 'Y';        % If 'N' then menu costs are zero

% Tolerances, iterations
options.Nbell       = 5;        % Number of Bellman (Contraction) iterations
options.Nnewt       = 15;       % Maximum number of Newton steps
options.tolc        = 1e-8;     % Tolerance on value functions
options.tolgolden   = 1e-6;     % Tolerance for golden search
options.itermaxL    = 5000;     % Maximum iterations to find stationary dist L
options.tolL        = 1e-11;    % Tolerance on L

% Print / plot 
options.print       = 'Y';      % Print out c-solution convergence
options.plotSD      = 'N';      % Plot stationary distribution while solving equilibrium

%% Statespace parameters
glob.n          = [7,5,3,3];   % Number of nodes in each dimension: [Np,Na,Nm,Ny]
glob.nf         = [100,100];   % Number of points for pP and a in histogram L
glob.curv       = 1;           % Curvature for pP (1 is no curvature, <1 is curvature)
glob.spliorder  = [3,1,1,1];   % Order of splines (use linear if exogenous vars are discrete (not AR1))
glob.pPmin       = 0.75;       % Lower bound on real price
glob.pPmax       = 1.25;        % Upper bound on real price

% % [ONLY USED FOR CONTINUOUS AR(1) APPROX]
% glob.Ne1        = 30;            % # of approx nodes of AR(1) iid shock in Expectation
% glob.Ne2        = 200;          % # of approx nodes of AR(1) iid shock in Approx of Q
% glob.plb        = 0.001;        % Lower bound on probability of iid shocks to AR(1) process, upper bound is 1-plb
% glob.palb       = 0.005;        % Lower bound on probability of z

% NOTE: Ne1 and Ne2 being very large only makes 'setup.m' run for longer,
% and not very much longer at that. They don't increase the dimensionality
% of *any* matrix


if glob.spliorder(2)>1 && strcmp(options.AR1,'N')
    fprintf('Error 1\n')
    return %     break
end

%% Model parameters

param.beta      = 0.99;     % discount factor
param.delta     = 0.3;      % relative weight of labor-to-consumption in utility
param.sigma     = 1;        % risk aversion coefficient
param.phielas   = 0.5;  % inveser labour supply elasticity
param.theta     = 5;        % elasticity of substitution
param.alpha     = 2/3;      % returns to labour
param.rhoa      = 0.35;     % persistence of productivity
param.sigzeta   = 0.225;      % stddev of productivity shocks
param.Phicost   = 0.156;    % menu cost in labour units
param.mu        = 0.006;    % s.s money growth
param.rhom      = 0.37;     % persistence of money growth
param.sigmaeps  = 0.0048;   % stddev of money growth shocks
param.tauc      = 0.005;    % tolerance for forecasting rule
param.n         = 5000;     % number of firms
param.T         = 96;       % simulation length
param.S         = 25;       % simulations for computing forecasting coeffs
param.s         = 100;      % simulations for moment computations

% Initial guess for state evolution equation
param.b0 = 0.001;
param.b1 = 0.5;
param.b2 = 0.1;

%% Setup problem

% NOTE: Don't forget to set up the fine grid inside 'setup'
fprintf('Setup\n');
[param,glob]    = setup(param,glob,options);      
fprintf('Setup complete\n');

%% Solve only c and L for a given price p

% UP TO HERE!!!!!!

switch options.solvecL
    case 'Y'
        Y                   = 1;  % Conjectured value of output, Y
        P                   = 1;  % Conjectured value of price level, P
        options.cresult     = [];   % Holds previous solution for c. Empty in this case.
        eq                  = solve_cL(Y,P,param,glob,options);  
        fprintf('Yin = %1.2f,\tYout = %1.2f\n',Y,eq.Y);
        fprintf('Pin = %1.2f,\tPout = %1.2f\n',P,eq.P);
        fprintf('--------------------------------------');
end

%% Solve equilibrium
switch options.solveeq
    case 'Y'
        options.tolp        = 0.0001;           % Tolerance on price
        options.Ylb         = 0.5;              % Output lower bound
        options.Yub         = 10;               % Output upper boud
        options.itermaxp    = 30;               % Max iterations of bisection
        options.eqplot      = 'Y';
        options.eqprint     = 'Y';
        options.print       = 'N';
        options.Loadc       = 'Y';              % For new guess of p use old c as starting guess
        options.plotSD      = 'Y';              % If Y plot steady state distribution
        options.fontsize    = 12;
        options.plotpolicyfun = 'N';            % If Y, plot policy functions
        eq                  = solve_eq(param,glob,options);
end

save TEMP

% % % %% Policy functions
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
