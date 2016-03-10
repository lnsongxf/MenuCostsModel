%------------------------------
% 
%
% Original: Simon Mongey
%--------------------------------
p = genpath('E:\Dropbox\Economics\Matlab_codes\CompEcon');
% p = genpath('C:\Users\James\Dropbox\economics\Matlab_codes\CompEcon');
addpath(p);
clc;clear;
cd('E:\Dropbox\Economics\2015_2016_material\QuantMacro_Violante\Simons_Stuff')
% cd('C:\Users\James\Dropbox\economics\2015_2016_material\QuantMacro_Violante\Simons_Stuff')
dbstop if error

%% Set all options

% Things to do
options.solvecL     = 'N';      % Solve only c and L 
options.solveeq     = 'N';      % Solve equilibrium
options.polfun      = 'N';
options.sim         = 'N';      % Solve simulation
options.solveIRF    = 'N';
options.solveKS     = 'N';      % Solve Krussel-Smith

% Model options
options.AR1         = 'N';      % Approx continuous AR1 process, if 'N' then Rouwenhurst discretization
options.ACdown      = 'N';      % If 'Y' then adjustment costs also for downwards movement
options.AC          = 'Y';      % If 'N' then adjustment costs are zero

% Tolerances, iterations
options.Nbell       = 2;        % Number of Bellman (Contraction) iterations
options.Nnewt       = 15;       % Maximum number of Newton steps
options.tolc        = 1e-8;     % Tolerance on value functions
options.tolgolden   = 1e-6;     % Tolerance for golden search
options.itermaxL    = 5000;     % Maximum iterations to find stationary dist L
options.tolL        = 1e-11;    % Tolerance on L

% Print / plot 
options.print       = 'N';      % Print out c-solution convergence
options.plotSD      = 'N';      % Plot stationary distribution while solving equilibrium

%% Statespace parameters
glob.n          = [4,4];        % Number of nodes in each dimension
glob.nf         = [100,100];    % Number of points for k and z in histogram L
glob.curv       = 0.4;          % Curvature for k (1 is no curvature)
glob.spliorder  = [3,3];        % Order of splines (always use linear if productivity is discrete (not AR1))
glob.kmin       = 0.001;        % Lower bound on capital
glob.kmax       = 20.0;         % Upper bound on capital
glob.pzlb       = 0.005;        % Lower bound on probability of z
glob.Ne1        = 30;            % # of approx nodes of AR(1) iid shock in Expectation
glob.plb        = 0.001;        % Lower bound on probability of iid shocks to AR(1) process, upper bound is 1-plb
glob.Ne2        = 200;          % # of approx nodes of AR(1) iid shock in Approx of Q

% NOTE: Ne1 and Ne2 being very large only makes 'setup.m' run for longer,
% and not very much longer at that. They don't increase the dimensionality
% of *any* matrix

if glob.spliorder(2)>1 && strcmp(options.AR1,'N')
    fprintf('Error 1\n')
    return %     break
end

%% Model parameters
% A. Outside parameters
glob.beta       = 0.99; 
glob.rhoz       = 0.95;         
glob.sige       = 0.05;
glob.A          = 1.00; 

% B. Inside parameters (those you might want to change in a calibration exercise)
param.psi       = 2.4;           
param.alpha     = 0.256;
param.nu        = 0.640;
param.delta     = 0.069;        
param.eta       = 0.5;

%% Setup problem
fprintf('Setup\n');
[param,glob]    = setup(param,glob,options);      
fprintf('Setup complete\n');

%% Solve only c and L for a given price p
switch options.solvecL
    case 'Y'
        p                   = 1.9;  % Conjectured value of p    
        options.cresult     = [];   % Holds previous solution for c. Empty in this case.
        eq                  = solve_cL(p,param,glob,options);  
        fprintf('pin = %1.2f,\tpout = %1.2f\n',p,eq.p);
end

%% Solve equilibrium
switch options.solveeq
    case 'Y'
        options.tolp        = 0.0001;             % Tolerance on price
        options.plb         = 0.5;              % Price lower bound
        options.pub         = 10;               % Price upper boud
        options.itermaxp    = 30;               % Max iterations of bisection
        options.eqplot      = 'N'; 
        options.eqprint     = 'Y'; 
        options.print       = 'N';
        options.Loadc       = 'Y';              % For new guess of p use old c as starting guess
        options.plotSD      = 'Y';              % If Y plot steady state distribution
        eq                  = solve_eq(param,glob,options); 
end

save TEMP

%% Policy functions
load TEMP
switch options.polfun
    case 'Y'
        %__________________________________________________________________
        % 1. Varying Z for 3 levels of K
        kbar    = eq.L'*glob.sf(:,1);
        kplot   = [.5*kbar,kbar,1.25*kbar]';
        zplot   = glob.zgridf;
        kpolZ   = [];
        ipolZ   = [];
        for i = 1:numel(kplot);
            s           = gridmake(kplot(i),zplot);
            glob.Phi_Z  = splibas(glob.zgrid0,0,glob.spliorder(2),s(:,2));
            v           = solve_valfunc(eq.c,s,1,param,glob,options,1);
            kpolZ(:,i)  = v.Kp;
            ipolZ(:,i)  = v.I;
        end
        %__________________________________________________________________
        % 2. Varying K for 3 levels of Z
        zbar    = eq.L'*glob.sf(:,2);
        zplot   = [.5*zbar,zbar,1.25*zbar]';
        kplot   = glob.kgridf;
        kpolK   = [];
        ipolK   = [];
        for i = 1:numel(zplot);
            s           = gridmake(kplot,zplot(i));
            glob.Phi_Z  = splibas(glob.zgrid0,0,glob.spliorder(2),s(:,2));
            v           = solve_valfunc(eq.c,s,1,param,glob,options,1);
            kpolK(:,i)  = v.Kp;
            ipolK(:,i)  = v.I;
        end
        %__________________________________________________________________
        figure(round(1000*rand));
        subplot(2,2,1);plot(glob.kgridf,kpolK);title('1A. Capital');xlabel('Capital');legend('Low Z','Med Z','High Z');grid on;
        subplot(2,2,2);plot(glob.kgridf,ipolK);title('1B. Investment');xlabel('Capital');legend('Low Z','Med Z','High Z');grid on;
        subplot(2,2,3);plot(glob.zgridf,kpolZ);title('2A. Capital');xlabel('Productivity');legend('Low K','Med K','High K');grid on;
        subplot(2,2,4);plot(glob.zgridf,ipolZ);title('2B. Investment');xlabel('Productivity');legend('Low K','Med K','High K');grid on;
end

% Interesting: Constprod + s0

%% IRF
% Note: An IRF is only possible from the following (k,z)
% z must be the unconditional mean of the productivity process
% k must be the long-run capital associated with that z under no-shocks
load TEMP
% options.solveIRF = 'Y';
switch options.solveIRF
    case 'Y'
        if ~strcmp(options.AR1,'Y')
            return
        end;
        %__________________________________________________________________
        Ez                  = glob.Psszf'*glob.zgridf; 
        %__________________________________________________________________
        % 1. Solve for long-run capital for z=E[z]=1 under no shocks
        options.simsolve    = 'N';
        options.constprod   = 'Y';
        options.IRF         = 'N';
        options.simplot     = 'N';
        T                   = 200;
        s0                  = gridmake(median(glob.kgridf),1);
        sim                 = solve_simulation(s0,T,eq,1,param,glob,options);  
        Kss                 = sim.Kt(end,:)'; 
        %__________________________________________________________________
        % 2. IRF - Start at [Kss,1] and solve with AR(1) shock in period 1
        options.simsolve    = 'Y';
        options.constprod   = 'N';
        options.IRF         = 'Y';
        options.simplot     = 'Y';
        T                   = 100;
        s0                  = [Kss,1];
        sim                 = solve_simulation(s0,T,eq,1,param,glob,options); 
end

%% Simulate model
switch options.sim
    case 'Y'
        %__________________________________________________________________
        % Initial vector of firms (can try different things)
        % 1. All k's average productivity
        s01                 = repmat(gridmake(glob.kgrid,glob.zgrid(round(glob.Nz/2))),500,1);
        % 2. Small k, all productivities - Growth
        s02                 = repmat(gridmake(glob.kgridf(10),glob.zgrid),1,1);
        % 3. Small k, all productivities - Growth
%         s04                 = repmat(gridmake(Kss,1*exp(-glob.sige)),5000,1);
        s04                 = repmat(gridmake(min(glob.kgrid),min(glob.zgridf)),10000,1);
        % 3. One individual point
        kbar                = eq.L'*glob.sf(:,1);
        zbar                = eq.L'*glob.sf(:,2);
        s03                 = [kbar,zbar];
        % Choose s0
        s0                  = s04; 
        %__________________________________________________________________
        % Number of periods
        T                   = 400;
        % Options
        options.simsolve    = 'N';          % Solve problem each period (EXACT)
                                            % If N, then linearly
                                            % interpolate policies (FAST) -
        options.constprod   = 'N';          % Keep productivity constant
        options.simplot     = 'Y';
        options.IRF         = 'N';          % Fix as 'N'
        %__________________________________________________________________
        % Solve
        options.AC          = 'Y';
        sim                 = solve_simulation(s0,T,eq,1,param,glob,options); 
end

% save TEMP3
%% Solve Krussell-Smith
% load TEMP3

options.solveKS = 'Y';

% return

switch options.solveKS
    case 'Y'
        NA              = 3;
        NK              = 3;
        glob.spliorder  = [3,1,3,1]; 
        %____________________________________________________
        % Nodes for k and z
        kgrid0          = glob.kgrid0;
        zgrid0          = glob.zgrid0;
        %____________________________________________________
        % Nodes for A
%         rhoA            = glob.rhoz;                     glob.rhoA   = rhoA;  
%         sigA            = glob.sige;                     glob.sigA   = sigA;  
        rhoA            = 0.999;                     glob.rhoA   = rhoA;     %%%% DEBUG
        sigA            = 0.005;                     glob.sigA   = sigA;     %%%% DEBUG
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





