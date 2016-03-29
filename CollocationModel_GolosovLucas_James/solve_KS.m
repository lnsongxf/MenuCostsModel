function [DMt,Yt,Pt] = solve_KS(cKS,eq,param,glob,options)
%SOLVE_KS Implements the Krussel Smith Algorithm for the model
%-------------------------------------------------
%   Simulates the aggregate states, solves the model at each step, 
%   
%   INPUTS
%   - cKS       = parameters for the law of motion (lnY = b0 + b1lnY-1 + b2 Dm)
%   - eq        = ???
%   - param     = model parameters
%   - glob      = global variables, including those used in the approximating functions
%   - options   = options, including continuous vs. discrete Markov, Tauchen vs. Rouwenhurst, 
%
%   OUTPUTS
%   - At        = Simulated path for aggregate productivity, A
%   - pt        = Simulated path for price level, p
%   - Kt        = Simulated path for aggregate capital, K
%-------------------------------------------------

%% Unpack
ns          = size(glob.s,1); 
T           = options.T;

%% Initialise guesses (if val.cresult has an old guess in it, use that)
cKold       = zeros(ns,1);
cCold       = zeros(ns,1);
cEold       = zeros(ns,1);
cold        = [cKold;cCold;cEold];


%% Solve value function problem
% Define equilibrium output and price
glob.P      = eq.P;     %%%% DEBUG
glob.Y      = eq.Y;     %%%% DEBUG
c           = solve_cKS(cold,cKS,param,glob,options);   % Get new val func approx coefficients

%% Simulate 

Pt  = zeros(T,1);   
Yt  = zeros(T,1);   
DMt = zeros(T,1);   
Mt  = zeros(T,1);

%% Initialize aggregate state

Lt          = eq.L;         % Set initial distribution at no agg. uncertainty stationary dist
Pt(1)       = eq.P;         % Initialize at eqm P under no agg. uncertainty 
Yt(1)       = eq.Y;         % Initialize at eqm Y under no agg. uncertainty 
DMt(1)      = median(glob.Mgrid);  % Initialize at median state (SS?)
Mt(1)       = 1;            % Set initial money stock to 1 (note this does not need to be on any grid)

%% Draw aggregate shocks from Markov chain (i.e. for money growth shocks)
% Simulate the true AR(1) and check to make sure it's within the gridspace, 
% adjusting to endpoints if simulation falls above or below the grid 

if strcmp(options.seed,'Y')
    rng(222);       % Set random number generator seed
end
DMshocks = randn(T,1);

for t = 2:T;
   DMt(t) = exp(param.mu*(1 - param.rhom) + param.rhom*log(DMt(t-1)) + param.sigmaeps*DMshocks(t));
   DMt(t) = max(min(DMt(t),max(glob.Mgrid)),min(glob.Mgrid));   % Adjust to keep within gridspace
   
   Mt(t)   = Mt(t-1)*DMt(t);         % Update level of money supply
end

% % % figure;
% % % subplot(1,2,1)
% % % plot(DMt)
% % % title('Money growth')
% % % set(gca,'fontsize',options.fontsize)
% % % subplot(1,2,2)
% % % plot(Mt)
% % % title('Money stock')
% % % set(gca,'fontsize',options.fontsize)


%% Simulation/KS updating

for t=1:T 
    
    %% 1. Construct current state vector [Np,Na,Nm,Ny]
    % Construct a state vector with a KNOWN value of output and money growth 
    % at a particular time, t. 
    
    st              = gridmake(glob.pPgridf,glob.agridf,DMt(t),Yt(t));

    %% 2. Compute basis matrix for continuation values
    
    % Computing new basis matrices based on the current state vector, st
    glob.Phi_pP          = splibas(glob.pPgrid0,0,glob.spliorder(1),st(:,1));       % Basis matrix for price grid
    glob.Phi_A           = splibas(glob.agrid0,0,glob.spliorder(2),st(:,2));        % Used in Bellman / Newton computing expected values
    glob.Phi_M           = splibas(glob.Mgrid0,0,glob.spliorder(3),st(:,3));        % Used in Bellman / Newton computing expected values
    glob.Phi_Y           = splibas(glob.Ygrid0,0,glob.spliorder(4),st(:,4));        % Used in Bellman / Newton computing expected values
    glob.Phi  = dprod(glob.Phi_Y,dprod(glob.Phi_M,dprod(glob.Phi_A,glob.Phi_pP)));
    
    
    %% 3. Solve partial equilibrium           

for iterY = 1:50;
    Yin         = Yt(t); 
    Pin         = Pt(t);
    Y           = ones(size(st(:,1)))*Yin;
    v           = solve_valfuncKS(c,st,Y,param,glob,options,1);    
    Pout        = ( Lt'*(v.pPdist*Pin).^(1 - param.theta) ).^(1/(1-param.theta));
    Yout        = Mt(t)/Pout;
    
%     Pt(t)       = glob.dampP*Pout + (1-glob.dampP)*Pin;
    Yt(t)       = glob.damp*Yout + (1-glob.damp)*Yin;    
    %         keyboard
    if strcmp(options.eqprint,'Y')
        fprintf('%2i. Yin:\t%2.4f\tYout:\t%2.4\n',iterY,Yin,Yout);
        fprintf('%2i. Pin:\t%2.4f\tPout:\t%2.4f\n',iterY,Pin,Pout);
        fprintf('---------------------------------\n');
    end
    if abs(Yin-Yout)<options.tolYks
        fprintf('Solved for output, t = %1i\n', t)
        break
    end    
end

fprintf('%2i. Yin:\t%2.4f\tYout:\t%2.4f\n',iterY,Yin,Yout);
fprintf('%2i. Pin:\t%2.4f\tPout:\t%2.4f\n',iterY,Pin,Pout);
fprintf('---------------------------------\n');
    Yt(t)           = Yout;
    Pt(t)           = Pout;
    
    %% 4. Compute next period states [Np,Na,Nm,Ny]
    % 1. Distribution of prices (Only need dist over idiosyncratic states
    % since the curren aggregate states are known. 
    
    pPdist = max( min(v.pPdist,max(glob.pPgridf)), min(glob.pPgridf)); 
    fspaceergpP      = fundef({'spli',glob.pPgridf,0,1});   % Linear interpolant
    QpP              = funbas(fspaceergpP,pPdist);
    Ltnew               = dprod(glob.QA,QpP)'*Lt;
    
    if strcmp(options.eqprint,'Y')
        fprintf('dL = %2.3e\n',norm(Ltnew-Lt)/norm(Lt));
    end
    Lt              = Ltnew;
    % 3. Productivity
    %%%% Already computed 
    
    % 4. Output: guess for t+1 is the same as last period eqm
    Yt(t+1)  = Yt(t);
    
    % 5. Price: guess for t+1 is the same as last period eqm
    Pt(t+1)  = Pt(t); 

    %% Plot paths
    if strcmp(options.simplot,'Y')
        figure(888);
        subplot(2,2,1);
        plot(DMt(1:t));   %,'bo','markersize',6,'markerfacecolor','b');
        title('Money growth, \Delta M_t');
        grid on;
        subplot(2,2,2);
        plot(Yt(1:t)); %,'bo','markersize',6,'markerfacecolor','b');
        title('Output, Y_t');
        grid on;
        subplot(2,2,3);
        plot(Pt(1:t),'color','b'); % ,'bo','markersize',6,'markerfacecolor','b');
        hold on
        plot(Mt(1:t),'color','r');
        title('Prices (blue) vs. Money (red), P_t, M_t');
        grid on;
        drawnow;
    end

    
end
    
   
    
    
%% Tidy output    
DMt = DMt(1:T);
Yt  = Yt(1:T);
Pt  = Pt(1:T); 

    
    
end














