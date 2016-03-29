function [DMt,Yt,Pt,Lt,ind,pPdist] = solve_KS(cKS,eq,param,glob,options)
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
%
%   NEED TO FIX THE TIMING OF THE MONEY SHOCKS!!!!
%
%
%% Run setup file again with new cKS params
glob = setup_agg(param,glob,cKS,options);


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
DMt = zeros(T+1,1);  % Dmt = (1-rho)mu + rho Dmt-1 + epst 
Mt  = zeros(T+1,1);
ind        = zeros(glob.nf(1)*glob.nf(2),T);
pPdist     = zeros(glob.nf(1)*glob.nf(2),T);
Lt         = zeros(glob.nf(1)*glob.nf(2),T);

%% Draw aggregate shocks from Markov chain (i.e. for money growth shocks)
% Simulate the true AR(1) and check to make sure it's within the gridspace, 
% adjusting to endpoints if simulation falls above or below the grid 

if strcmp(options.irf, 'N')
    if strcmp(options.seed,'Y')
        rng(222);       % Set random number generator seed
    end
    DMshocks = randn(T,1);
elseif strcmp(options.irf, 'Y')
    DMshocks = zeros(T,1);
    DMshocks(1) = 1;
end

% Initialize DM and M at t=0
DMt(1)    = exp(param.mu);  % Initialize at median state (SS?)
Mt(1)     = 1;            % Set initial money stock to 1 (note this does not need to be on any grid)
for t = 2:T+1;
    DMt(t) = exp(param.mu*(1 - param.rhom) + param.rhom*log(DMt(t-1)) + param.sigmaeps*DMshocks(t-1));
    DMt(t) = max(min(DMt(t),max(glob.Mgrid)),min(glob.Mgrid));   % Adjust to keep within gridspace
    Mt(t)  = Mt(t-1)*DMt(t);         % Update level of money supply
end


%% Initialize aggregate states

Lt(:,1)     = eq.L;         % Set initial distribution at no agg. uncertainty stationary dist
Lt(:,2)     = eq.L;         % Set initial distribution at no agg. uncertainty stationary dist
Pt(1)       = (1/eq.Y)*Mt(1);         % Initialize at eqm P under no agg. uncertainty 
Pt(2)       = Pt(1); 
Yt(1)       = Mt(1)/Pt(1);         % Initialize at eqm Y under no agg. uncertainty 
Yt(2)       = Yt(1);

%% Simulation/KS updating

for t=2:T 
    
    % 3. Solve partial equilibrium           
for iterY = 1:50;
    Yin         = Yt(t); 
    Pin         = Pt(t); 
    st          = gridmake(glob.pPgridf,glob.agridf,DMt(t),Yin);
    pi = Pin/Pt(t-1);
    
    % Computing new basis matrices based on the current state vector, st
    glob.Phi_A           = splibas(glob.agrid0,0,glob.spliorder(2),st(:,2));        % Used in Bellman / Newton computing expected values
    glob.Phi_M           = splibas(glob.Mgrid0,0,glob.spliorder(3),st(:,3));        % Used in Bellman / Newton computing expected values
    glob.Phi_Y           = splibas(glob.Ygrid0,0,glob.spliorder(4),st(:,4));        % Used in Bellman / Newton computing expected values
    glob.Phi             = funbas(glob.fspace,st);
    
    v           = solve_valfuncKS(c,[st(:,1)*(1/pi) st(:,2:end)],param,glob,options,1);    
    Pout        = ( Lt(:,t)'*(v.pPdist*Pin).^(1 - param.theta) ).^(1/(1-param.theta));
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
    
    % Save outputs
    Yt(t)           = Yout;
    Pt(t)           = Pout;    
    ind(:,t)        = v.ind;   
    
    %% 4. Compute next period states [Np,Na,Nm,Ny]
    % 1. Distribution of prices (Only need dist over idiosyncratic states
    % since the current aggregate states are known. 
    
    pPdist(:,t)      = max( min(v.pPdist,max(glob.pPgridf)), min(glob.pPgridf)); 
    fspaceergpP      = fundef({'spli',glob.pPgridf,0,1});   % Linear interpolant
    QpP              = funbas(fspaceergpP,pPdist(:,t));
    Ltnew            = dprod(glob.QA,QpP)'*Lt(:,t);
    Lt(:,t+1)        = Ltnew;
    
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
        plot(DMt(1:t+1)); 
        title('Money growth, \Delta M_t');
        grid on;
        subplot(2,2,2);
        plot(Yt(1:t));
        title('Output, Y_t');
        grid on;
        subplot(2,2,3);
        plot(Pt(1:t),'color','b');
        hold on
        plot(Mt(1:t+1),'color','r');
        title('Prices (blue) vs. Money (red), P_t, M_t');
        grid on;
        drawnow;
    end

    
end
    
    
%% Tidy output    
DMt = DMt(2:T+1);
Yt  = Yt(1:T);
Pt  = Pt(1:T); 

if nargout>3
    Lt         = Lt(:,1:end);
    ind        = ind(:,1:T);   
    pPdist     = pPdist(:,1:T);
end
    
%% Temp stuff
% 
% ones(1,length(Lt))*Lt
% v.pPstar'*Lt    
    
end














