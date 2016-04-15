function [cKS,R2,sim] = simulate_KS(c,v,eq,param,glob,options)

T           = options.T;

% Distribution over states
Lt          = zeros(glob.nf(1)*glob.nf(2),T);
Lt(:,1)     = eq.L;         % Set initial distribution at no agg. uncertainty stationary dist
Lt(:,2)     = eq.L;         % Set initial distribution at no agg. uncertainty stationary dist

% Policy functions 
pPdist      = zeros(glob.nf(1)*glob.nf(2),T);
ind         = zeros(glob.nf(1)*glob.nf(2),T);

% initial price states
p_state     = zeros(glob.nf(1)*glob.nf(2),options.T);

% Draw money growth shocks
DMt         = zeros(T+1,1);  % Dmt = (1-rho)mu + rho Dmt-1 + epst 
Mt          = zeros(T+1,1);
DMt(1)      = exp(param.mu);  % Initialize at median state (SS?)
Mt(1)       = 1;            % Set initial money stock to 1 (note this does not need to be on any grid)

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
%     DMshocks(2) = 1;
end

% Initialize DM and M at t=0
for t = 2:T;
    DMt(t)  = exp(param.mu*(1 - param.rhom) + param.rhom*log(DMt(t-1)) + param.sigmaeps*DMshocks(t));
    DMt(t)  = max(min(DMt(t),max(glob.Mgrid)),min(glob.Mgrid));   % Adjust to keep within gridspace
    Mt(t)   = Mt(t-1)*DMt(t);         % Update level of money supply
end


% Aggreagtes
Pt          = zeros(T,1);   
Yt          = zeros(T,1);  
Pt(1)       = (1/eq.Y)*Mt(1);         % Initialize at eqm P under no agg. uncertainty 
Yt(1)       = Mt(1)/Pt(1);         % Initialize at eqm Y under no agg. uncertainty 

%% Simulation/KS updating

for t=2:T 
    % Initial guess for period t is same as period t-1
    Yt(t)  = Yt(t-1);   
    Pt(t)  = Pt(t-1); 
    
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
    Pout        = ( Lt(:,t-1)'*(v.pPdist*Pin).^(1 - param.theta) ).^(1/(1-param.theta));   % WHY Lt-1?
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
    p_state(:,t)    = st(:,1).*(1/pi);
    pPdist(:,t)     = v.pPdist;
    
    %% 4. Compute next period states [Np,Na,Nm,Ny]
    % 1. Distribution of prices (Only need dist over idiosyncratic states
    % since the current aggregate states are known. 
    
    pP              = max( min(v.pPdist,max(glob.pPgridf)), min(glob.pPgridf)); 
    fspaceergpP     = fundef({'spli',glob.pPgridf,0,1});   % Linear interpolant
    QpP             = funbas(fspaceergpP,pP);
    Ltnew           = dprod(glob.QA,QpP)'*Lt(:,t-1);
    Lt(:,t)         = Ltnew;

    %% Plot paths
    if strcmp(options.simplot,'Y')
        figure(888);
        subplot(2,2,1);
        plot(DMt(1:t)); 
        title('Money growth, \Delta M_t');
        grid on;
        subplot(2,2,2);
        plot(Yt(1:t));
        title('Output, Y_t');
        grid on;
        subplot(2,2,3);
        plot(Pt(1:t),'color','b');
        hold on
        plot(Mt(1:t),'color','r');
        title('Prices (blue) vs. Money (red), P_t, M_t');
        grid on;
        drawnow;
    end

    
end
    
%% Conduct KS regression step
% HOW TO EXPLAIN THE TIMING OF THE DMt RELATIVE TO Yt? 
X = [ones(T-options.burn-1,1) log(Yt(options.burn:T-2)) log(DMt(options.burn+2:T))];
y = log(Yt(options.burn+1:T-1));
[b ,~,~,~,stats] = regress(y,X);
cKS = b;
R2     = stats(1);
    
%% Pack up output  

sim.Yt           = Yt;
sim.Pt           = Pt;   
sim.Lt           = Lt;
sim.ind          = ind;   
sim.p_state      = p_state;
sim.pPdist       = pPdist;
sim.DMt          = DMt;
sim.Mt           = Mt; 


end