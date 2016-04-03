function [coeffs,R2,paths] = simulate_KS(c,v,cKS,eq,param,glob,options)

% Distribution over states
Lt         = zeros(glob.nf(1)*glob.nf(2),T);
Lt(:,1)     = eq.L;         % Set initial distribution at no agg. uncertainty stationary dist
Lt(:,2)     = eq.L;         % Set initial distribution at no agg. uncertainty stationary dist

% Policy functions 
ind        = zeros(glob.nf(1)*glob.nf(2),T);
pPdist     = zeros(glob.nf(1)*glob.nf(2),T);

%% Simulate 

Pt  = zeros(T,1);   
Yt  = zeros(T,1);   
DMt = zeros(T+1,1);  % Dmt = (1-rho)mu + rho Dmt-1 + epst 
Mt  = zeros(T+1,1);


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