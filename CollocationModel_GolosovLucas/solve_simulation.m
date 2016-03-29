function sim = solve_simulation(s0,eq,val,param,glob,options)
%SOLVE_SIMULATION Simulates the model
%-------------------------------------------------
%   This file simulates the solved model. This version only solves with no
%   aggregate uncertainty.
%
%   INPUTS
%   - s0        = Initial vector/starting position of firms e.g. starting
%   at same point, all starting at average, starting all over the grid
%   - T         = Number of simulations
%   - eq        = solved approximating coefficients and equilibrium objects, 
%       including distribution, aggregate price, aggregate output 
%   - val       = Unused??
%   - param     = model parameters
%   - glob      = global variables, including those used in the approximating functions
%   - options   = options, including continuous vs. discrete Markov, Tauchen vs. Rouwenhurst, 
%-------------------------------------------------

%% UNPACK
agrid       = glob.agrid;
c           = eq.c;
fspace      = glob.fspace;
Na          = glob.Na;
aub         = max(agrid);
alb         = min(agrid);
rho         = param.rhoa;
sigma       = param.sigzeta; 
T           = options.T;
M           = 1;    % keep money stock constant at 1.

%% SIMULATION

% Initial state
st          = s0;

% Storage (firms)
ni          = size(st,1); 
pPt         = zeros(T,ni);
at          = zeros(T,ni);  
yt          = zeros(T-1,ni);
nt          = zeros(T-1,ni);
wPt          = zeros(T-1,ni);

% Storage (aggregate), not actually used again: just kept track of.
Pt         = zeros(T,1);                   
Yt         = zeros(T,1);                   
Nt         = zeros(T,1);                   

% Initial productivities
for i = 1:size(st,1);
    st(i,2) = find(glob.agrid==st(i,2));    % Convert into indexes
end
at(1,:)     = glob.agrid(st(:,2)');

% Initial prices/output
pPt(1,:)   = st(:,1)';
Pt(1,:)    = ( pPt(1,:).^(1-param.theta)*(1/ni)*ones(ni,1) )^(1/(1-param.theta));  % Agg price, not actually used just kept track of.
Yt(1,:)    = M/Pt(1,:);

% Policy interpolant
cpP         = funfitxy(fspace,glob.sf,eq.v.pPstar); 
cy          = funfitxy(fspace,glob.sf,eq.v.ystar);
cn          = funfitxy(fspace,glob.sf,eq.v.nstar);
cwP         = funfitxy(fspace,glob.sf,eq.v.wPstar);
% cv1         = funfitxy(fspace,glob.sf,eq.v.v1);

% Draw productivity shocks
rng(999);
PS          = rand(ni,T); 

% Option: Constant productivity simulation
P       = glob.P;        

% Simulation
for t = 2:T;
    % 1. Update state
    stm1                = st;
    
    % 2. Pull out state variables
     pP               = stm1(:,1);
     A               = glob.agrid(stm1(:,2));       % Put productivities onto the grid
    
    % 3. Compute values and policies given states stm1
    % Recompute basis matrix for continuation values
    glob.Phi_A      = splibas(glob.agrid0,0,glob.spliorder(2),A);
    Phi_pP          = splibas(glob.pPgrid0,0,glob.spliorder(1),pP);
    glob.Phi       = dprod(glob.Phi_A,Phi_pP); 
    
    % Solve problem (include 1 at end so don't compute E(v))
    v               = solve_valfunc_noagg(c,[pP,A],Yt(t,:),param,glob,options,1);

    % 4. Evolve A stochastically
    Ptrans              = P(st(:,2),:);             % Transition probabilities for each individual
    Pcumsum             = cumsum(Ptrans')';
    x                   = PS(:,t);                  % Draw productivity shocks
    x                   = repmat(x,1,Na);
    PPP                 = (Pcumsum>=x);
    PPP                 = cumsum(PPP')';            % Numbers elements greater than x
    [~,J]               = find(PPP==1);
    Jp                  = J;
    st                  = [v.pPstar,Jp];
    at(t,:)             = agrid(st(:,2)');
    % Record *next period state* for pP
    pPt(t,:)             = st(:,1)';   
    % Record period t objects, these depend on...
    yt(t-1,:)           = v.ystar';
    nt(t-1,:)           = v.nstar';
    wPt(t-1,:)          = v.wPstar';
    % Aggregates
    Pt(t)               = ( pPt(t,:).^(1-param.theta)*(1/ni)*ones(ni,1) )^(1/(1-param.theta)); 
    Yt(t)               = sum(yt(t-1,:));
    Nt(t)               = sum(nt(t-1,:));    
    
    fprintf('t: %3i\tAgg P: %1.3f\n',t,Pt(t)); 
end

%% PACK UP OUTPUT OBJECT {sim}
sim.pPt     = pPt;
sim.at      = at;
sim.nt      = nt;
sim.wPt     = wPt;
sim.Pt      = Pt;
sim.Yt      = Yt;
sim.Nt      = Nt;


%% Figures
if strcmp(options.simplot,'Y')

    varlist     = {'Kt','Zt','Nt','It'};
    titlelist   = { 'A. Capital','B. Productivity','C. Labor','D. Investment'}; 

%     h4 = figure(round(rand*1000));
%     tvec = 0:1:T-1;
%     for x = 1:4;
%         subplot(2,2,x);
%         eval(sprintf('X=%s;',char(varlist(x))));
%         if x<=2
%             plot(tvec,X,'-','LineWidth',2);grid on;hold on; 
%         else
%             plot(tvec(1:end-1),X,'-','LineWidth',2);grid on;hold on; 
%         end
%         xlabel('Periods','fontsize',16);
%         title(titlelist(x),'fontsize',16);
%         set(gca,'fontsize',16);
%         xlim([0,T-1]);
%     end
    
    varlist     = {'Kat','Cat','Nat','Iat'};
    titlelist   = { 'A. Capital','B. Consumption','C. Labor','D. Investment'}; 
    
    h5 = figure(round(rand*1000));
    tvec = 0:1:T-1;
    for x = 1:4;
        subplot(2,2,x);
        eval(sprintf('X=%s;',char(varlist(x))));
        if x<=1
            plot(tvec,X,'-','LineWidth',2);grid on;hold on; 
        else
            plot(tvec(1:end-1),X(2:end),'-','LineWidth',2);grid on;hold on;  
        end
        xlabel('Periods','fontsize',16);
        title(titlelist(x),'fontsize',16);
        set(gca,'fontsize',16);
        xlim([0,T-1]);
    end
    
end

