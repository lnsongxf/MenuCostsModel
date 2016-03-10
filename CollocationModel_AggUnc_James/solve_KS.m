function [At,pt,Kt] = solve_KS(cKS,eq,param,glob,options)

%% Unpack
ns          = size(glob.s,1); 

%% Initialise guesses (if val.cresult has an old guess in it, use that)
c1old       = zeros(ns,1);  
c2old       = zeros(ns,1); 
cold        = [c1old;c2old];

%% Solve value function problem
glob.p      = eq.p;     %%%% DEBUG
c           = solve_cKS(cold,cKS,param,glob,options);

%% Simulate 
T           = 100;
pt          = zeros(T,1);   
At          = zeros(T,1);
Kt          = zeros(T,1);

%% Initialize aggregate state
Lt          = eq.L;
Kt(1)       = eq.Ka; 
At(1)       = 1;

%% Draw aggregate shocks
% switch options.Krun
%     case 'Y'
%         At = ones(T,1);     %%%% DEBUG
%     case 'N'
T = 2000;
        rng(222);
        mu = (1/2)*(glob.sigA^2)/((1+glob.rhoA)*(1-glob.rhoA));
%         mu = 1+(1/2)*(glob.sigA^2)/((1+glob.rhoA)*(1-glob.rhoA));
        for t = 2:T;
            At(t) = exp(mu+glob.rhoA*log(At(t-1)) + glob.sigA*randn);
            At(t) = max(min(At(t),max(glob.Agrid)),min(glob.Agrid));
        end 
% end

% At = ones(T,1);     %%%% DEBUG
options.tolp = 1e-3;

for t=1:options.T;
    fprintf('t = %1i\n',t);
    %% 1. Construct state vector
    st              = gridmake(glob.kgridf,Kt(t),glob.zgridf,At(t));
    
    %% 2. Compute basis matrix for continuation values
    Phi_K           = splibas(glob.Kgrid0,0,glob.spliorder(2),st(:,2));
    Phi_z           = splibas(glob.zgrid0,0,glob.spliorder(3),st(:,3));
    Phi_A           = splibas(glob.Agrid0,0,glob.spliorder(4),st(:,4));
    Phi_zA          = dprod(Phi_A,Phi_z);
    glob.Phi_KzA    = dprod(Phi_zA,Phi_K);  
    
    %% 3. Solve partial equilibrium
    plb             = 0.50;
    pub             = 10.00;
    tictic          = tic; 
%     options.tolp    = Inf;
    for iterp = 1:50;
        pin         = (1/2)*(plb+pub);
%         pin         = eq.p;
        P           = ones(size(st(:,1)))*pin;
        v           = solve_valfuncKS(c,st,P,param,glob,options,1);
        Ca          = Lt'*v.C;
        pout        = 1/Ca; 
        down        = (pin>pout); 
        up          = (pin<pout);
        plb         = up*pin + down*plb;
        pub         = up*pub + down*pin;
%         keyboard 
        if strcmp(options.eqprint,'Y') 
            fprintf('%2i. pin:\t%2.6f\tpout:\t%2.6f\tt:%2.1f\n',iterp,pin,pout,toc(tictic));
        end
        if abs(pin-pout)<options.tolp;fprintf('Solved\n');break;end;
    end
    pt(t)           = pin;
    
    %% 4. Compute next period states
    % 1. Capital
    Kt(t)           = Lt'*v.kp;
    Kt(t)           = max(min(Kt(t),max(glob.Kgrid)),min(glob.Kgrid));
    % 2. Distribution
    kp              = max(min(v.kp,max(glob.kgridf)),min(glob.kgridf));
    fspaceergk      = fundef({'spli',glob.kgridf,0,1});
    Qk              = funbas(fspaceergk,kp);
    Lt              = dprod(glob.QZ,Qk)'*Lt;
%     fprintf('dL = %2.3e\n',norm(Ltnew-Lt)/norm(Lt));
%     Lt              = Ltnew;
    % 3. Productivity
    %%%% Already computed 
    
    %% Plot paths
    figure(888);      
    subplot(2,2,1);plot(At(1:t));title('At');grid on;
    subplot(2,2,2);plot(pt(1:t));title('p_t');grid on;
    subplot(2,2,3);plot(Kt(1:t));title('K_t');grid on;
%     strcmp(options.Krun,'Y') && (t>1)
%     if strcmp(options.Krun,'Y') && (t>1) 
%         fprintf('dK = %1.4e\n',abs(Kt(t)-Kt(t-1)));
%         if abs(Kt(t)-Kt(t-1))<0.0001;return;end;
%         keyboard
%     end
    
end
    
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end














