function eq  = solve_cL(Y,P,param,glob,options)
%SOLVE_CL Solve for value function coefficients and stationary distribution  
%-------------------------------------------------
%   Solves for the value function coeffivient vectors (cK,cC,cE) and the
%   stationary distribution matrix L. Solution is conditional on a
%   conjectured level of equilibrium output, Y. 
%
%   INPUTS
%   - Y         = conjectured value of output, Y
%   - P         = conjectured value of aggregate price level, P
%   - param     = parameters 
%   - glob      = includes state space, function space, approximating functions etc
%   - options   = 
%   OUTPUT
%   - eq        = 
%
%-------------------------------------------------


%% A. Globals 
s           = glob.s; 
sf          = glob.sf;
pPgrid      = glob.pPgrid;
% ns          = size(s,1);

%% B. Compute equilibrium objects that depend on p
% ----- None -----

%% Initialise guesses (if val.cresult has an old guess in it, use that)
cKold       = zeros(glob.Ns,1);
cCold       = zeros(glob.Ns,1);
cEold       = zeros(glob.Ns,1);

% c           = options.cresult;
% if ~isempty(c) && strcmp(options.Loadc,'Y')
%     cKold    = c(1:ns);
%     cCold    = c(ns+1:end);
% cEold       = zeros(ns,1);
% end
cold        = [cKold;cCold;cEold];
totaltic    = tic;

%% Bellman iteration
for citer = (1:options.Nbell)
    glob.citer  = citer;
    % 1. Compute values;
    
    v           = solve_valfunc(cold,s,Y,param,glob,options); 
    % 2. Update c
    cK          = glob.Phi\full(v.vK);      % Note: 'full' re-fills a sparse matrix for computations
    cC          = glob.Phi\full(v.vC);
    cE          = glob.Phi\full(v.vE);    
    c           = [cK;cC;cE];
    % 3. Compute distance and update
    dc          = norm(c-cold)/norm(cold); 
    cold        = c;
    if strcmp(options.print,'Y');
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n',citer,dc,toc(totaltic));
    end
end

%% Newton iterations
if strcmp(options.print,'Y');
    fprintf('~~~~~ Newton iterations ~~~~~\n');
end
eq.flag.cconv = false;
for citer = (1:options.Nnewt)
    % 1. Compute values
    [v,jac]     = solve_valfunc(cold,s,Y,param,glob,options);
    % 2. Update c 
    cKold       = cold(1:glob.Ns); 
    cCold       = cold(glob.Ns+1:2*glob.Ns);
    cEold       = cold(2*glob.Ns+1:end);
    c           = cold - jac\([glob.Phi*cKold - full(v.vK) ;
                               glob.Phi*cCold - full(v.vC) ;
                               glob.Phi*cEold - full(v.vE)]);  
    % 3. Compute distances and update
    dc          = norm(c-cold)/norm(cold);
    cold        = c;
    if strcmp(options.print,'Y');
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n',citer,dc,toc(totaltic));
    end
    % 4. Check convergence
    if (dc<options.tolc)
        eq.flag.cconv = true;
    end
    if eq.flag.cconv
        break
    end
end

%% Plot value functions

if strcmp(options.plotpolicyfun,'Y')
    cK       = c(1:glob.Ns);
    cC       = c(glob.Ns+1:2*glob.Ns);
    
    valK = funbas(glob.fspace,glob.sf)*cK;
    valC = funbas(glob.fspace,glob.sf)*cC;
    valtot = max(valK, valC);
    valtot = reshape(valtot, length(glob.pPgridf), length(glob.agridf), length(glob.Mgrid), length(glob.Ygrid));
    
    figure
    plot(glob.pPgridf, valtot(:,:, 2, 2))
    xlabel('Real price')
    ylabel('Value')
end


%% Solve again on a finer grid for pP
glob.Phi_A      = glob.Phi_Af; 
glob.Phi        = glob.Phif; 
v               = solve_valfunc(c,sf,Y,param,glob,options,1);


%% Compute stationary distribution
pPstar           = min(v.pPstar,max(pPgrid));
fspaceergpP      = fundef({'spli',glob.pPgridf,0,1});
QpP              = funbas(fspaceergpP,pPstar); 
QA              = glob.QA;
Q               = dprod(QA,QpP); 

% [vv,dd]         = eigs(Q');
% dd              = diag(dd);
% Lv              = vv(:,dd==max(dd));
% L               = Lv/sum(Lv);
L               = ones(size(Q,1),1);
L               = L/sum(L);

for itL = (1:options.itermaxL);
    Lnew    = Q'*L;  
    dL      = norm(Lnew-L)/norm(L);  
    if (dL<options.tolL)
        break
    end
    
    if mod(itL,100)==0 
        if strcmp(options.print,'Y')
            fprintf('dL:\t%1.3e\n',dL);
        end
    end
    L       = Lnew;
end

% Plot stationary distribution
if strcmp(options.plotSD,'Y');
    H = figure(888);
%     set(H,'Pos',[1          35        1920         964]);
    JpP  = numel(glob.pPgridf);
    Ja  = numel(glob.agridf);
    La  = kron(eye(Ja),ones(1,JpP))*L; 
    LpP  = kron(ones(1,Ja),eye(JpP))*L;
    % Marginal K
    subplot(2,2,1);
    plot(glob.pPgridf,LpP,'o-');title('Stationary Real Price Dist - LpP');
    grid on;
    % Marginal Z
    subplot(2,2,2);
    plot(exp(glob.agridf),La,'o-');title('Stationary Prod Dist - La');
    grid on;
    eq.LpP   = LpP; 
    eq.La   = La;
    % Joint (pP,A) - Surface plot
    subplot(2,2,4)
    Lmat    = reshape(L,JpP,Ja);
    Amat    = repmat(glob.agridf',JpP,1); 
    pPmat    = repmat(glob.pPgridf,1,Ja);
    if LpP(end)<0.001;
        pPub     = glob.pPgridf(find(cumsum(LpP)>0.98,1,'first'));
    else
        pPub     = max(glob.pPgridf);
    end
    Aub     = glob.agridf(find(cumsum(La)>0.98,1,'first'));
    mesh(Amat,pPmat,Lmat,'LineWidth',2);
    xlabel('Productivity - a');
    ylabel('Real Price - pP');
    title('Joint Distribution');
    xlim([min(glob.agridf),Aub]);
    ylim([min(glob.pPgridf),pPub]); 
    zlim([0,max(max(Lmat))]);
end

%% Compute aggregates and implied p

P = ( pPstar'.^(1-param.theta)*L )^(1/(1-param.theta));
Ynew = 1/P;

% Ya      = L'*v.Y;
% Ia      = L'*v.I;
% Ka      = L'*sf(:,1);
% ACa     = L'*v.AC;
% Ca      = Ya - Ia - ACa;     
% p       = 1/Ca;

%% Pack-up output
eq.v    = v;
eq.c    = c;
eq.P    = P;
eq.Y    = Ynew;
eq.L    = L;
eq.Q    = Q;

% eq.Ya   = Ya;
% eq.Ca   = Ca;
% eq.Ia   = Ia;
% eq.ACa  = ACa;
% eq.Ka   = Ka;


end

