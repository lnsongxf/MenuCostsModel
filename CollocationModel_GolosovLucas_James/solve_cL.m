function eq  = solve_cL(Y,param,glob,options)
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
%-------------------------------------------------


%% A. Globals 
s           = glob.s; 
sf          = glob.sf;
pPgrid      = glob.pPgrid;
ns          = size(s,1);

%% B. Compute equilibrium objects that depend on p
% ----- None -----

%% Initialise guesses
cKold       = zeros(ns,1);
cCold       = zeros(ns,1);
cEold       = zeros(ns,1);
cold        = [cKold;cCold;cEold];

% Check if previous solution exists
if exist('glob.c','var')
        cold = glob.c;
end

totaltic    = tic;
%% Bellman iteration
for citer = (1:options.Nbell)
    glob.citer  = citer;
    % 1. Compute values;
    v           = solve_valfunc_noagg(cold,s,Y,param,glob,options); 
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
    [v,jac]     = solve_valfunc_noagg(cold,s,Y,param,glob,options);
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


%% Compute stationary distribution
% Solve again on a finer grid for pP
glob.Phi_V      = glob.Phi_Vf; 
glob.Phi        = glob.Phif; 
glob.Phiprime   = glob.Phiprimef; 
v               = solve_valfunc_noagg(c,sf,Y,param,glob,options,1);

% Compute stationary distribution
pPdist           = min(v.pPdist,max(pPgrid)).*(1/exp(glob.piw));
fspaceergpP      = fundef({'spli',glob.pPgridf,0,1});
QpP              = funbas(fspaceergpP,pPdist);
QV              = glob.QV;
Q               = dprod(QV,QpP);

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
    H = figure(options.fignum);
    %     set(H,'Pos',[1          35        1920         964]);
    JpP  = numel(glob.pPgridf);
    Jv  = numel(glob.vgridf);
    Lv  = kron(eye(Jv),ones(1,JpP))*L;
    LpP  = kron(ones(1,Jv),eye(JpP))*L;
    % Marginal prices
    subplot(2,2,1);
    plot(glob.pPgridf,LpP,'o-');title('Stationary Real Price Dist - LpP');
    grid on;
    % Marginal productivity
    subplot(2,2,2);
    plot(exp(glob.vgridf),Lv,'o-');title('Stationary Prod Dist - Lv');
    grid on;
    eq.LpP   = LpP;
    eq.Lv   = Lv;
    % Joint (pP,A) - Surface plot
    subplot(2,2,4)
    Lmat    = reshape(L,JpP,Jv);
    Vmat    = repmat(glob.vgridf',JpP,1);
    pPmat    = repmat(glob.pPgridf,1,Jv);
    if LpP(end)<0.001;
        pPub     = glob.pPgridf(find(cumsum(LpP)>0.98,1,'first'));
    else
        pPub     = max(glob.pPgridf);
    end
    Vub     = glob.vgridf(find(cumsum(Lv)>0.98,1,'first'));
    mesh(Vmat,pPmat,Lmat,'LineWidth',2);
    xlabel('Productivity - v');
    ylabel('Real Price - pP');
    title('Joint Distribution');
    xlim([min(glob.vgridf),Vub]);
    ylim([min(glob.pPgridf),pPub]);
    zlim([0,max(max(Lmat))]);
    
end

% Compute implied Y (from the stationary agg consumption equation 17/21 in GL(2007))
Ynew = ( param.alpha^(1 - param.epsilon)*(pPdist'.^(1 - param.epsilon)*L) )^...
    (1/(param.gamma*(param.epsilon - 1)));

%% Pack-up output
eq.v    = v;
eq.c    = c;
% eq.P    = P;
eq.Y    = Ynew;
eq.L    = L;
eq.Q    = Q;

%% Plot value functions, other policy functions 

if strcmp(options.plotpolicyfun,'Y')
    cK       = c(1:glob.Ns);
    cC       = c(glob.Ns+1:2*glob.Ns);
    
    valK = funbas(glob.fspace,glob.sf)*cK;
    valC = funbas(glob.fspace,glob.sf)*cC;
    valtot = max(valK, valC);
    valtot = reshape(valtot, length(glob.pPgridf), length(glob.vgridf));
        
    pPdist = reshape(v.pPdist, length(glob.pPgridf), length(glob.vgridf));
    ind    = reshape(v.ind, length(glob.pPgridf), length(glob.vgridf));
%     ystar  = reshape(v.ystar, length(glob.pPgridf), length(glob.vgridf));
%     nstar  = reshape(v.nstar, length(glob.pPgridf), length(glob.vgridf)); 
%     wPstar = reshape(v.wPstar, length(glob.pPgridf), length(glob.vgridf));
    
    figure
    plot(glob.pPgridf, valtot)
    xlabel('Real price','fontsize',options.fontsize)
    ylabel('Value','fontsize',options.fontsize)
    set(gca, 'fontsize', options.fontsize)
    legend('v_1','v_2','v_3','v_4','v_5')

% 
%     figure
%     subplot(2,3,1)
%     plot(glob.pPgridf, valtot)
%     xlabel('Real price','fontsize',options.fontsize)
%     ylabel('Value','fontsize',options.fontsize)
%     set(gca, 'fontsize', options.fontsize)
%     legend('v_1','v_2','v_3','v_4','v_5')
%     subplot(2,3,2)
%     plot(glob.pPgridf, pPdist)
%     xlabel('Real price','fontsize',options.fontsize)
%     ylabel('Observed real price','fontsize',options.fontsize)
%     set(gca, 'fontsize', options.fontsize)
%     subplot(2,3,3)
%     plot(glob.pPgridf, ystar)
%     xlabel('Real price','fontsize',options.fontsize)
%     ylabel('Output','fontsize',options.fontsize)
%     set(gca, 'fontsize', options.fontsize)
%     subplot(2,3,4)
%     plot(glob.pPgridf, nstar)
%     xlabel('Real price','fontsize',options.fontsize)
%     ylabel('Labour demand','fontsize',options.fontsize)
%     set(gca, 'fontsize', options.fontsize)
%     subplot(2,3,5)
%     plot(glob.pPgridf, wPstar)
%     xlabel('Real price','fontsize',options.fontsize)
%     ylabel('Real wage','fontsize',options.fontsize)
%     set(gca, 'fontsize', options.fontsize)    
end



end

