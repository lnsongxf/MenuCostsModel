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
% ns          = size(s,1);

%% B. Compute equilibrium objects that depend on p
% ----- None -----

%% Initialise guesses
cKold       = zeros(glob.Ns,1);
cCold       = zeros(glob.Ns,1);
cEold       = zeros(glob.Ns,1);
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
glob.Phi_A      = glob.Phi_Af; 
glob.Phi        = glob.Phif; 
glob.Phiprime   = glob.Phiprimef; 
v               = solve_valfunc_noagg(c,sf,Y,param,glob,options,1);

% Compute stationary distribution
pPdist           = min(v.pPdist,max(pPgrid))*(1/exp(param.mu));
fspaceergpP      = fundef({'spli',glob.pPgridf,0,1});
QpK              = funbas(fspaceergpP,pPdist(1:end/2));
QpC              = funbas(fspaceergpP,pPdist(end/2+1:end));
QA              = glob.QA;
QK               = dprod(QA,QpK);
QC               = dprod(QA,QpC);
Q              = [param.lambda*QK, (1-param.lambda)*QK; (1-param.lambda)*QC, param.lambda*QC];


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

%% Compute aggregates and implied p
% P = ( L'*(pPdist).^(1-param.theta) )^(1/(1-param.theta));
% Ynew = 1/P;

gYfun = ( L'*(pPdist).^(1-param.theta) )^(1/(1-param.theta)) - 1;


%% Pack-up output
eq.v    = v;
eq.c    = c;
% eq.P    = P;
% eq.Y    = Ynew;
eq.L    = L;
eq.Q    = Q;
eq.gYfun = gYfun;

%% Plot stationary distribution
JpP  = numel(glob.pPgridf);
Ja  = numel(glob.agridf);
La  = kron(eye(Ja),ones(1,2*JpP))*L;
LpP  = kron(ones(1,Ja),eye(2*JpP))*L;
eq.LpP   = LpP;
eq.La   = La;
    
if strcmp(options.plotSD,'Y');
    H = figure(options.fignum);
    %     set(H,'Pos',[1          35        1920         964]);
    
    % Marginal prices
    subplot(2,2,1);
    plot(glob.pPgridf,LpP,'o-');title('Stationary Real Price Dist - LpP');
    grid on;
    % Marginal productivity
    subplot(2,2,2);
    plot(exp(glob.agridf),La,'o-');title('Stationary Prod Dist - La');
    grid on;

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


%% Plot value functions, other policy functions 

if strcmp(options.plotpolicyfun,'Y')
    cK       = c(1:glob.Ns);
    cC       = c(glob.Ns+1:2*glob.Ns);
    
    valK = funbas(glob.fspace,glob.sf)*cK;
    valC = funbas(glob.fspace,glob.sf)*cC;
    valtot = param.lambda*valK + (param.lambda-param.lambda)*valC;
    valtot = reshape(valtot, length(glob.pPgridf), length(glob.agridf));

    pPdistK = reshape(v.pPdist(1:end/2), length(glob.pPgridf), length(glob.agridf));
    pPdistC = reshape(v.pPdist(end/2+1:end), length(glob.pPgridf), length(glob.agridf));
%     ind    = reshape(v.ind, length(glob.pPgridf), length(glob.agridf));
    ystarK  = reshape(v.ystar(1:end/2), length(glob.pPgridf), length(glob.agridf));
    ystarC  = reshape(v.ystar(end/2+1:end), length(glob.pPgridf), length(glob.agridf));
    nstarK  = reshape(v.nstar(1:end/2), length(glob.pPgridf), length(glob.agridf)); 
    nstarC  = reshape(v.nstar(end/2+1:end), length(glob.pPgridf), length(glob.agridf)); 
    
    figure(111)
    subplot(2,2,1)
    plot(glob.pPgridf, valtot(:,3),'linewidth',3,'color','b')
%     grid on
%     ylim([51.4,51.8])
    legend('V = \lambda V^K + (1-\lambda)V^C','location','SouthEast')
    xlabel('Ex-ante real price','fontsize',options.fontsize)
    ylabel('Value function','fontsize',options.fontsize)
    set(gca, 'fontsize', options.fontsize)
    subplot(2,2,2)
    plot(glob.pPgridf, pPdistK(:,3),'linewidth',3,'color','b')
    hold on
    plot(glob.pPgridf, pPdistC(:,3),'linewidth',3,'color','r','linestyle','--')
    legend('Keep price','Change price','location','SouthEast')
%     grid on
    xlabel('Ex-ante real price','fontsize',options.fontsize)
    ylabel('Ex-post real price','fontsize',options.fontsize)
    set(gca, 'fontsize', options.fontsize)
    subplot(2,2,3)
    plot(glob.pPgridf, ystarK(:,3),'linewidth',3,'color','b')
    legend('Keep price','Change price','location','NorthEast')
    hold on
    plot(glob.pPgridf, ystarC(:,3),'linewidth',3,'color','r','linestyle','--')
%     grid on
    xlabel('Ex-ante real price','fontsize',options.fontsize)
    ylabel('Output','fontsize',options.fontsize)
    set(gca, 'fontsize', options.fontsize)
    subplot(2,2,4)
    plot(glob.pPgridf, nstarK(:,3),'linewidth',3,'color','b')
    hold on
    plot(glob.pPgridf, nstarC(:,3),'linewidth',3,'color','r','linestyle','--')
    legend('Keep price','Change price','location','NorthEast')
%     grid on
    xlabel('Ex-ante real price','fontsize',options.fontsize)
    ylabel('Labour demand','fontsize',options.fontsize)
    set(gca, 'fontsize', options.fontsize)
    
end



end

