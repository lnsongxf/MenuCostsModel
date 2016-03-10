%--------------------------------------------------------
%
% This file computes the firm pricing model with menu costs from Terry &
% Knotek (2008)
%
% James Graham
% 2/18/2016
%
%---------------------------------------------------------
% Add CompEcon package
% p = genpath('E:\Dropbox\Economics\Matlab_codes\CompEcon');
p = genpath('C:\Users\James\Dropbox\Economics\Matlab_codes\CompEcon');
addpath(p);

% cd('E:\Dropbox\Economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\MenuCostsModel\DiscretizedModel')
cd('C:\Users\James\Dropbox\economics\2015_2016_material\AdvMacro_Midrigan\TermPaper\MenuCostsModel\DiscretizedModel')
close all
clear
clc


%% parameters

parms.beta      = 0.99;     % discount factor
parms.delta     = 0.3;      % relative weight of labor-to-consumption in utility
parms.sigma     = 1;        % risk aversion coefficient
parms.phi       = 0.5;      % inveser labour supply elasticity
parms.theta     = 5;        % elasticity of substitution
parms.alpha     = 2/3;      % returns to labour
parms.rhoa      = 0.35;     % persistence of productivity
parms.sigmazeta = 0.225;    % stddev of productivity shocks
parms.Phi       = 0.156;    % menu cost in labour units
parms.mu        = 0.006;    % s.s money growth
parms.rhom       = 0.37;     % persistence of money growth
parms.sigmaeps  = 0.0048;   % stddev of money growth shocks
parms.tauc      = 0.005;    % tolerance for forecasting rule
parms.n         = 5000;     % number of firms
parms.T         = 96;       % simulation length
parms.S         = 25;       % simulations for computing forecasting coeffs
parms.s         = 100;      % simulations for moment computations

%% Create gridspace
% Number of grids for each variable
parms.Npp = 51;  %349;
parms.Na  = 5;
parms.Ny  = 5;
parms.Ndm = 5;
parms.Npi = 5;

% Real price grid
parms.pPmin = 0.75;
parms.pPmax = 1.25;
parms.pPgrid = linspace(parms.pPmin,parms.pPmax,parms.Npp);

% Output, money growth, and inflation grids
% NOTE: construct a VAR, then used VAR-Tauchen to construct grids. Note
% that for many parameterizations inflation is trending down...
b0 = 0.001;
b1 = 0.3;
b2 = 0.25;

% %Productivity grid
[parms.trans_a, lnagrid] = tauchen(parms.Na,parms.rhoa,parms.sigmazeta,0);
parms.agrid = exp(lnagrid);   % Take [ln(a_t)] back to non-logged vals 

Nvar = [parms.Ndm;...
    parms.Ny; ...
    parms.Npi];
muvar = [parms.mu*(1-parms.rhom); ...
    b0 + b2*parms.mu*(1-parms.rhom); ...
    -b0 + (1-b2)*parms.mu*(1-parms.rhom)];
Avar = [parms.rhom, 0, 0; ...
    b2*parms.rhom, b1, 0; ...
    (1-b2)*parms.rhom, (1-b1), 0];
Svar = [parms.sigmaeps; ...
    b2*parms.sigmaeps; ...
    (1-b2)*parms.sigmaeps];
[agggrid,parms.trans_agg,~]=tauchenvar(Nvar,muvar,Avar,Svar);
parms.agggrid = exp(agggrid);   % Take [D(ln(M_t)), ln(Y_t), ln(P_t/P_t-1)] back to non-logged vals  

parms.grid = [repmat(parms.agrid,1,parms.Ndm*parms.Ny*parms.Npi); kron(parms.agggrid,ones(1,parms.Na)) ];
Trans_admypi = kron(parms.trans_agg,parms.trans_a);

% Normalize so rows sum to one
s=size(Trans_admypi, 2);
Trans_admypi = spdiags(sum(Trans_admypi,2), 0, s, s)\Trans_admypi;
parms.trans = Trans_admypi;   % Make sure transiton probs are Prob(t->t+1)

%% Value function iteration
%--------------------------
% Need to speed up an insane amount if gonna get up to 349 grid points for
% real prices...
%
% Get rid of loops, vectorize as much as possible...
%
%--------------------------
tic 
[V_final, Vc_final, Vk_final] = val_iter(parms, 0.1);
toc

% Plot value function
idx = ((parms.Na*parms.Ny*parms.Ndm*parms.Npi) - 1)/2 + 1;

figure
plot( parms.pPgrid,squeeze(V_final(idx,:)) )
xlabel('Real price')
ylabel('Value')
title('Value function at mean value of state variables')


%% Compute stationary distribution
dmss = median(parms.grid(2,:));
Yss = median(parms.grid(3,:));
piss = median(parms.grid(4,:));

pPfun = @(pP,a) pricefunc(parms,pP,a,dmss,Yss,piss,Vk_final,Vc_final,V_final);
[stat_density, stat_dist_pPgrid] = statdist_eigen( parms, pPfun );


%% Krussel Smith Step...
% TO DO: simulate aggregate states, OLS step, find new coefficients, put
% into larger loop around the value function iteration step. Create a
% markov chain simulation stand-alone code.

% % % % First part of the Krussel Smith step


% 1. Guess P-1 = 1, M-1 = 1
% 2. Guess Y~
% 4. With dM and Y~ --> pi~
% 5. P~ = P-1*pi~

T = 50;         % simulation periods
F = 500;        % Number of firms
damp = 0.5;    % Dampening parameter on updating for Y

dM_sim = nan(T,1);
dM_sim(1) = dmss;   % initial money growth is steady state money growth

M_sim = nan(T+1,1);
M_sim(1) = 1;
M_sim(2) = dM_sim(1)*M_sim(1);

% To simulate dM, need to simulate indexes since dM is part of a VAR
agg_idx = nan(T+1,1);
agg_idx(1) = (parms.Ndm*parms.Ny*parms.Npi + 1)/2;
indexes = 1:1:parms.Ndm*parms.Ny*parms.Npi;

% 3. Simulate (dM, aj) for all t and j
a_sim = nan(T+1,F);
a_sim(1,:) = median(parms.agrid);       % In pre-history, all firms at steady state productivity
for t = 1: T
    agg_idx(t+1) = randsample(indexes,1,true,parms.trans_agg(agg_idx(t),:));
    dM_sim(t+1) = parms.agggrid(1,agg_idx(t+1));
    for f = 1:F
        a_idx = find(a_sim(t,f) == parms.agrid);
        a_sim(t+1,f) = randsample(parms.agrid,1,true,parms.trans_a(a_idx,:));
        %     pPout = newprice(parms,pP,statenum,Vk_final,Vc_final,V_final);
    end        
end

Ytilde_sim = nan(T+1,1);
Ytilde_sim(1) = Yss;        % Y in steady state in 'pre-history' period
Ytilde_sim(2) = Yss;          % Y guess for first period of history

Y_sim = nan(T+1,1);
Y_sim(1) = Yss;        % Y in steady state in 'pre-history' period

Ptilde_sim = nan(T+1,1);
Ptilde_sim(1) = 1;
pitilde_sim = nan(T,1);
P_sim = nan(T+1,1);
P_sim(1) = 1;              % P in the 'pre-history' period

pP_sim = nan(T,F);
pPdist = nan(T+1,F);
pPdist(1,:) = randsample(stat_dist_pPgrid,F,true,stat_density);
pdist = nan(T+1,F);

%% UPDATING FROM HERE
Y_sim(2) = Ytilde_sim(2) + 0.5;
while Y_sim(2) - Ytilde_sim(2) > 0.1


% Implied inflation rate, for some reason not on grids!
pitilde_sim(1) = dM_sim(1)*Ytilde_sim(2)/Ytilde_sim(1);     
% Put back onto grid... NOTE: THEN NOT EXACTLY IMPLIED BY DM, Y
tmp = abs( parms.grid(4,:) - pitilde_sim(1) );   % how far from grid?
[~, idx] = min(tmp);     % index of closest values   
pitilde_sim(1) = parms.grid(4,idx);       % Put inflation on the grid

Ptilde_sim(2) = Ptilde_sim(1)*pitilde_sim(1);               % Implied initial price level


% 6. Use pricing function to get the distribution for prices via pj =
% (pj/P)*P~  (Later on, use non-stochastic simulation to get the
% distribution on each loop exactly - much faster)

% pPdist = real price distribution, pdist = nominal price distribution

pdist(1,:) = pPdist(1,:)*P_sim(1);      % multiply by actual price level, known in "pre-history"

for f = 1:F
    pPdist(2,f) = pricefunc(parms,pPdist(1,f),a_sim(2,f),dM_sim(1),Ytilde_sim(2),pitilde_sim(1),Vk_final,Vc_final,V_final); 
end
pdist(2,:) = pPdist(2,:)*Ptilde_sim(2); % multiply by implied price level (given by implied inflation rate and pre-history price)

% 7. Get the actual price level: P = [ int_0^1 ((pj/P)P~)^(1-theta)
% dj]^(1/(1-theta))

P_sim(2) = (sum(pdist(2,:).^(1-parms.theta)))^(1/(1-parms.theta));  % Quite far guessed true P...


% 8. Note that Y = C = M/P. Check if Y~ = Y . If not, update Y~ using
% dampening.

Y_sim(2) = M_sim(2)/P_sim(2);

disp(['Diff in outputs is ' num2str(Y_sim(2)-Ytilde_sim(2))])
disp(['At computed price level ' num2str( P_sim(2) )])


% If Y is too low, increase Y and increase P
if Y_sim(2) - Ytilde_sim(2) > 0.1
   Ytilde_sim(2) =  damp*Ytilde_sim(2) + (1-damp)*Y_sim(2); 
   Ptilde_sim(1) = 1.1*Ptilde_sim(1);
   P_sim(1) = 1.1*P_sim(1);
   
   tmp = abs( parms.grid(3,:) - Ytilde_sim(2) );   % how far from grid?
   [~, idx] = min(tmp);     % index of closest values
   Ytilde_sim(2) = parms.grid(3,idx);       % Put inflation on the grid
end

disp(['New initial price level is ' num2str( P_sim(1) )])

   
end



%%
%{
% Compute values at each (pP,a) in the steady state
for i = 1:parms.Npp
    for j = 1:parms.Na
        pPss(i,j) = pPfun(parms.pPgrid(i),parms.grid(1,j));        
    end
end

tmp = Yss( (reshape(pPss,1,parms.Npp*parms.Na)).^(1-parms.theta)*stat_density')^(1/(1-parms.theta));


a_hat(a,y)
% Compute assets
% stat_density = W(:,index)'/sum(W(:,index));
A(rt) = 0;
A(rt) = reshape(a_hat,1,M*N)*stat_density';

% Compute capital
K(rt) = Lss*( (1/alpha)*(r(rt) + delta) )^(1/(alpha - 1));

% Compare aggregate assets to capital stock at r(rt)
AKdiff(rt+1) = A(rt) - K(rt);

% Update interest rate
weight = 0.99;  
weightextreme = 0.95;

if (alpha*A(rt)^(alpha-1)*Lss^(1-alpha) - delta) > 1/beta - 1
    r(rt+1) = weightextreme*r(rt) + (1-weightextreme)*(1/beta - 1);
    bound = 1;
else
    r(rt+1) = weight*r(rt) + ...
        (1-weight)*(alpha*A(rt)^(alpha-1)*Lss^(1-alpha) - delta);
    bound = 0;
end

% Update wage
w(rt+1) = (1-alpha)*(1/alpha*(r(rt+1) + delta))^(alpha/(alpha-1));

disp('%----------------------------------%')
disp('')
disp(['The current interest rate is: ' num2str(r(rt))])
disp(['The assets are: ' num2str(abs(A(rt)))])
disp(['The capital is: ' num2str(abs(K(rt)))])
if bound == 1
    disp(['UPPER BOUND VIOLATED AT NEW GUESS'])
end
disp(['The new interest rate is: ' num2str(r(rt+1))])
disp('')
disp('%----------------------------------%')


rt = rt+1;

% end
% toc
%}





